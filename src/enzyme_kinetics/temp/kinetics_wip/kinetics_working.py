from __future__ import annotations

import logging
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Callable, Final, Self

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from with_docs.kgdlibs.pathtools import ExportPathBuilder, to_absolute_path
from with_docs.plots.config import COUNT, filter_stats_params, MEAN, PEAK_ID, PROTEIN, REL_ACT, STD
from with_docs.process.analyze import COLUMN_NAME_MAP, normalize_peak_id
from with_docs.process.compounds import COMPOUND_PREFIXES
from with_docs.process.io import read_csv_peaks
from with_docs.process.normals import require_columns

_logger = logging.getLogger(__name__)

# Column labels for readability
SUBSTRATE_CONC: Final[str] = "substrate_conc"  # NEW!


@dataclass(frozen=True, slots=True)
class KineticsConfig:
    # E_t = protein_conc_uM
    prot_conc: float = 12.25  # protein concentration in units of uM
    rxn_time: float = 180.0  # in units of minutes
    # Import data
    peak_attr: str = "area"  # Peak attribute to plot. Options: "area", "height"
    detect_prefix: bool = False
    peak_prefix: frozenset[str] = field(default_factory=frozenset)
    # Export data
    overwrite: bool = False
    merge: bool = True
    verbose: bool = False

    def __post_init__(self) -> None:
        """Merges caller-supplied prefixes and include_cols with their library-wide defaults."""
        prefixes = self.peak_prefix | frozenset(COMPOUND_PREFIXES)
        object.__setattr__(self, "peak_prefix", prefixes)

    def load(self, path: Path, **kwargs: Any) -> pd.DataFrame:
        df = read_csv_peaks(path, convert_to_float32=False, verbose=self.verbose, **kwargs)
        df.rename(columns=COLUMN_NAME_MAP, inplace=True)

        # TODO: This may be unnecessary.
        # if PROTEIN in df.columns:
        #     # TODO: Verify this approach works. Unsure if it also creates a column "substrate conc"
        #     df = _clean_protein_column(df)
        # df = map_columns(df, [PROTEIN, PEAK_ID], reindex=False)

        # Column filtering and validation
        df = filter_stats_params(df, self.peak_attr)
        require_columns(df, (SUBSTRATE_CONC, PEAK_ID, MEAN, STD, COUNT, REL_ACT))

        # Memory optimization and normalization
        prefixes = tuple(self.peak_prefix) if self.detect_prefix else None
        df = normalize_peak_id(df, prefix=prefixes)
        return df

    def export(self, dest: Path, tag: str) -> Path:
        """Builds the export path and returns it.
        tag: string to prepend to the filename
        """
        base_name = f"{tag}_kinetics_{self.peak_attr}.png"
        build_path = ExportPathBuilder(
            base_name, overwrite=self.overwrite, merge=self.merge,
        ).interpolate()
        return build_path.path

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Creates an AnalyzerConfig from a dictionary, ignoring unknown keys."""
        valid_keys = {f.name for f in fields(cls)}
        filtered_data = {k: v for k, v in data.items() if k in valid_keys}
        return cls(**filtered_data)


def _clean_protein_column(df: pd.DataFrame) -> pd.DataFrame:
    """Clean the protein column by removing non-numeric characters and converting to float."""
    df[SUBSTRATE_CONC] = df[PROTEIN].astype(str).str.replace(r"[a-zA-Z]+", "", regex=True).astype(float)
    return df


class KineticsModels:
    """Mathematical models for enzyme kinetics."""

    @staticmethod
    def michaelis_menten(s: np.ndarray, vmax: float, km: float) -> np.ndarray:
        """Calculate the velocity term in the Michaelis–Menten kinetics equation."""
        return (vmax * s) / (km + s)
        # return (vmax * substrate_concentration) / (km + substrate_concentration)

    @staticmethod
    def hill_equation(s: np.ndarray, vmax: float, k_half: float, n: float) -> np.ndarray:
        """Hill equation for sigmoidal/cooperative kinetics."""
        return (vmax * s ** n) / (k_half ** n + s ** n)

    @staticmethod
    def threshold_michaelis_menten(s: np.ndarray, vmax: float, km: float, s0: float) -> np.ndarray:
        """Michaelis-Menten model with a dead-zone threshold (S0).
        The threshold model: v=0 if S <= S0, otherwise standard MM.
        """
        return np.where(s > s0, (vmax * (s - s0)) / (km + (s - s0)), 0)


@dataclass(frozen=True, slots=True)
class KineticParams:
    peak_id: str
    vmax: float  # area per minute // need a product standard curve to convert to uM per minute
    vmax_se: float
    km: float  # in millimolar (same as substrate concentration)
    km_se: float
    kcat: float  # kcat per minute
    kcat_se: float
    kcat_km: float  # catalytic efficiency
    popt: tuple[float, float, float]
    model_func: Callable
    extra: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        """Provides a flat dictionary suitable for pandas/CSV serialization."""
        base = {
            "peak_id": self.peak_id,
            "vmax": self.vmax, "vmax_se": self.vmax_se,
            "km": self.km, "km_se": self.km_se,
            "kcat": self.kcat, "kcat_se": self.kcat_se,
            "kcat_km": self.kcat_km,
        }
        if self.extra:
            # Only serialize scalar types for the CSV
            base.update({k: v for k, v in self.extra.items() if isinstance(v, (int, float, str))})
        return base

    """"
    vmax : area/ minutes ~ mM/min (velocity, concentration per time)
    km : mM (concentration)
    kcat : min-1 (events per time)
    kcat/km : area * uM-1 * min-1 * mM-1 ~ (mM-1 * min-1) or (M-1 * s-1)
    """
    """Calculate the kinetic parameters kcat and kcat/km.

                Variables:
                    vmax = units of mM/min
                    km = units of mM (per 3 hr)
                    km_um = units of uM (per 1 hr)
                    kcat = units of 1/min
                    kcat/km = units of 1/s * 1/M
            """


class EnzymeKineticsAnalysis:
    """Class to process, fit, and plot enzyme kinetics data."""

    def __init__(
        self,
        path: str | Path,
        target_dir: str | Path,
        *,
        config: KineticsConfig | None = None,
    ) -> None:
        self.path = to_absolute_path(path)
        self.dest = to_absolute_path(target_dir)
        self.config = config or KineticsConfig()

        self.df: pd.DataFrame = pd.DataFrame()
        self.stats_df: pd.DataFrame = pd.DataFrame()  # Replaces undefined stats_df in export
        self.fits: dict[str, KineticParams] = {}

    def _prepare_data(self, peak_id: str) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
        """Filters data and calculates velocity/SEM for a specific peak."""
        sub_df = self.df[self.df[PEAK_ID] == peak_id].dropna(subset=[SUBSTRATE_CONC, MEAN])
        if sub_df.empty or sub_df[MEAN].max() == 0:
            return None

        s = sub_df[SUBSTRATE_CONC].values
        v = sub_df[MEAN].values / self.config.rxn_time  # Velocity in Area/min
        # Calculate Standard Error of the Mean (SEM) for the points
        v_sem = (sub_df[STD].values / self.config.rxn_time) / np.sqrt(sub_df[COUNT].values)
        return s, v, v_sem

    def fit_all(self, **kwargs: Any) -> Self:
        """Executes data loading and math fitting logic cleanly without plotting."""
        self.df = self.config.load(self.path, **kwargs)
        peaks = self.df[PEAK_ID].unique()
        self.fits.clear()

        for peak in peaks:
            data = self._prepare_data(peak)
            if not data:
                continue

            s, v, _ = data
            try:
                if peak == 'OLV':
                    model_func = KineticsModels.threshold_michaelis_menten
                    p0 = [np.max(v) * 2, 1.0, 0.5]  # Guesses for vmax, km, s0
                    popt, pcov = curve_fit(model_func, s, v, p0=p0, bounds=(0, np.inf), maxfev=10000)
                    perr = np.sqrt(np.diag(pcov))
                    # Sum of squared residuals (SSR)
                    # Values closer to zero indicate a better fit.
                    ssr = np.sum((v - model_func(s, *popt)) ** 2)
                    _logger.info(f"Fit 1 (Best): S0 = {popt[2]:.3f} | SSR = {ssr:.4e}")

                    # vmax, km, s0 = popt
                    # extra = {'name': 's0', 'val': s0, 'se': np.sqrt(np.diag(pcov))[2]}

                    # Hill Equation (stored as extra for comparison)
                    hill_func = KineticsModels.hill_equation
                    p0_h = [np.max(v), 1.0, 2.0]  # Guesses for vmax, km, s0
                    popt_h, pcov_h = curve_fit(hill_func, s, v, p0=p0_h, bounds=(0, np.inf), maxfev=10000)
                    perr_h = np.sqrt(np.diag(pcov_h))

                    # vmax_h, km_h, s0_h = popt
                    # extra.update(
                    extra = {
                        'S0': popt[2], 'S0_SE': perr[2],
                        'Hill_Model': 'Hill',
                        'Hill_Vmax': popt_h[0], 'Hill_Vmax_SE': perr_h[0],
                        'Hill_n': popt_h[2], 'Hill_n_SE': perr_h[2],
                        'popt_h': tuple(popt_h),  # Retain for plotting
                    }
                else:
                    # Standard MM
                    model_func = KineticsModels.michaelis_menten
                    # 5. Fit the standard model to EVERY peak
                    # Initial guesses: Vmax is the maximum velocity, Km is the concentration at half Vmax
                    vmax_guess = np.max(v)
                    km_guess = s[np.argmin(np.abs(v - vmax_guess / 2))]
                    popt, pcov = curve_fit(
                        model_func, s, v,
                        p0=[vmax_guess, km_guess],
                        bounds=(0, np.inf),
                        maxfev=10000,
                    )
                    perr = np.sqrt(np.diag(pcov))  # Calculate Standard Errors for parameters
                    vmax, km = popt
                    extra = None

                kcat = vmax / self.config.prot_conc  # in units area * micromolar-1 * minutes-1
                kcat_se = perr[0] / self.config.prot_conc  # Propagation of error
                kcat_km = kcat / km if km > 0 else np.nan  # Catalytic efficiency

                self.fits[peak] = KineticParams(
                    peak_id=peak, vmax=popt[0], vmax_se=perr[0], km=popt[1], km_se=perr[1],
                    kcat=kcat, kcat_se=kcat_se, kcat_km=kcat_km, popt=tuple(popt),
                    model_func=model_func, extra=extra,
                )

            except Exception as exc:
                _logger.warning(f"Fit failed for peak {peak}: {exc}")

        # Construct analytical results dataframe
        self.stats_df = pd.DataFrame([fit.to_dict() for fit in self.fits.values()])
        return self

    def plot_results(self, show: bool = False) -> Self:
        """Generates and saves the matplotlib facet plot."""
        peaks = self.df[PEAK_ID].unique()

        # Create faceted plot
        cols = 3
        rows = int(np.ceil(len(peaks) / cols))

        # Memory explicit management
        fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
        axes = axes.flatten()

        for i, peak in enumerate(peaks):
            ax = axes[i]
            data = self._prepare_data(peak)

            if not data:
                ax.set_title(f"{peak} - No Data/Zero Velocity")
                ax.axis('off')
                continue

            s, v, v_sem = data
            ax.errorbar(s, v, yerr=v_sem, fmt='o', label='Data (Mean ± SEM)', color='blue', alpha=0.7, capsize=3)

            fit = self.fits.get(peak)
            if fit:
                s_fit = np.linspace(0, np.max(s) * 1.1, 200)

                if peak == 'OLV' and fit.extra:
                    label_text = (
                        f"Threshold Fit:\n"
                        f"$V_{{max}}$={fit.vmax:.2f}±{fit.vmax_se:.2f}\n"
                        f"$K_m$={fit.km:.2f}±{fit.km_se:.2f}\n"
                        f"$S_0$={fit.extra['S0']:.2f}\n"
                        f"$k_{{cat}}$={fit.kcat:.4f}"
                    )
                    ax.plot(
                        s_fit, KineticsModels.hill_equation(s_fit, *fit.extra['popt_h']),
                        label=f"Hill Fit (n={fit.extra['Hill_n']:.2f})", linestyle='-', color='blue',
                    )
                else:
                    label_text = (
                        f"MM Fit:\n"
                        f"$V_{{max}}$={fit.vmax:.2f}±{fit.vmax_se:.2f}\n"
                        f"$K_m$={fit.km:.2f}±{fit.km_se:.2f}\n"
                        f"$k_{{cat}}$={fit.kcat:.4f}"
                    )

                ax.plot(s_fit, fit.model_func(s_fit, *fit.popt), label=label_text, linestyle='--', color='red')
                ax.legend(fontsize=9, loc='lower right')

                # s_fit = np.linspace(0, np.max(s) * 1.1, 300)
                #
                # if fit.extra['Model'] == 'Hill':
                #     ax.plot(
                #         s_fit,
                #         fit.model_func(s_fit, *fit.extra['popt_h']),
                #         label=label_text,
                #         linestyle='-',
                #         color='blue',
                #         label=f'Hill Fit (n={popt_h[2]:.2f})',
                #     )

            else:
                ax.text(0.5, 0.5, 'Fit Failed', ha='center', va='center', transform=ax.transAxes)

            #     label_text = (
            #         f"Fit:\n"
            #         f"$V_{{max}}$: {fit.vmax:.2f}±{fit.vmax_se:.2f}\n"
            #         f"$K_m$: {fit.km:.2f}±{fit.km_se:.2f}"
            #         f"$k_{{cat}}$={fit.kcat:.4f}"
            #     )
            #
            #     s_fit = np.linspace(0, np.max(s) * 1.1, 300)
            #
            # ax.plot(
            #     s_fit, fit.model_func(s_fit, *fit.popt), label=label_text, linestyle='--', color='red',
            # )
            # ax.legend(fontsize=9, loc='lower right')
            #
            # if fit.vmax == fit.km == 0:
            #     ax.text(0.5, 0.5, 'Fit Failed', ha='center', va='center', transform=ax.transAxes)
            #
            # # Cleanup for CSV export
            # row = {k: v for k, v in fit.items() if k not in ['popt', 'model_func', 'extra']}
            # if fit.extra:
            #     row[fit.extra['name']] = fit.extra['val']
            #     row[f"{fit.extra['name']}_SE"] = fit.extra['se']
            # final_table_rows.append(row)

            ax.set_title(f'Michaelis-Menten Kinetics: {peak}', fontweight='bold')
            ax.set_xlabel('Substrate Concentration [HexCoA] (μM)')
            # ax.set_xlabel('Substrate Concentration [HexCoA] (μM)')
            ax.set_ylabel('Velocity (Area / min)')

        # Clean up empty subplots
        for j in range(len(peaks), len(axes)):
            axes[j].axis('off')

        plt.grid(True, linestyle=':', alpha=0.6)
        plt.tight_layout()

        # plot_name = 'kinetics_plot.png'
        # csv_name = 'kinetics_results.csv'
        # plt.savefig(plot_name, dpi=300)

        plot_path = self.dest / f"{Path(self.path).stem}_plot.png"
        plt.savefig(plot_path, dpi=300)

        # results_df = pd.DataFrame(final_table_rows)
        # results_df.to_csv(csv_name, index=False)
        # return results_df

        if show:
            plt.show()

        plt.close(fig)  # Prevent memory leaks
        _logger.info(f"Saved kinetics plot to: {plot_path}")
        return self

    def run(self, **kwargs: Any) -> Self:
        """Combined chain operation matching the old API."""
        self.fit_all(**kwargs)
        self.plot_results()
        self.export()
        return self

    def export(self) -> Self:
        # todo: export the original df after protein has been split into substrate conc and dna conc columns
        if self.stats_df.empty:
            _logger.warning("Export: No kinetics data to export.")
            return self

        csv_path = self.config.export(self.dest, tag=self.path.stem).with_suffix(".csv")
        self.stats_df.to_csv(csv_path, index=False, index_label="index")
        _logger.info(
            "Saved enzyme kinetics data (shape=%s) to: %s",
            self.stats_df.shape, csv_path,
        )

        return self


if __name__ == "__main__":
    # --- Usage Example ---
    df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")
    analyzer = EnzymeKineticsAnalysis(df, protein_conc_uM=12.25, time_min=180)
    results = analyzer.run()

    # 1. Load data
    df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")

    # 2. Define constants
    E_t = 12.25  # Enzyme concentration in microMolar
    time_min = 180  # 3 hour reaction time
