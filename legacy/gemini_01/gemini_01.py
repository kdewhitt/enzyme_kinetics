from __future__ import annotations

import logging
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Final, Self
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from kgdlibs.pathtools import to_absolute_path
from process.analyze import COLUMN_NAME_MAP, normalize_peak_id
from process.compounds import COMPOUND_PREFIXES
from process.io import read_csv_peaks
from process.normals import require_columns
from shimadzu.peaks import ExtraPeakColumns, PeakColumns

# Mocking ExportPathBuilder assuming it exists in your ecosystem
try:
    from process.io import ExportPathBuilder
except ImportError:
    pass

_logger = logging.getLogger(__name__)

# Column labels for readability
AREA: Final[str] = PeakColumns.AREA.snake_case
HEIGHT: Final[str] = PeakColumns.HEIGHT.snake_case
PEAK_ID: Final[str] = ExtraPeakColumns.PEAK_ID.value
PROTEIN: Final[str] = ExtraPeakColumns.PROTEIN.value
SUBSTRATE_CONC: Final[str] = "substrate_conc"

MEAN = "mean"
STD = "std"
COUNT = "count"
REL_ACTIVITY = "rel_activity"

SIMPLE_STATS_MAP = {
    "area_mean": MEAN,
    "area_std": STD,
    "area_count": COUNT,
    "area_rel": REL_ACTIVITY,
    "height_mean": MEAN,
    "height_std": STD,
    "height_count": COUNT,
    "height_rel": REL_ACTIVITY,
}


@dataclass(frozen=True, slots=True)
class KineticsConfig:
    prot_conc: float = 12.25  # uM
    rxn_time: float = 180.0  # minutes
    peak_attr: str = "area"
    detect_prefix: bool = False
    peak_prefix: frozenset[str] = field(default_factory=frozenset)
    overwrite: bool = False
    merge: bool = True
    verbose: bool = False
    scalar: float = 1.0
    convert_to_float32: bool = False
    include_cols: tuple[str, ...] = tuple()

    def __post_init__(self) -> None:
        prefixes = self.peak_prefix | frozenset(COMPOUND_PREFIXES)
        object.__setattr__(self, "peak_prefix", prefixes)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        valid_keys = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in valid_keys})

    @property
    def csv_args(self) -> dict[str, Any]:
        return {
            "case": "snake",
            "scalar": self.scalar,
            "convert_to_float32": self.convert_to_float32,
            "include_cols": self.include_cols,
            "verbose": self.verbose,
        }


def _filter_rename_columns(df: pd.DataFrame, keep_attr: str) -> pd.DataFrame:
    """Optimized to drop unused columns directly, saving memory."""
    search_str = HEIGHT if keep_attr == "area" else AREA
    cols_to_drop = [c for c in df.columns if search_str in c]
    df = df.drop(columns=cols_to_drop)
    df = df.rename(columns=SIMPLE_STATS_MAP)
    require_columns(df, (MEAN, STD, COUNT, REL_ACTIVITY))
    return df


def _clean_protein_column(df: pd.DataFrame) -> pd.DataFrame:
    df[SUBSTRATE_CONC] = (
        df[PROTEIN].astype(str).str.replace(r"[a-zA-Z]+", "", regex=True).astype(float)
    )
    return df


class KineticsModels:
    @staticmethod
    def michaelis_menten(s: np.ndarray, vmax: float, km: float) -> np.ndarray:
        return (vmax * s) / (km + s)

    @staticmethod
    def hill_equation(
        s: np.ndarray, vmax: float, k_half: float, n: float
    ) -> np.ndarray:
        return (vmax * s**n) / (k_half**n + s**n)

    @staticmethod
    def threshold_michaelis_menten(
        s: np.ndarray, vmax: float, km: float, s0: float
    ) -> np.ndarray:
        return np.where(s > s0, (vmax * (s - s0)) / (km + (s - s0)), 0.0)


@dataclass(frozen=True, slots=True)
class KineticParams:
    peak_id: str
    vmax: float
    vmax_se: float
    km: float
    km_se: float
    kcat: float
    kcat_se: float
    kcat_km: float
    popt: tuple[float, ...]
    model_func: Callable
    extra: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        """Provides a flat dictionary suitable for pandas/CSV serialization."""
        base = {
            "peak_id": self.peak_id,
            "Vmax": self.vmax,
            "Vmax_SE": self.vmax_se,
            "Km": self.km,
            "Km_SE": self.km_se,
            "kcat": self.kcat,
            "kcat_SE": self.kcat_se,
            "kcat_Km": self.kcat_km,
        }
        if self.extra:
            # Only serialize scalar types for the CSV
            base.update(
                {
                    k: v
                    for k, v in self.extra.items()
                    if isinstance(v, (int, float, str))
                }
            )
        return base


class EnzymeKineticsAnalysis:
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
        self.stats_df: pd.DataFrame = (
            pd.DataFrame()
        )  # Replaces undefined stats_df in export
        self.fits: dict[str, KineticParams] = {}

    def _load(self, **kwargs: Any) -> pd.DataFrame:
        df = read_csv_peaks(self.path, **self.config.csv_args, **kwargs)
        df = df.rename(columns=COLUMN_NAME_MAP)

        if PROTEIN in df.columns:
            df = _clean_protein_column(df)

        df = _filter_rename_columns(df, self.config.peak_attr)
        require_columns(df, (SUBSTRATE_CONC, PEAK_ID, MEAN, STD, COUNT, REL_ACTIVITY))

        prefixes = tuple(self.config.peak_prefix) if self.config.detect_prefix else None
        return normalize_peak_id(df, prefix=prefixes)

    def _prepare_data(
        self, peak_id: str
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
        sub_df = self.df[self.df[PEAK_ID] == peak_id].dropna(
            subset=[SUBSTRATE_CONC, MEAN]
        )
        if sub_df.empty or sub_df[MEAN].max() == 0:
            return None

        s = sub_df[SUBSTRATE_CONC].values
        v = sub_df[MEAN].values / self.config.rxn_time
        v_sem = (sub_df[STD].values / self.config.rxn_time) / np.sqrt(
            sub_df[COUNT].values
        )
        return s, v, v_sem

    def fit_all(self, **kwargs: Any) -> Self:
        """Executes data loading and math fitting logic cleanly without plotting."""
        self.df = self._load(**kwargs)
        peaks = self.df[PEAK_ID].unique()
        self.fits.clear()

        for peak in peaks:
            data = self._prepare_data(peak)
            if not data:
                continue

            s, v, _ = data
            try:
                if peak == "OLV":
                    # Threshold MM
                    model_func = KineticsModels.threshold_michaelis_menten
                    p0 = [np.max(v) * 2, 1.0, 0.5]
                    popt, pcov = curve_fit(
                        model_func, s, v, p0=p0, bounds=(0, np.inf), maxfev=10000
                    )
                    perr = np.sqrt(np.diag(pcov))
                    ssr = np.sum((v - model_func(s, *popt)) ** 2)
                    _logger.info(f"Fit 1 (Best): S0 = {popt[2]:.3f} | SSR = {ssr:.4e}")

                    # Hill Equation (stored as extra for comparison)
                    hill_func = KineticsModels.hill_equation
                    p0_h = [np.max(v), 1.0, 2.0]
                    popt_h, pcov_h = curve_fit(
                        hill_func, s, v, p0=p0_h, bounds=(0, np.inf), maxfev=10000
                    )
                    perr_h = np.sqrt(np.diag(pcov_h))

                    extra = {
                        "S0": popt[2],
                        "S0_SE": perr[2],
                        "Hill_Model": "Hill",
                        "Hill_Vmax": popt_h[0],
                        "Hill_Vmax_SE": perr_h[0],
                        "Hill_n": popt_h[2],
                        "Hill_n_SE": perr_h[2],
                        "popt_h": tuple(popt_h),  # Retain for plotting
                    }
                else:
                    # Standard MM
                    model_func = KineticsModels.michaelis_menten
                    vmax_guess = np.max(v)
                    km_guess = s[np.argmin(np.abs(v - vmax_guess / 2))]
                    popt, pcov = curve_fit(
                        model_func,
                        s,
                        v,
                        p0=[vmax_guess, km_guess],
                        bounds=(0, np.inf),
                        maxfev=10000,
                    )
                    perr = np.sqrt(np.diag(pcov))
                    extra = None

                kcat = popt[0] / self.config.prot_conc
                kcat_se = perr[0] / self.config.prot_conc
                kcat_km = kcat / popt[1] if popt[1] > 0 else np.nan

                self.fits[peak] = KineticParams(
                    peak_id=peak,
                    vmax=popt[0],
                    vmax_se=perr[0],
                    km=popt[1],
                    km_se=perr[1],
                    kcat=kcat,
                    kcat_se=kcat_se,
                    kcat_km=kcat_km,
                    popt=tuple(popt),
                    model_func=model_func,
                    extra=extra,
                )
            except Exception as e:
                _logger.warning(f"Fit failed for peak {peak}: {e}")

        # Construct analytical results dataframe
        self.stats_df = pd.DataFrame([fit.to_dict() for fit in self.fits.values()])
        return self

    def plot_results(self, show: bool = False) -> Self:
        """Generates and saves the matplotlib facet plot."""
        peaks = self.df[PEAK_ID].unique()
        cols = 3
        rows = int(np.ceil(len(peaks) / cols))

        # Memory explicit management
        fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
        axes = axes.flatten()

        for i, peak in enumerate(peaks):
            ax = axes[i]
            data = self._prepare_data(peak)

            if not data:
                ax.set_title(f"{peak} - No Data")
                ax.axis("off")
                continue

            s, v, v_sem = data
            ax.errorbar(
                s,
                v,
                yerr=v_sem,
                fmt="o",
                label="Data (Mean ± SEM)",
                color="blue",
                alpha=0.7,
                capsize=3,
            )

            fit = self.fits.get(peak)
            if fit:
                s_fit = np.linspace(0, np.max(s) * 1.1, 200)

                if peak == "OLV" and fit.extra:
                    label_text = (
                        f"Threshold Fit:\n$V_{{max}}$={fit.vmax:.2f}±{fit.vmax_se:.2f}\n"
                        f"$K_m$={fit.km:.2f}±{fit.km_se:.2f}\n$S_0$={fit.extra['S0']:.2f}\n$k_{{cat}}$={fit.kcat:.4f}"
                    )
                    ax.plot(
                        s_fit,
                        KineticsModels.hill_equation(s_fit, *fit.extra["popt_h"]),
                        label=f"Hill Fit (n={fit.extra['Hill_n']:.2f})",
                        linestyle="-",
                        color="blue",
                    )
                else:
                    label_text = (
                        f"MM Fit:\n$V_{{max}}$={fit.vmax:.2f}±{fit.vmax_se:.2f}\n"
                        f"$K_m$={fit.km:.2f}±{fit.km_se:.2f}\n$k_{{cat}}$={fit.kcat:.4f}"
                    )

                ax.plot(
                    s_fit,
                    fit.model_func(s_fit, *fit.popt),
                    label=label_text,
                    linestyle="--",
                    color="red",
                )
                ax.legend(fontsize=9, loc="lower right")
            else:
                ax.text(
                    0.5,
                    0.5,
                    "Fit Failed",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )

            ax.set_title(f"Kinetics: {peak}", fontweight="bold")
            ax.set_xlabel("Substrate Concentration [HexCoA] (μM)")
            ax.set_ylabel("Velocity (Area / min)")

        for j in range(len(peaks), len(axes)):
            axes[j].axis("off")

        plt.grid(True, linestyle=":", alpha=0.6)
        plt.tight_layout()

        plot_path = self.dest / f"{Path(self.path).stem}_plot.png"
        plt.savefig(plot_path, dpi=300)

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
        if self.stats_df.empty:
            _logger.warning("Export: No kinetics data to export.")
            return self

        try:
            # Using your system's builder
            build_path = ExportPathBuilder(
                f"{Path(self.path).stem}_kinetics.csv",
                overwrite=self.config.overwrite,
                merge=self.config.merge,
            ).build(self.dest)
            out_path = build_path.path
        except NameError:
            # Fallback if ExportPathBuilder isn't imported
            out_path = self.dest / f"{Path(self.path).stem}_kinetics.csv"

        self.stats_df.to_csv(out_path, index=False)
        _logger.info(
            f"Saved enzyme kinetics data (shape={self.stats_df.shape}) to: {out_path}"
        )
        return self
