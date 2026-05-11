from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final, Protocol

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Standardized Constants
SUBSTRATE_CONC: Final[str] = "substrate_conc"
VELOCITY: Final[str] = "velocity"
SEM: Final[str] = "v_sem"

_logger = logging.getLogger(__name__)


# --- Protocols & Models ---

class KineticModel(Protocol):

    def __call__(self, s: np.ndarray, *args: Any) -> np.ndarray: ...


# class KineticsModels:
#     """Vectorized mathematical models."""
#
#     @staticmethod
#     def michaelis_menten(s: np.ndarray, vmax: float, km: float) -> np.ndarray:
#         return (vmax * s) / (km + s)
#
#     @staticmethod
#     def hill_equation(s: np.ndarray, vmax: float, k_half: float, n: float) -> np.ndarray:
#         return (vmax * s ** n) / (k_half ** n + s ** n)
#
#     @staticmethod
#     def threshold_mm(s: np.ndarray, vmax: float, km: float, s0: float) -> np.ndarray:
#         return np.where(s > s0, (vmax * (s - s0)) / (km + (s - s0)), 0.0)


# --- Data Containers ---

@dataclass(frozen=True, slots=True)
class FitResult:
    peak_id: str
    vmax: float
    vmax_se: float
    km: float
    km_se: float
    kcat: float
    kcat_km: float
    popt: tuple[float, ...]
    model_name: str
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "peak_id": self.peak_id,
            "vmax": self.vmax, "vmax_se": self.vmax_se,
            "km": self.km, "km_se": self.km_se,
            "kcat": self.kcat, "kcat_km": self.kcat_km,
            **self.extra
        }


# --- Core Logic Classes ---

class KineticsFitter:
    """Handles the heavy lifting of curve fitting."""

    def __init__(self, prot_conc: float, rxn_time: float):
        self.prot_conc = prot_conc
        self.rxn_time = rxn_time

    def fit_peak(
        self,
        peak_id: str,
        s: np.ndarray,
        mean_area: np.ndarray,
        std_area: np.ndarray,
        counts: np.ndarray,
    ) -> FitResult | None:
        # Vectorized velocity calculation
        v = mean_area / self.rxn_time
        v_sem = (std_area / self.rxn_time) / np.sqrt(counts)

        try:
            # Model selection logic
            if peak_id == 'OLV':
                model = KineticsModels.threshold_mm
                p0, bounds = [np.max(v), 1.0, 0.5], (0, np.inf)
                name = "Threshold MM"
            else:
                model = KineticsModels.michaelis_menten
                p0, bounds = [np.max(v), np.median(s)], (0, np.inf)
                name = "Michaelis-Menten"

            popt, pcov = curve_fit(model, s, v, p0=p0, bounds=bounds, maxfev=10000)
            perr = np.sqrt(np.diag(pcov))

            vmax, km = popt[0], popt[1]
            kcat = vmax / self.prot_conc

            return FitResult(
                peak_id=peak_id,
                vmax=vmax, vmax_se=perr[0],
                km=km, km_se=perr[1],
                kcat=kcat,
                kcat_km=kcat / km if km > 0 else 0,
                popt=tuple(popt),
                model_name=name,
                extra={"s0": popt[2]} if len(popt) > 2 else {},
            )
        except Exception as e:
            _logger.error(f"Fit failed for {peak_id}: {e}")
            return None


class KineticsAnalyzer:
    """High-level API for orchestrating the analysis."""

    def __init__(self, config: dict[str, Any]):
        self.fitter = KineticsFitter(config.get('prot_conc', 12.25), config.get('rxn_time', 180.0))
        self.results: list[FitResult] = []
        self.data: pd.DataFrame = pd.DataFrame()

    def load_data(self, path: Path | str) -> pd.DataFrame:
        """Loads and optimizes memory footprint."""
        df = pd.read_csv(path)
        # Memory optimization: Downcast and Categoricals
        if 'peak_id' in df.columns:
            df['peak_id'] = df['peak_id'].astype('category')

        for col in [SUBSTRATE_CONC, 'mean', 'std', 'count']:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], downcast='float')

        self.data = df
        return df

    def analyze(self) -> pd.DataFrame:
        """Runs fitting for all unique peaks."""
        self.results = []
        for peak_id, group in self.data.groupby('peak_id', observed=True):
            res = self.fitter.fit_peak(
                str(peak_id),
                group[SUBSTRATE_CONC].values,
                group['mean'].values,
                group['std'].values,
                group['count'].values,
            )
            if res:
                self.results.append(res)

        return pd.DataFrame([r.to_dict() for r in self.results])

    def plot(self, output_path: Path):
        """Generates faceted plots."""
        n_peaks = len(self.results)
        cols = 3
        rows = (n_peaks + cols - 1) // cols

        fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows), squeeze=False)
        axes_flat = axes.flatten()

        for ax, res in zip(axes_flat, self.results):
            # Extract data for this peak
            group = self.data[self.data['peak_id'] == res.peak_id]
            s, v = group[SUBSTRATE_CONC].values, group['mean'].values / self.fitter.rxn_time

            # Plot raw data
            ax.errorbar(s, v, fmt='o', capsize=3, label='Data', alpha=0.6)

            # Plot fit
            s_fit = np.linspace(0, s.max() * 1.1, 100)
            model_func = (KineticsModels.threshold_mm if "Threshold" in res.model_name
                          else KineticsModels.michaelis_menten)
            ax.plot(s_fit, model_func(s_fit, *res.popt), 'r--', label=f'Fit ({res.model_name})')

            ax.set_title(f"Peak: {res.peak_id}")
            ax.legend(fontsize=8)

        # Cleanup empty plots
        for i in range(n_peaks, len(axes_flat)):
            axes_flat[i].axis('off')

        plt.tight_layout()
        plt.savefig(output_path, dpi=200)
        plt.close(fig)


# --- Updated API Usage ---

if __name__ == "__main__":
    config = {"prot_conc": 12.25, "rxn_time": 180.0}
    analyzer = KineticsAnalyzer(config)

    # 1. Load optimized data
    analyzer.load_data("kinetics_data.csv")

    # 2. Perform math fitting
    stats_df = analyzer.analyze()
    stats_df.to_csv("results.csv", index=False)

    # 3. Visualize
    analyzer.plot(Path("kinetics_plot.png"))
