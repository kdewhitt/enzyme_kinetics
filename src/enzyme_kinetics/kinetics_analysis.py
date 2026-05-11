from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Final, Self

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .calibration import Calibration
from .config import KineticsConfig
from .kinetics_core import (
    derive_constants,
    fit_hill,
    fit_lineweaver_burk,
    fit_michaelis_menten,
    FitResult,
    KineticConstants,
    lineweaver_burk_transform,
    prepare_velocity,
)

_logger = logging.getLogger(__name__)

SUBSTRATE_CONC: Final[str] = "substrate_conc"
COUNT = "count"
MEAN = "mean"
PEAK_ID = "peak_id"
# PROTEIN = "protein"
# REL_ACT = "rel_activity"
STD = "std"


# ---------------------------------------------------------------------------
# Per-peak data extraction
# FIX 3: unpack four values from prepare_velocity (now returns mean_signal too)
# ---------------------------------------------------------------------------

def _extract_peak_data(
    df: pd.DataFrame,
    peak_id: str,
    rxn_time: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Extract substrate concentrations, velocity, SEM, and raw area for a peak.

    FIX 3: Returns raw mean_signal as the fourth element, enabling
    apply_calibration to work directly on area rather than reconstructing
    it from velocity — avoiding the fragile v * rxn_time round-trip.

    Args:
        df: DataFrame containing all peaks.
        peak_id: The peak identifier to filter on.
        rxn_time: Reaction time in seconds used to convert area to velocity.

    Returns:
        Tuple of (s, v, v_sem, mean_signal), or None if no valid data.
    """
    sub = df.loc[df[PEAK_ID] == peak_id].dropna(subset=[SUBSTRATE_CONC, MEAN])
    if sub.empty or sub[MEAN].max() == 0:
        return None
    return prepare_velocity(
        sub[SUBSTRATE_CONC].values,
        sub[MEAN].values,  # velocity?
        sub[STD].values,
        sub[COUNT].values,
        rxn_time,
    )


# ---------------------------------------------------------------------------
# Model selection strategy
# ---------------------------------------------------------------------------

def _select_and_fit(
    peak_id: str,
    s: np.ndarray,
    v: np.ndarray,
    v_sem: np.ndarray,
    *,
    special_peaks: dict[str, str] | None = None,
) -> FitResult:
    special = (special_peaks or {}).get(peak_id)
    if special == "hill":
        return fit_hill(s, v, sigma=v_sem)
    return fit_michaelis_menten(s, v, sigma=v_sem)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

class EnzymeKineticsAnalysis:
    def __init__(
        self,
        path: str | Path,
        target_dir: str | Path,
        *,
        config: KineticsConfig | None = None,
        special_peaks: dict[str, str] | None = None,
    ) -> None:
        self.path = Path(path).expanduser().resolve()
        self.dest = Path(target_dir).expanduser().resolve()
        self.config = config or KineticsConfig()
        self.special_peaks = special_peaks or {}

        self.df: pd.DataFrame = pd.DataFrame()
        self.results: dict[str, KineticConstants] = {}
        self.calibration: Calibration | None = None

    # -- pipeline stages ---------------------------------------------------

    def load(self, **kwargs: Any) -> Self:
        self.df = self.config.load(self.path, **kwargs)
        self.df.sort_values(by=[PEAK_ID, SUBSTRATE_CONC], inplace=True)
        _logger.info("Loaded %d peaks from %s", len(self.df), self.path)
        return self

    def fit(self) -> Self:
        self.results.clear()
        for peak_id in self.df[PEAK_ID].unique():
            data = _extract_peak_data(self.df, peak_id, self.config.rxn_time)
            if data is None:
                continue
            s, v, v_sem, _ = data  # raw area not needed at fit stage
            try:
                mm_fit = _select_and_fit(
                    peak_id, s, v, v_sem, special_peaks=self.special_peaks,
                )
                lb_fit: FitResult | None = None
                try:
                    lb_fit = fit_lineweaver_burk(s, v)
                except Exception:
                    pass
                self.results[peak_id] = derive_constants(
                    peak_id, mm_fit, self.config.prot_conc, lb_fit=lb_fit,
                )
                print(f"✓ {peak_id}: {self.results[peak_id]}")

            except Exception as exc:
                _logger.warning("Fit failed for %s: %s", peak_id, exc)
        return self

    def apply_calibration(
        self,
        cal: Calibration,
        *,
        precise_sem: bool = True,  # FIX 4: user-selectable SEM scaling method
    ) -> Self:
        """Re-fit all peaks using calibrated concentration velocities.

        FIX 3: Calibration is now applied directly to the raw mean_signal
        (peak area) before dividing by rxn_time, avoiding the fragile
        reconstruction of area from velocity that existed in the original.

        FIX 4: Two v_sem scaling modes are available:
            precise_sem=True  (default): uses area_to_conc_with_error per
                data point for full delta-method propagation including
                slope_se and intercept_se contributions.
            precise_sem=False: simplified approximation v_sem / slope,
                appropriate when calibration uncertainty is negligible
                compared to replicate variance.

        This method rebuilds self.results from self.df on every call, making
        it safe to call more than once (e.g., after updating the calibration).

        FIX 6: Warns if the calibration does not pass is_valid().

        Args:
            cal: A Calibration object produced by fit_calibration().
            precise_sem: If True, use full delta-method SEM propagation
                (recommended). If False, use the simplified slope-only
                approximation.

        Returns:
            self, for method chaining.
        """
        # FIX 6: guard with calibration validity check
        if not cal.is_valid():
            _logger.warning(
                "Calibration does not meet quality thresholds "
                "(R²=%.4f, slope=%.4g). Results may be unreliable.",
                cal.r_squared, cal.slope,
            )

        self.calibration = cal
        recalculated: dict[str, KineticConstants] = {}

        for peak_id, kc in self.results.items():
            data = _extract_peak_data(self.df, peak_id, self.config.rxn_time)
            if data is None:
                continue

            # FIX 3: unpack raw area; apply calibration directly to area
            s, _v_raw, v_sem_raw, mean_signal = data

            # Convert area → µM concentration, then divide by rxn_time for velocity
            v_um = np.asarray(cal.area_to_conc(mean_signal)) / self.config.rxn_time

            # FIX 4: precise vs. simplified SEM scaling
            if precise_sem:
                # Full delta-method per data point: propagates slope_se and intercept_se
                v_sem_um = np.array(
                    [
                        cal.area_to_conc_with_error(float(area), float(area_se * self.config.rxn_time))[1]
                        / self.config.rxn_time
                        for area, area_se in zip(mean_signal, v_sem_raw * self.config.rxn_time)
                    ],
                )
            else:
                # Simplified approximation: ignores calibration parameter uncertainty
                v_sem_um = v_sem_raw / cal.slope

            try:
                mm_fit = fit_michaelis_menten(s, v_um, sigma=v_sem_um)
                recalculated[peak_id] = derive_constants(
                    peak_id, mm_fit, self.config.prot_conc, lb_fit=kc.lb_fit,
                )
            except Exception as exc:
                _logger.warning("Calibrated re-fit failed for %s: %s", peak_id, exc)

        self.results = recalculated
        return self

    # -- output ------------------------------------------------------------

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame([kc.to_dict() for kc in self.results.values()])

    def plot(self, *, show: bool = False, cols: int = 3) -> Self:
        """Plot MM (or Hill) fit curves for all peaks in a multi-panel grid.

        Each panel shows the raw data with error bars and the fitted curve,
        annotated with Vmax, Km (or k_half for Hill), and R².

        Args:
            show: If True, call plt.show() after saving.
            cols: Number of columns in the panel grid.

        Returns:
            self, for method chaining.
        """
        peaks = list(self.results)
        if not peaks:
            _logger.warning("No fits to plot.")
            return self

        rows = int(np.ceil(len(peaks) / cols))
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows))
        axes = np.atleast_1d(axes).flatten()

        for i, peak_id in enumerate(peaks):
            ax = axes[i]
            data = _extract_peak_data(self.df, peak_id, self.config.rxn_time)
            kc = self.results[peak_id]

            if data is None:
                ax.set_title(f"{peak_id} — no data")
                ax.axis("off")
                continue

            s, v, v_sem, _ = data
            ax.errorbar(
                s, v, yerr=v_sem,
                fmt="o", color="steelblue", alpha=0.7, capsize=3, label="Data",
            )

            s_fit = np.linspace(0, np.max(s) * 1.1, 200)
            v_fit = kc.fit.predict(s_fit)

            # FIX 7: use model-correct label for the affinity parameter
            is_hill = kc.fit.model_type == "hill"
            km_label = "k_half" if is_hill else "Km"
            label = (
                f"Vmax={kc.fit.vmax:.2f}±{kc.fit.vmax_se:.2f}\n"
                f"{km_label}={kc.fit.km:.2f}±{kc.fit.km_se:.2f}\n"
                f"R²={kc.fit.r_squared:.4f}"
            )
            ax.plot(s_fit, v_fit, "--", color="red", label=label)
            ax.legend(fontsize=8, loc="lower right")
            ax.set_title(peak_id, fontweight="bold")
            ax.set_xlabel("[S] (µM)")
            ax.set_ylabel("v (signal / min)")

        for j in range(len(peaks), len(axes)):
            axes[j].axis("off")

        fig.tight_layout()
        plot_path = self.dest / f"{self.path.stem}_kinetics.png"
        fig.savefig(plot_path, dpi=300)
        if show:
            plt.show()
        plt.close(fig)
        _logger.info("Saved kinetics plot to %s", plot_path)
        return self

    # FIX 5 — restore publication plots: LB, residuals, efficiency comparison
    def plot_lineweaver_burk(self, *, show: bool = False, cols: int = 3) -> Self:
        """Plot Lineweaver-Burk (double-reciprocal) panels for all peaks with LB fits.

        The fit line is drawn over the observed 1/[S] range. The y-intercept
        (= 1/Vmax) and x-intercept (= −1/Km) are annotated to aid visual
        inspection. Peaks without a valid LB fit are silently skipped.

        Args:
            show: If True, call plt.show() after saving.
            cols: Number of columns in the panel grid.

        Returns:
            self, for method chaining.
        """
        lb_peaks = [pid for pid, kc in self.results.items() if kc.lb_fit is not None]
        if not lb_peaks:
            _logger.warning("No Lineweaver-Burk fits available to plot.")
            return self

        rows = int(np.ceil(len(lb_peaks) / cols))
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
        axes = np.atleast_1d(axes).flatten()

        for i, peak_id in enumerate(lb_peaks):
            ax = axes[i]
            data = _extract_peak_data(self.df, peak_id, self.config.rxn_time)
            kc = self.results[peak_id]
            lb = kc.lb_fit  # guaranteed not None

            if data is None:
                ax.set_visible(False)
                continue

            s, v, _, _ = data
            s_inv, v_inv = lineweaver_burk_transform(s, v)

            ax.scatter(s_inv, v_inv, color="steelblue", zorder=3, label="Data (1/v vs 1/[S])")

            # Fit line over observed 1/[S] range
            x_line = np.linspace(s_inv.min(), s_inv.max(), 200)
            slope = lb.extra["slope"]
            intercept = lb.extra["intercept"]
            ax.plot(
                x_line, slope * x_line + intercept, "--", color="red",
                label=f"R²={lb.r_squared:.4f}",
            )

            # Annotate intercepts
            if not np.isnan(lb.vmax) and intercept != 0:
                ax.axhline(intercept, color="grey", lw=0.8, ls=":")
                ax.annotate(
                    f"1/Vmax={intercept:.3g}", xy=(s_inv.min(), intercept),
                    fontsize=7, color="grey", va="bottom",
                )
            if not np.isnan(lb.km) and slope != 0:
                x_int = -intercept / slope
                ax.axvline(x_int, color="orange", lw=0.8, ls=":")
                ax.annotate(
                    f"−1/Km={x_int:.3g}", xy=(x_int, v_inv.min()),
                    fontsize=7, color="orange", ha="right",
                )

            ax.set_title(peak_id, fontweight="bold")
            ax.set_xlabel("1/[S] (µM⁻¹)")
            ax.set_ylabel("1/v")
            ax.legend(fontsize=8)

        for j in range(len(lb_peaks), len(axes)):
            axes[j].axis("off")

        fig.tight_layout()
        plot_path = self.dest / f"{self.path.stem}_lineweaver_burk.png"
        fig.savefig(plot_path, dpi=300)
        if show:
            plt.show()
        plt.close(fig)
        _logger.info("Saved Lineweaver-Burk plot to %s", plot_path)
        return self

    def plot_residuals(self, *, show: bool = False, cols: int = 3) -> Self:
        """Plot MM fit residuals (v_observed − v_predicted) vs [S] for all peaks.

        A horizontal reference line at zero is drawn. Systematic curvature in
        residuals indicates model mis-specification (e.g., substrate inhibition
        or cooperativity not accounted for).

        Args:
            show: If True, call plt.show() after saving.
            cols: Number of columns in the panel grid.

        Returns:
            self, for method chaining.
        """
        peaks = list(self.results)
        if not peaks:
            _logger.warning("No fits to plot residuals for.")
            return self

        rows = int(np.ceil(len(peaks) / cols))
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
        axes = np.atleast_1d(axes).flatten()

        for i, peak_id in enumerate(peaks):
            ax = axes[i]
            data = _extract_peak_data(self.df, peak_id, self.config.rxn_time)
            kc = self.results[peak_id]

            if data is None:
                ax.set_visible(False)
                continue

            s, v, v_sem, _ = data
            v_pred = kc.fit.predict(s)
            residuals = v - v_pred

            ax.axhline(0, color="black", lw=0.8, ls="--")
            ax.errorbar(
                s, residuals, yerr=v_sem, fmt="o", color="steelblue",
                alpha=0.7, capsize=3,
            )
            ax.set_title(peak_id, fontweight="bold")
            ax.set_xlabel("[S] (µM)")
            ax.set_ylabel("Residual (v − v̂)")

        for j in range(len(peaks), len(axes)):
            axes[j].axis("off")

        fig.tight_layout()
        plot_path = self.dest / f"{self.path.stem}_residuals.png"
        fig.savefig(plot_path, dpi=300)
        if show:
            plt.show()
        plt.close(fig)
        _logger.info("Saved residuals plot to %s", plot_path)
        return self

    def plot_efficiency_comparison(self, *, show: bool = False) -> Self:
        """Plot a 4-panel bar chart comparing Km, Vmax, kcat, and kcat/Km across peaks.

        Panels with no valid data (e.g. kcat requires calibration) are labelled
        accordingly rather than silently omitting bars, preventing index mismatches
        between the peak list and the plotted values.

        Args:
            show: If True, call plt.show() after saving.

        Returns:
            self, for method chaining.
        """
        if not self.results:
            _logger.warning("No results available for efficiency comparison.")
            return self

        peaks = list(self.results)
        kcs = [self.results[p] for p in peaks]

        km_vals = [kc.fit.km for kc in kcs]
        km_errs = [kc.fit.km_se for kc in kcs]
        vmax_vals = [kc.fit.vmax for kc in kcs]
        vmax_errs = [kc.fit.vmax_se for kc in kcs]

        # kcat and kcat/Km may be nan if calibration not applied; handled explicitly
        kcat_vals = [kc.kcat for kc in kcs]
        kcat_errs = [kc.kcat_se for kc in kcs]
        kcat_km_vals = [kc.kcat_km for kc in kcs]
        kcat_km_errs = [kc.kcat_km_se for kc in kcs]

        x = np.arange(len(peaks))
        bar_kw: dict[str, Any] = dict(capsize=4, alpha=0.8, width=0.6)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        ax_km, ax_vmax, ax_kcat, ax_kcat_km = axes.flatten()

        def _bar(ax: plt.Axes, vals: list[float], errs: list[float], title: str, ylabel: str) -> None:
            colors = ["steelblue" if not np.isnan(v) else "lightgrey" for v in vals]
            safe_errs = [e if not np.isnan(e) else 0.0 for e in errs]
            safe_vals = [v if not np.isnan(v) else 0.0 for v in vals]
            ax.bar(x, safe_vals, yerr=safe_errs, color=colors, **bar_kw)
            ax.set_xticks(x)
            ax.set_xticklabels(peaks, rotation=45, ha="right", fontsize=8)
            ax.set_title(title, fontweight="bold")
            ax.set_ylabel(ylabel)
            if all(np.isnan(v) for v in vals):
                ax.text(
                    0.5, 0.5, "No data\n(calibration required?)",
                    transform=ax.transAxes, ha="center", va="center",
                    color="grey", fontsize=10,
                )

        _bar(ax_km, km_vals, km_errs, "Km", "Km (µM)")
        _bar(ax_vmax, vmax_vals, vmax_errs, "Vmax", "Vmax (signal/min)")
        _bar(ax_kcat, kcat_vals, kcat_errs, "kcat", "kcat (s⁻¹)")
        _bar(ax_kcat_km, kcat_km_vals, kcat_km_errs, "kcat / Km", "kcat/Km (µM⁻¹·s⁻¹)")

        fig.tight_layout()
        plot_path = self.dest / f"{self.path.stem}_efficiency.png"
        fig.savefig(plot_path, dpi=300)
        if show:
            plt.show()
        plt.close(fig)
        _logger.info("Saved efficiency comparison plot to %s", plot_path)
        return self

    def export_csv(self) -> Self:
        stats_df = self.to_dataframe()
        if stats_df.empty:
            _logger.warning("No kinetics data to export.")
            return self
        csv_path = self.config.build_export_path(
            self.dest, tag=self.path.stem,
        ).with_suffix(".csv")
        stats_df.to_csv(csv_path, index=False)
        _logger.info("Exported %s rows to %s", len(stats_df), csv_path)
        return self

    def run(self, **kwargs: Any) -> Self:
        """Full pipeline: load → fit → plot → Lineweaver-Burk → residuals → efficiency → CSV."""
        return (
            self.load(**kwargs)
            .fit()
            .plot()
            .plot_lineweaver_burk()
            .plot_residuals()
            .plot_efficiency_comparison()
            .export_csv()
        )
