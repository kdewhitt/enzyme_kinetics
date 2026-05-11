from __future__ import annotations

from pathlib import Path
from typing import Any, Self

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from kgdlibs.pathtools import PathBuilder
from logurich import RichLogAdapter

from .derive import KineticConstants
from .models import lineweaver_burk_transform
from .utils import extract_peak_data

_logger = RichLogAdapter(component=__name__)


class KineticPlots:
    def __init__(
        self,
        path: str | Path,
        data: pd.DataFrame,
        results: dict[str, KineticConstants],
        reaction_time_seconds: float,
        overwrite: bool = False,
    ) -> None:
        self.path = Path(path)
        self.df = data
        self.results = results
        self.reaction_time_seconds = reaction_time_seconds
        self.overwrite = overwrite

    def _make_builder(self) -> PathBuilder:
        """Returns a PathBuilder seeded with the acquisition date and destination path."""
        return PathBuilder.from_path(self.path, suffix=".png", overwrite=self.overwrite)

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
            data = extract_peak_data(self.df, peak_id, self.reaction_time_seconds)
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
        plot_path = self._make_builder().with_tag("kinetics")
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
            data = extract_peak_data(self.df, peak_id, self.reaction_time_seconds)
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
        plot_path = self._make_builder().with_tag("lineweaver_burk")
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
            data = extract_peak_data(self.df, peak_id, self.reaction_time_seconds)
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
        plot_path = self._make_builder().with_tag("residuals")
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
        kcat_km_vals = [kc.kcat_km_M for kc in kcs]
        kcat_km_errs = [kc.kcat_km_se_M for kc in kcs]

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
        _bar(ax_kcat_km, kcat_km_vals, kcat_km_errs, "kcat / Km", "kcat/Km (M⁻¹·s⁻¹)")

        fig.tight_layout()
        plot_path = self._make_builder().with_tag("efficiency")
        fig.savefig(plot_path, dpi=300)
        if show:
            plt.show()
        plt.close(fig)
        _logger.info("Saved efficiency comparison plot to %s", plot_path)
        return self
