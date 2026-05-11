from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from .models import r_squared

_logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Calibration — linear standard curve
# FIX 2, FIX 6
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class Calibration:
    """Linear HPLC calibration curve (peak area = slope × [CoA] + intercept).

    Assumes a linear detector response (Beer-Lambert regime). Fit via
    curve_fit with optional sigma weighting.

    Attributes:
        slope: Calibration slope (area / µM).
        slope_se: Standard error of slope.
        intercept: Calibration intercept (area at zero concentration).
        intercept_se: Standard error of intercept.
        r_squared: Coefficient of determination; ≥ 0.999 expected for HPLC.
        conc_range: (min, max) concentration of calibration standards in µM.
    """

    slope: float
    slope_se: float
    intercept: float
    intercept_se: float
    r_squared: float
    rmse: float
    n_points: int
    conc_range: tuple[float, float]

    def __str__(self) -> str:
        return (f"Slope={self.slope:.6f}±{self.slope_se:.6f} area/µM | "
                f"Intercept={self.intercept:.6f}±{self.intercept_se:.6f} | "
                f"R²={self.r_squared:.6f} | "
                f"RMSE={self.rmse:.4f}")

    # FIX 6: restore is_valid() guard from original CalibrationCurve
    def is_valid(self, r2_threshold: float = 0.999) -> bool:
        """Return True if the calibration meets quality thresholds.

        Requires slope > 0 (positive detector response) and R² ≥ r2_threshold.
        The default threshold is 0.999, consistent with typical CoA/HPLC
        calibrations; the original codebase used 0.95, which is retained as
        an acceptable lower bound via the keyword argument.

        Args:
            r2_threshold: Minimum acceptable R². Defaults to 0.999.

        Returns:
            True if slope > 0 and r_squared ≥ r2_threshold.
        """
        return self.slope > 0 and self.r_squared >= r2_threshold

    def area_to_conc(self, area: np.ndarray | float) -> np.ndarray | float:
        return (area - self.intercept) / self.slope

    def conc_to_area(self, conc: np.ndarray | float) -> np.ndarray | float:
        return self.slope * conc + self.intercept

    def area_to_conc_with_error(
        self,
        area: float,
        area_se: float,
    ) -> tuple[float, float]:
        """Convert a single area measurement to concentration with propagated error.

        Uses the full delta method, including intercept uncertainty (FIX 2):

            σ_C² = (σ_area / slope)²
                 + ((area − intercept) · σ_slope / slope²)²
                 + (σ_intercept / slope)²

        The intercept term was omitted in the original codebase. For well-behaved
        calibrations (intercept ≈ 0) the contribution is negligible, but it is
        included here for completeness.

        Args:
            area: Measured peak area.
            area_se: Standard error of the area measurement.

        Returns:
            Tuple of (concentration in µM, concentration SE in µM).
        """
        conc = float(self.area_to_conc(area))
        dc_da = 1.0 / self.slope
        dc_ds = -(area - self.intercept) / (self.slope ** 2)
        dc_di = -1.0 / self.slope  # FIX 2: intercept term
        conc_se = float(
            np.sqrt(
                (dc_da * area_se) ** 2
                + (dc_ds * self.slope_se) ** 2
                + (dc_di * self.intercept_se) ** 2,  # FIX 2
            ),
        )
        return conc, conc_se


# ---------------------------------------------------------------------------
# Fitting — thin wrappers around curve_fit
# ---------------------------------------------------------------------------


def fit_calibration(
    concentrations: np.ndarray,
    peak_areas: np.ndarray,
    *,
    sigma: np.ndarray | None = None,
) -> Calibration:
    def _linear(x: np.ndarray, slope: float, intercept: float) -> np.ndarray:
        return slope * x + intercept

    effective_sigma = sigma if sigma is not None and np.any(sigma > 0) else None
    popt, pcov = curve_fit(
        _linear, concentrations, peak_areas,
        sigma=effective_sigma, absolute_sigma=True,
    )
    perr = np.sqrt(np.diag(pcov))
    predicted = _linear(concentrations, *popt)
    return Calibration(
        slope=popt[0],
        slope_se=perr[0],
        intercept=popt[1],
        intercept_se=perr[1],
        r_squared=r_squared(peak_areas, predicted),
        rmse=float(np.sqrt(np.mean((peak_areas - predicted) ** 2))),
        n_points=len(concentrations),
        conc_range=(float(concentrations.min()), float(concentrations.max())),
    )


# ---------------------------------------------------------------------------
# Legacy
# ---------------------------------------------------------------------------


def save_calibration(cal: Calibration, path: Path, nice: bool = False) -> None:
    if not nice:
        data = {
            "slope": cal.slope,
            "slope_se": cal.slope_se,
            "intercept": cal.intercept,
            "intercept_se": cal.intercept_se,
            "r_squared": cal.r_squared,
            "conc_range": list(cal.conc_range),
        }
        path.write_text(json.dumps(data, indent=2))

    if nice:

        with open(path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("HPLC CALIBRATION PARAMETERS\n")
            f.write("=" * 80 + "\n\n")

            f.write("LINEAR MODEL:\n")
            f.write("  peak_area = slope × [concentration] + intercept\n\n")

            f.write("PARAMETERS:\n")
            f.write(f"  Slope:        {cal.slope:.10f} ± {cal.slope_se:.10f} area/µM\n")
            f.write(f"  Intercept:    {cal.intercept:.10f} ± {cal.intercept_se:.10f} area\n")
            f.write(f"  R²:           {cal.r_squared:.10f}\n")
            f.write(f"  RMSE:         {cal.rmse:.6f} area\n")
            f.write(f"  N points:     {cal.n_points}\n")
            f.write(f"  Conc range:   {cal.conc_range[0]:.2f}–{cal.conc_range[1]:.2f} µM\n\n")

            f.write("INVERSE FUNCTION:\n")
            f.write(f"  [concentration] = (peak_area - {cal.intercept:.10f}) / {cal.slope:.10f}\n\n")

            f.write("PYTHON CODE:\n")
            f.write("```python\n")
            f.write(f"CALIB_SLOPE = {cal.slope:.10f}\n")
            f.write(f"CALIB_INTERCEPT = {cal.intercept:.10f}\n")
            f.write(f"CALIB_SLOPE_SE = {cal.slope_se:.10f}\n")
            f.write(f"CALIB_INTERCEPT_SE = {cal.intercept_se:.10f}\n\n")
            f.write("concentration = (peak_area - CALIB_INTERCEPT) / CALIB_SLOPE\n")
            f.write("```\n")


def load_calibration(path: Path) -> Calibration:
    data = json.loads(path.read_text())
    return Calibration(
        slope=data["slope"],
        slope_se=data["slope_se"],
        intercept=data["intercept"],
        intercept_se=data["intercept_se"],
        r_squared=data["r_squared"],
        conc_range=tuple(data["conc_range"]),
    )


def plot_calibration(
    cal: Calibration,
    concentrations: np.ndarray,
    peak_areas: np.ndarray,
    *,
    output_path: Path | None = None,
    show: bool = False,
) -> None:
    predicted = cal.conc_to_area(concentrations)
    residuals = peak_areas - predicted
    std_res = residuals / np.std(residuals) if np.std(residuals) > 0 else residuals

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    ax_fit, ax_res, ax_std = axes

    # Panel 1 — calibration curve + fit
    ax_fit.scatter(concentrations, peak_areas, color="steelblue", zorder=3, label="Standards")
    c_line = np.linspace(concentrations.min(), concentrations.max(), 200)
    ax_fit.plot(
        c_line, cal.conc_to_area(c_line), "--", color="red",
        label=f"slope={cal.slope:.4g}, R²={cal.r_squared:.4f}",
    )
    ax_fit.set_xlabel("[CoA] (µM)")
    ax_fit.set_ylabel("Peak area")
    ax_fit.set_title("Calibration curve")
    ax_fit.legend(fontsize=8)

    # Panel 2 — absolute residuals
    ax_res.axhline(0, color="black", lw=0.8, ls="--")
    ax_res.scatter(concentrations, residuals, color="steelblue")
    ax_res.set_xlabel("[CoA] (µM)")
    ax_res.set_ylabel("Residual (area)")
    ax_res.set_title("Absolute residuals")

    # Panel 3 — standardized residuals with ±3σ bounds
    ax_std.axhline(0, color="black", lw=0.8, ls="--")
    ax_std.axhline(3, color="red", lw=0.8, ls=":", label="±3σ")
    ax_std.axhline(-3, color="red", lw=0.8, ls=":")
    ax_std.scatter(concentrations, std_res, color="steelblue")
    ax_std.set_xlabel("[CoA] (µM)")
    ax_std.set_ylabel("Standardized residual")
    ax_std.set_title("Standardized residuals")
    ax_std.legend(fontsize=8)

    fig.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=300)
        _logger.info("Saved calibration plot to %s", output_path)
    if show:
        plt.show()
    plt.close(fig)


"""
No, refactored v1 has no equivalent of CalibrationCurve.plot().
What the original did
The original produced a 4-panel figure attached to the CalibrationCurve object
itself:
- Calibration curve with the linear fit overlaid
- Absolute residuals vs. concentration
- Standardized residuals with ±3σ bounds — the most diagnostically useful panel,
as it flags outlier standards and heteroscedasticity
- A parameter summary table (slope, intercept, R², RMSE, n_points, conc_range)

This was a self-contained QC artifact. Because it lived on CalibrationCurve, it
could be called immediately after fitting and saved independently of any
kinetics analysis.

What refactored v1 has instead
Nothing directly equivalent. Calibration is a pure data object with no plot
method. The only calibration-related plotting in refactored v1 is implicit —
calibrated velocities appear in the MM fit panels produced by plot() on
EnzymeKineticsAnalysis, but there is no standalone figure that lets you inspect
the calibration curve itself before applying it.

Whether the omission is meaningful
Yes, more so than the missing save(). The calibration plot serves a verification
role that no scalar metric fully replaces. R² and RMSE tell you the fit is good
in aggregate; the residual panels tell you where it is not — a single outlier
standard, a nonlinear response at the extremes of the working range, or unequal
variance across the concentration range. Skipping this check and proceeding
directly to apply_calibration() means any of those problems propagate silently
into every kcat and kcat/Km value in the dataset.

A few notes on this relative to the original's 4-panel design. The parameter
summary table panel is dropped here because the information it contained —
slope, intercept, R², RMSE — is better accessed programmatically from the
Calibration dataclass fields directly, and a text-in-axes table is awkward to
maintain. If you want it back, matplotlib's ax.table() can render the dataclass
fields, but it adds layout complexity for marginal gain. The raw calibration
data (concentrations and peak_areas) must be passed in explicitly, which is the
correct design for a stateless function — unlike the original where
CalibrationCurve held the data internally as instance state.
"""

# Legacy plotting function
#     def plot(self, output_path: str | Path | None = None) -> None:
#         """
#         Generate calibration curve plot with validation panels.
#
#         Parameters
#         ----------
#         output_path : str or Path, optional
#             Save plot to this path. If None, displays interactively.
#         """
#
#         if self.params is None or self.fit_data is None:
#             raise ValueError("No calibration fitted to plot.")
#
#         conc = self.fit_data['concentrations']
#         area_measured = self.fit_data['peak_areas_measured']
#         area_std = self.fit_data['peak_areas_std']
#         area_pred = self.fit_data['peak_areas_predicted']
#         residuals = self.fit_data['residuals']
#         linear_fn = self.fit_data['linear_function']
#
#         fig, axes = plt.subplots(2, 2, figsize=(14, 10))
#
#         # Plot 1: Calibration curve
#         ax = axes[0, 0]
#         ax.errorbar(conc, area_measured, yerr=area_std, fmt='o', markersize=10,
#                     capsize=6, capthick=2.5, color='steelblue', elinewidth=2,
#                     markeredgewidth=2.5, markeredgecolor='navy', label='Measurements', zorder=3)
#
#         conc_smooth = np.linspace(conc.min() - conc.max()*0.1, conc.max() * 1.15, 200)
#         area_fit = linear_fn(conc_smooth, self.params.slope, self.params.intercept)
#         ax.plot(conc_smooth, area_fit, 'r-', linewidth=3, label='Linear fit', zorder=2)
#
#         ax.set_xlabel('[CoA] (µM)', fontsize=12, fontweight='bold')
#         ax.set_ylabel('Peak Area', fontsize=12, fontweight='bold')
#         ax.set_title(f'Calibration Curve (R² = {self.params.r2:.6f})', fontsize=13, fontweight='bold')
#         ax.legend(fontsize=11)
#         ax.grid(True, alpha=0.3)
#
#         # Plot 2: Residuals
#         ax = axes[0, 1]
#         ax.errorbar(conc, residuals, yerr=area_std if area_std is not None else None,
#                     fmt='o', markersize=10, capsize=6, capthick=2.5, color='coral',
#                     elinewidth=2, markeredgewidth=2.5, markeredgecolor='darkred', zorder=3)
#         ax.axhline(y=0, color='black', linestyle='-', linewidth=2, zorder=2)
#
#         ax.set_xlabel('[CoA] (µM)', fontsize=12, fontweight='bold')
#         ax.set_ylabel('Residuals (area)', fontsize=12, fontweight='bold')
#         ax.set_title('Residual Plot', fontsize=13, fontweight='bold')
#         ax.grid(True, alpha=0.3)
#
#         # Plot 3: Normalized residuals
#         ax = axes[1, 0]
#         if area_std is not None and np.any(area_std > 0):
#             norm_residuals = residuals / area_std
#         else:
#             norm_residuals = residuals / np.std(residuals)
#
#         ax.scatter(conc, norm_residuals, s=150, color='green', alpha=0.6,
#                    edgecolors='darkgreen', linewidth=2.5, zorder=3)
#         ax.axhline(y=0, color='black', linestyle='-', linewidth=2, zorder=2)
#         ax.axhline(y=3, color='red', linestyle='--', linewidth=2, alpha=0.5, label='±3σ')
#         ax.axhline(y=-3, color='red', linestyle='--', linewidth=2, alpha=0.5)
#
#         ax.set_xlabel('[CoA] (µM)', fontsize=12, fontweight='bold')
#         ax.set_ylabel('Standardized Residuals', fontsize=12, fontweight='bold')
#         ax.set_title('Normalized Residuals', fontsize=13, fontweight='bold')
#         ax.legend(fontsize=11)
#         ax.grid(True, alpha=0.3)
#         ax.set_ylim([-4, 4])
#
#         # Plot 4: Summary
#         ax = axes[1, 1]
#         ax.axis('off')
#
#         summary = f"""
# CALIBRATION SUMMARY
#
# Linear Model:
#   peak_area = slope × [CoA] + intercept
#
# Parameters:
#   Slope:       {self.params.slope:.6f} ± {self.params.slope_std:.6f} area/µM
#   Intercept:   {self.params.intercept:.6f} ± {self.params.intercept_std:.6f}
#
# Fit Quality:
#   R²:          {self.params.r2:.8f}
#   RMSE:        {self.params.rmse:.4f} area
#   N points:    {self.params.n_points}
#
# Concentration Range:
#   Min:         {self.params.conc_range[0]:.2f} µM
#   Max:         {self.params.conc_range[1]:.2f} µM
#         """
#
#         ax.text(0.05, 0.95, summary, transform=ax.transAxes,
#                 fontsize=11, verticalalignment='top', family='monospace',
#                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
#
#         fig.suptitle('HPLC Calibration Curve Analysis', fontsize=14, fontweight='bold')
#         fig.tight_layout()
#
#         if output_path is not None:
#             fig.savefig(output_path, dpi=300, bbox_inches='tight')
#             print(f"✓ Calibration plot saved: {output_path}")
#         else:
#             plt.show()
#
#         plt.close(fig)
