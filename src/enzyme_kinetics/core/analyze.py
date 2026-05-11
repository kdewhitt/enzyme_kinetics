from __future__ import annotations

from pathlib import Path
from typing import Self

import numpy as np
import pandas as pd
from kgdlibs.pathtools import PathBuilder
from logurich import RichLogAdapter

from .calibration import Calibration
from .derive import derive_constants, KineticConstants
from .models import fit_hill, fit_lineweaver_burk, fit_michaelis_menten, fit_threshold_michaelis_menten, FitResult
from .plots import KineticPlots
from .utils import extract_peak_data

_logger = RichLogAdapter(component=__name__)


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

    # if peak_id.lower() == "olv":
    #     return fit_threshold_michaelis_menten(s, v, sigma=v_sem)
    return fit_michaelis_menten(s, v, sigma=v_sem)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


class KineticAnalyzer:
    def __init__(
        self,
        enzyme_conc_um: float,
        reaction_time_seconds: float,
        substrate: str,
        data: pd.DataFrame,
        special_peaks: dict[str, str] | None = None,
    ) -> None:
        self.enzyme_conc_um = enzyme_conc_um
        self.reaction_time_seconds = reaction_time_seconds
        self.substrate = substrate
        self.df = data
        self.special_peaks = special_peaks or {}

        self.results: dict[str, KineticConstants] = {}
        self.calibration: Calibration | None = None

    def fit(self) -> Self:
        self.results.clear()
        for peak_id in self.df["peak_id"].unique():
            data = extract_peak_data(self.df, peak_id, self.reaction_time_seconds)
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
                    self.substrate,peak_id, mm_fit, self.enzyme_conc_um, lb_fit=lb_fit,
                )
                _logger.success(f"✓ {peak_id}: Km={self.results[peak_id].fit.km}")

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
            data = extract_peak_data(self.df, peak_id, self.reaction_time_seconds)
            if data is None:
                continue

            # FIX 3: unpack raw area; apply calibration directly to area
            s, _v_raw, v_sem_raw, mean_signal = data

            # Convert area → µM concentration, then divide by rxn_time for velocity
            v_um = np.asarray(cal.area_to_conc(mean_signal)) / self.reaction_time_seconds

            # FIX 4: precise vs. simplified SEM scaling
            if precise_sem:
                # Full delta-method per data point: propagates slope_se and intercept_se
                v_sem_um = np.array(
                    [
                        cal.area_to_conc_with_error(float(area), float(area_se * self.reaction_time_seconds))[1]
                        / self.reaction_time_seconds
                        for area, area_se in zip(mean_signal, v_sem_raw * self.reaction_time_seconds)
                    ],
                )
            else:
                # Simplified approximation: ignores calibration parameter uncertainty
                v_sem_um = v_sem_raw / cal.slope

            try:
                mm_fit = fit_michaelis_menten(s, v_um, sigma=v_sem_um)
                recalculated[peak_id] = derive_constants(
                    peak_id, mm_fit, self.enzyme_conc_um, lb_fit=kc.lb_fit,
                )
            except Exception as exc:
                _logger.warning("Calibrated re-fit failed for %s: %s", peak_id, exc)

        self.results = recalculated
        return self

    def plot(self, dest: Path, *, overwrite: bool = False) -> Self:
        plotter = KineticPlots(
            dest,
            data=self.df,
            results=self.results,
            reaction_time_seconds=self.reaction_time_seconds,
            overwrite=overwrite,
        )
        (
            plotter.plot()
            .plot_lineweaver_burk()
            .plot_residuals()
            .plot_efficiency_comparison()
        )
        return self

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame([kc.to_dict() for kc in self.results.values()])

    def export_csv(self, dest: Path, *, overwrite: bool = False) -> Self:
        stats_df = self.to_dataframe()
        if stats_df.empty:
            _logger.warning("No kinetics data to export.")
            return self

        out_path = PathBuilder.from_path(
            dest,
            suffix=".csv",
            overwrite=overwrite,
        ).with_tag("kinetics").path
        stats_df.to_csv(out_path, index=False)
        _logger.info("Exported %s rows to %s", len(stats_df), out_path)
        return self


# ---------------------------------------------------------------------------
# Batch analysis helper
# ---------------------------------------------------------------------------

# Unsure of the utility of this function

def analyze_peaks(
    substrate_conc: dict[str, np.ndarray],
    velocity: dict[str, np.ndarray],
    enzyme_conc_um: float,
    *,
    velocity_std: dict[str, np.ndarray] | None = None,
    fit_lb: bool = True,
    min_points: int = 3,
) -> dict[str, KineticConstants]:
    results: dict[str, KineticConstants] = {}
    velocity_std = velocity_std or {}

    for peak_id in substrate_conc:
        s = substrate_conc[peak_id]
        v = velocity[peak_id]
        if len(s) < min_points or np.all(v == 0):
            _logger.info("Skipping %s: insufficient data", peak_id)
            continue
        sigma = velocity_std.get(peak_id)
        try:
            mm_fit = fit_michaelis_menten(s, v, sigma=sigma)
            lb_fit = None
            if fit_lb:
                try:
                    lb_fit = fit_lineweaver_burk(s, v)
                except Exception:
                    pass
            results[peak_id] = derive_constants(
                peak_id, mm_fit, enzyme_conc_um, lb_fit=lb_fit,
            )
        except Exception as exc:
            _logger.warning("Fit failed for %s: %s", peak_id, exc)

    return results
