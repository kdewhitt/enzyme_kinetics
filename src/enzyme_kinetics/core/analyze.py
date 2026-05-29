"""Kinetic analysis orchestrator for fitting enzyme kinetic models to HPLC data.

Provides KineticAnalyzer, a stateful orchestrator that fits Michaelis-Menten (or
Hill) models to per-peak substrate-velocity data, optionally re-fits using a
calibrated concentration scale, and exports results to CSV. The module also exposes
analyze_peaks, a lower-level batch helper for fitting pre-extracted arrays directly.

Model selection is handled by _select_and_fit, which dispatches to fit_hill for
peaks listed in special_peaks and defaults to fit_michaelis_menten otherwise.
Lineweaver-Burk cross-validation fits are attempted for every peak and stored on
the KineticConstants result; failures are logged and silently skipped.

Typical usage example:
    >>> analyzer = KineticAnalyzer(
    ...     enzyme_conc_um=0.5,
    ...     reaction_time_seconds=3600.0,
    ...     substrate="HexCoA",
    ...     data=df,
    ... )
    >>> analyzer.fit().apply_calibration(cal).export_csv(dest)
"""

from __future__ import annotations

from pathlib import Path
from typing import Self

import numpy as np
import pandas as pd
from kgdlibs.pathtools import PathBuilder
from logurich import RichLogAdapter

from .calibration import Calibration
from .derive import derive_constants, KineticConstants
from .models import fit_hill, fit_lineweaver_burk, fit_michaelis_menten, FitResult
from .plots import KineticPlots
from .preprocess import extract_peak_data

__all__ = ["analyze_peaks", "KineticAnalyzer"]

_logger = RichLogAdapter(component=__name__)


# ---------------------------------------------------------------------
# Model selection strategy
# ---------------------------------------------------------------------


def _select_and_fit(
    peak_id: str,
    s: np.ndarray,
    v: np.ndarray,
    v_std: np.ndarray,
    *,
    special_peaks: dict[str, str] | None = None,
) -> FitResult:
    """Selects and fits the appropriate kinetic model for a single peak.

    Dispatches to fit_hill when peak_id appears in special_peaks with value
    "hill". All other peaks default to fit_michaelis_menten. Additional model
    types can be added to the dispatch block as needed.

    Args:
        peak_id: Identifier for the HPLC peak being fitted.
        s: Substrate concentration array in µM.
        v: Reaction velocity array (area/s before calibration; µM/s after).
        v_std: Per-point standard deviations of v in the same units as v.
        special_peaks: Optional mapping of peak_id to model type string. Only
            "hill" is currently dispatched; unrecognized values fall back to
            Michaelis-Menten. Defaults to None.

    Returns:
        FitResult from the selected model fit for the given peak.
    """
    special = (special_peaks or {}).get(peak_id)
    if special == "hill":
        return fit_hill(s, v, sigma=v_std)

    # if peak_id.lower() == "olv":
    #     return fit_threshold_michaelis_menten(s, v, sigma=v_std)
    return fit_michaelis_menten(s, v, sigma=v_std)


# ---------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------


class KineticAnalyzer:
    """Stateful orchestrator for per-peak enzyme kinetic analysis.

    Fits kinetic models to substrate-velocity data extracted from a DataFrame
    of HPLC peak measurements, optionally re-fits on calibrated concentration
    velocities, and exports results. The typical workflow is:
    fit() → apply_calibration() → plot() → export_csv().

    fit() populates self.results from raw area/s velocities. apply_calibration()
    rebuilds self.results using µM/s velocities derived from a Calibration object,
    making kcat and kcat/Km physically meaningful (s⁻¹ and s⁻¹·M⁻¹ respectively).
    Calling fit() again after apply_calibration() resets results to raw-area fits.

    Attributes:
        enzyme_conc_um: Total enzyme concentration used for kcat calculation in µM.
        reaction_time_seconds: Reaction duration used to convert peak area to
            velocity in seconds.
        substrate: Substrate name propagated to KineticConstants results.
        df: Source DataFrame containing peak_id, substrate concentration, and
            replicate signal columns.
        special_peaks: Mapping of peak_id to model type string for non-default
            model dispatch (e.g. {"peakA": "hill"}).
        results: Mapping of peak_id to KineticConstants populated by fit() or
            apply_calibration().
        calibration: Most recently applied Calibration object, or None if
            apply_calibration() has not been called.
    """

    def __init__(
        self,
        enzyme_conc_um: float,
        reaction_time_seconds: float,
        substrate: str,
        data: pd.DataFrame,
        special_peaks: dict[str, str] | None = None,
    ) -> None:
        """Initializes the KineticAnalyzer with experimental parameters and data.

        Args:
            enzyme_conc_um: Total enzyme concentration in µM used to compute kcat.
            reaction_time_seconds: Reaction duration in seconds used to convert
                integrated peak area to velocity (area/s).
            substrate: Substrate name string propagated to all KineticConstants
                results produced by this analyzer.
            data: DataFrame containing HPLC peak data with at minimum a "peak_id"
                column, substrate concentration values, and replicate signal columns
                as expected by extract_peak_data.
            special_peaks: Optional mapping of peak_id to model type string
                controlling model dispatch in _select_and_fit. Defaults to None
                (all peaks fitted with Michaelis-Menten).
        """
        self.enzyme_conc_um = enzyme_conc_um
        self.reaction_time_seconds = reaction_time_seconds
        self.substrate = substrate
        self.df = data
        self.special_peaks = special_peaks or {}

        self.results: dict[str, KineticConstants] = {}
        self.calibrations: dict[str, Calibration] = {}

    def fit(self) -> Self:
        """Fits kinetic models to all peaks in self.df using raw area/s velocities.

        Clears self.results before fitting, so calling fit() a second time replaces
        all previous results. For each unique peak_id, extracts substrate-velocity
        data via extract_peak_data, selects and fits the primary kinetic model via
        _select_and_fit, and attempts a Lineweaver-Burk cross-validation fit.
        Lineweaver-Burk failures are silently swallowed; primary model failures are
        logged as warnings and the peak is skipped.

        kcat and kcat/Km computed from raw-area fits carry non-physical units
        (area · µM⁻¹ · s⁻¹ and area · M⁻¹ · s⁻¹ respectively) until
        apply_calibration() is called.

        Returns:
            The KineticAnalyzer instance, enabling method chaining.
        """
        self.results.clear()
        for peak_id in self.df["peak_id"].unique():
            data = extract_peak_data(self.df, peak_id, self.reaction_time_seconds)
            if data is None:
                continue
            s, v, v_std, _ = data  # raw area not needed at fit stage
            try:
                mm_fit = _select_and_fit(
                    peak_id,
                    s,
                    v,
                    v_std,
                    special_peaks=self.special_peaks,
                )
                lb_fit: FitResult | None = None
                try:
                    lb_fit = fit_lineweaver_burk(s, v)
                except Exception:
                    _logger.warning("LB fit failed for %s", peak_id)
                self.results[peak_id] = derive_constants(
                    self.substrate,
                    peak_id,
                    mm_fit,
                    self.enzyme_conc_um,
                    lb_fit=lb_fit,
                )
                _logger.success(f"{peak_id}: Km={self.results[peak_id].fit.km}")

            except Exception as exc:
                _logger.warning("Fit failed for %s: %s", peak_id, exc)
        return self

    def apply_calibration(
        self,
        cal: Calibration,
        *,
        precise_std: bool = True,
        peak_ids: list[str] | None = None,
    ) -> Self:
        """Re-fits peaks using calibrated µM/s velocities derived from a Calibration.

        Converts raw mean peak area to µM concentration via the calibration curve,
        then divides by reaction_time_seconds to produce velocity in µM/s. This
        makes kcat (s⁻¹) and kcat/Km (s⁻¹·M⁻¹) physically meaningful for
        comparison with literature values.

        When peak_ids is None, all peaks are recalibrated and self.results is
        replaced entirely. When peak_ids is provided, only those peaks are
        recalibrated and self.results is updated in-place for those keys only;
        all other peaks retain their existing results. Model selection respects
        self.special_peaks via _select_and_fit in both modes.

        Two propagation modes are available. The precise mode (default) applies
        the full delta-method per data point, propagating both slope and intercept
        uncertainty from the calibration curve. The simplified mode divides v_std
        by the calibration slope only, which is appropriate when calibration
        parameter uncertainty is negligible relative to replicate variance.

        Rebuilds self.results from self.df on every call, making repeated calls safe
        when the calibration is updated between calls.

        Args:
            cal: Fitted Calibration object produced by fit_calibration(). Must have
                been fitted before passing; area_to_conc and area_to_conc_with_error
                are called on each peak's mean signal array.
            precise_std: If True, applies full delta-method propagation including
                slope and intercept uncertainty from the calibration curve (recommended).
                If False, uses the simplified approximation v_std / slope. Defaults
                to True.
            peak_ids: Optional list of peak_ids to recalibrate. If provided, only
                those peaks are re-fitted and merged into self.results; all other
                peaks retain their current results. All supplied peak_ids must already
                exist in self.results. Defaults to None.

        Returns:
            The KineticAnalyzer instance, enabling method chaining.

        Raises:
            ValueError: If peak_ids contains any peak_id not present in self.results.

        Note:
            A warning is logged if cal.is_valid() returns False (R² below threshold
            or non-positive slope), but fitting proceeds. Inspect calibration quality
            before interpreting kcat and kcat/Km values.
        """
        # Guard with calibration validity check
        if not cal.is_valid():
            _logger.warning(
                "Calibration does not meet quality thresholds "
                "(R²=%.4f, slope=%.4g). Results may be unreliable.",
                cal.r_squared,
                cal.slope,
            )

        if peak_ids is not None:
            unknown = set(peak_ids) - self.results.keys()
            if unknown:
                raise ValueError(
                    f"peak_ids not found in results (call fit() first): {sorted(unknown)}",
                )

        targets = (
            {k: v for k, v in self.results.items() if k in peak_ids}
            if peak_ids is not None
            else self.results
        )
        recalculated: dict[str, KineticConstants] = {}

        for peak_id, kc in targets.items():
            data = extract_peak_data(self.df, peak_id, self.reaction_time_seconds)
            if data is None:
                continue

            # Unpack raw area; apply calibration directly to area
            s, _v_raw, v_std_raw, mean_signal = data

            # Convert area → µM concentration, then divide by rxn_time for velocity
            v_um = np.asarray(cal.area_to_conc(mean_signal)) / self.reaction_time_seconds

            # Precise vs. simplified std propagation
            if precise_std:
                # Full delta-method per data point: propagates slope_std and intercept_std
                v_std_um = np.array(
                    [
                        cal.area_to_conc_with_error(
                            float(area),
                            float(area_std * self.reaction_time_seconds),
                        )[1]
                        / self.reaction_time_seconds
                        for area, area_std in zip(
                        mean_signal,
                        v_std_raw,
                    )
                    ],
                )
            else:
                # Simplified approximation: ignores calibration parameter uncertainty
                v_std_um = v_std_raw / cal.slope

            self.calibrations[peak_id] = cal
            try:
                mm_fit = _select_and_fit(
                    peak_id,
                    s,
                    v_um,
                    v_std_um,
                    special_peaks=self.special_peaks,
                )
                recalculated[peak_id] = derive_constants(
                    self.substrate,
                    peak_id,
                    mm_fit,
                    self.enzyme_conc_um,
                    lb_fit=kc.lb_fit,
                )
            except Exception as exc:
                _logger.warning("Calibrated re-fit failed for %s: %s", peak_id, exc)

        if peak_ids is not None:
            self.results.update(recalculated)
        else:
            self.results = recalculated

        return self

    def plot(self, dest: Path, *, is_calibrated: bool = False, overwrite: bool = False) -> Self:
        """Generates and saves all kinetic plots for the current results.

        Constructs a KineticPlots instance and renders the Michaelis-Menten curves,
        Lineweaver-Burk plots, residual plots, and efficiency comparison panel to
        dest. Plots are written to disk as a side effect.

        Args:
            dest: Base file path to use for writing plot files. A specific suffix
                will be appended to the path stem to identify each particular
                plot created.
            is_calibrated: If True, plots are generated assuming calibrated data. If False,
                plots are generated assuming uncorrected "raw" data. Defaults to False.
            overwrite: If True, overwrites existing plot files at dest. If False,
                existing files are preserved and new files receive a unique suffix.
                Defaults to False.

        Returns:
            The KineticAnalyzer instance, enabling method chaining.
        """
        plotter = KineticPlots(
            dest,
            data=self.df,
            results=self.results,
            reaction_time_seconds=self.reaction_time_seconds,
            overwrite=overwrite,
            calibrated=is_calibrated,
            calibrations=self.calibrations or None,
        )
        (
            plotter.plot()
            .plot_lineweaver_burk()
            .plot_residuals()
            .plot_efficiency_comparison(exclude="OLV")
        )
        return self

    def to_dataframe(self) -> pd.DataFrame:
        """Converts current results to a DataFrame with one row per peak."""
        return pd.DataFrame([kc.to_dict() for kc in self.results.values()])

    def export_csv(self, dest: Path, *, overwrite: bool = False) -> Self:
        """Exports current kinetics results to a tagged CSV file.

        Converts self.results to a DataFrame via to_dataframe() and writes it to
        dest with a "kinetics" tag appended to the filename via PathBuilder. If
        results are empty, logs a warning and returns without writing.

        Args:
            dest: Destination path for the CSV file. The final filename is
                constructed by PathBuilder with a "kinetics" tag and ".csv" suffix.
            overwrite: If True, overwrites an existing file at the resolved path.
                If False, a unique suffix is appended to avoid collisions.
                Defaults to False.

        Returns:
            The KineticAnalyzer instance, enabling method chaining.
        """
        stats_df = self.to_dataframe()
        if stats_df.empty:
            _logger.warning("No kinetics data to export.")
            return self

        out_path = (
            PathBuilder.from_path(
                dest,
                suffix=".csv",
                overwrite=overwrite,
            )
            .with_tag("kinetics")
            .path
        )
        stats_df.to_csv(out_path, index=False)
        _logger.info("Exported %s rows to %s", len(stats_df), out_path)
        return self


# ---------------------------------------------------------------------
# Batch analysis helper - unsure of its utility given KineticAnalyzer
# ---------------------------------------------------------------------


def analyze_peaks(
    substrate: str,
    substrate_conc: dict[str, np.ndarray],
    velocity: dict[str, np.ndarray],
    enzyme_conc_um: float,
    *,
    velocity_std: dict[str, np.ndarray] | None = None,
    fit_lb: bool = True,
    min_points: int = 3,
) -> dict[str, KineticConstants]:
    """Fits Michaelis-Menten models to pre-extracted substrate-velocity arrays.

    Iterates over peak_ids present in substrate_conc, skipping any peak with
    fewer than min_points data points or an all-zero velocity array. For each
    qualifying peak, fits the Michaelis-Menten model and optionally a
    Lineweaver-Burk cross-validation fit. Primary fit failures are logged as
    warnings; Lineweaver-Burk failures are silently skipped.

    kcat and kcat/Km computed from raw-area velocities carry non-physical units
    until a calibration has been applied upstream by the caller.

    Args:
        substrate: Substrate name propagated to KineticConstants results.
        substrate_conc: Mapping of peak_id to substrate concentration array in µM.
        velocity: Mapping of peak_id to reaction velocity array (area/s before
            calibration; µM/s after).
        enzyme_conc_um: Total enzyme concentration in µM used to compute kcat.
        velocity_std: Optional mapping of peak_id to per-point standard errors on
            velocity in the same units as velocity. If absent for a peak, unweighted
            fitting is used for that peak. Defaults to None.
        fit_lb: If True, attempts a Lineweaver-Burk cross-validation fit for each
            peak and stores it on the KineticConstants result. Defaults to True.
        min_points: Minimum number of data points required to attempt fitting.
            Peaks with fewer points are skipped. Defaults to 3.

    Returns:
        A mapping of peak_id to KineticConstants for all successfully fitted peaks.
        Peaks that fail fitting or are skipped are absent from the result.

    Note:
        This function operates on pre-extracted arrays and does not integrate with
        KineticAnalyzer's calibration workflow. For end-to-end analysis including
        calibration, prefer KineticAnalyzer.
    """
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
                substrate,
                peak_id,
                mm_fit,
                enzyme_conc_um,
                lb_fit=lb_fit,
            )
        except Exception as exc:
            _logger.warning("Fit failed for %s: %s", peak_id, exc)

    return results
