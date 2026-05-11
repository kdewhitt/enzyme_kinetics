from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any
from collections.abc import Callable

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress

_logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Pure model functions — no classes, no state
# ---------------------------------------------------------------------


def michaelis_menten(s: np.ndarray, vmax: float, km: float) -> np.ndarray:
    return (vmax * s) / (km + s)


def hill_equation(s: np.ndarray, vmax: float, k_half: float, n: float) -> np.ndarray:
    """Hill equation. k_half is the half-saturation constant, equal to Km only when n == 1."""
    return (vmax * s**n) / (k_half**n + s**n)


def substrate_inhibition(
    s: np.ndarray, vmax: float, km: float, ki: float
) -> np.ndarray:
    return (vmax * s) / (km + s + (s**2) / ki)


def lineweaver_burk_transform(
    s: np.ndarray,
    v: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mask = (s > 0) & (v > 0)
    return 1.0 / s[mask], 1.0 / v[mask]


# ---------------------------------------------------------------------
# Fit result — frozen, model-agnostic
# ---------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class FitResult:
    """Result of a non-linear or linear kinetic fit.

    Attributes:
        model_func: The callable used to generate predictions.
        popt: Optimal parameters as returned by curve_fit (vmax, km[, ...]).
        perr: Standard errors on popt (sqrt of pcov diagonal).
        pcov: Full covariance matrix of popt; used for correlated error propagation.
        r_squared: Coefficient of determination on the fitted data.
        n_points: Number of data points used in the fit.
        model_type: Short tag identifying the model, e.g. "mm", "hill", "si", "lb".
        extra: Supplementary scalar values forwarded to to_dict() output (e.g. Ki, n).
    """

    model_func: Callable[..., np.ndarray]
    popt: tuple[float, ...]
    perr: tuple[float, ...]
    pcov: np.ndarray  # FIX 1: carry full covariance for correlated propagation
    r_squared: float
    n_points: int
    model_type: str = "mm"  # FIX 7: track model identity for correct property semantics
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def vmax(self) -> float:
        return self.popt[0]

    @property
    def vmax_se(self) -> float:
        return self.perr[0]

    @property
    def km(self) -> float:
        """Km for MM/SI models; k_half for Hill (half-saturation, not equal to Km when n≠1)."""
        return self.popt[1]

    @property
    def km_se(self) -> float:
        return self.perr[1]

    # FIX 7 — Hill-specific accessors; raise clearly if called on wrong model type.
    @property
    def k_half(self) -> float:
        """Half-saturation constant from a Hill fit (popt[1]). Raises if model is not 'hill'."""
        if self.model_type != "hill":
            raise AttributeError(
                f"k_half is only defined for Hill fits; this is '{self.model_type}'"
            )
        return self.popt[1]

    @property
    def k_half_se(self) -> float:
        """SE of k_half from a Hill fit. Raises if model is not 'hill'."""
        if self.model_type != "hill":
            raise AttributeError(
                f"k_half_se is only defined for Hill fits; this is '{self.model_type}'"
            )
        return self.perr[1]

    @property
    def hill_n(self) -> float:
        """Hill cooperativity coefficient (popt[2]). Raises if model is not 'hill'."""
        if self.model_type != "hill":
            raise AttributeError(
                f"hill_n is only defined for Hill fits; this is '{self.model_type}'"
            )
        return self.popt[2]

    @property
    def hill_n_se(self) -> float:
        """SE of Hill n (perr[2]). Raises if model is not 'hill'."""
        if self.model_type != "hill":
            raise AttributeError(
                f"hill_n_se is only defined for Hill fits; this is '{self.model_type}'"
            )
        return self.perr[2]

    # FIX 4 — typed Ki accessor for substrate-inhibition fits.
    @property
    def ki(self) -> float:
        """Substrate inhibition constant Ki (popt[2]). Raises if model is not 'si'."""
        if self.model_type != "si":
            raise AttributeError(
                f"ki is only defined for substrate-inhibition fits; this is '{self.model_type}'"
            )
        return self.popt[2]

    @property
    def ki_se(self) -> float:
        """SE of Ki (perr[2]). Raises if model is not 'si'."""
        if self.model_type != "si":
            raise AttributeError(
                f"ki_se is only defined for substrate-inhibition fits; this is '{self.model_type}'"
            )
        return self.perr[2]

    def predict(self, s: np.ndarray) -> np.ndarray:
        return self.model_func(s, *self.popt)


def _r_squared(observed: np.ndarray, predicted: np.ndarray) -> float:
    ss_res = np.sum((observed - predicted) ** 2)
    ss_tot = np.sum((observed - np.mean(observed)) ** 2)
    if ss_tot == 0:
        return 0.0
    return 1.0 - ss_res / ss_tot


# ---------------------------------------------------------------------
# Fitting — thin wrappers around curve_fit
# ---------------------------------------------------------------------


def fit_model(
    model_func: Callable[..., np.ndarray],
    s: np.ndarray,
    v: np.ndarray,
    p0: list[float],
    *,
    sigma: np.ndarray | None = None,
    bounds: tuple = (0, np.inf),
    maxfev: int = 10_000,
    model_type: str = "mm",
) -> FitResult:
    effective_sigma = sigma if sigma is not None and np.any(sigma > 0) else None
    popt, pcov = curve_fit(
        model_func,
        s,
        v,
        p0=p0,
        sigma=effective_sigma,
        absolute_sigma=True,
        bounds=bounds,
        maxfev=maxfev,
    )
    perr = tuple(np.sqrt(np.diag(pcov)).tolist())
    predicted = model_func(s, *popt)
    return FitResult(
        model_func=model_func,
        popt=tuple(popt.tolist()),
        perr=perr,
        pcov=pcov,  # FIX 1: store full covariance matrix
        r_squared=_r_squared(v, predicted),
        n_points=len(s),
        model_type=model_type,  # FIX 7: propagate model tag
    )


def fit_michaelis_menten(
    s: np.ndarray,
    v: np.ndarray,
    *,
    sigma: np.ndarray | None = None,
    vmax_init: float | None = None,
    km_init: float | None = None,
) -> FitResult:
    """Fit the Michaelis-Menten equation to velocity data.

    The fitted Km is in the same units as s. Throughout this codebase s is
    expected in µM, so the returned Km is in µM. Passing s in any other unit
    (mM, nM) will produce a Km in that unit with no warning; the caller is
    responsible for unit consistency.

    Args:
        s: Substrate concentration array in µM.
        v: Reaction velocity array (area/s before calibration; µM/s after).
        sigma: Per-point standard errors on v for weighted fitting. If None
            or all-zero, unweighted fitting is used.
        vmax_init: Initial guess for Vmax. Defaults to max(v).
        km_init: Initial guess for Km in µM. Defaults to median(s).

    Returns:
        FitResult with model_type='mm', Km in µM, and full pcov matrix.
    """
    vmax_guess = vmax_init if vmax_init is not None else float(np.max(v))
    km_guess = km_init if km_init is not None else float(np.median(s))
    return fit_model(
        michaelis_menten, s, v, [vmax_guess, km_guess], sigma=sigma, model_type="mm"
    )


def fit_hill(
    s: np.ndarray,
    v: np.ndarray,
    *,
    sigma: np.ndarray | None = None,
) -> FitResult:
    """Fit the Hill equation. The second parameter is k_half, not Km."""
    p0 = [float(np.max(v)), float(np.median(s)), 2.0]
    return fit_model(hill_equation, s, v, p0, sigma=sigma, model_type="hill")


def fit_substrate_inhibition(
    s: np.ndarray,
    v: np.ndarray,
    *,
    sigma: np.ndarray | None = None,
) -> FitResult:
    """Fit the substrate-inhibition model. Ki is accessible via FitResult.ki."""
    p0 = [float(np.max(v)), float(np.median(s)), float(np.max(s))]
    return fit_model(substrate_inhibition, s, v, p0, sigma=sigma, model_type="si")


def fit_lineweaver_burk(
    s: np.ndarray,
    v: np.ndarray,
) -> FitResult:
    s_inv, v_inv = lineweaver_burk_transform(s, v)
    slope, intercept, r_value, _, std_err = linregress(s_inv, v_inv)
    r2 = r_value**2
    vmax = 1.0 / intercept if intercept != 0 else np.nan
    km = slope * vmax if not np.isnan(vmax) else np.nan
    # LB std errors live in reciprocal space and do not map cleanly to parameter
    # space without the delta method, so perr is set to nan. pcov is set to a
    # 2x2 nan matrix to satisfy the FitResult contract without implying false precision.
    return FitResult(
        model_func=michaelis_menten,
        popt=(vmax, km),
        perr=(np.nan, np.nan),
        pcov=np.full((2, 2), np.nan),
        r_squared=r2,
        n_points=len(s_inv),
        model_type="lb",
        extra={"slope": slope, "intercept": intercept, "std_err": std_err},
    )


# ---------------------------------------------------------------------
# FIX 1 — kcat/Km error propagation using full covariance matrix
# ---------------------------------------------------------------------

_UM_TO_M: float = 1e-6  # unit conversion factor: µM → M


def _kcat_km_with_covariance(
    vmax: float,
    km: float,
    enzyme_conc_um: float,
    pcov: np.ndarray,
) -> tuple[float, float]:
    """Propagate kcat/Km uncertainty using the full Vmax–Km covariance matrix.

    The standard relative-error formula assumes Vmax and Km are uncorrelated,
    which is generally false for MM fitting (they are anti-correlated). This
    function uses the exact first-order delta method:

        Var(kcat/Km) = (∂f/∂Vmax)² · Var(Vmax)
                     + (∂f/∂Km)²   · Var(Km)
                     + 2·(∂f/∂Vmax)·(∂f/∂Km) · Cov(Vmax, Km)

    where f = Vmax / (enzyme_conc · Km).

    For Hill and substrate-inhibition fits the pcov indices are the same
    (Vmax at [0,0], Km/k_half at [1,1], Cov at [0,1]) so the formula applies
    uniformly.

    Unit convention: Km is in µM throughout the codebase. kcat/Km is reported
    in the standard literature unit of M⁻¹·s⁻¹, so Km is converted to M
    (×1e-6) before computing the ratio. The partial derivatives are scaled by
    the same factor so that the returned SE is also in M⁻¹·s⁻¹.

    Args:
        vmax: Fitted Vmax in µM/s (after calibration) or area/s (before).
        km: Fitted Km in µM (or k_half in µM for Hill fits).
        enzyme_conc_um: Total enzyme concentration in µM.
        pcov: Full 2×2 (or larger) covariance matrix from curve_fit, with
            Vmax at index [0,0] and Km at [1,1], covariance at [0,1].
            Units of pcov entries must be consistent with vmax and km units.

    Returns:
        Tuple of (kcat_km_M, kcat_km_se_M) in M⁻¹·s⁻¹. Returns (nan, nan)
        if km ≤ 0, vmax ≤ 0, or pcov contains nan (e.g. for LB fits).
    """
    if km <= 0 or vmax <= 0:
        return np.nan, np.nan
    if np.any(np.isnan(pcov[:2, :2])):
        # Fallback for LB fits where pcov is undefined
        return np.nan, np.nan

    # Convert Km from µM to M for the standard M⁻¹·s⁻¹ unit of kcat/Km
    km_M = km * _UM_TO_M
    kcat_km_M = vmax / (enzyme_conc_um * km_M)

    # Partial derivatives of f = vmax / (enzyme_conc_um * km_M)
    # where km_M = km * 1e-6, so ∂f/∂km = ∂f/∂km_M · ∂km_M/∂km = (∂f/∂km_M) * 1e-6
    df_dvmax = 1.0 / (enzyme_conc_um * km_M)
    df_dkm = -vmax / (enzyme_conc_um * km_M**2) * _UM_TO_M

    var_vmax = pcov[0, 0]
    var_km = pcov[1, 1]
    cov_vmax_km = pcov[0, 1]

    var_kcat_km = (
        df_dvmax**2 * var_vmax
        + df_dkm**2 * var_km
        + 2.0 * df_dvmax * df_dkm * cov_vmax_km
    )
    kcat_km_se_M = float(np.sqrt(max(var_kcat_km, 0.0)))  # clamp numerical negatives
    return float(kcat_km_M), kcat_km_se_M


# ---------------------------------------------------------------------
# Derived kinetic constants
# ---------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class KineticConstants:
    """Fully-typed kinetic constants for one peak.

    Attributes:
        peak_id: Peak identifier string.
        fit: Primary FitResult (MM, Hill, or SI model).
        kcat: Turnover number (s⁻¹ after calibration; otherwise in raw area
            units per µM enzyme per second).
        kcat_se: Standard error of kcat, in the same units as kcat.
        kcat_km_M: Catalytic efficiency kcat/Km in M⁻¹·s⁻¹ (after
            calibration). Km is converted from µM to M before computing the
            ratio so that the result is in the standard literature unit.
            Before calibration the numerator (kcat) is in non-physical units,
            so kcat_km_M should not be compared to literature values until
            calibration has been applied.
        kcat_km_se_M: Standard error of kcat_km_M in M⁻¹·s⁻¹, propagated
            via the full Vmax–Km covariance matrix.
        ki: Substrate inhibition constant Ki (µM); None unless model_type == 'si'.
        ki_se: Standard error of Ki in µM; None unless model_type == 'si'.
        lb_fit: Optional Lineweaver-Burk FitResult for cross-validation.
    """

    peak_id: str
    fit: FitResult
    kcat: float
    kcat_se: float
    kcat_km_M: float  # M⁻¹·s⁻¹; Km converted µM→M before division
    kcat_km_se_M: float  # SE in M⁻¹·s⁻¹
    ki: float | None = None
    ki_se: float | None = None
    lb_fit: FitResult | None = None

    def to_dict(self) -> dict[str, Any]:
        # FIX 7: use model-correct labels for Hill output
        is_hill = self.fit.model_type == "hill"
        base: dict[str, Any] = {
            "peak_id": self.peak_id,
            "model": self.fit.model_type,
            "vmax": self.fit.vmax,
            "vmax_se": self.fit.vmax_se,
            # Label the second parameter correctly per model
            ("k_half" if is_hill else "km"): self.fit.km,
            ("k_half_se" if is_hill else "km_se"): self.fit.km_se,
            "r_squared": self.fit.r_squared,
            "kcat": self.kcat,
            "kcat_se": self.kcat_se,
            # Explicit M⁻¹·s⁻¹ suffix in the column name makes the unit
            # unambiguous in downstream CSV/DataFrame consumers
            "kcat_km_M": self.kcat_km_M,
            "kcat_km_se_M": self.kcat_km_se_M,
        }
        # FIX 7: emit hill_n only for Hill fits
        if is_hill:
            base["hill_n"] = self.fit.hill_n
            base["hill_n_se"] = self.fit.hill_n_se
        # FIX 4: emit Ki only for SI fits
        if self.ki is not None:
            base["ki"] = self.ki
            base["ki_se"] = self.ki_se
        if self.lb_fit is not None:
            base["km_lb"] = self.lb_fit.km
            base["vmax_lb"] = self.lb_fit.vmax
            base["r_squared_lb"] = self.lb_fit.r_squared
        # Forward any remaining extra scalars (e.g. LB slope/intercept)
        for k, v in self.fit.extra.items():
            if isinstance(v, (int, float, str)) and k not in base:
                base[k] = v
        return base


def derive_constants(
    peak_id: str,
    fit: FitResult,
    enzyme_conc_um: float,
    *,
    lb_fit: FitResult | None = None,
) -> KineticConstants:
    """Compute kcat and kcat/Km from a FitResult.

    kcat is computed as Vmax / enzyme_conc_um. Its unit is s⁻¹ only after a
    calibration curve has been applied (so that Vmax is in µM/s); before
    calibration Vmax is in area/s and kcat carries non-physical units.

    kcat/Km is reported in M⁻¹·s⁻¹ (the standard literature unit) by
    converting Km from µM to M before computing the ratio. The conversion
    factor _UM_TO_M (1e-6) is applied inside _kcat_km_with_covariance, which
    also propagates uncertainty via the full Vmax–Km covariance matrix,
    correctly accounting for the anti-correlation between Vmax and Km that
    arises from MM-family non-linear fitting.

    Substrate concentrations (s) passed to the fitting functions must be in
    µM for Km and kcat/Km to carry their stated units. No runtime check is
    performed; the caller is responsible for unit consistency.

    Ki and ki_se are extracted from the fit and promoted to typed fields when
    the model is a substrate-inhibition fit (model_type == 'si').

    Args:
        peak_id: Identifier for the peak being processed.
        fit: FitResult from any kinetic model fit. Km (popt[1]) must be in µM.
        enzyme_conc_um: Total enzyme concentration in µM.
        lb_fit: Optional Lineweaver-Burk cross-validation fit.

    Returns:
        KineticConstants with kcat (s⁻¹ post-calibration), kcat_km_M
        (M⁻¹·s⁻¹), their SEs, and optionally Ki/ki_se for SI fits.
    """
    kcat = fit.vmax / enzyme_conc_um
    kcat_se = fit.vmax_se / enzyme_conc_um

    # FIX 1 + units: full covariance propagation; Km converted µM→M inside
    kcat_km_M, kcat_km_se_M = _kcat_km_with_covariance(
        fit.vmax,
        fit.km,
        enzyme_conc_um,
        fit.pcov,
    )

    # FIX 4: promote Ki to typed fields for SI fits
    ki: float | None = None
    ki_se: float | None = None
    if fit.model_type == "si":
        ki = fit.ki
        ki_se = fit.ki_se

    return KineticConstants(
        peak_id=peak_id,
        fit=fit,
        kcat=kcat,
        kcat_se=kcat_se,
        kcat_km_M=kcat_km_M,
        kcat_km_se_M=kcat_km_se_M,
        ki=ki,
        ki_se=ki_se,
        lb_fit=lb_fit,
    )


# ---------------------------------------------------------------------
# Calibration — linear standard curve
# FIX 2, FIX 6
# ---------------------------------------------------------------------


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
    conc_range: tuple[float, float]

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
        dc_ds = -(area - self.intercept) / (self.slope**2)
        dc_di = -1.0 / self.slope  # FIX 2: intercept term
        conc_se = float(
            np.sqrt(
                (dc_da * area_se) ** 2
                + (dc_ds * self.slope_se) ** 2
                + (dc_di * self.intercept_se) ** 2  # FIX 2
            )
        )
        return conc, conc_se


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
        _linear,
        concentrations,
        peak_areas,
        sigma=effective_sigma,
        absolute_sigma=True,
    )
    perr = np.sqrt(np.diag(pcov))
    predicted = _linear(concentrations, *popt)
    return Calibration(
        slope=popt[0],
        slope_se=perr[0],
        intercept=popt[1],
        intercept_se=perr[1],
        r_squared=_r_squared(peak_areas, predicted),
        conc_range=(float(concentrations.min()), float(concentrations.max())),
    )


# ---------------------------------------------------------------------
# Velocity preparation from summary-statistics DataFrames
# FIX 3: return raw mean_signal alongside velocity so apply_calibration
#        can apply the calibration directly to area before dividing by rxn_time.
# ---------------------------------------------------------------------


def prepare_velocity(
    s: np.ndarray,
    mean_signal: np.ndarray,
    std_signal: np.ndarray,
    counts: np.ndarray,
    rxn_time: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute substrate concentrations, reaction velocities, and their SEMs.

    FIX 3: The raw mean_signal (area) is returned as the fourth element so that
    callers can apply a calibration curve directly to the area before converting
    to velocity, avoiding the fragile round-trip of reconstructing area from
    velocity inside apply_calibration.

    Args:
        s: Substrate concentration array (µM).
        mean_signal: Mean HPLC peak area per replicate group.
        std_signal: Standard deviation of peak area per replicate group.
        counts: Number of replicates per group.
        rxn_time: Reaction duration in seconds.

    Returns:
        Tuple of (s, v, v_sem, mean_signal) where:
            s        — substrate concentrations (µM, unchanged).
            v        — reaction velocity in area/time units (area / rxn_time).
            v_sem    — SEM of v = (std / rxn_time) / sqrt(n).
            mean_signal — raw peak area array, for direct calibration application.
    """
    v = mean_signal / rxn_time
    v_sem = (std_signal / rxn_time) / np.sqrt(counts)
    return s, v, v_sem, mean_signal


# ---------------------------------------------------------------------
# Batch analysis helper
# ---------------------------------------------------------------------


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
                peak_id,
                mm_fit,
                enzyme_conc_um,
                lb_fit=lb_fit,
            )
        except Exception as exc:
            _logger.warning("Fit failed for %s: %s", peak_id, exc)

    return results
