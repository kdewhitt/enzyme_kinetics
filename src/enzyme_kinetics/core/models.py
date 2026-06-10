"""Kinetic model functions and fitting utilities for enzyme kinetics analysis.

Provides pure mathematical functions for the Michaelis-Menten, Hill, threshold
Michaelis-Menten, substrate-inhibition, and Lineweaver-Burk models, together with
thin fitting wrappers around scipy.optimize.curve_fit. All fitting functions return
a FitResult dataclass that carries the full parameter covariance matrix for
downstream correlated error propagation.

Substrate concentrations (s) passed to any fitting function must be in µM for Km
and derived constants to carry their stated units. No runtime unit check is
performed; the caller is responsible for unit consistency.

Typical usage:
    import numpy as np
    s = np.array([1.0, 5.0, 10.0, 50.0, 100.0])
    v = np.array([0.1, 0.3, 0.5, 0.8, 0.9])
    result = fit_michaelis_menten(s, v)
    print(result.km, result.vmax)
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from logurich import DuoLogAdapter
from scipy.optimize import curve_fit
from scipy.stats import linregress

_logger = DuoLogAdapter.create(component=__name__)

__all__ = [
    "fit_hill",
    "fit_lineweaver_burk",
    "fit_michaelis_menten",
    "fit_model",
    "fit_substrate_inhibition",
    "fit_threshold_michaelis_menten",
    "FitResult",
    "hill_equation",
    "lineweaver_burk_transform",
    "michaelis_menten",
    "r_squared",
    "substrate_inhibition",
    "threshold_michaelis_menten",
]


def michaelis_menten(s: np.ndarray, vmax: float, km: float) -> np.ndarray:
    """Computes reaction velocity using the Michaelis-Menten equation.

    Equation: v = (Vmax · [S]) / (Km + [S])

    Args:
        s: Substrate concentration array in µM.
        vmax: Maximum reaction velocity (area/s before calibration; µM/s after).
        km: Michaelis constant in µM.

    Returns:
        Reaction velocity array in the same units as vmax.
    """
    return (vmax * s) / (km + s)


def hill_equation(s: np.ndarray, vmax: float, k_half: float, n: float) -> np.ndarray:
    """Computes reaction velocity using the Hill equation.

    k_half is the half-saturation constant — the substrate concentration at which
    v = Vmax / 2. k_half equals Km only when n == 1 (non-cooperative case); for
    n > 1 (positive cooperativity) k_half underestimates the intrinsic Km, and
    for n < 1 (negative cooperativity) it overestimates it.

    Equation: v = (Vmax · [S]ⁿ) / (k_half ⁿ + [S]ⁿ)

    Args:
        s: Substrate concentration array in µM.
        vmax: Maximum reaction velocity (area/s before calibration; µM/s after).
        k_half: Half-saturation constant in µM.
        n: Hill cooperativity coefficient (dimensionless). n > 1 indicates positive
            cooperativity; n < 1 indicates negative cooperativity.

    Returns:
        Reaction velocity array in the same units as vmax.
    """
    return (vmax * s ** n) / (k_half ** n + s ** n)


def threshold_michaelis_menten(
    s: np.ndarray,
    vmax: float,
    km: float,
    s0: float,
) -> np.ndarray:
    """Computes reaction velocity using a Michaelis-Menten model with a dead-zone threshold.

    Velocity is zero for all substrate concentrations at or below S0; above S0 the
    model follows standard Michaelis-Menten kinetics with the effective substrate
    concentration (s − S0).

    Equation: v = 0 if s ≤ S0, else (Vmax · (s − S0)) / (Km + (s − S0))

    Args:
        s: Substrate concentration array in µM.
        vmax: Maximum reaction velocity (area/s before calibration; µM/s after).
        km: Michaelis constant for the effective substrate concentration in µM.
        s0: Dead-zone threshold concentration in µM. No reaction occurs at or below
            this value.

    Returns:
        Reaction velocity array in the same units as vmax.
    """
    return np.where(s > s0, (vmax * (s - s0)) / (km + (s - s0)), 0.0)


def substrate_inhibition(
    s: np.ndarray,
    vmax: float,
    km: float,
    ki: float,
) -> np.ndarray:
    """Computes reaction velocity under substrate inhibition.

    Describes kinetics where excess substrate reduces velocity, producing a
    bell-shaped v vs. [S] curve. Inhibition becomes significant when [S] approaches
    Ki; at [S] >> Ki velocity declines toward zero.

    Equation: v = (Vmax · [S]) / (Km + [S] + [S]²/Ki)

    Args:
        s: Substrate concentration array in µM.
        vmax: Maximum reaction velocity in the absence of inhibition
            (area/s before calibration; µM/s after).
        km: Michaelis constant in µM.
        ki: Substrate inhibition constant in µM. Lower values indicate stronger
            inhibition at a given substrate concentration.

    Returns:
        Reaction velocity array in the same units as vmax.
    """
    return (vmax * s) / (km + s + (s ** 2) / ki)


def lineweaver_burk_transform(
    s: np.ndarray,
    v: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Transforms substrate concentration and velocity arrays into Lineweaver-Burk space.

    Computes the double-reciprocal transformation (1/[S], 1/v) required for a
    Lineweaver-Burk plot. Data points where s ≤ 0 or v ≤ 0 are excluded via a
    boolean mask before inversion to prevent division-by-zero and logarithmic
    artifacts. Error amplification is strongest at low [S], where small absolute
    errors in v produce large errors in 1/v.

    Args:
        s: Substrate concentration array in µM. Non-positive values are excluded.
        v: Reaction velocity array (area/s before calibration; µM/s after).
            Non-positive values are excluded.

    Returns:
        A tuple of (s_inv, v_inv) where s_inv contains 1/[S] values in µM⁻¹ and
        v_inv contains 1/v values, both filtered to the valid (positive) subset.
    """
    mask = (s > 0) & (v > 0)
    return 1.0 / s[mask], 1.0 / v[mask]


@dataclass(frozen=True, slots=True)
class FitResult:
    """Result of a non-linear or linear kinetic fit.

    Attributes:
        model_func: The callable used to generate predictions from popt.
        popt: Optimal parameters as returned by curve_fit (vmax, km[, ...]).
        perr: Standard errors on popt derived from the square root of the pcov diagonal.
        pcov: Full covariance matrix of popt in parameter-space units; used for
            correlated error propagation (e.g. kcat/Km via the delta method).
        r_squared: Coefficient of determination on the fitted data (dimensionless).
        n_points: Number of data points used in the fit.
        model_type: Short tag identifying the model ("mm", "hill", "si", "lb", "tmm").
        extra: Supplementary scalar values forwarded to downstream outputs (e.g. LB
            slope, intercept, std_err for Lineweaver-Burk fits).
    """

    model_func: Callable[..., np.ndarray]
    popt: tuple[float, ...]
    perr: tuple[float, ...]
    pcov: np.ndarray
    r_squared: float
    n_points: int
    model_type: str = "mm"
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def vmax(self) -> float:
        """Maximum reaction velocity Vmax from popt[0] (area/s before calibration; µM/s after)."""
        return self.popt[0]

    @property
    def vmax_std(self) -> float:
        """Standard deviation of Vmax in the same units as vmax."""
        return self.perr[0]

    @property
    def km(self) -> float:
        """Michaelis constant Km in µM (popt[1]).

        Note: Corresponds to k_half for Hill fits (half-saturation, not Km when n≠1).
        """
        return self.popt[1]

    @property
    def km_std(self) -> float:
        """Standard deviation of Km (or k_half for Hill fits) in µM."""
        return self.perr[1]

    @property
    def k_half(self) -> float:
        """Half-saturation constant k_half in µM from a Hill fit (popt[1]).

        Raises if model is not "hill".
        """
        if self.model_type != "hill":
            raise AttributeError(
                f"k_half is only defined for Hill fits; this is '{self.model_type}'",
            )
        return self.popt[1]

    @property
    def k_half_std(self) -> float:
        """Standard deviation of k_half in µM from a Hill fit.

        Raises if model is not "hill".
        """
        if self.model_type != "hill":
            raise AttributeError(
                f"k_half_std is only defined for Hill fits; this is '{self.model_type}'",
            )
        return self.perr[1]

    @property
    def hill_n(self) -> float:
        """Hill cooperativity coefficient n (dimensionless) from popt[2].

        Raises if model is not "hill".
        """
        if self.model_type != "hill":
            raise AttributeError(
                f"hill_n is only defined for Hill fits; this is '{self.model_type}'",
            )
        return self.popt[2]

    @property
    def hill_n_std(self) -> float:
        """Standard deviation of the Hill cooperativity coefficient n (dimensionless).

        Raises if model is not "hill".
        """
        if self.model_type != "hill":
            raise AttributeError(
                f"hill_n_std is only defined for Hill fits; this is '{self.model_type}'",
            )
        return self.perr[2]

    @property
    def ki(self) -> float:
        """Substrate inhibition constant Ki in µM from popt[2].

        Raises if model is not "si".
        """
        if self.model_type != "si":
            raise AttributeError(
                f"ki is only defined for substrate-inhibition fits; this is '{self.model_type}'",
            )
        return self.popt[2]

    @property
    def ki_std(self) -> float:
        """Standard deviation of Ki in µM.

        Raises if model is not "si".
        """
        if self.model_type != "si":
            raise AttributeError(
                f"ki_std is only defined for substrate-inhibition fits; "
                f"this is '{self.model_type}'",
            )
        return self.perr[2]

    def predict(self, s: np.ndarray) -> np.ndarray:
        """Generates model predictions for the given substrate concentrations in µM."""
        return self.model_func(s, *self.popt)


def r_squared(observed: np.ndarray, predicted: np.ndarray) -> float:
    """Computes the coefficient of determination (R²) between observed and predicted values."""
    ss_res = np.sum((observed - predicted) ** 2)
    ss_tot = np.sum((observed - np.mean(observed)) ** 2)
    if ss_tot == 0:
        return 0.0
    return 1.0 - (ss_res / ss_tot)


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
    """Fits a kinetic model function to substrate-velocity data via Levenberg-Marquardt.

    Delegates to scipy.optimize.curve_fit with absolute_sigma=True so that
    supplied per-point standard errors are treated as absolute measurement
    uncertainties rather than relative weights. The full parameter covariance
    matrix is stored in FitResult.pcov for downstream correlated error propagation.

    Args:
        model_func: Callable with signature f(s, *params) -> np.ndarray representing
            the kinetic model to fit.
        s: Substrate concentration array in µM.
        v: Reaction velocity array (area/s before calibration; µM/s after).
        p0: Initial parameter guesses in the same order as model_func expects.
        sigma: Per-point standard errors on v in the same units as v. If None or
            all-zero, unweighted (OLS) fitting is used. Defaults to None.
        bounds: Lower and upper bounds for parameters passed directly to curve_fit.
            Defaults to (0, np.inf).
        maxfev: Maximum number of function evaluations for the optimizer. Defaults
            to 10_000.
        model_type: Short tag stored on the returned FitResult to identify the model
            (e.g. "mm", "hill", "si"). Defaults to "mm".

    Returns:
        FitResult with optimal parameters in popt, standard errors in perr, full
        covariance matrix in pcov, R², and n_points equal to len(s).

    Raises:
        RuntimeError: If curve_fit fails to converge within maxfev evaluations.
        ValueError: If s and v have incompatible shapes or p0 has the wrong length
            for model_func.
    """
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
        r_squared=r_squared(v, predicted),
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
    """Fits the Michaelis-Menten equation to substrate-velocity data.

    The fitted Km is in the same units as s. Throughout this codebase s is
    expected in µM, so the returned Km is in µM. Passing s in any other unit
    (mM, nM) will produce a Km in that unit with no warning; the caller is
    responsible for unit consistency.

    Args:
        s: Substrate concentration array in µM.
        v: Reaction velocity array (area/s before calibration; µM/s after).
        sigma: Per-point standard errors on v for weighted fitting. If None
            or all-zero, unweighted fitting is used. Defaults to None.
        vmax_init: Initial guess for Vmax in the same units as v. Defaults to max(v).
        km_init: Initial guess for Km in µM. Defaults to median(s).

    Returns:
        FitResult with model_type="mm", Km in µM, and full pcov matrix.

    Raises:
        RuntimeError: If curve_fit fails to converge within the default maxfev.
    """
    vmax_guess = vmax_init if vmax_init is not None else float(np.max(v))
    km_guess = km_init if km_init is not None else float(np.median(s))
    return fit_model(
        michaelis_menten,
        s,
        v,
        [vmax_guess, km_guess],
        sigma=sigma,
        model_type="mm",
    )


def fit_threshold_michaelis_menten(
    s: np.ndarray,
    v: np.ndarray,
    *,
    sigma: np.ndarray | None = None,
    vmax_init: float | None = None,
    km_init: float | None = None,
) -> FitResult:
    """Fits the threshold Michaelis-Menten equation to substrate-velocity data.

    The threshold model includes a dead-zone parameter S0: velocity is zero for
    s ≤ S0 and follows standard MM kinetics for the effective concentration (s − S0)
    above the threshold. S0 is initialized to 0.5 µM; vmax_init defaults to
    2 · max(v) to account for the compressed dynamic range introduced by the dead zone.

    The fitted Km reflects the apparent Michaelis constant for the effective substrate
    concentration (s − S0), not the total substrate concentration. Units of Km and S0
    match the units of s; throughout this codebase s is expected in µM.

    Args:
        s: Substrate concentration array in µM.
        v: Reaction velocity array (area/s before calibration; µM/s after).
        sigma: Per-point standard errors on v for weighted fitting. If None
            or all-zero, unweighted fitting is used. Defaults to None.
        vmax_init: Initial guess for Vmax in the same units as v. Defaults to 2 · max(v).
        km_init: Initial guess for Km in µM. Defaults to median(s).

    Returns:
        FitResult with model_type="mm", Km in µM, S0 in popt[2], and full pcov matrix.

    Raises:
        RuntimeError: If curve_fit fails to converge within the default maxfev.
    """
    vmax_guess = vmax_init if vmax_init is not None else float(np.max(v) * 2)
    km_guess = km_init if km_init is not None else float(np.median(s))
    return fit_model(
        threshold_michaelis_menten,
        s,
        v,
        [vmax_guess, km_guess, 0.5],
        sigma=sigma,
        model_type="mm",
    )


def fit_hill(
    s: np.ndarray,
    v: np.ndarray,
    *,
    sigma: np.ndarray | None = None,
) -> FitResult:
    """Fits the Hill equation to substrate-velocity data.

    The second fitted parameter is k_half (the half-saturation constant in µM),
    not Km. k_half equals Km only when the Hill coefficient n == 1. Initial guesses
    are max(v) for Vmax, median(s) for k_half, and 2.0 for n.

    Args:
        s: Substrate concentration array in µM.
        v: Reaction velocity array (area/s before calibration; µM/s after).
        sigma: Per-point standard errors on v for weighted fitting. If None
            or all-zero, unweighted fitting is used. Defaults to None.

    Returns:
        FitResult with model_type="hill", k_half in µM (popt[1]), Hill coefficient n
        in popt[2], and full pcov matrix. Use FitResult.k_half and FitResult.hill_n
        rather than FitResult.km to avoid semantic ambiguity.

    Raises:
        RuntimeError: If curve_fit fails to converge within the default maxfev.
    """
    p0 = [float(np.max(v)), float(np.median(s)), 2.0]
    return fit_model(hill_equation, s, v, p0, sigma=sigma, model_type="hill")


def fit_substrate_inhibition(
    s: np.ndarray,
    v: np.ndarray,
    *,
    sigma: np.ndarray | None = None,
) -> FitResult:
    """Fits the substrate-inhibition model to substrate-velocity data.

    Appropriate when enzyme activity decreases at high substrate concentrations,
    producing a bell-shaped v vs. [S] curve. Ki is initialized to max(s) as a
    heuristic starting point; Vmax to max(v) and Km to median(s).

    Equation: v = (Vmax · [S]) / (Km + [S] + [S]²/Ki)

    Args:
        s: Substrate concentration array in µM.
        v: Reaction velocity array (area/s before calibration; µM/s after).
        sigma: Per-point standard errors on v for weighted fitting. If None
            or all-zero, unweighted fitting is used. Defaults to None.

    Returns:
        FitResult with model_type="si", Km in µM (popt[1]), Ki in µM accessible
        via FitResult.ki (popt[2]), and full pcov matrix.

    Raises:
        RuntimeError: If curve_fit fails to converge within the default maxfev.
    """
    p0 = [float(np.max(v)), float(np.median(s)), float(np.max(s))]
    return fit_model(substrate_inhibition, s, v, p0, sigma=sigma, model_type="si")


def fit_lineweaver_burk(s: np.ndarray, v: np.ndarray) -> FitResult:
    """Fits a Lineweaver-Burk (double-reciprocal) model to substrate-velocity data.

    Transforms (s, v) to reciprocal space via lineweaver_burk_transform, then fits
    a linear regression (OLS) to 1/v vs. 1/[S]. Vmax and Km are back-calculated
    from the regression intercept and slope respectively.

    Standard errors in reciprocal space do not map cleanly to Vmax and Km parameter
    space without the delta method, so perr is set to (nan, nan) and pcov is a 2×2
    nan matrix. These values satisfy the FitResult contract without implying false
    precision. Downstream consumers should treat LB parameter uncertainties as
    unavailable and use the MM fit for error propagation.

    Args:
        s: Substrate concentration array in µM. Non-positive values are excluded
            by lineweaver_burk_transform before fitting.
        v: Reaction velocity array (area/s before calibration; µM/s after).
            Non-positive values are excluded before fitting.

    Returns:
        FitResult with model_type="lb", Vmax in the same units as v (popt[0]),
        Km in µM (popt[1]), perr of (nan, nan), pcov of shape (2, 2) filled with
        nan, and extra dict containing the raw regression slope, intercept, and
        std_err from scipy.stats.linregress.
    """
    s_inv, v_inv = lineweaver_burk_transform(s, v)
    slope, intercept, r_value, _, std_err = linregress(s_inv, v_inv)
    r2 = r_value ** 2
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
