from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
from logurich import RichLogAdapter
from scipy.optimize import curve_fit
from scipy.stats import linregress

_logger = RichLogAdapter(component=__name__)


# ---------------------------------------------------------------------------
# Pure Mathematical Functions — no classes, no state
# ---------------------------------------------------------------------------

def michaelis_menten(s: np.ndarray, vmax: float, km: float) -> np.ndarray:
    """
    Michaelis-Menten enzyme kinetics model.

    Equation: v = (Vmax * [S]) / (Km + [S])

    substrate_conc : np.ndarray
            Substrate concentrations (µM)
        Vmax : float
            Maximum reaction velocity
        Km : float
            Michaelis constant

        Returns
        -------
        velocity : np.ndarray
            Reaction velocities
    """
    return (vmax * s) / (km + s)


def hill_equation(s: np.ndarray, vmax: float, k_half: float, n: float) -> np.ndarray:
    """Hill equation. k_half is the half-saturation constant, equal to Km only when n == 1."""
    return (vmax * s ** n) / (k_half ** n + s ** n)


def threshold_michaelis_menten(s: np.ndarray, vmax: float, km: float, s0: float) -> np.ndarray:
    return np.where(s > s0, (vmax * (s - s0)) / (km + (s - s0)), 0.0)


def substrate_inhibition(s: np.ndarray, vmax: float, km: float, ki: float) -> np.ndarray:
    """
    Substrate inhibition equation.

    Parameters
    ----------
    substrate_conc : np.ndarray
        Substrate concentrations (µM)
    Vmax : float
        Maximum reaction velocity
    Km : float
        Michaelis constant
    Ki : float
        Inhibition constant

    Equation: v = (Vmax * [S]) / (Km + [S] + [S]²/Ki)

    Useful when enzyme activity decreases at high substrate concentrations.

    Returns
    -------
    velocity : np.ndarray
        Reaction velocities
    """
    return (vmax * s) / (km + s + (s ** 2) / ki)


def lineweaver_burk_transform(
    s: np.ndarray,
    v: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Lineweaver-Burk equation.

    Parameters
    ----------
    substrate_conc_inv : np.ndarray
        Inverse substrate concentrations (1/[S])
    Vmax : float
        Maximum reaction velocity
    Km : float
        Michaelis constant

    Returns
    -------
    velocity_inv : np.ndarray
        Inverse velocities (1/v)

    (Km / Vmax) * substrate_conc_inv + 1 / Vmax
    """
    mask = (s > 0) & (v > 0)
    return 1.0 / s[mask], 1.0 / v[mask]


# ---------------------------------------------------------------------------
# Fit Result — frozen, model-agnostic
# ---------------------------------------------------------------------------

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
            raise AttributeError(f"k_half is only defined for Hill fits; this is '{self.model_type}'")
        return self.popt[1]

    @property
    def k_half_se(self) -> float:
        """SE of k_half from a Hill fit. Raises if model is not 'hill'."""
        if self.model_type != "hill":
            raise AttributeError(f"k_half_se is only defined for Hill fits; this is '{self.model_type}'")
        return self.perr[1]

    @property
    def hill_n(self) -> float:
        """Hill cooperativity coefficient (popt[2]). Raises if model is not 'hill'."""
        if self.model_type != "hill":
            raise AttributeError(f"hill_n is only defined for Hill fits; this is '{self.model_type}'")
        return self.popt[2]

    @property
    def hill_n_se(self) -> float:
        """SE of Hill n (perr[2]). Raises if model is not 'hill'."""
        if self.model_type != "hill":
            raise AttributeError(f"hill_n_se is only defined for Hill fits; this is '{self.model_type}'")
        return self.perr[2]

    # FIX 4 — typed Ki accessor for substrate-inhibition fits.
    @property
    def ki(self) -> float:
        """Substrate inhibition constant Ki (popt[2]). Raises if model is not 'si'."""
        if self.model_type != "si":
            raise AttributeError(f"ki is only defined for substrate-inhibition fits; this is '{self.model_type}'")
        return self.popt[2]

    @property
    def ki_se(self) -> float:
        """SE of Ki (perr[2]). Raises if model is not 'si'."""
        if self.model_type != "si":
            raise AttributeError(f"ki_se is only defined for substrate-inhibition fits; this is '{self.model_type}'")
        return self.perr[2]

    def predict(self, s: np.ndarray) -> np.ndarray:
        return self.model_func(s, *self.popt)


def r_squared(observed: np.ndarray, predicted: np.ndarray) -> float:
    ss_res = np.sum((observed - predicted) ** 2)
    ss_tot = np.sum((observed - np.mean(observed)) ** 2)
    if ss_tot == 0:
        return 0.0
    return 1.0 - (ss_res / ss_tot)


# ---------------------------------------------------------------------------
# Fitting — thin wrappers around curve_fit
# ---------------------------------------------------------------------------

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
        model_func, s, v,
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
    return fit_model(michaelis_menten, s, v, [vmax_guess, km_guess], sigma=sigma, model_type="mm")

def fit_threshold_michaelis_menten(
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
    vmax_guess = vmax_init if vmax_init is not None else float(np.max(v)*2)
    km_guess = km_init if km_init is not None else float(np.median(s))
    return fit_model(threshold_michaelis_menten, s, v, [vmax_guess, km_guess, 0.5], sigma=sigma, model_type="mm")

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
    """Fit the substrate-inhibition model. Ki is accessible via FitResult.ki.

    p0 = [Vmax_init, Km_init, Ki_init]

    Equation: v = (Vmax * [S]) / (Km + [S] + [S]²/Ki)

    Useful when enzyme activity decreases at high substrate concentrations.
    """
    p0 = [float(np.max(v)), float(np.median(s)), float(np.max(s))]
    return fit_model(substrate_inhibition, s, v, p0, sigma=sigma, model_type="si")


def fit_lineweaver_burk(s: np.ndarray, v: np.ndarray) -> FitResult:
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
