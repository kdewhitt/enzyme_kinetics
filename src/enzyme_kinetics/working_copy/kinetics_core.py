from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
from scipy.optimize import curve_fit

_logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pure model functions — no classes, no state
# ---------------------------------------------------------------------------

def michaelis_menten(s: np.ndarray, vmax: float, km: float) -> np.ndarray:
    return (vmax * s) / (km + s)


def hill_equation(s: np.ndarray, vmax: float, k_half: float, n: float) -> np.ndarray:
    return (vmax * s ** n) / (k_half ** n + s ** n)


def substrate_inhibition(s: np.ndarray, vmax: float, km: float, ki: float) -> np.ndarray:
    return (vmax * s) / (km + s + (s ** 2) / ki)


def lineweaver_burk_transform(
        s: np.ndarray,
        v: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mask = (s > 0) & (v > 0)
    return 1.0 / s[mask], 1.0 / v[mask]


# ---------------------------------------------------------------------------
# Fit result — frozen, model-agnostic
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class FitResult:
    model_func: Callable[..., np.ndarray]
    popt: tuple[float, ...]
    perr: tuple[float, ...]
    r_squared: float
    n_points: int
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def vmax(self) -> float:
        return self.popt[0]

    @property
    def vmax_se(self) -> float:
        return self.perr[0]

    @property
    def km(self) -> float:
        return self.popt[1]

    @property
    def km_se(self) -> float:
        return self.perr[1]

    def predict(self, s: np.ndarray) -> np.ndarray:
        return self.model_func(s, *self.popt)


def _r_squared(observed: np.ndarray, predicted: np.ndarray) -> float:
    ss_res = np.sum((observed - predicted) ** 2)
    ss_tot = np.sum((observed - np.mean(observed)) ** 2)
    if ss_tot == 0:
        return 0.0
    return 1.0 - ss_res / ss_tot


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
        r_squared=_r_squared(v, predicted),
        n_points=len(s),
    )


def fit_michaelis_menten(
        s: np.ndarray,
        v: np.ndarray,
        *,
        sigma: np.ndarray | None = None,
        vmax_init: float | None = None,
        km_init: float | None = None,
) -> FitResult:
    vmax_guess = vmax_init if vmax_init is not None else float(np.max(v))
    km_guess = km_init if km_init is not None else float(np.median(s))
    return fit_model(michaelis_menten, s, v, [vmax_guess, km_guess], sigma=sigma)


def fit_hill(
        s: np.ndarray,
        v: np.ndarray,
        *,
        sigma: np.ndarray | None = None,
) -> FitResult:
    p0 = [float(np.max(v)), float(np.median(s)), 2.0]
    return fit_model(hill_equation, s, v, p0, sigma=sigma)


def fit_substrate_inhibition(
        s: np.ndarray,
        v: np.ndarray,
        *,
        sigma: np.ndarray | None = None,
) -> FitResult:
    p0 = [float(np.max(v)), float(np.median(s)), float(np.max(s))]
    return fit_model(substrate_inhibition, s, v, p0, sigma=sigma)


def fit_lineweaver_burk(
        s: np.ndarray,
        v: np.ndarray,
) -> FitResult:
    from scipy.stats import linregress

    s_inv, v_inv = lineweaver_burk_transform(s, v)
    slope, intercept, r_value, _, std_err = linregress(s_inv, v_inv)
    r2 = r_value ** 2
    vmax = 1.0 / intercept if intercept != 0 else np.nan
    km = slope * vmax if not np.isnan(vmax) else np.nan
    return FitResult(
        model_func=michaelis_menten,
        popt=(vmax, km),
        perr=(np.nan, np.nan),
        r_squared=r2,
        n_points=len(s_inv),
        extra={"slope": slope, "intercept": intercept, "std_err": std_err},
    )


# ---------------------------------------------------------------------------
# Derived kinetic constants
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class KineticConstants:
    peak_id: str
    fit: FitResult
    kcat: float
    kcat_se: float
    kcat_km: float
    kcat_km_se: float
    lb_fit: FitResult | None = None

    def to_dict(self) -> dict[str, Any]:
        base: dict[str, Any] = {
            "peak_id": self.peak_id,
            "vmax": self.fit.vmax,
            "vmax_se": self.fit.vmax_se,
            "km": self.fit.km,
            "km_se": self.fit.km_se,
            "r_squared": self.fit.r_squared,
            "kcat": self.kcat,
            "kcat_se": self.kcat_se,
            "kcat_km": self.kcat_km,
            "kcat_km_se": self.kcat_km_se,
        }
        if self.lb_fit is not None:
            base["km_lb"] = self.lb_fit.km
            base["vmax_lb"] = self.lb_fit.vmax
            base["r_squared_lb"] = self.lb_fit.r_squared
        for k, v in self.fit.extra.items():
            if isinstance(v, (int, float, str)):
                base[k] = v
        return base


def derive_constants(
        peak_id: str,
        fit: FitResult,
        enzyme_conc_um: float,
        *,
        lb_fit: FitResult | None = None,
) -> KineticConstants:
    kcat = fit.vmax / enzyme_conc_um
    kcat_se = fit.vmax_se / enzyme_conc_um
    if fit.km > 0 and kcat > 0:
        kcat_km = kcat / fit.km
        rel_kcat = kcat_se / kcat
        rel_km = fit.km_se / fit.km
        kcat_km_se = kcat_km * np.sqrt(rel_kcat ** 2 + rel_km ** 2)
    else:
        kcat_km = np.nan
        kcat_km_se = np.nan
    return KineticConstants(
        peak_id=peak_id,
        fit=fit,
        kcat=kcat,
        kcat_se=kcat_se,
        kcat_km=kcat_km,
        kcat_km_se=kcat_km_se,
        lb_fit=lb_fit,
    )


# ---------------------------------------------------------------------------
# Calibration — linear standard curve
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class Calibration:
    slope: float
    slope_se: float
    intercept: float
    intercept_se: float
    r_squared: float
    conc_range: tuple[float, float]

    def area_to_conc(self, area: np.ndarray | float) -> np.ndarray | float:
        return (area - self.intercept) / self.slope

    def conc_to_area(self, conc: np.ndarray | float) -> np.ndarray | float:
        return self.slope * conc + self.intercept

    def area_to_conc_with_error(
            self, area: float, area_se: float,
    ) -> tuple[float, float]:
        conc = float(self.area_to_conc(area))
        dc_da = 1.0 / self.slope
        dc_ds = -(area - self.intercept) / (self.slope ** 2)
        conc_se = float(np.sqrt((dc_da * area_se) ** 2 + (dc_ds * self.slope_se) ** 2))
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
        r_squared=_r_squared(peak_areas, predicted),
        conc_range=(float(concentrations.min()), float(concentrations.max())),
    )


# ---------------------------------------------------------------------------
# Velocity preparation from summary-statistics DataFrames
# ---------------------------------------------------------------------------

def prepare_velocity(
        s: np.ndarray,
        mean_signal: np.ndarray,
        std_signal: np.ndarray,
        counts: np.ndarray,
        rxn_time: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    v = mean_signal / rxn_time
    v_sem = (std_signal / rxn_time) / np.sqrt(counts)
    return s, v, v_sem


# ---------------------------------------------------------------------------
# Batch analysis helper
# ---------------------------------------------------------------------------

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
