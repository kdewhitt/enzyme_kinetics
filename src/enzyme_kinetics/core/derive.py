"""Utilities for computing kcat and kcat/Km from kinetic fit results.

Provides derive_constants, the primary entry point for computing turnover number
(kcat) and catalytic efficiency (kcat/Km) from a FitResult produced by any
MM-family fitting function. Uncertainty propagation uses the full Vmax–Km
covariance matrix via the delta method, correctly accounting for the
anti-correlation between Vmax and Km that arises from non-linear fitting.
Results are returned as a frozen KineticConstants dataclass.

kcat and kcat/Km carry physically meaningful units (s⁻¹ and s⁻¹·M⁻¹
respectively) only after a calibration curve has been applied upstream so that
Vmax is in µM/s. Before calibration both quantities are in non-physical units
and should not be compared to literature values.

Typical usage example:
    >>> from enzyme_kinetics.core.models import fit_michaelis_menten
    >>> fit = fit_michaelis_menten(s, v)
    >>> constants = derive_constants("HexCoA", "peak_1", fit, enzyme_conc_um=0.5)
    >>> print(constants.kcat, constants.kcat_km_M)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from logurich import DuoLogAdapter

from .models import FitResult

__all__ = ["derive_constants", "KineticConstants"]

_logger = DuoLogAdapter.create(component=__name__)

# ---------------------------------------------------------------------
# kcat/Km error propagation using full covariance matrix
# ---------------------------------------------------------------------

# Unit conversion factor: µM → M.
_UM_TO_M: float = 1e-6


def _kcat_km_with_covariance(
    vmax: float,
    km: float,
    enzyme_conc_um: float,
    pcov: np.ndarray,
) -> tuple[float, float]:
    """Propagates kcat/Km uncertainty using the full Vmax–Km covariance matrix.

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
    in the standard literature unit of s⁻¹·M⁻¹, so Km is converted to M
    (×1e-6) before computing the ratio. The partial derivatives are scaled by
    the same factor so that the returned SE is also in s⁻¹·M⁻¹.

    Args:
        vmax: Fitted Vmax in µM/s (after calibration) or area/s (before).
        km: Fitted Km in µM (or k_half in µM for Hill fits).
        enzyme_conc_um: Total enzyme concentration in µM.
        pcov: Full 2×2 (or larger) covariance matrix from curve_fit, with
            Vmax at index [0,0] and Km at [1,1], covariance at [0,1].
            Units of pcov entries must be consistent with vmax and km units.

    Returns:
        A tuple of (kcat_km_M, kcat_km_std_M) in s⁻¹·M⁻¹. Returns (nan, nan)
        if km ≤ 0, vmax ≤ 0, or pcov contains nan (e.g. for LB fits).
    """
    if km <= 0 or vmax <= 0:
        return np.nan, np.nan
    if not np.all(np.isfinite(pcov[:2, :2])):
        # Fallback for LB fits or singular covariance (pcov contains nan or inf)
        return np.nan, np.nan

    # Convert Km from µM to M for the standard s⁻¹·M⁻¹ unit of kcat/Km
    km_M = km * _UM_TO_M
    kcat_km_M = vmax / (enzyme_conc_um * km_M)

    # Partial derivatives of f = vmax / (enzyme_conc_um * km_M)
    # where km_M = km * 1e-6, so ∂f/∂km = ∂f/∂km_M · ∂km_M/∂km = (∂f/∂km_M) * 1e-6
    df_dvmax = 1.0 / (enzyme_conc_um * km_M)
    df_dkm = -vmax / (enzyme_conc_um * km_M ** 2) * _UM_TO_M

    var_vmax = pcov[0, 0]
    var_km = pcov[1, 1]
    cov_vmax_km = pcov[0, 1]

    var_kcat_km = (
            df_dvmax ** 2 * var_vmax
            + df_dkm ** 2 * var_km
            + 2.0 * df_dvmax * df_dkm * cov_vmax_km
    )
    kcat_km_std_M = float(np.sqrt(max(var_kcat_km, 0.0)))  # clamp numerical negatives
    return float(kcat_km_M), kcat_km_std_M


# ---------------------------------------------------------------------
# Derived kinetic constants
# ---------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class KineticConstants:
    """Fully-typed kinetic constants for one peak.

    Attributes:
        substrate: Substrate name string identifying the acyl-CoA donor used
            in the assay.
        peak_id: HPLC peak identifier string.
        fit: Primary FitResult from the MM, Hill, or SI model fit.
        kcat: Turnover number in s⁻¹ after calibration; in non-physical units
            (area · µM⁻¹ · s⁻¹) before calibration.
        kcat_std: Standard deviation of kcat in the same units as kcat.
        kcat_km_M: Catalytic efficiency kcat/Km in s⁻¹·M⁻¹ after calibration.
            Km is converted from µM to M before computing the ratio so that the
            result is in the standard literature unit. Before calibration the
            numerator carries non-physical units and kcat_km_M should not be
            compared to literature values.
        kcat_km_std_M: Standard deviation of kcat_km_M in s⁻¹·M⁻¹, propagated via
            the full Vmax–Km covariance matrix.
        ki: Substrate inhibition constant Ki in µM; None unless model_type == "si".
        ki_std: Standard deviation of Ki in µM; None unless model_type == "si".
        lb_fit: Optional Lineweaver-Burk FitResult stored for cross-validation.
    """

    substrate: str
    peak_id: str
    fit: FitResult
    kcat: float
    kcat_std: float
    kcat_km_M: float  # s⁻¹·M⁻¹; Km converted µM→M before division
    kcat_km_std_M: float  # std in s⁻¹·M⁻¹
    ki: float | None = None
    ki_std: float | None = None
    lb_fit: FitResult | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serializes the KineticConstants instance to a flat dictionary.

        Constructs a dictionary suitable for DataFrame or CSV export. The key
        for the second kinetic parameter is model-dependent: "k_half" for Hill
        fits and "km" for all other models. hill_n and hill_n_std are emitted
        only for Hill fits; ki and ki_std are emitted only for SI fits; lb_fit
        fields (km_lb, vmax_lb, r_squared_lb) are emitted only when lb_fit is
        not None. Any extra scalar values on fit.extra are forwarded unless
        their key already exists in the base dictionary.

        Returns:
            A flat dictionary mapping column name strings to scalar values
            (float, int, or str). Keys vary by model type as described above.
        """
        # FIX 7: use model-correct labels for Hill output
        is_hill = self.fit.model_type == "hill"
        base: dict[str, Any] = {
            "substrate": self.substrate,
            "peak_id": self.peak_id,
            "model": self.fit.model_type,
            "vmax": self.fit.vmax,
            "vmax_std": self.fit.vmax_std,
            # Label the second parameter correctly per model
            ("k_half" if is_hill else "km"): self.fit.km,
            ("k_half_std" if is_hill else "km_std"): self.fit.km_std,
            "r_squared": self.fit.r_squared,
            "kcat": self.kcat,
            "kcat_std": self.kcat_std,
            # Explicit s⁻¹·M⁻¹ suffix in the column name makes the unit
            # unambiguous in downstream CSV/DataFrame consumers
            "kcat_km_M": self.kcat_km_M,
            "kcat_km_std_M": self.kcat_km_std_M,
        }
        # FIX 7: emit hill_n only for Hill fits
        if is_hill:
            base["hill_n"] = self.fit.hill_n
            base["hill_n_std"] = self.fit.hill_n_std
        # FIX 4: emit Ki only for SI fits
        if self.ki is not None:
            base["ki"] = self.ki
            base["ki_std"] = self.ki_std
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
    substrate: str,
    peak_id: str,
    fit: FitResult,
    enzyme_conc_um: float,
    *,
    lb_fit: FitResult | None = None,
) -> KineticConstants:
    """Computes kcat and kcat/Km from a FitResult.

    kcat is computed as Vmax / enzyme_conc_um. Its unit is s⁻¹ only after a
    calibration curve has been applied (so that Vmax is in µM/s); before
    calibration Vmax is in area/s and kcat carries non-physical units.

    kcat/Km is reported in s⁻¹·M⁻¹ (the standard literature unit) by
    converting Km from µM to M before computing the ratio. The conversion
    factor _UM_TO_M (1e-6) is applied inside _kcat_km_with_covariance, which
    also propagates uncertainty via the full Vmax–Km covariance matrix,
    correctly accounting for the anti-correlation between Vmax and Km that
    arises from MM-family non-linear fitting.

    Substrate concentrations (s) passed to the fitting functions must be in
    µM for Km and kcat/Km to carry their stated units. No runtime check is
    performed; the caller is responsible for unit consistency.

    Ki and ki_std are extracted from the fit and promoted to typed fields when
    the model is a substrate-inhibition fit (model_type == "si").

    Args:
        substrate: Substrate name string identifying the acyl-CoA donor used
            in the assay (e.g. "HexCoA").
        peak_id: Identifier for the HPLC peak being processed.
        fit: FitResult from any kinetic model fit. Km (popt[1]) must be in µM.
        enzyme_conc_um: Total enzyme concentration in µM.
        lb_fit: Optional Lineweaver-Burk cross-validation FitResult. Defaults
            to None.

    Returns:
        KineticConstants with substrate, peak_id, kcat in s⁻¹ (post-calibration),
        kcat_km_M in s⁻¹·M⁻¹, their stds, and optionally ki/ki_std in µM for
        SI fits.
    """
    # kcat = Vmax / [E]
    kcat = fit.vmax / enzyme_conc_um
    kcat_std = fit.vmax_std / enzyme_conc_um

    # FIX 1 + units: full covariance propagation; Km converted µM→M inside
    kcat_km_M, kcat_km_std_M = _kcat_km_with_covariance(
        fit.vmax,
        fit.km,
        enzyme_conc_um,
        fit.pcov,
    )

    # FIX 4: promote Ki to typed fields for SI fits
    ki: float | None = None
    ki_std: float | None = None
    if fit.model_type == "si":
        ki = fit.ki
        ki_std = fit.ki_std

    return KineticConstants(
        substrate=substrate,
        peak_id=peak_id,
        fit=fit,
        kcat=kcat,
        kcat_std=kcat_std,
        kcat_km_M=kcat_km_M,
        kcat_km_std_M=kcat_km_std_M,
        ki=ki,
        ki_std=ki_std,
        lb_fit=lb_fit,
    )
