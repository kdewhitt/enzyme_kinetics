from __future__ import annotations

from typing import Final

import numpy as np
import pandas as pd

SUBSTRATE_CONC = "substrate_conc"
PEAK_ID: Final = "peak_id"
MEAN: Final = "mean"
STD: Final = "std"
COUNT: Final = "count"


# ---------------------------------------------------------------------------
# Velocity preparation from summary-statistics DataFrames
# FIX 3: return raw mean_signal alongside velocity so apply_calibration
#        can apply the calibration directly to area before dividing by rxn_time.
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Per-peak data extraction
# FIX 3: unpack four values from prepare_velocity (now returns mean_signal too)
# ---------------------------------------------------------------------------

def extract_peak_data(
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
        sub[MEAN].values,
        sub[STD].values,
        sub[COUNT].values,
        rxn_time,
    )
