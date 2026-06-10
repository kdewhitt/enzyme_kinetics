"""Velocity preparation and per-peak data extraction utilities for kinetic analysis."""

from __future__ import annotations

from typing import Final

import numpy as np
import pandas as pd

__all__ = [
    "COUNT",
    "extract_peak_data",
    "MEAN",
    "PEAK_ID",
    "prepare_velocity",
    "STD",
    "SUBSTRATE_CONC",
]

SUBSTRATE_CONC = "substrate_conc"
"""DataFrame column name for substrate concentration values."""

PEAK_ID: Final = "peak_id"
"""DataFrame column name for HPLC peak identifiers."""

MEAN: Final = "mean"
"""DataFrame column name for mean peak area per replicate group."""

STD: Final = "std"
"""DataFrame column name for standard deviation of peak area per replicate group."""

COUNT: Final = "count"
"""DataFrame column name for replicate count per group."""


def prepare_velocity(
    s: np.ndarray,
    mean_signal: np.ndarray,
    std_signal: np.ndarray,
    counts: np.ndarray,
    rxn_time: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Computes reaction velocities and their SEMs from replicate summary statistics.

    Returns raw mean_signal as the fourth element so that callers can apply a
    calibration curve directly to area before converting to velocity, avoiding
    the fragile round-trip of reconstructing area from velocity.

    Args:
        s: Substrate concentration array in µM.
        mean_signal: Mean HPLC peak area per replicate group in area units.
        std_signal: Standard deviation of peak area per replicate group in area units.
        counts: Number of replicates per group.
        rxn_time: Reaction duration in seconds.

    Returns:
        A tuple of (s, v, v_std, mean_signal) where s is the substrate
        concentration array in µM (unchanged), v is reaction velocity in area/s,
        v_std is the standard deviation of v in area/s, and mean_signal is the raw
        peak area array in area units for direct calibration application.
    """
    v = mean_signal / rxn_time
    v_std = std_signal / rxn_time
    return s, v, v_std, mean_signal


def extract_peak_data(
    df: pd.DataFrame,
    peak_id: str,
    rxn_time: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Extracts substrate concentrations, velocity, SEM, and raw area for a single peak.

    Filters df to rows matching peak_id, drops rows with missing substrate
    concentration or mean signal, and returns None if the filtered subset is
    empty or all-zero. Raw mean_signal is returned as the fourth element to
    enable apply_calibration to work directly on area rather than reconstructing
    it from velocity.

    Args:
        df: DataFrame containing summary statistics for all peaks, with columns
            matching PEAK_ID, SUBSTRATE_CONC, MEAN, STD, and COUNT.
        peak_id: HPLC peak identifier to filter on.
        rxn_time: Reaction duration in seconds used to convert area to velocity.

    Returns:
        A tuple of (s, v, v_std, mean_signal) as returned by prepare_velocity,
        or None if no valid data exists for the given peak_id.
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
