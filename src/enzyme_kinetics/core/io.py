"""Data loading utility."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd
from kgdlibs.datatools import clean_names

__all__ = ["drop_indexlike_columns", "load_plottable_data"]

# Matches bare "index", "idx", or pandas auto-generated "Unnamed: N" column names.
_INDEXLIKE_RE = re.compile(r"^(?:index|idx|unnamed:\s*\d+(?:\.\d+)?)$", re.IGNORECASE)


def drop_indexlike_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Removes columns whose names match accidental index artifact patterns.

    Args:
        df: The DataFrame from which index-like columns are to be removed.
            Column names are matched case-insensitively against _INDEXLIKE_RE.

    Returns:
        A view of the DataFrame with all index-like columns excluded.
    """
    return df.loc[:, ~df.columns.str.contains(_INDEXLIKE_RE)]


def load_plottable_data(
    path: Path,
    *,
    sanitize: bool = False,
    drop_indexlike: bool = True,
    **kwargs: Any,
) -> pd.DataFrame:
    """Loads a CSV file into a DataFrame and applies optional pre-processing.

    Args:
        path: Path to the CSV file to load.
        sanitize: If True, applies clean_names to normalize column headers to
            snake_case after loading. Defaults to False.
        drop_indexlike: If True, removes columns whose names match index
            artifact patterns via drop_indexlike_columns. Defaults to True.
        **kwargs: Additional keyword arguments forwarded to pd.read_csv.

    Returns:
        A non-empty DataFrame ready for plotting, with index-like columns and
        optional column name normalization applied.

    Raises:
        ValueError: If the loaded DataFrame is empty.
    """
    df = pd.read_csv(path, **kwargs)

    if df.empty:
        raise ValueError(f"No data to plot: {path!r}")

    if drop_indexlike:
        df = drop_indexlike_columns(df)

    if sanitize:
        df = clean_names(df)

    return df
