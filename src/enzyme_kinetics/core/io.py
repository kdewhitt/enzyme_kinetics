"""Data loading and pre-processing utilities for plot modules."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd
from kgdlibs.datatools import clean_names

from enzyme_kinetics.compounds import canonicalize_peak_ids, categorize_proteins

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


def categorize_columns(
    df: pd.DataFrame,
    prefixes: str | Iterable[str] | None,
    ref_protein: str | None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Canonicalizes peak ID prefixes and categorizes protein labels in a DataFrame.

    Delegates peak ID normalization to canonicalize_peak_ids, then applies
    protein categorization via categorize_proteins using ref_protein as the
    reference baseline.

    Args:
        df: The DataFrame containing at least "protein" and "peak_id" columns.
        prefixes: One or more peak ID prefix strings used by
            canonicalize_peak_ids to normalize peak identifiers. When None,
            no prefix filtering is applied.
        ref_protein: The reference protein label passed to categorize_proteins
            as the baseline category. When None, no reference is set.
        **kwargs: Additional keyword arguments forwarded to categorize_proteins.

    Returns:
        A new DataFrame with canonicalized peak IDs and categorized protein
        labels.
    """
    df = canonicalize_peak_ids(df, prefixes=prefixes)
    return categorize_proteins(df, reference=ref_protein, **kwargs)


def wide_to_long(
    df: pd.DataFrame,
    *,
    id_vars: str | Sequence[str] | None = None,
    var_name: str = "variable",
    value_name: str = "value",
) -> pd.DataFrame:
    """Converts a wide-format DataFrame to long format via reset_index and melt.

    Args:
        df: A wide-format DataFrame, typically indexed by retention time, with
            one column per sample or variable to unpivot.
        id_vars: Column name or sequence of column names to retain as
            identifier variables after reset_index. When None, all columns are
            treated as value variables. Defaults to None.
        var_name: Name assigned to the new column containing the original
            column headers. Defaults to "variable".
        value_name: Name assigned to the new column containing the melted
            values. Defaults to "value".

    Returns:
        A long-format DataFrame with one row per (id_var, variable) combination,
        with the original index moved into a regular column by reset_index.
    """
    return df.reset_index().melt(
        id_vars=id_vars,
        var_name=var_name,
        value_name=value_name,
    )
