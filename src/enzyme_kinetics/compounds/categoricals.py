"""Categorical normalization utilities for peak identifiers and protein ordering."""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from typing import Any

import pandas as pd

from .compounds import CATEGORICAL_ORDER, Compounds

__all__ = ["canonicalize_peak_ids", "categorize_proteins", "CUSTOM_PROTEIN_ORDER"]

CUSTOM_PROTEIN_ORDER = [
    "TKS", "TKS_OAC", "TKS (Red.)", "Red.", "TKS (Ox.)", "Ox.",
    "K51R", "K51R/", "K51R/K301R", "A125S", "A125S/", "A125S/D185E",
    "A125T", "S126A", "S126C", "M130A", "M130K", "M130R", "D185E",
    "M187A", "M187A_OAC", "C189A", "C189S", "L190F", "L190G", "L190T",
    "L257A", "L257G", "L261A", "L261N", "K263A", "K263R", "K301R",
    "S332A", "S332C", "S332T",
]
"""Default ordered list of protein variant labels for categorical axis sorting."""


def _ensure_column_exists(df: pd.DataFrame, column: str) -> None:
    """Raises KeyError if column is not present in the DataFrame."""
    if column not in df.columns:
        raise KeyError(f"Column {column!r} not found in DataFrame")


def _coerce_to_list(value: str | Iterable[str] | None) -> list[str]:
    """Coerces None, a single string, or an iterable of strings to a list."""
    if value is None:
        return []

    if isinstance(value, str):
        return [value]

    return list(value)


def _canonicalize_prefixes(prefixes: str | Iterable[str] | None) -> list[str]:
    """Normalizes prefix strings to a canonical trailing-hyphen format.

    Strips leading and trailing whitespace and hyphens from each prefix, then
    appends a single trailing hyphen. This ensures that prefixes like "C5",
    "C5-", or " C5- " all produce the canonical form "C5-".

    Args:
        prefixes: A single prefix string, an iterable of prefix strings, or
            None. Passed through _coerce_to_list before normalization.

    Returns:
        A list of normalized prefix strings, each ending with a single hyphen.
        Returns an empty list if prefixes is None or empty.
    """
    return [f"{p.strip().strip('-')}-" for p in _coerce_to_list(prefixes)]


def _compile_prefix_pattern(prefixes: Sequence[str]) -> re.Pattern[str] | None:
    """Compiles a regex pattern matching any of the supplied prefixes at the start of a string.

    Args:
        prefixes: A sequence of normalized prefix strings (as produced by
            _canonicalize_prefixes) to include in the alternation.

    Returns:
        A compiled regex Pattern matching any of the prefixes at position zero,
        or None if prefixes is empty.
    """
    if not prefixes:
        return None

    escaped = "|".join(re.escape(p) for p in prefixes)
    return re.compile(rf"^({escaped})")


def _ordered_unique(values: Iterable[Any]) -> list[Any]:
    """Returns a list of unique values from an iterable, preserving first-seen order."""
    return list(dict.fromkeys(values))


def _build_peak_category_order(
    *,
    base_order: Sequence[str],
    prefixes: Sequence[str],
    observed: Iterable[Any],
) -> list[str]:
    """Builds an ordered list of peak category labels for use in a pandas Categorical.

    Interleaves prefixed compound names into the base ordering so that each
    base compound is immediately followed by all of its prefixed variants. Any
    observed values not covered by the base order or prefix combinations are
    appended in sorted order at the end.

    Args:
        base_order: The canonical sequence of base compound label strings,
            typically CATEGORICAL_ORDER.
        prefixes: The sequence of normalized prefix strings (trailing-hyphen
            form) to interleave after each base compound.
        observed: The full set of normalized peak ID values observed in the
            DataFrame column, used to collect any extra labels not covered by
            base_order and prefixes.

    Returns:
        An ordered list of category label strings with base compounds and
        their prefixed variants interleaved, followed by any additional
        observed values in sorted order.
    """
    ordered: list[str] = []

    for name in base_order:
        ordered.append(name)
        ordered.extend(f"{prefix}{name}" for prefix in prefixes)

    priority = set(ordered)
    extras = sorted(str(value) for value in observed if value not in priority)

    return ordered + extras


def canonicalize_peak_ids(
    df: pd.DataFrame,
    *,
    column: str = "peak_id",
    prefixes: str | Iterable[str] | None = None,
) -> pd.DataFrame:
    """Normalizes a peak ID column to canonical Compounds values and converts it to a Categorical.

    Strips and uppercases raw peak ID strings, optionally extracts and
    reattaches chain-length prefixes (e.g. "C5-"), resolves the base
    identifier to a canonical Compounds value via Compounds.resolve, and
    assigns an ordered pandas Categorical using the interleaved compound
    ordering produced by _build_peak_category_order.

    Args:
        df: The peak table DataFrame containing the column to normalize.
        column: The name of the column holding raw peak identifier strings.
            Defaults to "peak_id".
        prefixes: Optional prefix or prefixes to extract before compound
            resolution and reattach after normalization. Accepts a single
            string (e.g. "C5"), an iterable of strings, or None to skip
            prefix handling entirely. Defaults to None.

    Returns:
        A copy of df with the peak ID column replaced by an ordered
        Categorical whose categories follow the interleaved compound ordering
        defined by CATEGORICAL_ORDER and the supplied prefixes.

    Raises:
        KeyError: If column is not present in df.

    Examples:
        >>> canonicalize_peak_ids(df,column="peak_id")
        >>> canonicalize_peak_ids(df,column="peak_id",prefixes=["C5", "C6"])
    """
    _ensure_column_exists(df, column)

    df = df.copy()
    prefix_list = _canonicalize_prefixes(prefixes)
    prefix_pattern = _compile_prefix_pattern(prefix_list)

    raw = df[column].astype(str).str.strip()

    if prefix_pattern is None:
        extracted = ""
        peaks = raw.str.upper()
    else:
        pattern = prefix_pattern.pattern
        extracted = raw.str.extract(pattern, expand=False).fillna("")
        peaks = raw.str.replace(pattern, "", regex=True).str.upper()

    unique_peaks = pd.unique(peaks.dropna())
    label_map = {peak: str(Compounds.resolve(peak)) for peak in unique_peaks}

    normalized = peaks.map(label_map)

    if prefix_pattern is not None:
        normalized = extracted + normalized

    categories = _build_peak_category_order(
        base_order=CATEGORICAL_ORDER,
        prefixes=prefix_list,
        observed=pd.unique(normalized.dropna()),
    )

    df[column] = pd.Categorical(normalized, categories=categories, ordered=True)
    return df


def categorize_proteins(
    df: pd.DataFrame,
    *,
    column: str = "protein",
    reference: str | None = "BSA",
    category_order: Sequence[str] | None = None,
    include_missing_categories: bool = False,
) -> pd.DataFrame:
    """Normalizes a protein column to an ordered Categorical with a pinned reference entry.

    Derives the category list from category_order if supplied, appending any
    observed values not present in the specified order. If reference is
    provided, that value is moved to the front of the category list. By
    default, only categories observed in the column are retained; setting
    include_missing_categories to True preserves all specified categories even
    if absent from the data.

    Args:
        df: The DataFrame containing the column to normalize.
        column: The name of the column holding protein identifier strings.
            Defaults to "protein".
        reference: A protein label to pin as the first category, typically a
            standard or control sample (e.g. "BSA"). If None, no pinning is
            applied. Defaults to "BSA".
        category_order: An optional sequence of protein label strings
            specifying the desired sort order. Observed values not in this
            sequence are appended in sorted order. If None, all observed
            values are sorted alphabetically. Defaults to None.
        include_missing_categories: If True, retains all labels from
            category_order in the Categorical even if they are not observed
            in the column. If False, only observed labels are kept as
            categories. Defaults to False.

    Returns:
        A copy of df with the protein column replaced by an ordered
        Categorical whose category sequence respects category_order, appends
        unrecognized observed values in sorted order, and pins reference at
        position zero when present.

    Raises:
        KeyError: If column is not present in df.

    Examples:
        >>> categorize_proteins(df,column="protein",reference="BSA")
        >>> categorize_proteins(
        ...     df,
        ...     column="protein",
        ...     category_order=CUSTOM_PROTEIN_ORDER,
        ...     include_missing_categories=True
        ... )
    """
    _ensure_column_exists(df, column)

    df = df.copy()

    observed = _ordered_unique(df[column].dropna())
    observed_set = set(observed)

    if category_order is None:
        categories = sorted(observed)
    else:
        ordered = _ordered_unique(category_order)
        ordered_set = set(ordered)
        extras = sorted(p for p in observed if p not in ordered_set)
        categories = ordered + extras

    if not include_missing_categories:
        categories = [p for p in categories if p in observed_set]

    if reference is not None:
        categories = [p for p in categories if p != reference]

        if reference in observed_set or include_missing_categories:
            categories.insert(0, reference)

    df[column] = pd.Categorical(df[column], categories=categories, ordered=True)
    return df
