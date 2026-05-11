"""Categorical normalization utilities for HPLC peak identifiers and compound ordering.

Defines the canonical set of analyte compounds tracked across the enzyme kinetics
pipeline, together with flexible alias resolution and label normalization
infrastructure. The Compounds enum is the authoritative identifier for each
analyte; canonicalize_peak_ids resolves raw peak-ID strings to canonical
Compounds values and converts the result to an ordered pandas.Categorical
for downstream grouping and plotting.

The module also exposes coerce_peak_prefixes and the FlexPeakPrefixes
annotated alias for Pydantic models, and pre-built ordering constants
(CATEGORICAL_ORDER, CATEGORICAL_ORDER_EXTENDED) for use in axis sorting
and categorical assignment workflows.

Typical usage example:
    >>> compound = Compounds.resolve("olivetolic acid")
    >>> df = canonicalize_peak_ids(df, column="peak_id", prefixes=["C5", "C6"])
    >>> df["peak_id"].cat.categories.tolist()
    ['HTAL', 'C5-HTAL', 'C6-HTAL', 'PDAL', 'C5-PDAL', 'C6-PDAL', ...]
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from functools import lru_cache
from typing import Annotated, Any, Final

import pandas as pd
from pydantic import BeforeValidator

__all__ = [
    "canonicalize_label",
    "canonicalize_peak_ids",
    "CATEGORICAL_ORDER",
    "CATEGORICAL_ORDER_EXTENDED",
    "coerce_peak_prefixes",
    "COMPOUND_PREFIXES",
    "Compounds",
    "FlexPeakPrefixes",
    "LabelProfile",
]


# ---------------------------------------------------------------------
# Validators
# ---------------------------------------------------------------------


def coerce_peak_prefixes(value: Any) -> frozenset[str] | None:
    """Normalizes and merges peak prefix input with the default COMPOUND_PREFIXES set.

    Args:
        value: The raw field value supplied by the caller. Accepts None,
            a string, a Mapping, an iterable of strings, or any scalar.

    Returns:
        None if value is None, otherwise a frozenset combining
        COMPOUND_PREFIXES with the extracted string elements from value.
    """
    if value is None:
        return None

    if isinstance(value, str):
        extra = {value}
    elif isinstance(value, Mapping):
        extra = set(value)  # intentionally keys-only; values are display labels
    elif isinstance(value, Iterable):
        extra = {str(v) for v in value if v is not None}
    else:
        extra = {str(value)}

    return frozenset(COMPOUND_PREFIXES) | frozenset(extra)


FlexPeakPrefixes = Annotated[
    frozenset[str] | None,
    BeforeValidator(coerce_peak_prefixes),
]
"""Annotated type alias applying coerce_peak_prefixes as a Pydantic BeforeValidator."""

# ---------------------------------------------------------------------
# Compounds - parsing, categorization, and label formatting
# ---------------------------------------------------------------------

# Pre-compiled pattern matching any run of non-alphanumeric characters for label normalization.
_NORM_RE: Final = re.compile(r"[^A-Z0-9]+")

# Ordered chain-length prefixes used to build extended categorical compound labels.
_COMPOUND_PREFIX_VALUES: Final = ("C3", "C4", "C5", "C6", "C7")


@lru_cache(maxsize=256)
def canonicalize_label(value: object) -> str:
    """Normalizes a compound label to a canonical uppercase underscore-delimited string.

    Args:
        value: Any object whose string representation is to be normalized.
            Non-string values are coerced via str().

    Returns:
        The normalized label string with non-alphanumeric runs replaced by
        underscores and leading or trailing underscores stripped.

    Examples:
        >>> canonicalize_label("olivetolic acid")
        'OLIVETOLIC_ACID'
        >>> canonicalize_label("  OLA  ")
        'OLA'
        >>> canonicalize_label("hexanoyl-triacetic acid")
        'HEXANOYL_TRIACETIC_ACID'
    """
    return _NORM_RE.sub("_", str(value).strip().upper()).strip("_")


class Compounds(StrEnum):
    """Enumeration of canonical analyte compound identifiers.

    Serves as the authoritative identity for each tracked compound across the
    HPLC analysis pipeline. Each member's value is the canonical string label
    used in column names, categorical ordering, and retention time registry keys.

    Attributes:
        HTAL: Hexanoyl triacetic acid lactone.
        PDAL: Pentyl diacetic acid lactone.
        OLA: Olivetolic acid.
        OLV: Olivetol.
    """

    HTAL = "HTAL"
    PDAL = "PDAL"
    OLA = "OLA"
    OLV = "OLV"

    # VOID_PEAK = "VOID_PEAK"

    @classmethod
    def resolve(cls, label: str | Compounds) -> Compounds:
        """Resolves a string or Compounds instance to a canonical Compounds member.

        If label is already a Compounds instance it is returned directly.
        Otherwise the label is normalized via canonicalize_label and looked up in
        the internal alias registry, which covers member names, values, and
        entries in _COMPOUND_ALIASES.

        Args:
            label: The compound identifier to resolve. Accepts a Compounds
                member, a canonical name or value string, or any alias string
                registered in _COMPOUND_ALIASES.

        Returns:
            The matching Compounds member.

        Raises:
            ValueError: If label does not match any known compound name, value,
                or alias after normalization.

        Examples:
            >>> Compounds.resolve("ola")
            <Compounds.OLA: 'OLA'>
            >>> Compounds.resolve("olivetolic acid")
            <Compounds.OLA: 'OLA'>
        """
        if isinstance(label, cls):
            return label

        try:
            return _COMPOUNDS_LOOKUP[canonicalize_label(label)]
        except KeyError as exc:
            raise ValueError(f"Unknown or invalid compound label: {label!r}") from exc

    @classmethod
    def values(cls) -> tuple[str, ...]:
        """Returns a tuple of all member values in definition order."""
        return tuple(member.value for member in cls)

    @classmethod
    def labels(cls, profile: LabelProfile) -> tuple[str, ...]:
        """Returns formatted label strings for all members under a given profile.

        Args:
            profile: The LabelProfile controlling prefix, separator, and any
                per-compound value overrides.

        Returns:
            A tuple of formatted label strings in member definition order,
            each produced by format_label.
        """
        return tuple(member.format_label(profile) for member in cls)

    def format_label(self, profile: LabelProfile) -> str:
        """Formats a compound label string using the supplied profile.

        Args:
            profile: The LabelProfile controlling prefix, separator, and any
                per-compound value overrides.

        Returns:
            A formatted label string of the form "<prefix><sep><base>".

        Examples:
            >>> p = LabelProfile(prefix="C5")
            >>> Compounds.OLA.format_label(p)
            'C5-OLA'
        """
        base = profile.overrides.get(self, self.value)
        return f"{profile.prefix}{profile.sep}{base}"


@dataclass(frozen=True, slots=True)
class LabelProfile:
    """Configuration for compound label formatting.

    Controls the prefix, separator, and optional per-compound value overrides
    applied by Compounds.format_label and Compounds.labels.

    Attributes:
        prefix: The string prepended to every formatted label (e.g. "C5").
        sep: The separator inserted between prefix and compound base label.
            Defaults to "-".
        overrides: Optional mapping from Compounds members to replacement base
            label strings, applied before the prefix is prepended. Defaults
            to an empty dict.
    """

    prefix: str
    sep: str = "-"
    overrides: Mapping[Compounds, str] = field(default_factory=dict)


# Mapping from each Compounds member to its accepted alias strings for label resolution.
_COMPOUND_ALIASES: Final[Mapping[Compounds, tuple[str, ...]]] = {
    Compounds.HTAL: ("hexanoyl triacetic acid", "hexanoyl-triacetic acid"),
    Compounds.PDAL: ("pentyl diacetic acid", "pentyl-diacetic acid"),
    Compounds.OLA: ("olivetolic acid", "oa"),
    Compounds.OLV: ("olivetol", "olivetolate"),
}


def _iter_lookup_keys(member: Compounds) -> Iterable[str]:
    """Yields all lookup key strings for a Compounds member, including aliases."""
    yield member.name
    yield member.value
    yield from _COMPOUND_ALIASES.get(member, ())


def _build_compound_lookup() -> dict[str, Compounds]:
    """Builds the normalized label-to-Compounds lookup dict from all members and their aliases.

    Iterates every Compounds member, normalizes each of its lookup keys via
    canonicalize_label, and inserts them into a shared dict. Raises on any
    collision where two distinct members map to the same normalized key, which
    would indicate an ambiguous alias definition.

    Returns:
        A dict mapping normalized label strings to their corresponding
        Compounds members, covering member names, values, and all aliases.

    Raises:
        ValueError: If two distinct Compounds members produce the same
            normalized key, indicating an alias collision.
    """
    lookup: dict[str, Compounds] = {}

    for member in Compounds:
        for key in _iter_lookup_keys(member):
            normalized = canonicalize_label(key)
            existing = lookup.get(normalized)

            if existing is not None and existing is not member:
                raise ValueError(
                    f"Collision: {key!r} maps to both "
                    f"{existing.value!r} and {member.value!r}",
                )

            lookup[normalized] = member

    return lookup


# Normalized label-to-Compounds mapping built at import time from all members and aliases.
_COMPOUNDS_LOOKUP: Final = _build_compound_lookup()

COMPOUND_PREFIXES: Final[tuple[str, ...]] = _COMPOUND_PREFIX_VALUES
"""Ordered chain-length prefix strings used to construct extended compound labels."""

CATEGORICAL_ORDER: Final[tuple[str, ...]] = Compounds.values()
"""Canonical ordering of base compound labels for categorical axis sorting."""

CATEGORICAL_ORDER_EXTENDED: Final[tuple[str, ...]] = (
    *CATEGORICAL_ORDER,
    *(
        f"{prefix}-{compound}"
        for prefix in COMPOUND_PREFIXES
        for compound in CATEGORICAL_ORDER
    ),
)
"""Extended categorical ordering including all prefix-compound combinations after base compounds."""


# ---------------------------------------------------------------------
# Canonicalization
# ---------------------------------------------------------------------


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
