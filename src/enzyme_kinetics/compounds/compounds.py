"""Compound identity, retention time windows, and label formatting for HPLC analysis.

Defines the canonical set of analyte compounds tracked across the enzyme kinetics
pipeline, together with their expected chromatographic retention time windows and
flexible label formatting infrastructure. The Compounds enum is the authoritative
identifier for each analyte; PeakWindowRegistry manages per-compound RTWindow
objects and supports runtime overrides for method-specific retention time tuning.

The module also exposes pre-built ordering constants and a default registry
instance for use by downstream assignment and quantitation workflows.

Typical usage example:
    >>> compound = Compounds.resolve("ola")
    >>> registry = PeakWindowRegistry.from_overrides([["OLA", "13.5", "15.5"]])
    >>> window = registry.get(compound)
    >>> window.as_tuple()
    (13.5, 15.5)
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from functools import lru_cache
from typing import Final, Self

__all__ = [
    "CATEGORICAL_ORDER",
    "CATEGORICAL_ORDER_EXTENDED",
    "COMPOUND_PREFIXES",
    "COMPOUND_RT_MAP",
    "Compounds",
    "DEFAULT_RT_REGISTRY",
    "LabelProfile",
    "canonicalize_label",
    "PeakWindowRegistry",
    "profile_labels",
    "RTWindow",
]

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
    def resolve(cls, label: str | Compounds) -> Self:
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


@dataclass(frozen=True, slots=True)
class RTWindow:
    """A half-open retention time window defining the expected elution range for a compound.

    Attributes:
        min: The inclusive lower bound of the retention time window in minutes.
        max: The exclusive upper bound of the retention time window in minutes.
    """

    min: float
    max: float

    def __post_init__(self) -> None:
        """Validates that min is strictly less than max.

        Raises:
            ValueError: If min is greater than or equal to max.
        """
        if self.min >= self.max:
            raise ValueError(f"Invalid RT window: {self.min} must be < {self.max}")

    def as_tuple(self) -> tuple[float, float]:
        """Returns the retention time window as a (min, max) float tuple."""
        return self.min, self.max


class PeakWindowRegistry:
    """Registry mapping compound identifiers to their expected retention time windows.

    Stores RTWindow objects keyed by canonical compound strings. Supports
    programmatic updates via set() and apply(), bulk construction via classmethods,
    and serialization back to a plain dict via as_dict().
    """

    __slots__ = ("_windows",)

    def __init__(
        self,
        initial: Mapping[str | Compounds, tuple[float, float]] | None = None,
    ) -> None:
        """Initializes the registry, optionally pre-populating it from a mapping.

        Args:
            initial: Optional mapping from compound identifiers to (min, max)
                float tuples. Each entry is passed through set() so keys are
                normalized and values are validated as RTWindow objects. If
                None, the registry is initialized empty. Defaults to None.
        """
        self._windows: dict[str, RTWindow] = {}

        if initial:
            for compound, window in initial.items():
                self.set(compound, *window)

    @staticmethod
    def key(compound: str | Compounds) -> str:
        """Returns the canonical registry key for a compound identifier.

        Attempts to resolve compound via Compounds.resolve to obtain the
        canonical member value; falls back to canonicalize_label if resolution
        fails, accommodating unknown or ad-hoc compound strings.

        Args:
            compound: The compound identifier to convert to a registry key.

        Returns:
            The canonical Compounds value string if compound is resolvable,
            otherwise the canonicalize_label result.
        """
        try:
            return Compounds.resolve(compound).value
        except ValueError:
            return canonicalize_label(compound)

    def set(self, compound: str | Compounds, rt_min: float, rt_max: float) -> None:
        """Sets the retention time window for a compound.

        Args:
            compound: The compound identifier to register. Resolved to a
                canonical key via key().
            rt_min: The inclusive lower bound of the retention time window
                in minutes.
            rt_max: The exclusive upper bound of the retention time window
                in minutes.

        Raises:
            ValueError: If rt_min is greater than or equal to rt_max.
        """
        self._windows[self.key(compound)] = RTWindow(float(rt_min), float(rt_max))

    def get(self, compound: str | Compounds) -> RTWindow:
        """Returns the RTWindow registered for a compound.

        Args:
            compound: The compound identifier to look up. Resolved to a
                canonical key via key().

        Returns:
            The RTWindow registered under the compound's canonical key.

        Raises:
            KeyError: If no window has been registered for the compound.
        """
        return self._windows[self.key(compound)]

    def apply(self, overrides: Iterable[Sequence[str]] | None) -> None:
        """Applies a sequence of retention time window overrides to the registry.

        Each override entry must be a three-element sequence of strings in the
        form ["COMPOUND", "MIN", "MAX"]. Entries are applied in iteration order,
        so later entries overwrite earlier ones for the same compound.

        Args:
            overrides: An iterable of three-element string sequences, each
                representing a compound identifier, a minimum retention time,
                and a maximum retention time. If None, the method returns
                without modifying the registry.

        Raises:
            ValueError: If any entry does not contain exactly three elements.
            ValueError: If the min or max value of any entry cannot be
                converted to float, or if the resulting RTWindow is invalid.
        """
        if overrides is None:
            return

        for entry in overrides:
            if len(entry) != 3:
                raise ValueError(
                    f"Invalid override format ['COMPOUND', 'MIN', 'MAX']: {entry!r}",
                )

            label, rt_min, rt_max = entry

            try:
                self.set(label, float(rt_min), float(rt_max))
            except ValueError as exc:
                raise ValueError(f"Error updating window for {label!r}: {exc}") from exc

    def as_dict(self) -> dict[str, tuple[float, float]]:
        """Returns the registry contents as a plain dict of (min, max) tuples."""
        return {key: window.as_tuple() for key, window in self._windows.items()}

    @classmethod
    def default(cls) -> Self:
        """Returns a new registry pre-populated from COMPOUND_RT_MAP."""
        return cls(COMPOUND_RT_MAP)

    @classmethod
    def from_dict(cls, data: Mapping[str, Sequence[float]]) -> Self:
        """Returns a new registry constructed from a mapping of compound keys to window bounds.

        Args:
            data: A mapping from compound identifier strings to two-element
                sequences of (min, max) float values.

        Returns:
            A new PeakWindowRegistry with all entries from data registered
            as RTWindow objects.
        """
        return cls({key: (values[0], values[1]) for key, values in data.items()})

    @classmethod
    def from_overrides(cls, overrides: Iterable[Sequence[str]] | None) -> Self:
        """Returns a default registry with the supplied overrides applied.

        Constructs a registry from COMPOUND_RT_MAP via default(), then calls
        apply() with overrides to update any specified windows.

        Args:
            overrides: An iterable of three-element string sequences accepted
                by apply(), or None to return the unmodified default registry.

        Returns:
            A new PeakWindowRegistry seeded from COMPOUND_RT_MAP with all
            overrides applied.

        Raises:
            ValueError: If any override entry is malformed or contains invalid
                numeric values, as raised by apply().
        """
        registry = cls.default()
        registry.apply(overrides)
        return registry


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


def profile_labels(profile: LabelProfile) -> tuple[str, ...]:
    """Returns formatted label strings for all Compounds members under a given profile."""
    return Compounds.labels(profile)


# Normalized label-to-Compounds mapping built at import time from all members and aliases.
_COMPOUNDS_LOOKUP: Final = _build_compound_lookup()

COMPOUND_PREFIXES: Final[tuple[str, ...]] = _COMPOUND_PREFIX_VALUES
"""Ordered chain-length prefix strings used to construct extended compound labels."""

CATEGORICAL_ORDER: Final[tuple[str, ...]] = Compounds.values()
"""Canonical ordering of base compound labels for categorical axis sorting."""

CATEGORICAL_ORDER_EXTENDED: Final[tuple[str, ...]] = (
    *CATEGORICAL_ORDER,
    *(f"{prefix}-{compound}" for prefix in COMPOUND_PREFIXES for compound in CATEGORICAL_ORDER),
)
"""Extended categorical ordering including all prefix-compound combinations after base compounds."""

COMPOUND_RT_MAP: Final[Mapping[str, tuple[float, float]]] = {
    Compounds.HTAL.value: (9.0, 10.0),
    Compounds.PDAL.value: (11.0, 12.0),
    Compounds.OLA.value: (14.0, 16.0),
    Compounds.OLV.value: (16.5, 18.0),
}
"""Default retention time windows (min, max) in minutes for each compound."""

DEFAULT_RT_REGISTRY: Final = PeakWindowRegistry.default()
"""Default PeakWindowRegistry instance pre-populated from COMPOUND_RT_MAP."""
