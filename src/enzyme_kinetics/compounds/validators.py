"""Peak prefix coercion validator and annotated type alias for Shimadzu HPLC peak identification."""

from __future__ import annotations

from typing import Annotated, Any, Iterable, Mapping

from pydantic import BeforeValidator

from .compounds import COMPOUND_PREFIXES

__all__ = ["coerce_peak_prefixes", "FlexPeakPrefixes"]


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


FlexPeakPrefixes = Annotated[frozenset[str] | None, BeforeValidator(coerce_peak_prefixes)]
"""Annotated type alias applying coerce_peak_prefixes as a Pydantic BeforeValidator."""
