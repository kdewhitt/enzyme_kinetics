"""Label normalization helpers for column names, slugs, and filesystem-safe tokens."""

from __future__ import annotations

import math
import re
from collections.abc import Iterable
from functools import lru_cache
from typing import Final

import unicodedata

__all__ = ["asciify", "to_snake_case", "to_slug"]

# ---------------------------------------------------------------------------
# Tokenization
# ---------------------------------------------------------------------------

# Word boundaries inside a compound token, matched in one pass:
#   1. an acronym running into a capitalized word ("HTTPServer" -> "HTTP Server")
#   2. a lowercase letter or digit running into a capital ("fooBar" -> "foo Bar",
#      "v2Api" -> "v2 Api")
# Both alternatives are pure lookaround, so the match is zero-width and the
# replacement is a literal space instead of a backreference expansion.
_WORD_BOUNDARY_RE: Final[re.Pattern[str]] = re.compile(
    r"(?<=[A-Z])(?=[A-Z][a-z])|(?<=[a-z0-9])(?=[A-Z])",
)

# Runs of everything that is not an ASCII letter or digit.
_NON_ALNUM_RE: Final[re.Pattern[str]] = re.compile(r"[^a-zA-Z0-9]+")

# Combining marks that can outlive NFKC composition, either because the base
# letter has no precomposed form or because lowercasing introduced them. They
# belong to the word they sit on, so the Unicode filter below has to keep them.
_COMBINING_MARKS: Final[str] = r"\u0300-\u036f\u1ab0-\u1aff\u1dc0-\u1dff\u20d0-\u20f0\ufe20-\ufe2f"

# Runs of everything that is not a Unicode letter, digit, or combining mark. \W
# is Unicode-aware but treats the underscore as a word character, so the
# underscore is named as a separator explicitly.
_NON_ALNUM_UNICODE_RE: Final[re.Pattern[str]] = re.compile(
    rf"(?:[^\w{_COMBINING_MARKS}]|_)+",
)

# Delimiter placed between words in a normalized label.
_SEPARATOR: Final[str] = "_"

# Upper bound on the label vocabulary held by the conversion worker. Column
# labels and schema keys repeat heavily across frames and detectors, so a small
# cache absorbs almost every repeat call.
_CACHE_SIZE: Final[int] = 2048


def _insert_word_boundaries(text: str) -> str:
    """Insert spaces at case-change boundaries within a compound token."""
    return _WORD_BOUNDARY_RE.sub(" ", text)


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

# Characters NFKD cannot reduce to an ASCII base, which asciify would otherwise
# drop: ligatures and barred or slashed Latin letters, plus the micro sign that
# instrument exports use for microlitre and microabsorbance units.
_TRANSLITERATIONS: Final[dict[int, str]] = str.maketrans(
    {
        "ß": "ss",
        "ẞ": "SS",
        "æ": "ae",
        "Æ": "AE",
        "œ": "oe",
        "Œ": "OE",
        "ø": "o",
        "Ø": "O",
        "đ": "d",
        "Đ": "D",
        "ð": "d",
        "Ð": "D",
        "þ": "th",
        "Þ": "TH",
        "ł": "l",
        "Ł": "L",
        "ħ": "h",
        "Ħ": "H",
        "ı": "i",
        "İ": "I",  # noqa: RUF001 - the table exists to resolve look-alikes
        "ŋ": "ng",
        "Ŋ": "NG",
        "µ": "u",
        "μ": "u",
    },
)


def asciify(text: str) -> str:
    """Reduce text to ASCII, transliterating what decomposition cannot reach.

    NFKD decomposition separates base letters from their combining marks, so an
    accented Latin letter survives as its unaccented base. Characters that have
    no ASCII base are transliterated first when they appear in _TRANSLITERATIONS
    and dropped otherwise, which is the fate of every non-Latin script.

    Args:
        text: The string to reduce.

    Returns:
        An ASCII-only string, returned unchanged when text is already ASCII.
    """
    if text.isascii():
        return text

    return (
        unicodedata.normalize("NFKD", text.translate(_TRANSLITERATIONS))
        .encode("ascii", "ignore")
        .decode("ascii")
    )


# ---------------------------------------------------------------------------
# Transformation
# ---------------------------------------------------------------------------


@lru_cache(maxsize=_CACHE_SIZE)
def _snake_case(text: str, preserve_unicode: bool) -> str:
    """Convert one non-empty, already-stringified label to snake_case.

    Ordering carries the correctness here. Folding to ASCII precedes boundary
    detection so that accented capitals are still recognized as boundaries, and
    lowercasing comes last so that the combining marks it can introduce -- the
    dot a dotted capital I leaves behind -- are kept as part of the word instead
    of being read as separators. Results are cached because callers convert the
    same small vocabulary of labels repeatedly; see _CACHE_SIZE.

    Args:
        text: A stripped, non-empty label.
        preserve_unicode: Whether to keep Unicode letters and digits.

    Returns:
        A snake_case string, empty when the label holds no letters or digits.
    """
    if preserve_unicode:
        # NFKC yields composed, canonically equivalent forms, so labels that
        # differ only in how an accent is encoded normalize to the same name.
        text = _insert_word_boundaries(unicodedata.normalize("NFKC", text))
        text = _NON_ALNUM_UNICODE_RE.sub(_SEPARATOR, text)
    else:
        text = _insert_word_boundaries(asciify(text))
        text = _NON_ALNUM_RE.sub(_SEPARATOR, text)

    # Both patterns match greedily, so repeated separators cannot survive and
    # only the outer padding is left to trim.
    return text.strip(_SEPARATOR).lower()


def to_snake_case(
    value: object,
    *,
    preserve_unicode: bool = False,
    fallback: str = "",
) -> str:
    """Convert a label of any type to snake_case.

    Compound tokens are split at case-change boundaries before delimiting, so
    acronyms and camelCase survive as separate words.

    Args:
        value: Object to convert; coerced with str unless it is null.
        preserve_unicode: Whether to keep Unicode letters and digits instead of
            reducing the label to ASCII.
        fallback: Returned when the label is null, blank, or holds no letters or
            digits to keep.

    Returns:
        A snake_case string, or fallback when nothing survives conversion.
    """
    # Labels arriving from a pandas Index may be None or NaN; neither is a name,
    # and str would silently turn them into the columns "none" and "nan".
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return fallback

    text = str(value).strip()
    if not text:
        return fallback

    return _snake_case(text, preserve_unicode) or fallback


def to_slug(
    value: object,
    *,
    separator: str = _SEPARATOR,
    lowercase: bool = True,
) -> str:
    """Join one or more labels into a single filesystem-safe slug.

    Unlike to_snake_case this does not retokenize: a label is filtered and
    joined as written, so an existing compound such as DluHKS stays intact
    rather than splitting at its internal capital. Each label is reduced to
    ASCII before non-alphanumeric runs collapse to the separator, and labels
    that contribute no characters are dropped rather than leaving empty
    segments behind.

    Args:
        value: A single label, or an iterable of labels to join.
        separator: Single character placed between labels and in place of
            non-alphanumeric runs within them.
        lowercase: Whether to fold the result to lowercase, which keeps slugs
            stable on case-insensitive filesystems.

    Returns:
        A slug built from the surviving labels, empty when none survive.
    """
    items: Iterable[object]
    if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
        items = (value,)
    else:
        items = value

    parts: list[str] = []
    for item in items:
        if item is None:
            continue
        cleaned = _NON_ALNUM_RE.sub(separator, asciify(str(item))).strip(separator)
        if cleaned:
            parts.append(cleaned.lower() if lowercase else cleaned)

    return separator.join(parts)
