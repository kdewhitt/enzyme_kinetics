"""Column-label cleaning for instrument tables read from CSV or LabSolutions exports.

Instrument software names its columns for people, not for programs: R.Time,
Area (%), Unnamed: 3, and the same detector heading repeated once per channel.
Attribute access, groupby keys, and schema lookups all need the opposite --
stable, unique, snake_case identifiers -- so tables entering the pipeline pass
through here before anything else touches them.

Cleaning is two independent steps that the options record turns on separately.
Artifact removal discards the columns a spreadsheet round trip leaves behind,
with _INDEXLIKE_RE deciding which those are. Relabelling routes each surviving
label through to_snake_case and then resolves collisions by numeric suffix, so
it never drops a column: a frame with three identical headings keeps all three.
Removal comes first when both run, because relabelling rewrites the very
punctuation the artifact pattern looks for. Both steps accept flat and
MultiIndex columns, both ignore labels that are not strings, and neither
disturbs the frame's values, row index, or attrs -- relabelling is applied with
set_axis, so no data is copied unless the caller writes to the result.

Example:
    frame = read_clean_csv(
        "peaks.csv",
        options=CleaningOptions(drop_indexlike=True, verbose=True),
    )
"""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

import pandas as pd
from rich.console import Console
from rich.table import Table

from .case import to_snake_case

__all__ = [
    "DEFAULT_CLEANING_OPTIONS",
    "CleaningOptions",
    "clean_index",
    "column_renames",
    "drop_index_artifacts",
    "read_clean_csv",
    "render_rename_table",
    "snake_case_columns",
]

type Label = str | tuple[str, ...]

# Bare "index" or "idx", and the "Unnamed: N" headings pandas writes for an
# unlabelled index column, in both their raw and already-relabelled spellings.
# Matched with fullmatch, so the real columns "index_value" and "unnamed 0"
# survive.
_INDEXLIKE_RE: Final[re.Pattern[str]] = re.compile(
    r"(?:index|idx|unnamed[:_]\s*\d+(?:\.\d+)?)",
    flags=re.IGNORECASE,
)

# Delimiter between a repeated label and the counter that disambiguates it.
_SUFFIX_SEPARATOR: Final[str] = "_"

# Separator between MultiIndex levels when a label is rendered for a terminal.
_LEVEL_SEPARATOR: Final[str] = " / "


@dataclass(frozen=True, slots=True)
class CleaningOptions:
    """Label-cleaning policy shared by every table-reading entry point.

    Attributes:
        snake_case: Whether column labels are rewritten to unique snake_case
            identifiers. False leaves labels exactly as read, which suits frames
            that are plotted rather than queried by name.
        preserve_unicode: Whether Unicode letters and digits survive relabelling
            instead of the label being reduced to ASCII.
        drop_indexlike: Whether columns matching _INDEXLIKE_RE are discarded.
        verbose: Whether a summary of the changed labels is rendered.
        console: Rich console receiving that summary; a default console is
            created when None.
    """

    snake_case: bool = True
    preserve_unicode: bool = False
    drop_indexlike: bool = False
    verbose: bool = False
    console: Console | None = None


DEFAULT_CLEANING_OPTIONS: Final[CleaningOptions] = CleaningOptions()


# ---------------------------------------------------------------------------
# Label helpers
# ---------------------------------------------------------------------------


def _with_suffix[LabelT: (str, tuple[str, ...])](label: LabelT, count: int) -> LabelT:
    """Return the label with a counter appended, on its last level if it is a tuple."""
    if isinstance(label, tuple):
        if not label:
            return label
        *prefix, final = label
        return (*prefix, f"{final}{_SUFFIX_SEPARATOR}{count}")  # type: ignore[return-value]
    return f"{label}{_SUFFIX_SEPARATOR}{count}"  # type: ignore[return-value]


def _make_unique[LabelT: (str, tuple[str, ...])](labels: Sequence[LabelT]) -> list[LabelT]:
    """Disambiguate repeated labels by appending a counter, keeping the original order.

    Nothing is dropped: the first occurrence of a label keeps it and later ones
    take the next free counter. A candidate suffix that already exists elsewhere
    in the sequence is skipped, so a frame whose headings already include both
    "a" and "a_2" cannot be given a duplicate by the disambiguation itself.

    Args:
        labels: Column labels in their original order.

    Returns:
        A list of unique labels, one per input label.
    """
    counts: defaultdict[LabelT, int] = defaultdict(int)
    used: set[LabelT] = set()
    result: list[LabelT] = []

    for label in labels:
        count = counts[label] + 1
        candidate = label if count == 1 else _with_suffix(label, count)

        while candidate in used:
            count += 1
            candidate = _with_suffix(label, count)

        counts[label] = count
        used.add(candidate)
        result.append(candidate)

    return result


def _join_label(label: object) -> str:
    """Render a flat or MultiIndex label as one string for terminal display."""
    if isinstance(label, tuple):
        return _LEVEL_SEPARATOR.join(map(str, label))
    return str(label)


def _is_indexlike(label: object) -> bool:
    """Report whether a label is an index artifact rather than a real column.

    A MultiIndex label counts as an artifact when any of its levels does, which
    covers the ("Unnamed: 0", "") heading a pandas round trip writes for an
    unlabelled index. Non-string labels are never artifacts, so integer and
    datetime columns pass through untouched instead of raising.

    Args:
        label: A column label of any type.

    Returns:
        True when the label matches _INDEXLIKE_RE.
    """
    if isinstance(label, tuple):
        return any(_is_indexlike(part) for part in label)
    return isinstance(label, str) and _INDEXLIKE_RE.fullmatch(label) is not None


# ---------------------------------------------------------------------------
# Index-level cleaning
# ---------------------------------------------------------------------------


def clean_index(columns: pd.Index, *, preserve_unicode: bool = False) -> pd.Index:
    """Return a column index of unique snake_case labels.

    Every level of a MultiIndex is converted independently and the resulting
    tuples are disambiguated as a whole, so two columns differing in one level
    only are left alone. Level names and the flat index name are carried over.

    Args:
        columns: The column index to convert.
        preserve_unicode: Whether Unicode letters and digits survive conversion.

    Returns:
        A new index of the same shape, with unique labels.
    """
    if isinstance(columns, pd.MultiIndex):
        cleaned_tuples = _make_unique(
            [
                tuple(to_snake_case(part, preserve_unicode=preserve_unicode) for part in label)
                for label in columns
            ],
        )
        return pd.MultiIndex.from_tuples(cleaned_tuples, names=columns.names)

    cleaned_labels = _make_unique(
        [to_snake_case(label, preserve_unicode=preserve_unicode) for label in columns],
    )
    return pd.Index(cleaned_labels, name=columns.name)


def column_renames(
    columns: pd.Index,
    *,
    preserve_unicode: bool = False,
    changed_only: bool = False,
) -> list[tuple[Label, Label]]:
    """Pair each original column label with the label clean_index gives it.

    Pairs are returned in column order rather than as a mapping because the
    originals need not be unique -- disambiguating them is the point -- and a
    dict would silently keep only the last of each repeated label. Callers that
    know their labels are unique can build one with dict(column_renames(...))
    and hand it to DataFrame.rename; callers auditing or logging a conversion
    want the full list.

    Args:
        columns: The column index to convert.
        preserve_unicode: Whether Unicode letters and digits survive conversion.
        changed_only: Whether pairs whose label is unchanged are omitted.

    Returns:
        A list of (original, cleaned) pairs.
    """
    cleaned = clean_index(columns, preserve_unicode=preserve_unicode)
    pairs = list(zip(columns, cleaned, strict=True))
    return [(old, new) for old, new in pairs if old != new] if changed_only else pairs


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def render_rename_table(
    changes: Sequence[tuple[Label, Label]],
    *,
    console: Console | None = None,
) -> None:
    """Print a Rich table of renamed columns, or a note when nothing changed.

    Args:
        changes: The (original, cleaned) pairs to show, already filtered to the
            labels that differ.
        console: Rich console receiving the output; a default console is created
            when None.
    """
    console = console or Console()

    if not changes:
        console.print("[dim]No column labels changed.[/dim]")
        return

    table = Table(title="Column Label Changes", show_header=True, header_style="bold")
    table.add_column("#", justify="right", style="dim")
    table.add_column("Original")
    table.add_column("Cleaned")

    for position, (old, new) in enumerate(changes, start=1):
        table.add_row(str(position), _join_label(old), _join_label(new))

    console.print(table)


# ---------------------------------------------------------------------------
# Frame-level cleaning
# ---------------------------------------------------------------------------


def snake_case_columns(
    frame: pd.DataFrame,
    *,
    options: CleaningOptions = DEFAULT_CLEANING_OPTIONS,
    copy: bool = True,
) -> pd.DataFrame:
    """Give a frame unique snake_case column labels.

    The default relabels through set_axis, which under copy-on-write shares the
    frame's blocks instead of duplicating them: the result is still isolated
    from the source, but a chromatogram table costs a new index rather than a
    second copy of its values. Pass copy=False to relabel the caller's frame in
    place, which is only worth it when no one else holds a reference to it.

    Args:
        frame: The frame to relabel.
        options: Cleaning policy; only preserve_unicode, verbose, and console
            are consulted here.
        copy: Whether the caller's frame is left untouched.

    Returns:
        The relabelled frame -- a new object unless copy is False.
    """
    original = frame.columns
    cleaned = clean_index(original, preserve_unicode=options.preserve_unicode)

    if options.verbose:
        changes = [(old, new) for old, new in zip(original, cleaned, strict=True) if old != new]
        render_rename_table(changes, console=options.console)

    if not copy:
        frame.columns = cleaned
        return frame

    return frame.set_axis(cleaned, axis="columns")


def drop_index_artifacts(frame: pd.DataFrame) -> pd.DataFrame:
    """Drop the leftover index columns a spreadsheet round trip writes into a table.

    Selection is done label by label against _INDEXLIKE_RE rather than through
    the string accessor, so frames with integer, null, or MultiIndex columns are
    filtered instead of raising. The frame is returned unchanged, and not copied,
    when it holds no artifacts -- the common case for every caller.

    Args:
        frame: The frame to filter.

    Returns:
        The frame without its index-like columns.
    """
    keep = [not _is_indexlike(label) for label in frame.columns]
    return frame if all(keep) else frame.loc[:, keep]


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def read_clean_csv(
    path: str | Path,
    *,
    options: CleaningOptions = DEFAULT_CLEANING_OPTIONS,
    **read_csv_kwargs: Any,
) -> pd.DataFrame:
    """Read a CSV and apply the cleaning steps its options ask for.

    Artifact removal runs first, against the headings as written: relabelling
    rewrites the colon in Unnamed: 0 and would otherwise hand _INDEXLIKE_RE a
    label it no longer recognizes. Dropping first also keeps the discarded
    columns out of the disambiguation, so a real column is never pushed onto a
    numeric suffix by an artifact that is about to disappear. Keyword arguments
    are passed straight to pandas.read_csv, which is where dtype, separator, and
    index_col belong.

    Args:
        path: The CSV to read.
        options: Cleaning policy applied to the parsed frame.
        **read_csv_kwargs: Additional arguments for pandas.read_csv.

    Returns:
        The parsed and cleaned frame.
    """
    frame = pd.read_csv(path, **read_csv_kwargs)

    if options.drop_indexlike:
        frame = drop_index_artifacts(frame)

    return snake_case_columns(frame, options=options, copy=False) if options.snake_case else frame
