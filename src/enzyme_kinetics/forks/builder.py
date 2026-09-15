"""Fluent filesystem path construction utilities with versioning and tag support.

Provides PathBuilder, an immutable, dataclass-based builder for assembling output
file paths with prefix, suffix, and end tags, collision-driven version numbering,
and byte-safe filename truncation. Resolution is explicit: resolve reads the
filesystem once and caches its answer, while reserve additionally claims the name
atomically so concurrent writers targeting one directory cannot agree on the same
path. Designed for batch processing pipelines where output paths must be
deterministic, collision-safe, and length-constrained.

Typical usage:
    builder = PathBuilder.for_file("/data/output/report.csv", create=True)
    path = builder.with_timestamp().with_tag("final").reserve()
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterator
from dataclasses import dataclass, field, replace
from datetime import datetime, tzinfo
from pathlib import Path
from typing import Final, Literal, Self

from loguru import logger

__all__ = ["PathBuilder", "TagPosition", "WINDOWS_MAX_PATH", "sanitize_component"]

# Matches a stem ending in _vN to extract the base stem and version number.
_VERSION_RE: Final = re.compile(r"^(.*)_v(\d+)$")

# Matches characters outside the portable cross-platform filename set.
_UNSAFE_ASCII_RE: Final = re.compile(r"[^a-zA-Z0-9\-_.]")

# Matches unsafe characters once unicode letters and digits are permitted.
_UNSAFE_UNICODE_RE: Final = re.compile(r"[^\w\-.]", re.UNICODE)

# Matches runs of underscores left behind by unsafe character replacement.
_UNDERSCORE_RUN_RE: Final = re.compile(r"_{2,}")

# Device names Windows reserves in every directory, with or without an extension.
_WINDOWS_RESERVED: Final[frozenset[str]] = frozenset(
    ("CON", "PRN", "AUX", "NUL")
    + tuple(f"COM{digit}" for digit in range(1, 10))
    + tuple(f"LPT{digit}" for digit in range(1, 10)),
)

# Highest version the collision sequence emits before it reports exhaustion.
_MAX_VERSION: Final = 999

# Stem bytes held in reserve so truncation is identical across every version.
_VERSION_RESERVE_BYTES: Final = len(f"_v{_MAX_VERSION}".encode())

# Maximum characters in a full path on Windows without long path support enabled.
WINDOWS_MAX_PATH: Final = 260

type TagPosition = Literal["prefix", "suffix", "end"]

# Maps each tag position onto the field holding the tags placed there.
_TAG_FIELDS: Final[dict[str, str]] = {
    "prefix": "prefix_tags",
    "suffix": "suffix_tags",
    "end": "end_tags",
}


def sanitize_component(text: str, *, allow_unicode: bool = False) -> str:
    """Reduce text to characters that are safe in a filename on every platform.

    Unsafe characters become underscores, runs of underscores collapse to one,
    and leading or trailing underscores and trailing dots are dropped, since
    Windows rejects names ending in a dot. A leading dot survives so dotfile
    stems stay intact.

    Args:
        text: Text to reduce to a single safe filename component.
        allow_unicode: If True, unicode letters and digits survive instead of
            being replaced, leaving multi-byte characters for byte-safe
            truncation to handle.

    Returns:
        The sanitized component, or "" if no safe characters remain.
    """
    pattern = _UNSAFE_UNICODE_RE if allow_unicode else _UNSAFE_ASCII_RE
    collapsed = _UNDERSCORE_RUN_RE.sub("_", pattern.sub("_", text))
    cleaned = collapsed.strip("_").rstrip(".")
    if cleaned != text:
        logger.debug("Sanitized {!r} to {!r}.", text, cleaned)
    return cleaned


@dataclass(frozen=True, slots=True)
class PathBuilder:
    """An immutable, chainable builder for versioned, collision-safe file paths.

    Every with_ method returns a new instance through dataclasses.replace, so the
    original is never altered and one configured builder can safely seed many
    output names. Resolution touches the filesystem, so it is explicit: resolve
    reports the path a write would land on, reserve claims that path atomically,
    and the path property caches the first answer so one instance always reports
    the same location.

    Filenames are assembled as prefix tags, stem, suffix tags, end tags, version,
    and extension. The unversioned filename counts as version 1, so the first
    collision produces _v2 and the _v1 segment is never emitted.

    Attributes:
        directory: Destination directory the built path resolves into.
        stem: Base filename stem, before any prefix or suffix tags.
        extension: File extension, including the leading dot.
        overwrite: Whether an existing file at the resolved path is overwritten
            rather than triggering version numbering.
        max_name_bytes: Maximum length of the resulting filename, in UTF-8 bytes.
            Defaults to 255.
        max_path_chars: Upper bound on the length of the full resolved path. If
            None, path length goes unchecked; WINDOWS_MAX_PATH enforces the
            default Windows limit.
        source_version: Version parsed off a trailing _vN stem by for_file, kept
            so continue_versioning can restore it.
        version_floor: Lowest version the collision sequence may emit. If None,
            the unversioned filename is tried first.
        prefix_tags: Tags prepended before stem in the base filename stem.
        suffix_tags: Tags appended after stem in the base filename stem.
        end_tags: Tags appended after the suffix tags, before the version.
        preserve_unicode: Whether unicode letters and digits survive sanitization
            instead of being replaced with underscores.
    """

    directory: Path
    stem: str = ""
    extension: str = ""
    overwrite: bool = False
    max_name_bytes: int = 255
    max_path_chars: int | None = None
    source_version: int = 0
    version_floor: int | None = None
    prefix_tags: tuple[str, ...] = ()
    suffix_tags: tuple[str, ...] = ()
    end_tags: tuple[str, ...] = ()
    preserve_unicode: bool = False
    _resolved: Path | None = field(default=None, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Coerce directory to a Path and reject out-of-range numeric fields.

        Raises:
            ValueError: If max_name_bytes, max_path_chars, source_version, or
                version_floor falls outside its permitted range.
        """
        object.__setattr__(self, "directory", Path(self.directory))

        if self.max_name_bytes < 1:
            raise ValueError(f"max_name_bytes must be at least 1, got {self.max_name_bytes}.")
        if self.max_path_chars is not None and self.max_path_chars < 1:
            raise ValueError(f"max_path_chars must be at least 1, got {self.max_path_chars}.")
        if self.source_version < 0:
            raise ValueError(f"source_version cannot be negative, got {self.source_version}.")
        if self.version_floor is not None and not 1 <= self.version_floor <= _MAX_VERSION:
            raise ValueError(
                f"version_floor must fall between 1 and {_MAX_VERSION}, "
                f"got {self.version_floor}.",
            )

    @classmethod
    def for_file(
        cls,
        destination: str | Path,
        *,
        extension: str | None = None,
        overwrite: bool = False,
        max_name_bytes: int = 255,
        max_path_chars: int | None = None,
        preserve_unicode: bool = False,
        expanduser: bool = True,
        resolve: bool = False,
        create: bool = False,
    ) -> Self:
        """Construct a builder for a file path, inferring its directory and parts.

        A trailing _vN version segment is split off the stem and recorded as
        source_version instead of being carried into the output name, so a
        versioned input does not compound into report_v3_v4. Call
        continue_versioning to put it back.

        Args:
            destination: File path to derive the initial state from. Its parent
                becomes the destination directory.
            extension: File extension to use instead of the one carried by
                destination, normalized to a single leading dot. If None, the
                extension on destination is used.
            overwrite: If True, an existing file at the resolved path is
                overwritten rather than triggering version numbering.
            max_name_bytes: Maximum length of the resulting filename, in UTF-8
                bytes. Defaults to 255.
            max_path_chars: Upper bound on the length of the full resolved path.
                If None, path length goes unchecked.
            preserve_unicode: If True, unicode letters and digits survive
                sanitization instead of being replaced with underscores.
            expanduser: If True, a leading tilde in destination is expanded.
                Defaults to True.
            resolve: If True, destination is resolved to an absolute path with
                symlinks followed.
            create: If True, the destination directory is created on disk.

        Returns:
            A new PathBuilder seeded from destination.

        Raises:
            ValueError: If extension holds no usable characters, or if a numeric
                bound falls outside its permitted range.
        """
        path = cls._normalize_destination(destination, expanduser=expanduser, resolve=resolve)
        raw_stem, source_version = cls._split_version(path.stem)

        return cls._seed(
            directory=path.parent,
            stem=raw_stem,
            extension=path.suffix if extension is None else extension,
            overwrite=overwrite,
            max_name_bytes=max_name_bytes,
            max_path_chars=max_path_chars,
            preserve_unicode=preserve_unicode,
            source_version=source_version,
            create=create,
        )

    @classmethod
    def for_directory(
        cls,
        directory: str | Path,
        *,
        extension: str | None = None,
        overwrite: bool = False,
        max_name_bytes: int = 255,
        max_path_chars: int | None = None,
        preserve_unicode: bool = False,
        expanduser: bool = True,
        resolve: bool = False,
        create: bool = False,
    ) -> Self:
        """Construct a builder that writes into directory, with no stem yet.

        The directory name never leaks into the stem, so the caller supplies one
        through with_stem or relies on tags alone.

        Args:
            directory: Destination directory the built path resolves into.
            extension: Default file extension, normalized to a single leading
                dot. If None, the builder starts with no extension.
            overwrite: If True, an existing file at the resolved path is
                overwritten rather than triggering version numbering.
            max_name_bytes: Maximum length of the resulting filename, in UTF-8
                bytes. Defaults to 255.
            max_path_chars: Upper bound on the length of the full resolved path.
                If None, path length goes unchecked.
            preserve_unicode: If True, unicode letters and digits survive
                sanitization instead of being replaced with underscores.
            expanduser: If True, a leading tilde in directory is expanded.
                Defaults to True.
            resolve: If True, directory is resolved to an absolute path with
                symlinks followed.
            create: If True, the directory is created on disk.

        Returns:
            A new PathBuilder rooted at directory.

        Raises:
            ValueError: If extension holds no usable characters, or if a numeric
                bound falls outside its permitted range.
        """
        path = cls._normalize_destination(directory, expanduser=expanduser, resolve=resolve)

        return cls._seed(
            directory=path,
            stem="",
            extension="" if extension is None else extension,
            overwrite=overwrite,
            max_name_bytes=max_name_bytes,
            max_path_chars=max_path_chars,
            preserve_unicode=preserve_unicode,
            source_version=0,
            create=create,
        )

    # @classmethod
    # def from_path(
    #     cls,
    #     destination: str | Path,
    #     *,
    #     extension: str | None = None,
    #     overwrite: bool = False,
    #     max_name_bytes: int = 255,
    #     max_path_chars: int | None = None,
    #     preserve_unicode: bool = False,
    #     expanduser: bool = True,
    #     resolve: bool = False,
    #     create: bool = False,
    # ) -> Self:
    #     """Construct a builder from a path of unknown kind, probing disk first.
    #
    #     An existing directory dispatches to for_directory and anything else to
    #     for_file. Prefer the explicit constructors when the kind is known, and
    #     for_target for a user-supplied destination that should fall back to a
    #     default filename, since the guess for a path that does not exist yet
    #     rests only on the presence of an extension.
    #
    #     Args:
    #         destination: Path to derive the initial state from. It may name
    #             either a file or a directory.
    #         extension: File extension to use instead of the one carried by
    #             destination, normalized to a single leading dot. If None, the
    #             extension on destination is used.
    #         overwrite: If True, an existing file at the resolved path is
    #             overwritten rather than triggering version numbering.
    #         max_name_bytes: Maximum length of the resulting filename, in UTF-8
    #             bytes. Defaults to 255.
    #         max_path_chars: Upper bound on the length of the full resolved path.
    #             If None, path length goes unchecked.
    #         preserve_unicode: If True, unicode letters and digits survive
    #             sanitization instead of being replaced with underscores.
    #         expanduser: If True, a leading tilde in destination is expanded.
    #             Defaults to True.
    #         resolve: If True, destination is resolved to an absolute path with
    #             symlinks followed.
    #         create: If True, the destination directory is created on disk.
    #
    #     Returns:
    #         A new PathBuilder seeded from destination.
    #
    #     Raises:
    #         ValueError: If extension holds no usable characters, or if a numeric
    #             bound falls outside its permitted range.
    #     """
    #     path = cls._normalize_destination(destination, expanduser=expanduser, resolve=resolve)
    #     constructor = cls.for_directory if cls._names_directory(path) else cls.for_file
    #
    #     return constructor(
    #         path,
    #         extension=extension,
    #         overwrite=overwrite,
    #         max_name_bytes=max_name_bytes,
    #         max_path_chars=max_path_chars,
    #         preserve_unicode=preserve_unicode,
    #         expanduser=False,
    #         resolve=False,
    #         create=create,
    #     )

    @classmethod
    def for_target(
        cls,
        target: str | Path,
        *,
        default_name: str | Path | None = None,
        extension: str | None = None,
        overwrite: bool = False,
        max_name_bytes: int = 255,
        max_path_chars: int | None = None,
        preserve_unicode: bool = False,
        expanduser: bool = True,
        resolve: bool = True,
        create: bool = False,
    ) -> Self:
        """Construct a builder from a command-line target of either kind.

        A target naming a directory takes its stem and extension from
        default_name, so one argument can accept both an output directory and an
        explicit output file. A target naming a file keeps its own stem and
        extension and ignores default_name. Tags added afterwards reach the
        filename in either case, which a caller folding a tag into default_name
        cannot achieve.

        Args:
            target: Destination the user supplied, naming either a directory to
                write into or the output file itself.
            default_name: Filename to fall back on when target names a
                directory, usually the input file's name. Only its final
                component is used, and a trailing version segment on it becomes
                source_version. If None, the builder starts with no stem.
            extension: File extension to use instead of the one carried by
                target or default_name, normalized to a single leading dot. If
                None, the extension is inferred.
            overwrite: If True, an existing file at the resolved path is
                overwritten rather than triggering version numbering.
            max_name_bytes: Maximum length of the resulting filename, in UTF-8
                bytes. Defaults to 255.
            max_path_chars: Upper bound on the length of the full resolved path.
                If None, path length goes unchecked.
            preserve_unicode: If True, unicode letters and digits survive
                sanitization instead of being replaced with underscores.
            expanduser: If True, a leading tilde in target is expanded.
                Defaults to True.
            resolve: If True, target is resolved to an absolute path with
                symlinks followed. Defaults to True, since a target arrives as
                user input relative to the working directory.
            create: If True, the destination directory is created on disk.

        Returns:
            A new PathBuilder seeded from target, falling back to default_name.

        Raises:
            ValueError: If extension holds no usable characters, or if a numeric
                bound falls outside its permitted range.
        """
        path = cls._normalize_destination(target, expanduser=expanduser, resolve=resolve)

        if not cls._names_directory(path):
            return cls.for_file(
                path,
                extension=extension,
                overwrite=overwrite,
                max_name_bytes=max_name_bytes,
                max_path_chars=max_path_chars,
                preserve_unicode=preserve_unicode,
                expanduser=False,
                resolve=False,
                create=create,
            )

        fallback = Path(Path(default_name).name) if default_name else Path()
        raw_stem, source_version = cls._split_version(fallback.stem)

        return cls._seed(
            directory=path,
            stem=raw_stem,
            extension=fallback.suffix if extension is None else extension,
            overwrite=overwrite,
            max_name_bytes=max_name_bytes,
            max_path_chars=max_path_chars,
            preserve_unicode=preserve_unicode,
            source_version=source_version,
            create=create,
        )

    def with_tag(self, text: str, *, position: TagPosition = "suffix") -> Self:
        """Return a new builder with one sanitized tag added at the given position.

        Args:
            text: Tag text to sanitize and add. A tag left empty by
                sanitization is dropped rather than contributing a separator.
            position: Where the tag belongs: "prefix" before the stem, "suffix"
                after the stem, or "end" after the suffix tags and before the
                version. Defaults to "suffix".

        Returns:
            A new PathBuilder carrying the added tag.

        Raises:
            ValueError: If position is not "prefix", "suffix", or "end".
        """
        return self.with_tags(text, position=position)

    def with_tags(self, *texts: str, position: TagPosition = "suffix") -> Self:
        """Return a new builder with several sanitized tags added in order.

        Args:
            *texts: Tag texts to sanitize and add. Tags left empty by
                sanitization are dropped rather than contributing separators.
            position: Where the tags belong: "prefix" before the stem, "suffix"
                after the stem, or "end" after the suffix tags and before the
                version. Defaults to "suffix".

        Returns:
            A new PathBuilder carrying the added tags.

        Raises:
            ValueError: If position is not "prefix", "suffix", or "end".
        """
        attribute = self._tag_attribute(position)
        tags = tuple(tag for tag in (self._sanitize(text) for text in texts) if tag)
        current: tuple[str, ...] = getattr(self, attribute)

        return replace(self, **{attribute: (*current, *tags)})

    def without_tags(self, *, position: TagPosition | None = None) -> Self:
        """Return a new builder with tags cleared at one position or everywhere.

        Args:
            position: Position to clear. If None, prefix, suffix, and end tags
                are all cleared.

        Returns:
            A new PathBuilder with the selected tags removed.

        Raises:
            ValueError: If position is neither None nor one of "prefix",
                "suffix", or "end".
        """
        if position is None:
            return replace(self, prefix_tags=(), suffix_tags=(), end_tags=())
        return replace(self, **{self._tag_attribute(position): ()})

    def with_stem(self, stem: str) -> Self:
        """Return a new builder with the stem replaced by a sanitized value.

        Args:
            stem: New base filename stem.

        Returns:
            A new PathBuilder carrying the given stem.
        """
        return replace(self, stem=self._sanitize(stem))

    def with_extension(self, extension: str) -> Self:
        """Return a new builder with the file extension replaced.

        Args:
            extension: New file extension, with or without a leading dot.
                Normalized to exactly one leading dot.

        Returns:
            A new PathBuilder carrying the given extension.

        Raises:
            ValueError: If extension holds no usable characters once its leading
                dots and whitespace are removed.
        """
        return replace(self, extension=self._normalize_extension(extension))

    def with_timestamp(
        self,
        when: str | datetime | None = None,
        *,
        fmt: str = "%Y%m%d",
        tz: tzinfo | None = None,
        position: TagPosition = "prefix",
    ) -> Self:
        """Return a new builder with a timestamp tag added at the given position.

        Args:
            when: A datetime is formatted with fmt; None uses the current time
                formatted with fmt; a string is added directly, after the same
                sanitization every tag receives.
            fmt: strftime format applied to a datetime or to the current time.
                Defaults to "%Y%m%d".
            tz: Time zone for the current time. If None, local time is used.
                Ignored when when supplies its own value.
            position: Where the tag belongs: "prefix" before the stem, "suffix"
                after the stem, or "end" after the suffix tags and before the
                version. Defaults to "prefix".

        Returns:
            A new PathBuilder carrying the formatted timestamp as a tag.
        """
        moment = datetime.now(tz) if when is None else when
        text = moment if isinstance(moment, str) else moment.strftime(fmt)

        return self.with_tag(text, position=position)

    def in_directory(self, directory: str | Path) -> Self:
        """Return a new builder pointed at a different directory.

        Args:
            directory: Directory that replaces the current destination.

        Returns:
            A new PathBuilder rooted at the given directory.
        """
        return replace(self, directory=Path(directory))

    def nested_under(self, base: str | Path, *, flatten: bool = False) -> Self:
        """Return a new builder whose directory is re-rooted beneath base.

        Args:
            base: Directory to nest the current destination beneath.
            flatten: If True, only the final component of the current directory
                is kept; otherwise its full relative structure is preserved.

        Returns:
            A new PathBuilder rooted at base joined with the current directory.
        """
        subdirectory = Path(self.directory.name) if flatten else self._relative_directory()
        return replace(self, directory=Path(base) / subdirectory)

    def with_overwrite(self, overwrite: bool = True) -> Self:
        """Return a new builder with collision handling switched on or off.

        Args:
            overwrite: If True, the resolved path is returned even when a file
                already occupies it. Defaults to True.

        Returns:
            A new PathBuilder with the given overwrite behavior.
        """
        return replace(self, overwrite=overwrite)

    def with_max_name_bytes(self, max_name_bytes: int) -> Self:
        """Return a new builder with a different filename byte budget.

        Args:
            max_name_bytes: Maximum length of the resulting filename, in UTF-8 bytes.

        Returns:
            A new PathBuilder with the given filename budget.

        Raises:
            ValueError: If max_name_bytes is less than 1.
        """
        return replace(self, max_name_bytes=max_name_bytes)

    def with_version(self, version: int | None) -> Self:
        """Return a new builder whose version sequence starts at version.

        Args:
            version: Lowest version the collision sequence may emit. If None,
                the unversioned filename is tried first.

        Returns:
            A new PathBuilder with the given version floor.

        Raises:
            ValueError: If version falls outside the supported version range.
        """
        return replace(self, version_floor=version)

    def continue_versioning(self) -> Self:
        """Return a new builder resuming from the version parsed off the input.

        Returns:
            A new PathBuilder whose version floor is source_version, or an
            unchanged copy when no version was parsed.
        """
        return replace(self, version_floor=self.source_version or None)

    def ensure_directory(self) -> Self:
        """Create the builder's directory on disk, including missing parents.

        Returns:
            This same PathBuilder instance, for chaining.
        """
        self.directory.mkdir(parents=True, exist_ok=True)
        return self

    def sibling(self, name: str, *, create: bool = False) -> Path:
        """Return the path of another file in this builder's directory.

        Stem, tags, version, and extension play no part: only the sanitized bare
        filename from name is joined onto the directory, so any directory
        components in name are discarded.

        Args:
            name: Name of the neighboring file. Only its final component is used.
            create: If True, the directory is created on disk before the path is
                returned.

        Returns:
            The directory joined with the sanitized bare filename from name.
        """
        if create:
            self.ensure_directory()
        return self.directory / self._sanitize(Path(name).name)

    def resolve(self) -> Path:
        """Return the path a write would land on, reading the directory once.

        The answer is cached on this instance, so repeated calls stay consistent
        even as the directory changes underneath. Nothing is created, so two
        processes resolving concurrently can pick the same name; use reserve
        when that matters.

        Returns:
            The full path for this state, version-numbered unless overwrite is
            enabled or no file occupies the unversioned name.

        Raises:
            ValueError: If the filename would be empty, if the extension and
                tags leave no room for a stem, or if the resolved path exceeds
                max_path_chars.
        """
        if self._resolved is not None:
            return self._resolved

        stem = self._require_stem()
        version = self.version_floor if self.overwrite else self._select_version(stem)

        return self._remember(self.directory / self._name_for(stem, version))

    def reserve(self) -> Path:
        """Claim an unused path atomically and return it, creating the directory.

        Candidates are tried in version order and each is created with an
        exclusive flag, so a name is only returned once it belongs to this
        caller. That makes concurrent writers in one directory safe, at the cost
        of leaving an empty file behind if the caller never writes to it.
        Resolution is cached, so the path property afterwards reports the
        reserved path.

        Returns:
            The reserved path, already present on disk as an empty file, or the
            resolved path when overwrite is enabled.

        Raises:
            ValueError: If the filename would be empty, if the extension and
                tags leave no room for a stem, if the resolved path exceeds
                max_path_chars, or if every version up to the internal ceiling
                is taken.
        """
        if self.overwrite:
            return self.ensure_directory().resolve()

        stem = self._require_stem()
        self.ensure_directory()

        for version in self._version_candidates(self._select_version(stem)):
            candidate = self.directory / self._name_for(stem, version)
            self._check_path_length(candidate)
            try:
                candidate.touch(exist_ok=False)
            except FileExistsError:
                continue
            return self._remember(candidate, force=True)

        raise ValueError(
            f"Every version up to {_MAX_VERSION} is taken for {stem!r} in {self.directory}.",
        )

    @property
    def filename(self) -> str:
        """The resolved filename for this state, version-numbered to avoid collisions."""
        return self.resolve().name

    @property
    def path(self) -> Path:
        """The resolved full path for this state, cached after the first read."""
        return self.resolve()

    def __repr__(self) -> str:
        """Return a string showing the directory, stem, extension, and versions."""
        return (
            f"{self.__class__.__name__}("
            f"directory={self.directory!r}, "
            f"stem={self.stem!r}, "
            f"extension={self.extension!r}, "
            f"source_version={self.source_version!r}, "
            f"version_floor={self.version_floor!r})"
        )

    def __str__(self) -> str:
        """Return the string form of the resolved path."""
        return str(self.resolve())

    def __fspath__(self) -> str:
        """Return the string form of the resolved path, per os.PathLike."""
        return str(self.resolve())

    @staticmethod
    def _normalize_destination(
        destination: str | Path,
        *,
        expanduser: bool,
        resolve: bool,
    ) -> Path:
        """Apply tilde expansion and absolute resolution to a constructor argument.

        Args:
            destination: Path supplied to a constructor.
            expanduser: If True, a leading tilde is expanded.
            resolve: If True, the result is made absolute with symlinks followed.

        Returns:
            The normalized path.
        """
        path = Path(destination)
        if expanduser:
            path = path.expanduser()
        if resolve:
            path = path.resolve()
        return path

    @staticmethod
    def _names_directory(path: Path) -> bool:
        """Judge whether a path names a directory rather than a file.

        What is already on disk decides it when the path exists. For a path that
        does not exist yet, only the absence of an extension and of a leading dot
        are available as signals, which is why the explicit constructors are
        preferred whenever the kind is known.

        Args:
            path: Normalized path to classify.

        Returns:
            True if path names a directory.
        """
        if path.exists():
            return path.is_dir()
        return not (path.suffix or path.name.startswith("."))

    @classmethod
    def _seed(
        cls,
        *,
        directory: Path,
        stem: str,
        extension: str,
        overwrite: bool,
        max_name_bytes: int,
        max_path_chars: int | None,
        preserve_unicode: bool,
        source_version: int,
        create: bool,
    ) -> Self:
        """Build the initial instance shared by every public constructor.

        Args:
            directory: Destination directory the built path resolves into.
            stem: Base filename stem, sanitized here so every entry point agrees.
            extension: File extension, normalized here to a single leading dot.
            overwrite: If True, an existing file at the resolved path is
                overwritten rather than triggering version numbering.
            max_name_bytes: Maximum length of the resulting filename, in UTF-8 bytes.
            max_path_chars: Upper bound on the length of the full resolved path.
            preserve_unicode: If True, unicode letters and digits survive sanitization.
            source_version: Version parsed off a trailing _vN stem.
            create: If True, the destination directory is created on disk.

        Returns:
            A new PathBuilder with normalized parts.
        """
        builder = cls(
            directory=directory,
            stem=sanitize_component(stem, allow_unicode=preserve_unicode),
            extension=cls._normalize_extension(extension),
            overwrite=overwrite,
            max_name_bytes=max_name_bytes,
            max_path_chars=max_path_chars,
            source_version=source_version,
            preserve_unicode=preserve_unicode,
        )

        return builder.ensure_directory() if create else builder

    @staticmethod
    def _tag_attribute(position: TagPosition) -> str:
        """Return the field name storing tags for a position.

        Args:
            position: Tag position to look up.

        Returns:
            Name of the field holding tags for that position.

        Raises:
            ValueError: If position is not "prefix", "suffix", or "end".
        """
        if (attribute := _TAG_FIELDS.get(position)) is None:
            raise ValueError(
                f"Unsupported tag position: {position!r}. "
                f"Expected one of {', '.join(map(repr, _TAG_FIELDS))}.",
            )
        return attribute

    @staticmethod
    def _normalize_extension(extension: str | None) -> str:
        """Normalize an extension to a single leading dot and safe characters.

        Args:
            extension: Extension to normalize, with or without a leading dot. If
                None or empty, an empty string is returned.

        Returns:
            The normalized extension, or "" if extension is empty.

        Raises:
            ValueError: If extension holds no usable characters once its leading
                dots and whitespace are removed.
        """
        if not extension:
            return ""

        cleaned = sanitize_component(extension.strip().lstrip("."))
        if not cleaned:
            raise ValueError(f"Extension {extension!r} contains no usable characters.")

        return f".{cleaned}"

    @staticmethod
    def _split_version(stem: str) -> tuple[str, int]:
        """Split a trailing _vN version segment off a filename stem.

        Args:
            stem: Filename stem to split.

        Returns:
            A (base_stem, version) tuple, where version is 0 when stem carries
            no version segment.
        """
        if match := _VERSION_RE.match(stem):
            return match.group(1), int(match.group(2))
        return stem, 0

    def _sanitize(self, text: str) -> str:
        """Sanitize text under this builder's unicode policy.

        Args:
            text: Text to reduce to a safe filename component.

        Returns:
            The sanitized component, or "" if no safe characters remain.
        """
        return sanitize_component(text, allow_unicode=self.preserve_unicode)

    def _base_stem(self) -> str:
        """Join prefix tags, the stem, and suffix tags into the base stem.

        Returns:
            The underscore-joined base stem, with empty parts omitted.
        """
        parts = (*self.prefix_tags, self.stem, *self.suffix_tags)
        return "_".join(part for part in parts if part)

    def _end_segment(self) -> str:
        """Join end tags into a single underscore-prefixed segment.

        Returns:
            The underscore-joined end tags prefixed with "_", or "" when there are none.
        """
        return f"_{'_'.join(self.end_tags)}" if self.end_tags else ""

    def _require_stem(self) -> str:
        """Return the on-disk stem every version candidate shares.

        The stem is truncated against a budget that always reserves room for a
        version segment, so the same characters survive whether or not a version
        is appended. Version scanning therefore recognizes the files that
        earlier writes actually produced.

        Returns:
            The truncated base stem with the end segment appended.

        Raises:
            ValueError: If the resulting filename would be empty, or if the
                extension and tags leave no room for a stem.
        """
        end_segment = self._end_segment()
        stem = f"{self._truncate_base(self._base_stem(), end_segment)}{end_segment}"

        if stem.upper() in _WINDOWS_RESERVED:
            # Windows reserves these device names in every directory, extension
            # or not, so a trailing underscore keeps the name usable.
            stem = f"{stem}_"
        if not stem and not self.extension:
            raise ValueError("Cannot resolve a path: stem, tags, and extension are all empty.")

        return stem

    def _truncate_base(self, base_stem: str, end_segment: str) -> str:
        """Shorten a base stem to fit the filename byte budget.

        Filesystem limits are measured in bytes rather than characters, so the
        stem is cut in UTF-8 byte space to avoid splitting a multi-byte
        character. The budget subtracts a fixed reserve for the version segment
        whether or not one is in use, which keeps the truncated stem stable
        across versions.

        Args:
            base_stem: Stem portion of the filename, before the end segment.
            end_segment: End-tag segment that follows the base stem.

        Returns:
            base_stem unchanged when it fits, otherwise its byte-safe prefix.

        Raises:
            ValueError: If the extension, end segment, and version reserve leave
                no room for a stem, or if byte-safe truncation empties the stem.
        """
        fixed = len(f"{end_segment}{self.extension}".encode()) + _VERSION_RESERVE_BYTES
        available = self.max_name_bytes - fixed

        if available < 1:
            raise ValueError(
                f"The extension, end tags, and version reserve occupy {fixed} bytes, "
                f"which exceeds max_name_bytes={self.max_name_bytes}. Shorten the tags "
                f"or raise max_name_bytes.",
            )

        encoded = base_stem.encode()
        if len(encoded) <= available:
            return base_stem

        truncated = encoded[:available].decode("utf-8", "ignore")
        if not truncated:
            # Every byte in the window fell inside one multi-byte character, so
            # nothing decodable remains and a stem-less name would be wrong.
            raise ValueError(
                f"Byte-safe truncation of {base_stem!r} to {available} bytes produced an "
                f"empty stem. Raise max_name_bytes or shorten the stem.",
            )

        logger.warning("Filename too long. Truncated stem to {!r}.", truncated)
        return truncated

    def _name_for(self, stem: str, version: int | None) -> str:
        """Assemble a filename from a resolved stem and a version.

        Args:
            stem: On-disk stem shared by every version candidate.
            version: Version to append. If None, no version segment is added.

        Returns:
            The assembled filename.
        """
        version_segment = "" if version is None else f"_v{version}"
        return f"{stem}{version_segment}{self.extension}"

    def _version_candidates(self, start: int | None) -> Iterator[int | None]:
        """Yield versions to try, in the order a reservation should attempt them.

        Args:
            start: First version to try. If None or 1, the unversioned filename
                is yielded first, since it counts as version 1.

        Yields:
            None for the unversioned filename, then ascending version numbers up
            to the internal ceiling.
        """
        if start is None or start <= 1:
            yield None
        yield from range(max(start or 2, 2), _MAX_VERSION + 1)

    def _select_version(self, stem: str) -> int | None:
        """Choose the version a write should use, scanning only on collision.

        Args:
            stem: On-disk stem shared by every version candidate.

        Returns:
            The version floor when its filename is free, otherwise the next
            version after the highest already present.
        """
        floor = self.version_floor
        if not (self.directory / self._name_for(stem, floor)).exists():
            return floor

        scanned = self._scan_for_next_version(stem)

        return scanned if floor is None else max(scanned, floor + 1)

    def _scan_for_next_version(self, stem: str) -> int:
        """Scan the directory for existing versions of a stem.

        The unversioned filename counts as version 1, so the lowest value this
        can return is 2 and a _v1 segment is never produced.

        Args:
            stem: On-disk stem shared by every version candidate.

        Returns:
            One past the highest version present, or 1 when the directory cannot
            be read.
        """
        pattern = re.compile(
            rf"^{re.escape(stem)}(?:_v(\d+))?{re.escape(self.extension)}$",
        )
        highest = 0

        try:
            with os.scandir(self.directory) as entries:
                for entry in entries:
                    # The cheap name test runs first: an empty extension makes
                    # endswith always true, so it is skipped in that case.
                    if self.extension and not entry.name.endswith(self.extension):
                        continue
                    if (match := pattern.match(entry.name)) and entry.is_file():
                        highest = max(highest, int(match.group(1) or 1))
        except OSError as error:
            logger.error("Cannot scan directory {}: {}", self.directory, error)
            return 1

        return highest + 1

    def _relative_directory(self) -> Path:
        """Return the directory with any root, drive, or share anchor removed.

        The anchor decides this rather than is_absolute, which on Windows is
        False for a rooted path that carries no drive letter. Such a path would
        otherwise survive intact and, once joined onto a base, replace
        everything after that base's drive.

        Returns:
            directory unchanged when it carries no anchor, otherwise the same
            path with its anchor stripped.
        """
        anchor = self.directory.anchor
        return self.directory.relative_to(anchor) if anchor else self.directory

    def _check_path_length(self, path: Path) -> None:
        """Reject a path longer than the configured character budget.

        Args:
            path: Candidate path to measure.

        Raises:
            ValueError: If max_path_chars is set and path exceeds it.
        """
        if self.max_path_chars is None:
            return

        length = len(str(path))
        if length > self.max_path_chars:
            raise ValueError(
                f"Resolved path is {length} characters, which exceeds "
                f"max_path_chars={self.max_path_chars}: {path}",
            )

    def _remember(self, path: Path, *, force: bool = False) -> Path:
        """Validate a resolved path, cache it on this instance, and return it.

        Args:
            path: Resolved path to cache.
            force: If True, an existing cached path is replaced.

        Returns:
            The cached path.

        Raises:
            ValueError: If path exceeds max_path_chars.
        """
        self._check_path_length(path)
        if force or self._resolved is None:
            object.__setattr__(self, "_resolved", path)
        return path
