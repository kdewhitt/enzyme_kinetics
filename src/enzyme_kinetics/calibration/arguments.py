"""Validated command-line argument model for the per-peak HPLC calibration pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import tyro
from pydantic import BaseModel, ConfigDict, Field, field_validator

from enzyme_kinetics.core.compounds import FlexPeakPrefixes

__all__ = ["CalibrationArgs"]


class CalibrationArgs(BaseModel):
    """CLI configuration model for the per-peak HPLC calibration pipeline.

    Accepts a pre-aggregated CSV produced by PeakAnalyzer (the
    "_[area|height]_analysis.csv" output) containing one row per
    peak_id × substrate_conc combination with pre-computed replicate
    statistics. Column names are normalized to snake_case before
    validation.

    In the calibration context, substrate_conc holds the known CoA
    standard concentrations (µM) used as the x-axis of the calibration
    curve, and mean holds the corresponding measured peak areas. This is
    the inverse of the roles these columns play in the enzyme kinetics
    pipeline. A separate calibration curve is fitted and exported for
    each unique peak_id in the input file.

    Required columns:
        substrate_conc: Known CoA concentrations of calibration standards
            in µM. Values must be in µM; passing values in any other unit
            will silently produce an incorrectly scaled calibration curve.
        peak_id: HPLC peak identifier string.
        mean: Mean integrated peak area across replicates.
        std: Standard deviation of peak area across replicates.
        count: Number of replicates per group.

    Attributes:
        path: Path to the pre-aggregated PeakAnalyzer CSV file to process.
        outfile: Base output file path. A peak-specific suffix is appended
            to the stem for each peak_id, producing one report (.txt when
            pretty, .json otherwise) and one .png diagnostic plot per peak.
        peak_prefixes: Optional FlexPeakPrefixes configuration controlling
            which peak ID prefixes are recognized during canonicalization.
            If None, default canonicalization rules apply. Defaults to None.
        pretty: Whether to write calibration parameters as a human-readable
            plain-text report (nice=True) rather than compact JSON.
            Defaults to True.
        overwrite: Whether to overwrite existing output files. Defaults to
            False.
        verbose: Whether to enable verbose logging output. Excluded from
            model serialization. Defaults to False.
    """

    model_config = ConfigDict(
        extra="ignore",
        frozen=True,
        arbitrary_types_allowed=True,
        str_strip_whitespace=True,
        validate_default=True,
    )

    path: Annotated[Path, tyro.conf.Positional]
    outfile: Annotated[Path, tyro.conf.Positional]

    peak_prefixes: FlexPeakPrefixes = None

    pretty: bool = True

    overwrite: bool = False
    verbose: bool = Field(default=False, exclude=True)

    @field_validator("path", "outfile", mode="before")
    @classmethod
    def _resolve_paths(cls, value) -> Path:
        """Resolves paths to absolute paths and expands user-home."""
        return Path(value).expanduser().resolve()
