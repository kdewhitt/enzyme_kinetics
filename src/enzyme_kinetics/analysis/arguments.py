"""Validated command-line argument model for the enzyme kinetics analysis pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import tyro
from pydantic import BaseModel, computed_field, ConfigDict, Field, field_validator

from enzyme_kinetics.core.compounds import FlexPeakPrefixes


class KineticArgs(BaseModel):
    """CLI configuration model for the enzyme kinetics analysis pipeline.

    Accepts a pre-aggregated CSV produced by **PeakAnalyzer** (the
    "_[area|height]_analysis.csv" output) containing one row per
    peak_id × substrate_conc combination with pre-computed replicate
    statistics. Column names are normalized to snake_case before
    validation.

    Required columns:
        substrate_conc: Substrate concentrations in µM. All kinetic
            parameters (Km, kcat, kcat/Km) assume µM as the concentration
            unit; passing values in any other unit will silently produce
            incorrect results.
        peak_id: HPLC peak identifier string.
        mean: Mean integrated peak area across replicates.
        std: Standard deviation of peak area across replicates.
        count: Number of replicates per group.

    Attributes:
        path: Path to the pre-aggregated PeakAnalyzer CSV file to process.
        target_dir: Destination directory where output CSV and plot files
            are written.
        calibration_path: Optional path to a calibration standards CSV in the
            same PeakAnalyzer format as path (not a saved Calibration JSON).
            If provided, a calibration curve is fitted per peak and applied
            after fit() to convert velocities to µM/s and produce
            physically meaningful kcat (s⁻¹) and kcat/Km (M⁻¹·s⁻¹) values.
            If None, kcat and kcat/Km are computed from raw area velocities
            and carry non-physical units. Defaults to None.
        enzyme_conc_um: Total enzyme concentration in µM used to compute
            kcat = Vmax / [E]. Defaults to 12.25.
        rxn_time: Reaction duration in minutes. Converted to seconds
            internally via the reaction_time_seconds computed property.
            Defaults to 180.0 (3 hours).
        substrate: Substrate name string propagated to KineticConstants
            results and plot axis labels (e.g. "HexCoA"). Defaults to
            "HexCoA".
        peak_prefixes: Optional FlexPeakPrefixes configuration controlling
            which peak ID prefixes are recognized during canonicalization.
            If None, default canonicalization rules apply. Defaults to None.
        is_calibrated: Whether the input mean integrated peak area values have already
            had calibration corrections applied. Defaults to False.
        overwrite: Whether to overwrite existing output files. If False,
            PathBuilder appends a unique suffix to avoid collisions.
            Defaults to False.
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
    target_dir: Annotated[Path, tyro.conf.Positional]
    calibration_path: Path | None = None

    enzyme_conc_um: float = 12.25
    rxn_time: float = 180.0
    substrate: str = "HexCoA"

    peak_prefixes: FlexPeakPrefixes = None

    is_calibrated: bool = False

    overwrite: bool = False
    verbose: bool = Field(default=False, exclude=True)

    @field_validator("path", "target_dir", mode="before")
    @classmethod
    def _resolve_paths(cls, value) -> Path:
        """Resolves paths to absolute paths and expands user-home."""
        return Path(value).expanduser().resolve()

    @field_validator("calibration_path", mode="before")
    @classmethod
    def _resolve_calibration_path(cls, value) -> Path | None:
        """Resolves calibration path to absolute path if provided."""
        return Path(value).expanduser().resolve() if value is not None else None

    @computed_field
    @property
    def reaction_time_seconds(self) -> float:
        """Reaction duration in seconds."""
        return self.rxn_time * 60
