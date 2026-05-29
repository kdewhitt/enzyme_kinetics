"""CLI entrypoint for the enzyme kinetics analysis pipeline."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Annotated

import tyro
from logurich import configure_richloguru, RichLogAdapter
from pydantic import BaseModel, computed_field, ConfigDict, Field, field_validator
from rich.console import Console

from enzyme_kinetics.core import (
    canonicalize_peak_ids,
    FlexPeakPrefixes,
    KineticAnalyzer,
    load_plottable_data,
)
from enzyme_kinetics.core.calibration import fit_calibration

configure_richloguru(level="INFO")
_logger = RichLogAdapter(component=__name__)

# Environment setup for Rich
os.environ.setdefault("FORCE_COLOR", "1")
os.environ.setdefault("TERM", "xterm-256color")
console = Console(force_terminal=True, color_system="truecolor")


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
        calibration_path: Optional path to a fitted Calibration JSON file
            produced by fit_calibration(). If provided, apply_calibration()
            is called after fit() to convert velocities to µM/s and produce
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
        return self.rxn_time * 60


def run_enzyme_kinetic_analysis_pipeline(args: KineticArgs) -> None:
    """Executes the full enzyme kinetics analysis pipeline from CLI arguments.

    Loads and validates the input DataFrame, canonicalizes peak identifiers,
    fits Michaelis-Menten models to each peak, and exports results as a
    tagged CSV and a set of kinetic plots. Steps are logged at INFO level.

    kcat and kcat/Km values in the output carry non-physical units unless
    a calibration curve is applied upstream. The calibration step is
    currently disabled; enable analyzer.apply_calibration() with a fitted
    Calibration object before calling analyzer.fit() to produce physically
    meaningful s⁻¹ and s⁻¹·M⁻¹ values.

    Args:
        args: Fully validated KineticArgs instance produced by tyro.cli.

    Raises:
        KeyError: If the input DataFrame is missing any of the required
            columns: "substrate_conc", "peak_id", "mean", "std", "count".

    Note:
        substrate_conc values in the input file must be in µM. Passing
        values in any other unit will silently produce incorrect Km,
        kcat, and kcat/Km results.
    """
    # 1. Load data
    df = load_plottable_data(args.path, sanitize=True, drop_indexlike=True)
    required_cols = {"substrate_conc", "peak_id", "mean", "std", "count"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    # 2. Canonicalize peak IDs
    df = canonicalize_peak_ids(df, prefixes=args.peak_prefixes)

    # 3. Sort data (defensive guard)
    df.sort_values(by=["peak_id", "substrate_conc"], inplace=True)
    _logger.info("Loaded %d peaks from %s", len(df), args.path)

    # 4. Initialize pipeline
    analyzer = KineticAnalyzer(
        args.enzyme_conc_um,
        args.reaction_time_seconds,
        args.substrate,
        data=df,
        special_peaks=None,
    )

    # 5. Fit models to raw area/s velocities
    analyzer.fit()

    is_calibrated = args.is_calibrated

    # 6. Apply calibration if a calibration file was provided
    if args.calibration_path is not None:
        apply_calibration(args.calibration_path, args.peak_prefixes, analyzer)
        is_calibrated = True
    elif not args.is_calibrated:
        _logger.warning(
            "No calibration path provided. kcat and kcat/Km carry non-physical "
            "units until a Calibration is applied.",
        )

    # 7. Plot and save results
    base_dest_path = args.target_dir / args.path.name
    analyzer.export_csv(base_dest_path, overwrite=args.overwrite)
    analyzer.plot(base_dest_path, is_calibrated=is_calibrated, overwrite=args.overwrite)

    _logger.info("Analysis complete.")

    _logger.info(
        "Assuming input concentrations were in units µM, "
        "kinetic constants are in:"
        "\n **µM** (Km)"
        "\n **s⁻¹** (kcat)"
        "\n **s⁻¹·M⁻¹** (kcat/Km)",
    )


def apply_calibration(path: Path, peak_prefixes: frozenset[str] | None, analyzer: KineticAnalyzer) -> KineticAnalyzer:
    """Applies a fitted Calibration object to the KineticAnalyzer."""
    # 1. Load data
    df = load_plottable_data(path, sanitize=True, drop_indexlike=True)
    required_cols = {"substrate_conc", "peak_id", "mean", "std", "count"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    # 2. Canonicalize peak IDs
    df = canonicalize_peak_ids(df, prefixes=peak_prefixes)

    # 3. Sort data (defensive guard)
    df.sort_values(by=["peak_id", "substrate_conc"], inplace=True)
    _logger.info("Loaded %d peaks from %s", len(df), path)

    for peak_id in df["peak_id"].unique():
        if peak_id.lower() == "olv":
            continue
        sub_df = df[df["peak_id"] == peak_id]
        # 4. Fit calibration model
        model = fit_calibration(sub_df["substrate_conc"], sub_df["mean"], sigma=sub_df["std"])
        analyzer.apply_calibration(model, peak_ids=[peak_id])
        _logger.info("Applied calibration to %s", peak_id)

    return analyzer


def main() -> None:
    """Parses CLI arguments and runs the enzyme kinetics analysis pipeline."""
    args = tyro.cli(KineticArgs)
    run_enzyme_kinetic_analysis_pipeline(args)


if __name__ == "__main__":
    main()
