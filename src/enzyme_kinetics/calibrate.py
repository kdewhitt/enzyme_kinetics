"""CLI entrypoint for the calibration pipeline."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Annotated

import tyro
from logurich import DuoLogAdapter
from pydantic import BaseModel, ConfigDict, Field, field_validator
from rich.console import Console

from .core.calibration import fit_calibration, plot_calibration, save_calibration
from .core.compounds import canonicalize_peak_ids, FlexPeakPrefixes
from .core.io import load_plottable_data

_logger = DuoLogAdapter.create(component=__name__)

# Environment setup for Rich
os.environ.setdefault("FORCE_COLOR", "1")
os.environ.setdefault("TERM", "xterm-256color")
console = Console(force_terminal=True, color_system="truecolor")


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
            to the stem for each peak_id, producing one .txt report and one
            .png diagnostic plot per peak.
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


def run_calibration_pipeline(args: CalibrationArgs) -> None:
    """Executes the per-peak HPLC calibration pipeline from CLI arguments.

    Loads and validates the input DataFrame, canonicalizes peak identifiers,
    and for each unique peak_id fits a linear calibration curve to the
    substrate_conc (µM) vs. mean peak area data. Two output files are produced
    per peak: a plain-text parameter report (or JSON if pretty=False) and a
    three-panel diagnostic plot (.png).

    Args:
        args: Fully validated CalibrationArgs instance produced by tyro.cli.

    Raises:
        KeyError: If the input DataFrame is missing any of the required
            columns: "substrate_conc", "peak_id", "mean", "std", "count".

    Note:
        substrate_conc values must be in µM. Passing values in any other
        unit will silently produce an incorrectly scaled calibration curve
        and corrupt all downstream kcat and kcat/Km calculations.
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

    for peak_id in df["peak_id"].unique():
        sub_df = df[df["peak_id"] == peak_id]
        # 4. Fit calibration model
        model = fit_calibration(sub_df["substrate_conc"], sub_df["mean"], sigma=sub_df["std"])

        out_path = args.outfile / args.outfile.with_name(f"{args.outfile.stem}_{peak_id}.txt")
        # out_path = PathBuilder(args.outfile).with_tag(peak_id).path

        # 7. Plot and save results
        save_calibration(model, peak_id, out_path.with_suffix(".txt"), nice=args.pretty)
        plot_calibration(
            model,
            sub_df["substrate_conc"],
            sub_df["mean"],
            output_path=out_path.with_suffix(".png"),
        )

    _logger.info("Analysis complete.")


def main() -> None:
    """Parses CLI arguments and runs the per-peak HPLC calibration pipeline."""
    args = tyro.cli(CalibrationArgs)
    run_calibration_pipeline(args)


if __name__ == "__main__":
    main()
