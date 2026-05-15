from __future__ import annotations

import os
from pathlib import Path

import tyro
from logurich import configure_richloguru, RichLogAdapter
from pydantic import BaseModel, ConfigDict, Field, field_validator
from rich.console import Console
from tyro.conf import Positional

from enzyme_kinetics.core import (
    canonicalize_peak_ids,
    FlexPeakPrefixes,
    load_plottable_data,
)
from enzyme_kinetics.core.calibration import fit_calibration, plot_calibration, save_calibration

configure_richloguru(level="INFO")
_logger = RichLogAdapter(component=__name__)

# TODO: Update docstrings


# Environment setup for Rich
os.environ.setdefault("FORCE_COLOR", "1")
os.environ.setdefault("TERM", "xterm-256color")

console = Console(force_terminal=True, color_system="truecolor")


class CalibrationArgs(BaseModel):
    model_config = ConfigDict(extra="ignore", frozen=True, arbitrary_types_allowed=True)

    path: Positional[Path]
    outfile: Positional[Path]

    peak_prefixes: FlexPeakPrefixes = None

    pretty: bool = True

    overwrite: bool = False
    verbose: bool = Field(default=False, exclude=True)

    @field_validator("path", "outfile", mode="before")
    @classmethod
    def _resolve_paths(cls, value) -> Path:
        """Resolve paths to absolute paths and expand user-home."""
        return Path(value).expanduser().resolve()


def run_calibration_pipeline(args: CalibrationArgs) -> None:
    """Run the enzyme kinetics calibration pipeline with provided arguments."""
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
        plot_calibration(model, sub_df["substrate_conc"], sub_df["mean"], output_path=out_path.with_suffix(".png"))

    _logger.info("Analysis complete.")


def main() -> None:
    """Parses CLI arguments and runs the enzyme kinetics analysis pipeline."""
    args = tyro.cli(CalibrationArgs)
    run_calibration_pipeline(args)


if __name__ == "__main__":
    main()
