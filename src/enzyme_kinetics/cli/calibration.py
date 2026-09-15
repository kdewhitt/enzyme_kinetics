"""CLI entrypoint for the calibration pipeline."""

from __future__ import annotations

import os
import sys

import tyro
from loguru import logger

from enzyme_kinetics.calibration.arguments import CalibrationArgs
from enzyme_kinetics.calibration.calibrate import (
    fit_calibration,
    plot_calibration,
    save_calibration,
)
from enzyme_kinetics.core.compounds import canonicalize_peak_ids
from enzyme_kinetics.forks import CleaningOptions, read_clean_csv

# Maximize system reliability and output consistency
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

# Configure stdout and stderr encoding and handling
if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# Environment setup for tyro
os.environ.setdefault("FORCE_COLOR", "1")
os.environ.setdefault("TERM", "xterm-256color")


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
    df = read_clean_csv(args.path, options=CleaningOptions(snake_case=True, drop_indexlike=True))
    required_cols = {"substrate_conc", "peak_id", "mean", "std", "count"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    # 2. Canonicalize peak IDs
    df = canonicalize_peak_ids(df, prefixes=args.peak_prefixes)

    # 3. Sort data (defensive guard)
    df.sort_values(by=["peak_id", "substrate_conc"], inplace=True)
    logger.info("Loaded {} peaks from {}", len(df), args.path)

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

    logger.success("Analysis complete.")


def main():
    """Fit per-peak HPLC calibration curves and write a report and plot for each peak.

    Reads a PeakAnalyzer CSV of calibration standards from PATH and writes a
    parameter report and diagnostic plot per peak alongside OUTFILE. Pass
    --no-pretty to write the parameters as JSON instead of a text report.
    """
    args = tyro.cli(
        CalibrationArgs,
        description="Fit calibration curves to enzyme kinetics data.",
        compact_help=True,
    )
    run_calibration_pipeline(args)


if __name__ == "__main__":
    main()
