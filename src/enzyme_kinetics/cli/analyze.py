"""CLI entrypoint for the enzyme kinetics analysis pipeline."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import tyro
from loguru import logger

from enzyme_kinetics.analysis.analyze import KineticAnalyzer
from enzyme_kinetics.analysis.arguments import KineticArgs
from enzyme_kinetics.calibration.calibrate import fit_calibration
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


def apply_calibration(
    path: Path,
    peak_prefixes: frozenset[str] | None,
    analyzer: KineticAnalyzer,
) -> KineticAnalyzer:
    """Applies a fitted Calibration object to the KineticAnalyzer."""
    # 1. Load data
    df = read_clean_csv(path, options=CleaningOptions(snake_case=True, drop_indexlike=True))
    required_cols = {"substrate_conc", "peak_id", "mean", "std", "count"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    # 2. Canonicalize peak IDs
    df = canonicalize_peak_ids(df, prefixes=peak_prefixes)

    # 3. Sort data (defensive guard)
    df.sort_values(by=["peak_id", "substrate_conc"], inplace=True)
    logger.info("Loaded %d peaks from %s", len(df), path)

    for peak_id in df["peak_id"].unique():
        # if peak_id.lower() == "olv":
        #     continue
        sub_df = df[df["peak_id"] == peak_id]
        # 4. Fit calibration model
        model = fit_calibration(sub_df["substrate_conc"], sub_df["mean"], sigma=sub_df["std"])
        analyzer.apply_calibration(model, peak_ids=[peak_id])
        logger.info("Applied calibration to %s", peak_id)

    return analyzer


def main(args: KineticArgs) -> None:
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
    df = read_clean_csv(args.path, options=CleaningOptions(snake_case=True, drop_indexlike=True))
    required_cols = {"substrate_conc", "peak_id", "mean", "std", "count"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    # 2. Canonicalize peak IDs
    df = canonicalize_peak_ids(df, prefixes=args.peak_prefixes)

    # 3. Sort data (defensive guard)
    df.sort_values(by=["peak_id", "substrate_conc"], inplace=True)
    logger.info("Loaded %d peaks from %s", len(df), args.path)

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
        logger.warning(
            "No calibration path provided. kcat and kcat/Km carry non-physical "
            "units until a Calibration is applied.",
        )

    # 7. Plot and save results
    base_dest_path = args.target_dir / args.path.name
    analyzer.export_csv(base_dest_path, overwrite=args.overwrite)
    analyzer.plot(base_dest_path, is_calibrated=is_calibrated, overwrite=args.overwrite)

    logger.info("Analysis complete.")

    logger.info(
        "Assuming input concentrations were in units µM, "
        "kinetic constants are in:"
        "\n **µM** (Km)"
        "\n **s⁻¹** (kcat)"
        "\n **s⁻¹·M⁻¹** (kcat/Km)",
    )


if __name__ == "__main__":
    main(tyro.cli(KineticArgs))
