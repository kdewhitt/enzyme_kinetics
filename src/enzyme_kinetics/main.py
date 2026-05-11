import os
from pathlib import Path

import tyro
from logurich import configure_richloguru, RichLogAdapter
from pydantic import BaseModel, computed_field, ConfigDict, Field, field_validator
from rich.console import Console
from tyro.conf import Positional

from enzyme_kinetics.compounds.categoricals import canonicalize_peak_ids
from enzyme_kinetics.compounds.validators import FlexPeakPrefixes
from enzyme_kinetics.core.analyze import KineticAnalyzer
from enzyme_kinetics.core.io import load_plottable_data

configure_richloguru(level="INFO")
_logger = RichLogAdapter(component=__name__)

# Environment setup for Rich
os.environ.setdefault("FORCE_COLOR", "1")
os.environ.setdefault("TERM", "xterm-256color")

console = Console(force_terminal=True, color_system="truecolor")


class KineticArgs(BaseModel):
    """Orchestrates batch processing of Shimadzu HPLC exported data files.

    Attributes:
        path_dir: Path to the directory containing Shimadzu text files to
            process.
        target_dir: Destination directory where exported CSV files are
            written.
        config: Pipeline configuration controlling scope, cleaning, rescaling,
            and export behavior.
        debug: Whether to enable debug-level logging output.
    """
    model_config = ConfigDict(extra="ignore", frozen=True, arbitrary_types_allowed=True)

    path: Positional[Path]
    target_dir: Positional[Path]

    enzyme_conc_um: float = 12.25  # in units µM
    rxn_time: float = 180.0  # in units seconds?

    peak_prefixes: FlexPeakPrefixes = None

    overwrite: bool = False
    verbose: bool = Field(default=False, exclude=True)

    @field_validator("path", "target_dir", mode="before")
    @classmethod
    def _resolve_paths(cls, value) -> Path:
        """Resolve paths to absolute paths and expand user-home."""
        return Path(value).expanduser().resolve()

    @computed_field
    @property
    def reaction_time_seconds(self) -> float:
        return self.rxn_time * 60


def run_enzyme_kinetic_analysis_pipeline(args: KineticArgs):
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
        data=df,
        special_peaks=None,
    )

    # 5. Apply calibration
    # analyzer.apply_calibration()

    # 6. Perform enzyme kinetics analysis
    analyzer.fit()

    # 7. Plot and save results
    base_dest_path = args.target_dir / args.path.name
    analyzer.plot(base_dest_path, overwrite=args.overwrite)
    analyzer.export_csv(base_dest_path, overwrite=args.overwrite)

    _logger.info("Analysis complete.")


def main():
    args = tyro.cli(KineticArgs)
    run_enzyme_kinetic_analysis_pipeline(args)


if __name__ == "__main__":
    main()
