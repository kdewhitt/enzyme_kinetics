"""Command-line interface for enzyme_kinetics."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Annotated

import tyro
from loguru import logger
from rich.console import Console

from enzyme_kinetics.analysis.arguments import KineticArgs
from enzyme_kinetics.calibration.arguments import CalibrationArgs

__all__ = ["main"]

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

DEFAULT_LOG_DIR = Path("./logs")
DEFAULT_LOG_FILE = "enzyme_kinetics.log"

CONSOLE_LOG_FORMAT = (
    "<yellow>{time:HH:mm:ss}</yellow> | "
    "<level>{level: <8}</level> | "
    "<cyan>{name}</cyan>:<green>{function}</green> | "
    "<level>{message}</level>"
)

FILE_LOG_FORMAT = (
    "{time:YYYY-MM-DD HH:mm:ss.SSS} | "
    "{level: <8} | "
    "{process.id} | "
    "{thread.name} | "
    "{name}:{function}:{line} | "
    "{message}"
)


def configure_logging(
    console_level: str = "INFO",
    *,
    file_level: str = "DEBUG",
    log_dir: Path = DEFAULT_LOG_DIR,
    log_file: str = DEFAULT_LOG_FILE,
) -> Path:
    """Configure console and file logging for the CLI."""
    console = Console(stderr=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / log_file

    logger.remove()

    # Console sink
    logger.add(
        lambda message: console.print(
            str(message),
            end="",
            markup=False,
            highlight=False,
            soft_wrap=True,
        ),
        level=console_level,
        format=CONSOLE_LOG_FORMAT,
        colorize=True,
        backtrace=False,
        diagnose=False,
    )

    # File sink
    logger.add(
        log_path,
        level=file_level,
        format=FILE_LOG_FORMAT,
        encoding="utf-8",
        rotation="10 MB",
        compression="zip",
        enqueue=True,
        backtrace=False,
        diagnose=False,
    )

    logger.enable("enzyme_kinetics")
    return log_path


Command = Annotated[
    Annotated[
        KineticArgs,
        tyro.conf.subcommand("analyze", description="Run enzyme kinetics analysis pipeline."),
    ]
    | Annotated[
        CalibrationArgs,
        tyro.conf.subcommand(
            "calibrate",
            description="Run per-peak HPLC calibration curve fitting pipeline.",
        ),
    ],
    tyro.conf.OmitArgPrefixes,
]


def main() -> int:
    """Run the enzyme_kinetics command-line interface."""
    cmd = tyro.cli(
        Command,
        description="Compute Michaelis-Menten enzyme kinetics and fit HPLC calibration curves.",
        compact_help=True,
    )

    log_level = "DEBUG" if getattr(cmd, "verbose", False) else "INFO"
    configure_logging(console_level=log_level)

    try:
        if isinstance(cmd, KineticArgs):
            from enzyme_kinetics.cli.analyze import main as run_enzyme_kinetic_analysis_pipeline

            run_enzyme_kinetic_analysis_pipeline(cmd)
        elif isinstance(cmd, CalibrationArgs):
            from enzyme_kinetics.cli.calibration import main as run_calibration_pipeline

            run_calibration_pipeline(cmd)
    except (OSError, ValueError, KeyError) as exc:
        logger.error("Command failed: {}", exc)
        return 2
    except KeyboardInterrupt:
        logger.warning("Command interrupted")
        return 130
    except Exception:
        logger.exception("Command failed unexpectedly")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
