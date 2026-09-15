# enzyme-kinetics

A Python CLI and library for computing Michaelis-Menten enzyme kinetics from HPLC peak-area data.

## Overview

`enzyme-kinetics` takes pre-aggregated HPLC data (produced by PeakAnalyzer) and:

- fits a linear calibration curve (area = slope × [CoA] + intercept) per peak from a standards run,
- fits Michaelis-Menten (or Hill) kinetic models per peak to assay data,
- optionally converts raw peak areas to µM-based velocities via those calibration curves, so kcat
  and kcat/Km carry physical units, and
- exports results as a CSV plus diagnostic plots.

## Installation

Requires Python ≥ 3.14. Install with [uv](https://docs.astral.sh/uv/):

```bash
uv sync
```

## Usage

Both pipelines are subcommands of a single CLI:

```bash
uv run python -m enzyme_kinetics --help
uv run python -m enzyme_kinetics analyze --help
uv run python -m enzyme_kinetics calibrate --help
```

### Input format

Both subcommands read a pre-aggregated PeakAnalyzer CSV (`*_area_analysis.csv` or
`*_height_analysis.csv`) with one row per `peak_id` × `substrate_conc` combination. Column names
are normalized to snake_case on load, and index-like artifact columns (e.g. `Unnamed: 0`) are
dropped. Peak IDs are canonicalized to known compound names (HTAL, PDAL, OA, OLV) before fitting.

**Required columns:** `substrate_conc` (µM), `peak_id`, `mean`, `std`, `count`

---

### `analyze` — Kinetics pipeline

Fits kinetic models to each peak in the input CSV and writes a results CSV plus plots.

```
python -m enzyme_kinetics analyze PATH TARGET_DIR [options]
```

**Positional arguments**

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `PATH`       | Pre-aggregated PeakAnalyzer CSV of assay data |
| `TARGET_DIR` | Directory to write output files into       |

**Options**

| Flag                                 | Default  | Description                                                                                                                                          |
|--------------------------------------|----------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--calibration-path PATH`            | `None`   | Calibration **standards CSV** (same format as the input). A calibration curve is fitted per peak and applied, converting velocities to µM/s.          |
| `--enzyme-conc-um FLOAT`             | `12.25`  | Total enzyme concentration in µM (kcat = Vmax / [E]).                                                                                                |
| `--rxn-time FLOAT`                   | `180.0`  | Reaction duration in **minutes**.                                                                                                                    |
| `--substrate STR`                    | `HexCoA` | Substrate name propagated to results and plot labels.                                                                                                |
| `--peak-prefixes [STR ...]`          | `None`   | Chain-length prefixes (e.g. `C5`) to split off before peak ID canonicalization; merged with the built-in set (`C3`–`C7`). Omit to skip prefix handling. |
| `--is-calibrated / --no-is-calibrated` | `False`  | Mark input areas as already calibrated (affects plot units; suppresses the non-physical-units warning).                                              |
| `--overwrite / --no-overwrite`       | `False`  | Overwrite existing output files instead of writing a versioned copy.                                                                                 |
| `--verbose / --no-verbose`           | `False`  | Log at DEBUG level on the console.                                                                                                                   |

When `--calibration-path` is given, every peak in the assay CSV that also appears in the standards
CSV is re-fitted on calibrated velocities. A standards peak with no matching assay peak causes the
command to fail.

**Outputs** (in `TARGET_DIR`, named from the input file stem)

| File                         | Contents                                                        |
|------------------------------|-----------------------------------------------------------------|
| `<stem>_kinetics.csv`        | One row per peak: Km, Vmax, kcat, kcat/Km, and uncertainties    |
| `<stem>_kinetics.png`        | Michaelis-Menten (or Hill) fit curves                           |
| `<stem>_lineweaver_burk.png` | Lineweaver-Burk double-reciprocal plots                         |
| `<stem>_residuals.png`       | Fit residuals                                                   |
| `<stem>_efficiency.png`      | Km, Vmax, kcat, and kcat/Km comparison across peaks (OLV excluded) |

---

### `calibrate` — Calibration pipeline

Fits a linear HPLC calibration curve per peak and writes a parameter report and diagnostic plot
for each.

```
python -m enzyme_kinetics calibrate PATH OUTFILE [options]
```

**Positional arguments**

| Argument  | Description                                                              |
|-----------|--------------------------------------------------------------------------|
| `PATH`    | Pre-aggregated PeakAnalyzer CSV of calibration standards                 |
| `OUTFILE` | Base output path; `_<peak_id>` is appended to its stem for each peak     |

**Options**

| Flag                           | Default | Description                                                                 |
|--------------------------------|---------|-----------------------------------------------------------------------------|
| `--peak-prefixes [STR ...]`    | `None`  | Chain-length prefixes for peak ID canonicalization (see `analyze`).         |
| `--pretty / --no-pretty`       | `True`  | Write a human-readable text report; `--no-pretty` writes compact JSON.      |
| `--overwrite / --no-overwrite` | `False` | Accepted but currently unused; reports and plots are always overwritten.   |
| `--verbose / --no-verbose`     | `False` | Log at DEBUG level on the console.                                          |

In the calibration CSV, `substrate_conc` holds the known CoA standard concentrations (µM) and
`mean` holds the measured peak areas — the inverse of their roles in `analyze`.

**Outputs per peak:** `<outfile_stem>_<peak_id>.txt` (parameters: slope, intercept, R², RMSE,
inverse function) and `<outfile_stem>_<peak_id>.png` (diagnostic plot). The `.txt` extension is
used for both the pretty report and JSON output.

---

## Typical Workflow

```bash
# 1. (Optional) Inspect calibration quality for a standards run
uv run python -m enzyme_kinetics calibrate standards_area_analysis.csv results/calibration.txt

# 2. Run kinetics, fitting and applying calibration from the same standards CSV
uv run python -m enzyme_kinetics analyze experiment_area_analysis.csv results/ \
    --calibration-path standards_area_analysis.csv \
    --enzyme-conc-um 12.25 \
    --rxn-time 180
```

Without `--calibration-path` (and without `--is-calibrated`), kcat and kcat/Km are computed from
raw area/s velocities and carry non-physical units — suitable for relative comparisons only.

A calibration that fails the quality check (low R² or non-positive slope) is logged as a warning,
but fitting still proceeds; check the `calibrate` report before trusting the constants.

## Units

All concentration inputs must be in **µM**. Outputs assume µM input:

| Parameter | Unit    |
|-----------|---------|
| Km        | µM      |
| kcat      | s⁻¹     |
| kcat/Km   | s⁻¹·M⁻¹ |

Passing concentrations in any other unit will silently produce incorrect results.

## Logging and exit codes

Console logs go to stderr at INFO (DEBUG with `--verbose`). A DEBUG-level log file is written to
`./logs/enzyme_kinetics.log` relative to the working directory (rotated at 10 MB).

| Exit code | Meaning                                                      |
|-----------|--------------------------------------------------------------|
| `0`       | Success                                                      |
| `1`       | Unexpected error (traceback in log)                          |
| `2`       | Expected failure: missing file, missing columns, bad values  |
| `130`     | Interrupted                                                  |

## Python API

```python
from enzyme_kinetics.analysis.analyze import KineticAnalyzer
from enzyme_kinetics.calibration import fit_calibration

analyzer = KineticAnalyzer(
    enzyme_conc_um=12.25,
    reaction_time_seconds=180 * 60,
    substrate="HexCoA",
    data=df,
    special_peaks={"OLV": "hill"},  # per-peak model override
)
analyzer.fit()

cal = fit_calibration(std_df["substrate_conc"], std_df["mean"], sigma=std_df["std"])
analyzer.apply_calibration(cal, peak_ids=["OA"])

analyzer.export_csv(dest).plot(dest, is_calibrated=True)
```

**Model selection**

- **Default:** Michaelis-Menten (`v = Vmax·s / (Km + s)`)
- **Hill:** per peak via `special_peaks={peak_id: "hill"}` (Python API only; the CLI uses
  Michaelis-Menten for every peak)
- **Lineweaver-Burk** cross-validation is attempted for every peak; failures are logged and skipped

## Project Layout

```
src/enzyme_kinetics/
├── __main__.py        — `python -m enzyme_kinetics` entry point
├── cli/               — tyro CLI: subcommand dispatch, logging, analyze/calibrate runners
├── analysis/          — KineticArgs, KineticAnalyzer orchestrator, KineticPlots
├── calibration/       — CalibrationArgs, Calibration dataclass, fit/save/load/plot
├── core/              — kinetic models, kcat derivation, preprocessing, compound IDs
└── forks/             — PathBuilder, snake_case labels, CSV cleaning
```

## Development

```bash
uv sync --group dev
ruff check src/
ruff format src/
ty check src/
```

There is no test suite yet.
