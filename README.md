# enzyme-kinetics

A Python CLI for computing Michaelis-Menten enzyme kinetics from HPLC peak-area data.

## Overview

`enzyme-kinetics` processes pre-aggregated HPLC data (produced by **PeakAnalyzer**) to fit
Michaelis-Menten (or Hill) kinetic models per peak, apply optional HPLC calibration curves to
convert raw peak areas to physically meaningful concentration-based velocities, and export
results as tagged CSV files and diagnostic plots.

## Installation

Requires Python ≥ 3.14. Install with [uv](https://docs.astral.sh/uv/):

```bash
uv sync
```

The two private dependencies (`kgdlibs`, `logurich`) are pulled directly from GitHub and
require network access during installation.

## CLI Commands

### `analyze` — Kinetics pipeline

Fits kinetic models to each HPLC peak in the input CSV and writes a results CSV plus plots.

```
analyze <path> <target_dir> [options]
```

**Positional arguments**

| Argument | Description |
|---|---|
| `path` | Pre-aggregated PeakAnalyzer CSV (`_area_analysis.csv` or `_height_analysis.csv`) |
| `target_dir` | Directory to write output files into |

**Options**

| Flag | Default | Description |
|---|---|---|
| `--calibration-path PATH` | `None` | Path to a calibration `.txt` file. Converts velocities to µM/s and makes kcat (s⁻¹) and kcat/Km (M⁻¹·s⁻¹) physically meaningful. |
| `--enzyme-conc-um FLOAT` | `12.25` | Total enzyme concentration in µM (used to compute kcat = Vmax / [E]). |
| `--rxn-time FLOAT` | `180.0` | Reaction duration in minutes. |
| `--substrate STR` | `"HexCoA"` | Substrate name propagated to result labels. |
| `--is-calibrated` | `False` | Pass if input mean areas have already had calibration applied. |
| `--overwrite` | `False` | Overwrite existing output files instead of appending a unique suffix. |
| `--verbose` | `False` | Enable verbose logging. |

**Required CSV columns:** `substrate_conc` (µM), `peak_id`, `mean`, `std`, `count`

**Outputs**
- `<target_dir>/<input_stem>_kinetics.csv` — one row per peak with Km, Vmax, kcat, kcat/Km
- Michaelis-Menten curve plots, Lineweaver-Burk plots, residual plots, efficiency comparison panel

---

### `calibrate` — Calibration pipeline

Fits a linear HPLC calibration curve (area = slope × [CoA] + intercept) per peak and writes
parameter reports and diagnostic plots.

```
calibrate <path> <outfile> [options]
```

**Positional arguments**

| Argument | Description |
|---|---|
| `path` | Pre-aggregated PeakAnalyzer CSV containing calibration standards |
| `outfile` | Base output path; a `_{peak_id}` suffix is appended per peak |

**Options**

| Flag | Default | Description |
|---|---|---|
| `--peak-prefixes ...` | `None` | Custom peak-ID prefix configuration for canonicalization. |
| `--pretty / --no-pretty` | `True` | Write human-readable plain-text reports (vs. compact JSON). |
| `--overwrite` | `False` | Overwrite existing output files. |
| `--verbose` | `False` | Enable verbose logging. |

In the calibration CSV, `substrate_conc` holds known CoA standard concentrations (µM) and `mean`
holds measured peak areas — the inverse of their roles in `analyze`.

**Outputs per peak:** `<outfile>_{peak_id}.txt` (parameters) + `<outfile>_{peak_id}.png` (3-panel diagnostic plot)

---

## Typical Workflow

```bash
# 1. Fit calibration curves from a standards run
calibrate standards_area_analysis.csv results/calibration

# 2. Run kinetics with calibration applied
analyze experiment_area_analysis.csv results/ \
    --calibration-path results/calibration_HexCoA.txt \
    --enzyme-conc-um 12.25 \
    --rxn-time 180
```

Without `--calibration-path`, kcat and kcat/Km are computed from raw area/s velocities and
carry non-physical units — suitable for relative comparisons only.

## Units

All concentration inputs must be in **µM**. Outputs assume µM input:

| Parameter | Unit |
|---|---|
| Km | µM |
| kcat | s⁻¹ |
| kcat/Km | M⁻¹·s⁻¹ |

Passing concentrations in any other unit will silently produce incorrect results.

## Model Selection

- **Default:** Michaelis-Menten (`v = Vmax·s / (Km + s)`)
- **Hill model:** enabled per peak via `special_peaks` in the Python API
- **Lineweaver-Burk** cross-validation is attempted for every peak; failures are logged and skipped

## Development

```bash
uv sync --group dev
ruff check src/
ruff format src/
pytest
```
