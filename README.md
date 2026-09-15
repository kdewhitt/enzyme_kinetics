# enzyme-kinetics

A Python CLI and library for computing Michaelis-Menten enzyme kinetics from HPLC peak-area data.

## Overview

`enzyme-kinetics` takes pre-aggregated HPLC data (produced by PeakAnalyzer) and:

- fits a linear calibration curve (area = slope × [CoA] + intercept) per peak from a standards run,
- fits a Michaelis-Menten model per peak to assay data (Hill is available through the Python API),
  with a Lineweaver-Burk fit alongside for cross-checking,
- optionally converts peak areas to µM-based velocities via those calibration curves, so kcat and
  kcat/Km carry physical units, and
- exports results as a CSV plus diagnostic plots.

## Installation

Requires Python ≥ 3.14. Install with [uv](https://docs.astral.sh/uv/):

```bash
uv sync
```

## Usage

Installing the package provides three console scripts:

| Command     | Equivalent                             | Description                               |
|-------------|----------------------------------------|-------------------------------------------|
| `pipeline`  | `python -m enzyme_kinetics`            | Single CLI with `analyze` and `calibrate` subcommands |
| `analyze`   | `python -m enzyme_kinetics analyze`    | Kinetics pipeline only                    |
| `calibrate` | `python -m enzyme_kinetics calibrate`  | Calibration pipeline only                 |

```bash
uv run pipeline --help
uv run pipeline analyze --help
uv run pipeline calibrate --help
```

The options are the same either way. The only difference is error handling: `pipeline` /
`python -m enzyme_kinetics` catches errors, logs them and returns an exit code (see
[Logging and exit codes](#logging-and-exit-codes)). The standalone `analyze` and `calibrate`
scripts let exceptions propagate with a traceback.

### Input format

Both subcommands read a pre-aggregated PeakAnalyzer CSV (`*_area_analysis.csv` or
`*_height_analysis.csv`) with one row per `peak_id` × `substrate_conc` combination.

**Required columns:** `substrate_conc` (µM), `peak_id`, `mean`, `std`, `count`

On load:

- Column names are normalized to snake_case (e.g. `Substrate Conc` → `substrate_conc`), and
  index-like artifact columns (`index`, `idx`, `Unnamed: 0`) are dropped.
- Peak IDs are canonicalized, case-insensitively, to the known compounds **HTAL**, **PDAL**,
  **OA** and **OLV**. Accepted aliases include `hexanoyl triacetic acid`, `pentyl diacetic acid`,
  `olivetolic acid`, `OLA`, `olivetol` and `olivetolate`. An unrecognized peak ID is an error.
- Rows are sorted by `peak_id` then `substrate_conc`.

Rows missing `substrate_conc` or `mean` are ignored, and so is any peak whose `mean` values are all
zero.

---

### `analyze`: kinetics pipeline

Fits kinetic models to each peak in the input CSV and writes a results CSV plus plots.

```
pipeline analyze PATH TARGET_DIR [options]
```

**Positional arguments**

| Argument     | Description                                    |
|--------------|------------------------------------------------|
| `PATH`       | Pre-aggregated PeakAnalyzer CSV of assay data  |
| `TARGET_DIR` | Directory to write output files into (created if missing) |

**Options**

| Flag                                  | Default  | Description |
|---------------------------------------|----------|-------------|
| `--calibration-path PATH`             | `None`   | Calibration **standards CSV** (same format as the input, not a saved calibration JSON). A calibration curve is fitted per peak and applied, converting velocities to µM/s. |
| `--enzyme-conc-um FLOAT`              | `12.25`  | Total enzyme concentration in µM (kcat = Vmax / [E]). |
| `--rxn-time FLOAT`                    | `180.0`  | Reaction duration in **minutes**. |
| `--substrate STR`                     | `HexCoA` | Substrate name propagated to results and plot labels. |
| `--peak-prefixes [STR ...]`           | `None`   | Chain-length prefixes (e.g. `C5`) to split off before peak ID canonicalization and reattach afterwards (`c5-ola` → `C5-OA`). Values are merged with the built-in set `C3`–`C7`. Omit to skip prefix handling. |
| `--is-calibrated / --no-is-calibrated` | `False` | Mark input areas as already calibrated. This switches plot units to µM/s and suppresses the non-physical-units warning; it does not change any calculation. |
| `--overwrite / --no-overwrite`        | `False`  | Overwrite existing output files instead of writing a versioned copy (`_v2`, `_v3`, …). |
| `--verbose / --no-verbose`            | `False`  | Log at DEBUG level on the console. |

#### How velocities are computed

Velocity is an endpoint rate: `v = mean peak area / reaction time (s)`, with `v_std = std / reaction
time`. Each peak is fitted first on these raw area/s velocities.

When `--calibration-path` is given, a calibration is fitted for each peak in the standards CSV.
Each assay peak that already has a raw fit is then re-fitted on calibrated velocities:
`v = area_to_conc(mean area) / reaction time`, in µM/s. The velocity uncertainty comes from the
delta method and combines replicate std, slope and intercept uncertainty, and their covariance.
A standards peak with no matching assay fit is skipped with a warning. The Lineweaver-Burk fit is
also recomputed on the calibrated velocities.

A calibration that fails the quality check (R² < 0.999 or slope ≤ 0) is logged as a warning, but
fitting still proceeds. Check the `calibrate` report before trusting the constants.

**Outputs** (in `TARGET_DIR`, named from the input file stem)

| File                         | Contents |
|------------------------------|----------|
| `<stem>_kinetics.csv`        | One row per peak (columns below) |
| `<stem>_kinetics.png`        | Observed velocities with error bars and the fitted curve, annotated with Km (K½ for Hill), Vmax and R² |
| `<stem>_lineweaver_burk.png` | 1/v vs 1/[S] with the regression line and intercepts |
| `<stem>_residuals.png`       | Fit residuals vs [S] |
| `<stem>_efficiency.png`      | Bar charts of Km, Vmax, kcat and kcat/Km across peaks (**OLV is excluded**) |

**Results CSV columns**

| Column | Description |
|--------|-------------|
| `substrate`, `peak_id` | Identifiers |
| `model` | Model tag: `mm` (Michaelis-Menten) or `hill` |
| `vmax`, `vmax_std` | Vmax (area/s uncalibrated, µM/s calibrated) |
| `km`, `km_std` | Km in µM (`k_half`, `k_half_std` for Hill fits) |
| `r_squared` | R² of the primary fit |
| `kcat`, `kcat_std` | Vmax / [E] |
| `kcat_km_M`, `kcat_km_std_M` | kcat/Km in s⁻¹·M⁻¹, with the error propagated through the full Vmax–Km covariance |
| `hill_n`, `hill_n_std` | Hill coefficient (Hill fits only) |
| `km_lb`, `vmax_lb`, `r_squared_lb` | Lineweaver-Burk estimates (when that fit succeeds). Here R² is measured in reciprocal space. |

---

### `calibrate`: calibration pipeline

Fits a linear HPLC calibration curve per peak and writes a parameter report and diagnostic plot
for each. Use it to inspect calibration quality; `analyze --calibration-path` refits from the
standards CSV itself.

```
pipeline calibrate PATH OUTFILE [options]
```

**Positional arguments**

| Argument  | Description |
|-----------|-------------|
| `PATH`    | Pre-aggregated PeakAnalyzer CSV of calibration standards |
| `OUTFILE` | Base output path; `_<peak_id>` is appended to its stem for each peak. The parent directory is created if missing. |

**Options**

| Flag                           | Default | Description |
|--------------------------------|---------|-------------|
| `--peak-prefixes [STR ...]`    | `None`  | Chain-length prefixes for peak ID canonicalization (see `analyze`). |
| `--pretty / --no-pretty`       | `True`  | Write a human-readable text report; `--no-pretty` writes compact JSON. |
| `--overwrite / --no-overwrite` | `False` | Accepted but currently unused; reports and plots are always overwritten. |
| `--verbose / --no-verbose`     | `False` | Log at DEBUG level on the console. |

In the calibration CSV, `substrate_conc` holds the known CoA standard concentrations (µM) and
`mean` holds the measured peak areas, the reverse of their roles in `analyze`.

**Outputs per peak**

| File | Contents |
|------|----------|
| `<outfile_stem>_<peak_id>.txt` (`--pretty`) | Slope, intercept, their std devs, R², RMSE, point count, concentration range, inverse function, and a Python snippet with the constants |
| `<outfile_stem>_<peak_id>.json` (`--no-pretty`) | The same parameters plus the slope–intercept covariance, readable with `load_calibration` |
| `<outfile_stem>_<peak_id>.png` | Three panels: calibration curve, absolute residuals, standardized residuals (±3σ) |

---

## Typical workflow

```bash
# 1. (Optional) Inspect calibration quality for a standards run
uv run pipeline calibrate standards_area_analysis.csv results/calibration.txt

# 2. Run kinetics, fitting and applying calibration from the same standards CSV
uv run pipeline analyze experiment_area_analysis.csv results/ \
    --calibration-path standards_area_analysis.csv \
    --enzyme-conc-um 12.25 \
    --rxn-time 180
```

Without `--calibration-path` (and without `--is-calibrated`), kcat and kcat/Km are computed from
raw area/s velocities and carry non-physical units. Use them for relative comparisons only.

## Fitting details

- **Models:** kinetic fits use `scipy.optimize.curve_fit` with parameters bounded to be
  non-negative; the calibration line is fitted unbounded.
  The CLI fits Michaelis-Menten (`v = Vmax·[S] / (Km + [S])`) for every peak.
- **Weighting:** the per-point `std` is used as the fit sigma, and uncertainties are treated as
  absolute. Zero or non-finite std values (e.g. blanks, single replicates) are raised to the
  smallest positive std in the set. If no positive std exists, the fit is unweighted and
  parameter errors are scaled by the residual variance.
- **Note:** `std` is the replicate standard deviation as supplied, not the standard error of the
  mean; `count` is not currently used.
- **Calibration:** rows with a non-numeric or missing concentration, area or std are excluded
  before fitting, and at least two valid standards are required.
- **Lineweaver-Burk:** an ordinary least-squares regression of 1/v on 1/[S], using only positive
  [S] and v. It is a cross-check only; its parameter uncertainties are not reported.

## Units

All concentration inputs must be in **µM**. Outputs assume µM input:

| Parameter | Unit (calibrated) |
|-----------|-------------------|
| Km        | µM                |
| Vmax      | µM/s              |
| kcat      | s⁻¹               |
| kcat/Km   | s⁻¹·M⁻¹           |

Passing concentrations in any other unit will silently produce incorrect results.

## Logging and exit codes

All three commands log to the console (stderr) at INFO, or DEBUG with `--verbose`. They also write
a DEBUG-level log file to `./logs/enzyme_kinetics.log`, relative to the working directory. The
file rotates at 10 MB, and rotated files are zip-compressed.

Exit codes for `pipeline` / `python -m enzyme_kinetics`:

| Exit code | Meaning |
|-----------|---------|
| `0`       | Success |
| `1`       | Unexpected error (traceback in log) |
| `2`       | Expected failure (`OSError`, `ValueError`, `KeyError`): missing file, missing columns, unknown peak ID, too few calibration points |
| `130`     | Interrupted |

A fit that fails for a single peak does not stop the run: it is logged as a warning and that peak
is left out of the results.

## Python API

```python
from enzyme_kinetics.analysis.analyze import KineticAnalyzer
from enzyme_kinetics.calibration import fit_calibration

analyzer = KineticAnalyzer(
    enzyme_conc_um=12.25,
    reaction_time_seconds=180 * 60,
    substrate="HexCoA",
    data=df,                         # cleaned + canonicalized DataFrame
    special_peaks={"OLV": "hill"},   # per-peak model override
)
analyzer.fit()                       # raw area/s fits; must run before apply_calibration

cal = fit_calibration(std_df["substrate_conc"], std_df["mean"], sigma=std_df["std"])
analyzer.apply_calibration(cal, peak_ids=["OA"])   # omit peak_ids to recalibrate every peak

analyzer.export_csv(dest).plot(dest, is_calibrated=True)
results_df = analyzer.to_dataframe()
```

- `apply_calibration(..., peak_ids=[...])` raises `ValueError` for peak IDs without an existing
  fit. `precise_std=False` switches to the simpler `v_std / slope` error estimate.
- `Calibration` provides `area_to_conc`, `conc_to_area`, `area_to_conc_with_error` and
  `is_valid(r2_threshold=0.999)`. `save_calibration` and `load_calibration` live in
  `enzyme_kinetics.calibration.calibrate`.
- `enzyme_kinetics.core.models` also provides substrate-inhibition (`fit_substrate_inhibition`)
  and threshold Michaelis-Menten (`fit_threshold_michaelis_menten`) fits. `KineticAnalyzer` does
  not use them yet.
- `enzyme_kinetics.core.derive.derive_constants` turns any `FitResult` into `KineticConstants`
  (kcat, kcat/Km and, for substrate-inhibition fits, Ki).

**Model selection in `KineticAnalyzer`**

- **Default:** Michaelis-Menten
- **Hill:** per peak via `special_peaks={peak_id: "hill"}`. This is available only in the Python
  API; the CLI uses Michaelis-Menten for every peak.
- **Lineweaver-Burk** is fitted for every peak; failures are logged and skipped.

## Project layout

```
src/enzyme_kinetics/
├── __main__.py          `python -m enzyme_kinetics` entry point
├── cli/
│   ├── main.py          `pipeline` subcommand dispatch, logging setup, exit codes
│   ├── analyze.py       `analyze` runner (loads data, applies calibration per peak)
│   └── calibration.py   `calibrate` runner
├── analysis/
│   ├── arguments.py     KineticArgs (CLI options)
│   ├── analyze.py       KineticAnalyzer orchestrator, analyze_peaks batch helper
│   └── plots.py         KineticPlots: fit curves, LB, residuals, efficiency charts
├── calibration/
│   ├── arguments.py     CalibrationArgs (CLI options)
│   └── calibrate.py     Calibration dataclass; fit/save/load/plot
├── core/
│   ├── models.py        kinetic model functions, curve_fit wrappers, FitResult
│   ├── derive.py        kcat and kcat/Km with covariance propagation, KineticConstants
│   ├── preprocess.py    per-peak extraction and velocity computation
│   └── compounds.py     Compounds enum, peak ID aliases and canonicalization
└── forks/
    ├── builder.py       PathBuilder: tagged, versioned, collision-safe output paths
    ├── case.py          snake_case / slug normalization
    └── clean.py         CSV column cleaning (read_clean_csv)
```

## Development

```bash
uv sync --group dev
uv run ruff check src/
uv run ruff format src/
uv run ty check src/
```

Pre-commit hooks (`.pre-commit-config.yaml`) run `uv sync`, `ruff check --fix` and `ruff format`.
`pytest` is included in the dev group, but there is no test suite yet.
