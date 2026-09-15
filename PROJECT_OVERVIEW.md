# Project Overview: enzyme-kinetics

## Purpose

Computes Michaelis-Menten enzyme kinetics (Km, Vmax, kcat, kcat/Km) from pre-aggregated HPLC
peak-area data produced by PeakAnalyzer, with optional per-peak linear calibration so velocities
and derived constants carry physical units.

The package has two pipelines, both reachable from one CLI:

| Pipeline    | Runner                                             | Purpose |
|-------------|----------------------------------------------------|---------|
| `analyze`   | `cli/analyze.py::run_enzyme_kinetic_analysis_pipeline` | Fit kinetic models per peak, optionally calibrate, export CSV and plots |
| `calibrate` | `cli/calibration.py::run_calibration_pipeline`     | Fit and report a linear calibration curve per peak |

Entry points (`pyproject.toml` `[project.scripts]`):

| Script      | Target                               | Error handling |
|-------------|--------------------------------------|----------------|
| `pipeline`  | `enzyme_kinetics.cli.main:main`      | Subcommands `analyze` / `calibrate`; catches errors and returns exit codes 0/1/2/130 |
| `analyze`   | `enzyme_kinetics.cli.analyze:main`   | Exceptions propagate |
| `calibrate` | `enzyme_kinetics.cli.calibration:main` | Exceptions propagate |

`python -m enzyme_kinetics` is equivalent to `pipeline`. All three call
`cli.main.configure_logging`, which re-enables the `enzyme_kinetics` logger (the package disables
it on import in `__init__.py`) and adds a console sink plus a rotating DEBUG file sink at
`./logs/enzyme_kinetics.log`.

Typical workflow:

1. `pipeline calibrate standards.csv out/cal.txt` to inspect calibration quality per peak.
2. `pipeline analyze assay.csv out/ --calibration-path standards.csv`. The standards CSV is
   refitted inside `analyze`; the saved report or JSON from step 1 is not consumed by the CLI.

User-facing usage, options and output formats are documented in `README.md`.

---

## Module Map

```
src/enzyme_kinetics/
├── __init__.py            logger.disable("enzyme_kinetics")
├── __main__.py            python -m entry point → cli.main.main
├── cli/
│   ├── main.py            tyro subcommand union (KineticArgs | CalibrationArgs), configure_logging,
│   │                      exit-code mapping
│   ├── analyze.py         run_enzyme_kinetic_analysis_pipeline; apply_calibration (fits a
│   │                      Calibration per standards peak and applies it to the analyzer)
│   └── calibration.py     run_calibration_pipeline (fit → save report/JSON → diagnostic plot)
├── analysis/
│   ├── arguments.py       KineticArgs (Pydantic, frozen); reaction_time_seconds computed field
│   ├── analyze.py         _select_and_fit, KineticAnalyzer, analyze_peaks (batch helper)
│   └── plots.py           KineticPlots: fit curves, Lineweaver-Burk, residuals, efficiency bars
├── calibration/
│   ├── __init__.py        re-exports Calibration, CalibrationArgs, fit/plot/save_calibration
│   ├── arguments.py       CalibrationArgs (Pydantic, frozen)
│   └── calibrate.py       Calibration dataclass, fit_calibration, save/load_calibration,
│                          plot_calibration
├── core/
│   ├── models.py          model equations, effective_sigma, fit_model + per-model wrappers,
│   │                      FitResult, r_squared
│   ├── derive.py          derive_constants, KineticConstants, covariance-aware kcat/Km
│   ├── preprocess.py      column-name constants, extract_peak_data, prepare_velocity
│   └── compounds.py       Compounds enum, aliases, canonicalize_peak_ids, prefix handling
└── forks/                 general-purpose utilities vendored from other projects
    ├── __init__.py        re-exports PathBuilder, to_snake_case, CleaningOptions, read_clean_csv
    ├── builder.py         PathBuilder: immutable, tagged, versioned, collision-safe output paths
    ├── case.py            to_snake_case, to_slug, asciify
    └── clean.py           read_clean_csv, index-artifact removal, unique snake_case columns
```

Outside `src/`: `excluded/fitting.py` holds the former `core/calibration.py`. It is not packaged
or imported.

### Import dependencies

```
cli ──► analysis ──► calibration ──► core.models
 │         │
 │         ├──► core.derive ──► core.models
 │         ├──► core.preprocess
 │         └──► forks (PathBuilder)
 ├──► core.compounds
 └──► forks (read_clean_csv)
```

`core` has no dependency on `analysis`, `calibration`, `cli` or `forks`. `cli.analyze` and
`cli.calibration` import `configure_logging` from `cli.main`. `cli.main` imports only the two
argument models at module level and imports the runners lazily inside `main()`.

---

## Data Flow

### `analyze`

```
assay CSV (PeakAnalyzer output)
  └─ read_clean_csv(drop_indexlike=True, snake_case=True)   forks/clean.py
       └─ required-column check: substrate_conc, peak_id, mean, std, count   cli/analyze.py
            └─ canonicalize_peak_ids(prefixes=args.peak_prefixes)          core/compounds.py
                 → peak_id becomes an ordered pandas.Categorical
                 └─ sort by peak_id, substrate_conc
                      └─ KineticAnalyzer(enzyme_conc_um, reaction_time_seconds, substrate, df)
                           │
                           ├─ fit()                              raw pass, area/s
                           │    for each peak_id:
                           │      extract_peak_data → prepare_velocity   v = mean / t, v_std = std / t
                           │      _select_and_fit → fit_michaelis_menten | fit_hill   (sigma = v_std)
                           │      fit_lineweaver_burk(s, v)             failure → warning, lb_fit=None
                           │      derive_constants → KineticConstants
                           │    primary-fit failure → warning, peak omitted from results
                           │
                           ├─ [--calibration-path] cli.analyze.apply_calibration
                           │    read + clean + canonicalize standards CSV
                           │    for each standards peak_id present in analyzer.results:
                           │      fit_calibration(substrate_conc, mean, sigma=std)
                           │      KineticAnalyzer.apply_calibration(cal, peak_ids=[peak_id])
                           │        is_valid() false → warning only
                           │        v_um     = area_to_conc(mean) / t
                           │        v_std_um = area_to_conc_with_error(mean, std)[1] / t
                           │        refit primary model and LB on (s, v_um)
                           │        derive_constants → replaces results[peak_id]
                           │    standards peaks with no raw fit → warning, skipped
                           │
                           ├─ export_csv(target_dir / input.name)
                           │    KineticConstants.to_dict() rows → <stem>_kinetics.csv
                           │    (PathBuilder.reserve(): atomic, versioned unless --overwrite)
                           │
                           └─ plot(target_dir / input.name, is_calibrated)
                                KineticPlots.plot → plot_lineweaver_burk → plot_residuals
                                  → plot_efficiency_comparison(exclude={"OLV"})
                                  PNGs: <stem>_kinetics / _lineweaver_burk / _residuals / _efficiency
```

`is_calibrated` passed to `plot()` is `True` when `--calibration-path` is given or when
`--is-calibrated` is set. It changes only labels and whether plotted points are converted; no
calculation depends on it.

### `calibrate`

```
standards CSV
  └─ read_clean_csv → required-column check → canonicalize_peak_ids → sort
       └─ for each peak_id:
            fit_calibration(substrate_conc, mean, sigma=std)
            save_calibration(nice=args.pretty)
              → <outfile_stem>_<peak_id>.txt   (pretty text report)
              → <outfile_stem>_<peak_id>.json  (--no-pretty; loadable by load_calibration)
            plot_calibration → <outfile_stem>_<peak_id>.png  (fit, residuals, standardized residuals)
```

Output paths are built with `Path.with_name`, not `PathBuilder`, so they are always overwritten.
`CalibrationArgs.overwrite` is accepted but unused.

---

## Units of Measure

### Input columns

| Column | Unit | Notes |
|---|---|---|
| `substrate_conc` | µM | Must be µM throughout; no runtime check. In `calibrate` this is the known standard concentration |
| `mean` | area | Integrated peak area (instrument counts) |
| `std` | area | Replicate standard deviation; used directly as fit sigma |
| `count` | replicates | Required by the column check and passed to `prepare_velocity`, but not used in any calculation |

### Time and enzyme

| Variable | Unit |
|---|---|
| `KineticArgs.rxn_time` (CLI `--rxn-time`) | minutes |
| `KineticArgs.reaction_time_seconds`, `KineticAnalyzer.reaction_time_seconds`, `rxn_time` args in `core.preprocess` | seconds |
| `enzyme_conc_um` | µM |

### Calibration (`Calibration` dataclass)

| Field | Unit | Notes |
|---|---|---|
| `slope`, `slope_std` | area/µM | |
| `intercept`, `intercept_std` | area | |
| `slope_intercept_cov` | area²/µM | Defaults to 0.0 when loaded from JSON saved without it |
| `r_squared` | — | `is_valid()` requires ≥ 0.999 and slope > 0 |
| `rmse` | area | |
| `n_points` | count | Valid standards used |
| `conc_range` | µM | (min, max) of standards |

### Velocity and fit parameters

| Quantity | Uncalibrated | Calibrated |
|---|---|---|
| `v`, `v_std` | area/s | µM/s |
| `FitResult.vmax`, `vmax_std` | area/s | µM/s |
| `FitResult.km` / `k_half`, `ki`, S0 | µM | µM |
| `FitResult.hill_n` | — | — |
| `FitResult.pcov` | mixed: each entry is the product of its two parameters' units | |
| LB `s_inv` | µM⁻¹ | µM⁻¹ |
| LB `v_inv` | s/area | s/µM |

### Derived constants (`KineticConstants`)

| Field | Calibrated unit | Uncalibrated unit |
|---|---|---|
| `kcat`, `kcat_std` | s⁻¹ | area·µM⁻¹·s⁻¹ (non-physical) |
| `kcat_km_M`, `kcat_km_std_M` | s⁻¹·M⁻¹ | area·µM⁻¹·s⁻¹·M⁻¹ (non-physical) |
| `ki`, `ki_std` | µM | µM (substrate-inhibition fits only, else `None`) |

---

## Key Design Decisions

### Raw-then-calibrated fitting

`fit()` always runs first on area/s velocities; `apply_calibration()` rebuilds results for the
selected peaks from the original `mean` areas. `extract_peak_data` / `prepare_velocity` return raw
`mean_signal` as a fourth element so calibration is applied to area directly rather than to
velocity × time. With `peak_ids=None` the results dict is replaced wholesale, so peaks whose
calibrated refit fails disappear. With `peak_ids=[...]` (the CLI path) results are updated in
place, and IDs not already in `results` raise `ValueError`. Each applied `Calibration` is stored
in `KineticAnalyzer.calibrations[peak_id]` so plots can convert observed points to the same units.

### Weighting and sigma sanitization

`core.models.effective_sigma` is shared by `fit_model` and `fit_calibration`. Zero or non-finite
sigma entries are floored at the smallest positive finite sigma; if none exists, sigma is `None`.
`absolute_sigma` is `True` only when a sigma array is actually used. Unweighted fits therefore
get their covariance scaled by the residual variance instead of assuming unit errors.
`fit_calibration` additionally drops rows whose concentration, area or sigma is non-finite before
fitting, and it requires at least two points.

### Covariance-aware error propagation

- **kcat/Km** (`derive._kcat_km_with_covariance`): first-order delta method using `pcov[0,0]`,
  `pcov[1,1]` and `pcov[0,1]`. Km is converted µM→M (×1e-6) so the result is in s⁻¹·M⁻¹, with the
  partial derivatives scaled to match. Returns `(nan, nan)` when Km ≤ 0, Vmax ≤ 0 or the
  covariance is non-finite, as with LB fits.
- **kcat**: `vmax / [E]` and `vmax_std / [E]`; enzyme concentration is treated as exact.
- **Area → concentration** (`Calibration.area_to_conc_with_error`): delta method over area std,
  slope std, intercept std and the slope–intercept covariance, clamped at zero variance.
- **Calibrated velocity std** in `apply_calibration`: the precise delta method by default;
  `precise_std=False` uses `v_std / slope`.
- **Plot error bars** (`KineticPlots._observed`) always use the simple `v_std / slope`, so after
  calibration they can differ slightly from the sigmas used in the fit.

### Lineweaver-Burk as a cross-check only

`fit_lineweaver_burk` is an OLS regression of 1/v on 1/[S] over positive points
(`scipy.stats.linregress`). It back-calculates Vmax = 1/intercept and Km = slope·Vmax, and stores
slope, intercept and std_err in `extra`. Its `perr` and `pcov` are NaN by design. It never feeds
kcat or kcat/Km; it only supplies the `km_lb`, `vmax_lb` and `r_squared_lb` columns and the LB
plot. After calibration it is refitted on µM/s velocities so its units match the primary fit.

### Model tagging and parameter access

`FitResult` stores `popt` as a tuple, with index 0 = Vmax and index 1 = Km / k_half, so
`derive_constants` works uniformly across models. Model-specific properties (`k_half`, `hill_n`,
`ki` and their `_std` variants) raise `AttributeError` for the wrong `model_type`.
`KineticConstants.to_dict()` names the second parameter `k_half` for Hill fits and `km`
otherwise, and emits `hill_n` / `ki` columns only for the relevant models.

### Peak identity and ordering

`Compounds` (`StrEnum`: HTAL, PDAL, OA, OLV) is the canonical identity. Aliases are normalized by
`canonicalize_label`, which uppercases and joins non-alphanumeric runs with `_`, and are checked
for collisions at import time. `canonicalize_peak_ids`:

- uppercases raw IDs;
- strips optional chain-length prefixes, matching them case-insensitively;
- resolves the remainder with `Compounds.resolve`, raising `ValueError` on unknown labels;
- reattaches the prefix;
- stores the column as an ordered `Categorical` in which each base compound is followed by its
  prefixed variants.

Prefix handling is off when `peak_prefixes` is `None`. When the CLI supplies any value,
`coerce_peak_prefixes` merges it with the built-in `C3`–`C7`.

### Output path safety

`forks.builder.PathBuilder` is a frozen dataclass whose `with_*` methods return new instances. It
sanitizes tags and stems, truncates filenames by byte count with a fixed reserve for the version
suffix, avoids Windows reserved device names, and versions on collision (`_v2`, `_v3`, … up to
`_v999`). `reserve()` claims a name atomically with an exclusive create; `export_csv` uses it.
`KineticPlots` passes the builder straight to `savefig` through `__fspath__`, which resolves
without reserving, and builds with `create=True`.

### Column cleaning

`forks.clean.read_clean_csv` drops index artifacts (`index`, `idx`, `Unnamed: N`) before
relabelling, because snake_casing would rewrite the colon that pattern relies on. It then converts
labels with `to_snake_case`, which splits camelCase and acronyms and transliterates to ASCII,
and resolves duplicate labels with numeric suffixes instead of dropping columns.

---

## Kinetic Models

| `model_type` | Equation function | Fit wrapper | Parameters (`popt` order) | Initial guess | Used by `KineticAnalyzer` |
|---|---|---|---|---|---|
| `"mm"` | `michaelis_menten` | `fit_michaelis_menten` | Vmax, Km | max(v), median(s) | Default for every peak |
| `"hill"` | `hill_equation` | `fit_hill` | Vmax, k_half, n | max(v), median(s), 2.0 | `special_peaks={id: "hill"}` (Python API only) |
| `"tmm"` | `threshold_michaelis_menten` | `fit_threshold_michaelis_menten` | Vmax, Km, S0 | 2·max(v), median(s), 0.5 | No (dispatch commented out) |
| `"si"` | `substrate_inhibition` | `fit_substrate_inhibition` | Vmax, Km, Ki | max(v), median(s), max(s) | No |
| `"lb"` | `michaelis_menten` (for `predict`) | `fit_lineweaver_burk` | Vmax, Km (back-calculated) | — (linear) | Always, as cross-check |

Non-linear fits go through `fit_model`: `curve_fit` with bounds `(0, inf)` on every parameter and
`maxfev=10_000`. The CLI passes `special_peaks=None`, so it fits Michaelis-Menten everywhere.
`analyze_peaks` is a standalone batch helper: it runs Michaelis-Menten plus LB fits on
pre-extracted arrays, skips peaks with fewer than 3 points, and does not calibrate.

## Analytes Tracked

| Enum member | Full name | Accepted aliases |
|---|---|---|
| `HTAL` | Hexanoyl triacetic acid lactone | hexanoyl triacetic acid, hexanoyl-triacetic acid |
| `PDAL` | Pentyl diacetic acid lactone | pentyl diacetic acid, pentyl-diacetic acid |
| `OA` | Olivetolic acid | olivetolic acid, ola |
| `OLV` | Olivetol | olivetol, olivetolate |

`plots._LABEL_MAP` gives display names for `OA` and `OLV` only.

---

## Known Gaps and Inconsistencies

Observed in the current code; not yet addressed.

- **Std vs SEM:** `std` is used as the absolute fit sigma for points that are replicate means, and
  `count` is unused, so parameter uncertainties may be inflated by √n.
- **Hard-coded choices:** the efficiency plot always excludes `OLV`, and Hill selection cannot be
  set from the CLI.
- **Unused options:** `CalibrationArgs.overwrite`.
- **Unwired models:** the threshold Michaelis-Menten and substrate-inhibition fits exist but
  `_select_and_fit` never dispatches to them.
- **Repo hygiene:** `logs/` is not in `.gitignore`.
- **Tests:** there is no test suite (`pytest` is in the `test` dependency group).

---

## Tooling and Conventions

- **Python:** ≥ 3.14 (`.python-version` = 3.14). Code uses PEP 695 generics (`def f[T: ...]`,
  `type X = ...`) and `typing.Self`.
- **Build:** `uv_build`; dependencies are managed with uv (`uv.lock`).
- **Lint and format:** ruff with line length 99 and Google pydocstyle convention. The extended rule
  sets are A, B, C4, COM, D, E, EM, F, FA, FURB, SIM, UP and W.
- **Type checking:** `ty`, with `pandas-stubs` and `scipy-stubs`.
- **Pre-commit:** `uv-sync`, `ruff-check --fix`, `ruff-format`.
- **Docstrings:** Google style throughout.
  - Public functions and classes have `Args`, `Returns` and `Raises` sections, dataclasses and
    Pydantic models have `Attributes`, and physical quantities state their units.
  - Most modules include a `Typical usage:` example.
  - Private helpers usually have a one-line summary.
- **Comments:** runners use numbered step comments (`# 1. Load data`, …). Some fixes are tagged
  inline (`# FIX 1`, `# FIX 2`, `# FIX 4`, `# FIX 7`), referring to an earlier refactor.
- **Logging:** loguru with `{}`-style placeholders, disabled at package import and enabled by the
  CLI.
