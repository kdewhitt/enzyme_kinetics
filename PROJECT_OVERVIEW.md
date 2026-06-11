# Project Overview: enzyme-kinetics

## Purpose

Computes Michaelis-Menten enzyme kinetics from pre-aggregated HPLC peak area data.
Two CLI entry points are provided: `analyze` (main kinetic pipeline) and `calibrate`
(per-peak linear calibration curve fitting).

The typical workflow is:
1. Run `calibrate` on CoA standard data to produce a `Calibration` per peak.
2. Run `analyze` on assay data, optionally supplying the calibration file from step 1.

---

## Module Map

```
src/enzyme_kinetics/
├── main.py          — `analyze` CLI; KineticArgs Pydantic model
├── calibrate.py     — `calibrate` CLI; CalibrationArgs Pydantic model
└── core/
    ├── io.py        — CSV loading, column sanitization
    ├── preprocess.py — velocity computation from area/time, per-peak extraction
    ├── compounds.py  — Compounds enum, peak ID canonicalization, categorical ordering
    ├── models.py     — kinetic model functions + curve_fit wrappers, FitResult
    ├── calibration.py — linear calibration fit, Calibration dataclass, plot/save/load
    ├── derive.py     — kcat / kcat/Km computation, KineticConstants dataclass
    ├── analyze.py    — KineticAnalyzer orchestrator, analyze_peaks batch helper
    └── plots.py      — KineticPlots: MM curves, LB plots, residuals, efficiency bar chart
```

---

## Data Flow

```
CSV (PeakAnalyzer output)
  └─ load_plottable_data()          # io.py  — loads, sanitizes columns to snake_case
       └─ canonicalize_peak_ids()   # compounds.py  — resolves peak IDs to Compounds enum,
            |                         produces ordered pandas.Categorical
            └─ KineticAnalyzer.fit()
                 └─ extract_peak_data()      # preprocess.py
                      └─ prepare_velocity()  # area / rxn_time → v (area/s)
                 └─ _select_and_fit()        # models.py  — MM or Hill fit → FitResult
                 └─ fit_lineweaver_burk()    # models.py  — cross-validation fit
                 └─ derive_constants()       # derive.py  — kcat, kcat/Km → KineticConstants

  [optional] KineticAnalyzer.apply_calibration(Calibration)
                 └─ area_to_conc()           # calibration.py  — area → µM
                 └─ v = µM / rxn_time        # velocity now µM/s
                 └─ re-fit → KineticConstants with physical units

  └─ KineticAnalyzer.export_csv()   # KineticConstants.to_dict() → DataFrame → CSV
  └─ KineticAnalyzer.plot()         # plots.py  — saves PNG figures
```

---

## Units of Measure — Comprehensive Reference

### Input / DataFrame columns

| Column | Unit | Notes |
|---|---|---|
| `substrate_conc` | µM | **Critical**: must be µM throughout; no runtime check performed |
| `mean` | area (dimensionless HPLC counts) | Integrated peak area; instrument-specific |
| `std` | area | Standard deviation of peak area across replicates |
| `count` | integer (replicates) | Number of replicates per group |

### Calibration (`Calibration` dataclass)

| Field | Unit | Notes |
|---|---|---|
| `slope` | area / µM | Calibration sensitivity |
| `slope_std` | area / µM | Uncertainty on slope |
| `intercept` | area | y-intercept at [CoA] = 0 |
| `intercept_std` | area | Uncertainty on intercept |
| `r_squared` | dimensionless | ≥ 0.999 expected for HPLC CoA |
| `rmse` | area | Root mean square error of the fit |
| `conc_range` | µM, µM | (min, max) of calibration standards |

### Velocity

| Variable | Before calibration | After `apply_calibration()` |
|---|---|---|
| `v` (velocity) | area/s | µM/s |
| `v_std` | area/s | µM/s |
| `vmax` (`FitResult.vmax`) | area/s | µM/s |

### Kinetic model parameters (`FitResult`)

| Parameter | Unit | Notes |
|---|---|---|
| `km` / `k_half` | µM | Michaelis constant or Hill half-saturation constant |
| `km_std` / `k_half_std` | µM | |
| `ki` | µM | Substrate inhibition constant (SI model only) |
| `ki_std` | µM | |
| `hill_n` | dimensionless | Hill cooperativity coefficient |
| `r_squared` | dimensionless | Coefficient of determination |
| `pcov` | (vmax units)² / µM² | Full covariance matrix; used for delta-method error propagation |

### Derived constants (`KineticConstants`)

| Field | Unit | Notes |
|---|---|---|
| `kcat` | s⁻¹ (post-calibration) | Non-physical before calibration (area · µM⁻¹ · s⁻¹) |
| `kcat_std` | s⁻¹ (post-calibration) | |
| `kcat_km_M` | s⁻¹ · M⁻¹ | Km is converted µM → M (×1e-6) before computing ratio |
| `kcat_km_std_M` | s⁻¹ · M⁻¹ | Propagated via full Vmax–Km covariance matrix |
| `ki` | µM | Forwarded from FitResult; None unless model is "si" |
| `ki_std` | µM | |

### Time

| Variable | Unit |
|---|---|
| `rxn_time` (CLI arg) | minutes (converted to seconds via `reaction_time_seconds`) |
| `reaction_time_seconds` / `rxn_time` (internal) | seconds |
| `enzyme_conc_um` | µM |

### Lineweaver-Burk reciprocal space

| Variable | Unit |
|---|---|
| `s_inv` (1/[S]) | µM⁻¹ |
| `v_inv` (1/V) | s/area (before cal) or µM⁻¹ · s (after cal) |

---

## Key Design Decisions

### Two-phase velocity: area/s → µM/s

Before calibration, all velocities are in **area/s** (dimensionless HPLC counts per
second). `apply_calibration()` converts using:

```
v_µM_per_s = Calibration.area_to_conc(mean_area) / reaction_time_seconds
```

This avoids a lossy area ↔ velocity round-trip; raw `mean_signal` (area) is explicitly
returned as the fourth element from `extract_peak_data()` / `prepare_velocity()`.

### kcat/Km unit conversion

`_kcat_km_with_covariance()` in `derive.py` converts Km from µM to M (×1e-6) before
computing the ratio, so that `kcat_km_M` is always in the standard literature unit of
s⁻¹ · M⁻¹ — even though all internal concentration math uses µM. The partial
derivatives are scaled by the same factor.

### Covariance-aware error propagation

`pcov` (the full parameter covariance matrix from `scipy.optimize.curve_fit`) is stored
on every `FitResult`. Vmax and Km are anti-correlated in MM fitting; summing variances
separately underestimates the true uncertainty in kcat/Km. The delta method with the
off-diagonal covariance term is applied in `_kcat_km_with_covariance`.

### Calibration validity guard

`Calibration.is_valid(r2_threshold=0.999)` checks `slope > 0` and `R² ≥ threshold`.
`apply_calibration()` logs a warning if the calibration fails but proceeds; the caller
must inspect kcat/Km quality after a failed-validity calibration.

### Compounds enum and categorical ordering

`Compounds` (a `StrEnum`) is the authoritative identity for each tracked analyte
(HTAL, PDAL, OA, OLV). `canonicalize_peak_ids()` resolves raw HPLC string labels to
`Compounds` values and produces an **ordered** `pandas.Categorical` so that downstream
grouping and plotting respect the defined compound order. Chain-length prefixes (C3–C7)
are optionally interleaved.

---

## Docstring Style

**Convention**: Google style (`ruff.lint.pydocstyle.convention = "google"`), enforced
by ruff `D` rules.

**Level of detail applied**: Comprehensive. All public functions and classes carry:

- A one-line summary sentence.
- An extended description when the function has non-obvious behavior (e.g. unit
  transformations, numeric edge cases, covariance propagation, dead-zone thresholds,
  anti-correlation effects, calibration validity semantics).
- `Args:` section with per-parameter descriptions that **always include units** for
  every physically meaningful quantity (µM, s, area, area/µM, s⁻¹, s⁻¹·M⁻¹).
- `Returns:` section specifying units of returned quantities.
- `Raises:` section for `ValueError` / `RuntimeError` paths.
- `Note:` sections for unit-consistency warnings (e.g. "substrate_conc must be in µM;
  passing other units will silently produce incorrect results").
- `Attributes:` sections on all dataclasses and Pydantic models, each attribute
  annotated with units where applicable.
- Module-level docstrings include a `Typical usage:` code snippet.

**Private helpers** (`_func`) have minimal docstrings — one summary line is sufficient
unless the logic is non-obvious.

**No narrative comments** describing what the code does; inline comments are used only
for the WHY (e.g. `# clamp numerical negatives`, `# FIX 2: intercept term`).

---

## Analytes Tracked

| Enum member | Full name |
|---|---|
| `HTAL` | Hexanoyl triacetic acid lactone |
| `PDAL` | Pentyl diacetic acid lactone |
| `OA` | Olivetolic acid |
| `OLV` | Olivetol |

---

## Kinetic Models Available

| Tag | Function | Parameters |
|---|---|---|
| `"mm"` | `michaelis_menten` | Vmax, Km |
| `"hill"` | `hill_equation` | Vmax, k_half, n |
| `"tmm"` | `threshold_michaelis_menten` | Vmax, Km, S0 |
| `"si"` | `substrate_inhibition` | Vmax, Km, Ki |
| `"lb"` | `fit_lineweaver_burk` (linear) | Vmax, Km (back-calc) |

Default model for all peaks is `"mm"`. Pass `special_peaks={"peak_id": "hill"}` to
`KineticAnalyzer` to override per-peak.
