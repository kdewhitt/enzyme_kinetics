# Enzyme Kinetics Project — Overview & Analysis

Generated: 2026-05-10  
Scope: `src/enzyme_kinetics/previous_02/` — four core modules plus examples  

---

## 1. Package Architecture

### Module Dependency Graph

```
models.py
    ↑ imported by
calibration.py      (standalone; no dependency on models.py)
    ↑
    ├─── analyzer.py     (imports MichaelisMentenModel, LineweauerBurkModel, KineticParameters
    │                     from models; CalibrationCurve, CalibrationParameters from calibration)
    └─── plotting.py     (imports EnzymeKineticsResult, EnzymeKineticsAnalyzer from analyzer)

__init__.py          (empty — no public re-exports defined yet)
```

**Key observation:** `plotting.py` imports from `analyzer.py`, which means `KineticsPlotter`
is tightly coupled to `EnzymeKineticsResult`. This is intentional (plotting operates on
result objects), but it means `plotting` cannot be imported without pulling in the full
`analyzer` + `models` + `calibration` chain.

**Key observation:** `__init__.py` is empty. The API documentation describes public exports
but they are not implemented. Any `from enzyme_kinetics import EnzymeKineticsAnalyzer`
call will fail until `__init__.py` is populated.

---

## 2. Module-by-Module Summary

### 2.1 `models.py` (~400 lines)

**Purpose:** Core kinetic equations and non-linear fitting via `scipy.optimize.curve_fit`.

#### Classes

| Class | Role |
|---|---|
| `KineticParameters` | Frozen dataclass holding Km, Vmax, uncertainties, R², n_points |
| `EnzymeKineticModel` | Abstract base — defines `__call__`, `fit`, `calculate_r2` interface |
| `MichaelisMentenModel` | Non-linear MM fit; Levenberg-Marquardt via `curve_fit` |
| `LineweauerBurkModel` | Double-reciprocal linearization via `scipy.stats.linregress` |
| `SubstrateInhibitionModel` | Extended MM with competitive substrate inhibition term (Ki) |

#### Scientific details

- **MM fitting:** Uses `curve_fit` with `sigma=velocity_std` and `absolute_sigma=True`,
  so uncertainties are propagated correctly into `pcov` when replicate standard deviations
  are supplied. Initial guesses default to `max(velocity)` for Vmax and
  `median(substrate_conc)` for Km — sensible heuristics for typical substrate ranges.

- **Lineweaver-Burk:** Transforms to 1/v vs 1/[S] space, fits by OLS. Km and Vmax are
  back-calculated from slope and intercept. `Km_std` and `Vmax_std` are set to `nan`
  (standard errors in reciprocal space do not map cleanly to parameter space without the
  delta method). This is correct behavior but should be documented explicitly.

- **Substrate inhibition:** Three-parameter model `v = (Vmax·[S]) / (Km + [S] + [S]²/Ki)`.
  Ki (inhibition constant) is fitted as a third free parameter. The `KineticParameters`
  dataclass does not have a Ki field — Ki is stored only in `fit_data`. This is a design
  gap: there is no type-safe way to retrieve Ki from a `SubstrateInhibitionModel` result
  without inspecting the dict.

#### Docstring quality issues (pre-pipeline)

- `EnzymeKineticModel.__call__` and `fit` are abstract stubs that raise
  `NotImplementedError` but are documented as if they have real behavior.
- `LineweauerBurkModel.fit` does not document that `Km_std` and `Vmax_std` are always
  `nan` in the returned `KineticParameters`.
- `SubstrateInhibitionModel.fit` documents "note: Ki stored in n_points field" in a
  misleading way — Ki is in `fit_data`, not in `n_points`.
- `calculate_r2` is a simple 4-line method and should have a single-line docstring only.
- Module docstring is domain-accurate but does not mention the base class pattern or
  that all subclasses share the `fit → (KineticParameters, dict)` contract.

---

### 2.2 `calibration.py` (~450 lines)

**Purpose:** Fit a linear HPLC calibration curve (peak area vs. CoA concentration) and
provide area→concentration conversion with analytical error propagation.

#### Classes

| Class | Role |
|---|---|
| `CalibrationParameters` | Dataclass: slope, intercept, uncertainties, R², RMSE, n_points, conc_range |
| `CalibrationCurve` | Fits linear model; exposes conversion, saving, and validation methods |

#### Scientific details

- **Model:** `peak_area = slope × [CoA] + intercept`. Linear response of UV/fluorescence
  detector to CoA-SH concentration, as expected for HPLC with electrochemical or UV
  detection in the working range.

- **Fitting:** Uses `curve_fit` with optional `sigma` weighting. Equivalent to weighted
  least squares when standard deviations are provided; reduces to OLS otherwise.

- **Error propagation in `area_to_concentration_with_error`:** Uses the delta method:
  ```
  σ_C² = (σ_area / slope)² + ((area − intercept) · σ_slope / slope²)²
  ```
  Intercept uncertainty is not included in the propagation, which is a minor omission
  (intercept_std is available in `CalibrationParameters` but unused in the formula).
  For typical calibration curves where intercept ≈ 0 and intercept_std is small, the
  error is negligible.

- **Validation:** `is_valid()` requires R² > 0.95 and slope > 0. This is a reasonable
  threshold for HPLC calibration; typical CoA calibrations achieve R² > 0.999.

- **Plot:** 4-panel figure: (1) calibration curve + fit, (2) absolute residuals,
  (3) standardized residuals with ±3σ bounds, (4) parameter summary table.

#### Docstring quality issues (pre-pipeline)

- `CalibrationParameters.__str__` should be a single-line docstring.
- `CalibrationCurve.fit` correctly documents `sigma` / `absolute_sigma` behavior but
  does not mention that `curve_fit` is used internally (the method name `fit` is
  consistent with the models API).
- `area_to_concentration` raises `ValueError` if not fitted — documented in the Raises
  section, which is correct.
- `area_to_concentration_with_error` documents the error formula as a comment in the
  body but omits it from the docstring itself, where it would be most useful.
- `is_valid` has no docstring.
- `save` and `plot` docstrings are accurate and appropriately detailed.
- Module docstring should explicitly state that the calibration model assumes a linear
  detector response (Beer-Lambert regime).

---

### 2.3 `analyzer.py` (~350 lines)

**Purpose:** Orchestration layer. Wraps the fitting models and calibration into a
high-level API for per-peak and batch DataFrame analysis.

#### Classes

| Class | Role |
|---|---|
| `EnzymeKineticsResult` | Dataclass aggregating MM params, LB params, kcat, kcat/Km, raw arrays |
| `EnzymeKineticsAnalyzer` | Stateful analyzer: fits peaks, holds results list, applies calibration |

#### Scientific details

- **kcat calculation:** `kcat = (Vmax / t_reaction) / [E]` where Vmax is in raw area
  units. This is dimensionally consistent only if velocity is already in concentration/time
  units (i.e., after calibration). When using raw peak areas, kcat has units of
  `area · µM⁻¹ · s⁻¹` — not the standard `s⁻¹`. The docstrings and API documentation
  describe kcat as having units `s⁻¹` without noting this caveat, which is a significant
  accuracy issue for users who have not applied a calibration curve.

- **kcat/Km propagation:** Uses the standard relative-error formula:
  ```
  σ(kcat/Km) / (kcat/Km) = sqrt((σ_kcat/kcat)² + (σ_Km/Km)²)
  ```
  This is correct under the assumption that Km and Vmax are uncorrelated, which is
  typically not true for MM fitting (they are anti-correlated). For high-quality data
  (R² > 0.95) the correlation is small, but for noisy data it can lead to underestimated
  uncertainties. A more rigorous approach would use the full covariance matrix.

- **`apply_calibration`:** Re-fits MM model on calibrated velocities in-place, mutating
  `EnzymeKineticsResult` objects stored in `self.results`. This is a side-effecting
  design — calling it twice would re-calibrate already-calibrated data, leading to
  incorrect results. There is no guard against this.

- **LB error suppression:** `analyze_peak` silently catches all exceptions during LB
  fitting (`except Exception: pass`). This means fitting failures are invisible.

- **Batch analysis:** `analyze_dataframe` groups by `peak_id_column`, sorts by substrate
  concentration, and calls `analyze_peak` per group. Results are stored in both the
  returned list and `self.results` — callers must be aware that the analyzer accumulates
  state across multiple `analyze_dataframe` calls.

#### Docstring quality issues (pre-pipeline)

- `EnzymeKineticsResult` is a dataclass — attributes should be documented in an
  `Attributes:` section in the class docstring, not repeated in each method.
- `EnzymeKineticsResult.to_dict` docstring is accurate.
- `EnzymeKineticsAnalyzer.__init__` attributes (`mm_model`, `lb_model`, `results`)
  are not documented.
- `analyze_peak` does not document that results are appended to `self.results` as a
  side effect.
- `apply_calibration` does not warn that it is not idempotent.
- `analyze_dataframe` correctly documents all parameters and the `min_points` guard.
- `results_to_dataframe` and `save_results` are trivially short and should have
  single-line docstrings only.
- Module docstring mentions "different models" but only MM and LB are actually
  dispatched; SI model is accessible through `models.py` directly, not through the
  analyzer.

---

### 2.4 `plotting.py` (~350 lines)

**Purpose:** Static methods for publication-quality figures using `matplotlib`.

#### Classes

| Class | Role |
|---|---|
| `KineticsPlotter` | Namespace of static plotting methods — not intended to be instantiated |

#### Method inventory

| Method | Output |
|---|---|
| `plot_michaelis_menten_curve` | Single MM curve with error bars, fit line, Km/Vmax annotations |
| `plot_lineweaver_burk` | Double-reciprocal plot with regression line |
| `plot_residuals` | Residuals from MM fit vs. substrate concentration |
| `plot_comparison_panel` | Multi-panel grid of MM curves + LB plots for all peaks |
| `plot_efficiency_comparison` | 4-panel bar chart: Km, Vmax, kcat, kcat/Km |

#### Scientific details

- MM curve plots annotate Km with a vertical dashed line and Vmax/2 with a horizontal
  dashed line — the standard graphical definition of Km. This is correct.
- LB plot correctly uses 1/[S] on x-axis and 1/v on y-axis. Fit line is drawn only over
  the observed 1/[S] range (no extrapolation to y-intercept visible on plot), which is
  a minor presentational limitation — the y-intercept (= 1/Vmax) and x-intercept
  (= -1/Km) are informative and commonly shown in publications.
- `plot_efficiency_comparison` silently skips results where `kcat_over_km is None`
  (i.e., when enzyme concentration was not provided or calibration not applied).
  The kcat panel also silently skips None values, but the bar index is not adjusted,
  which could cause a length mismatch between `peaks` and the kcat list if some results
  have None kcat.

#### Docstring quality issues (pre-pipeline)

- `KineticsPlotter` class docstring is missing entirely.
- All static methods have docstrings but none include `Examples:`.
- `plot_michaelis_menten_curve` return annotation `→ plt.Axes` is documented but the
  method also creates a figure when `ax=None` — this side effect (figure creation) is
  not mentioned.
- `plot_comparison_panel` and `plot_efficiency_comparison` return `None` but their
  docstrings have no `Returns:` section (acceptable for `None`-returning functions, but
  the `output_path` save side effect should be documented).

---

## 3. Cross-Cutting Issues

### 3.1 Typo in class name
`LineweauerBurkModel` (note: "Lineweauer") is a consistent misspelling of
"Lineweaver-Burk" across all modules and documentation. This is a public API name and
should be corrected to `LineweaverBurkModel` with a deprecation alias.

### 3.2 `__init__.py` is empty
All `from enzyme_kinetics import X` calls in examples and documentation will fail.
The `__init__.py` must be populated with the exports listed in `API_DOCUMENTATION.md`.

### 3.3 Units and calibration dependency of kcat

`kcat` and `kcat/Km` only have physically meaningful units (s⁻¹ and s⁻¹·M⁻¹
respectively) after a calibration curve has been applied. When operating on raw peak
areas, these values are in non-standard units. This distinction is critical for
comparing results with literature values and should be prominently documented
in `EnzymeKineticsAnalyzer` and `EnzymeKineticsResult`.

### 3.4 `apply_calibration` is not idempotent
Calling `apply_calibration` a second time will re-calibrate already-calibrated velocity
values, producing incorrect kinetic parameters. A `_calibration_applied: bool` guard
should be added.

### 3.5 `SubstrateInhibitionModel` Ki not in `KineticParameters`
Ki and Ki_std live only in `fit_data`, a plain dict. There is no type-safe way to
access them. Either `KineticParameters` needs optional Ki/Ki_std fields, or a
`SubstrateInhibitionParameters` subclass should be introduced.

### 3.6 LB fitting failures are silently swallowed
The bare `except Exception: pass` in `analyze_peak` makes debugging difficult.
At minimum, the exception should be logged.

### 3.7 State accumulation in `EnzymeKineticsAnalyzer`
`self.results` grows with every call to `analyze_peak` or `analyze_dataframe`.
Multiple calls on different datasets will mix results. Users who reuse an analyzer
instance across datasets will get incorrect `results_to_dataframe()` output.

---

## 4. Terminology Reference (for docstring authoring)

The following domain-specific terms should be used consistently in docstrings:

| Term                           | Definition / Context                                                                                                                 |
|--------------------------------|--------------------------------------------------------------------------------------------------------------------------------------|
| Km (Michaelis constant)        | Substrate concentration at half-maximal velocity; units µM; reflects apparent substrate affinity                                     |
| Vmax                           | Maximum reaction velocity at saturating substrate; units depend on whether calibration applied                                       |
| kcat (turnover number)         | Catalytic rate constant; Vmax / [E]total; units s⁻¹ (after calibration)                                                              |
| kcat/Km (catalytic efficiency) | Second-order rate constant for substrate capture; units s⁻¹·M⁻¹; the standard metric for comparing enzyme activity across substrates |
| Ki (inhibition constant)       | Substrate concentration at which inhibition halves the rate in the SI model                                                          |
| Lineweaver-Burk plot           | Double-reciprocal linearization of MM equation; susceptible to error amplification at low [S]                                        |
| CoA / CoA-SH                   | Coenzyme A (free thiol form); the product measured by HPLC in acyltransferase assays                                                 |
| HexCoA                         | Hexanoyl-CoA; acyl-CoA substrate used in the assay represented in the example data                                                   |
| HPLC peak area                 | Integrated chromatographic signal; proportional to analyte mass injected                                                             |
| Calibration curve              | Linear regression of peak area vs. known CoA concentration; converts raw HPLC signal to µM                                           |
| R²                             | Coefficient of determination; ≥ 0.95 considered acceptable for MM fitting; ≥ 0.999 expected for HPLC calibration                     |
| Levenberg-Marquardt            | Non-linear least-squares algorithm used by `scipy.optimize.curve_fit`; robust to poor initial guesses                                |
| Error propagation              | Delta-method propagation of measurement uncertainty through derived quantities (kcat, kcat/Km, concentration)                        |
| Reaction velocity (v)          | Rate of product formation per unit time; in raw area units unless calibrated                                                         |
| Substrate inhibition           | Kinetic phenomenon where excess [S] reduces v, observable as a bell-shaped v vs. [S] curve                                           |

---

## 5. Docstring Pipeline Pre-Assessment

Summary of issues across all four modules, classified by the pipeline audit codes:

| Code | Description | Count (approx.) |
|---|---|---|
| MISSING | No docstring at all | 4 (`CalibrationCurve.is_valid`, `KineticsPlotter` class, `__init__` attrs) |
| STALE | Signature/body has changed since docstring written | 3 (`SubstrateInhibitionModel.fit` Ki note; `analyze_peak` side-effect; LB std=nan) |
| INACCURATE | Contradicts actual behavior | 2 (kcat units claim; `apply_calibration` idempotency) |
| STYLE | Violates Google-style rules (wrong mood, wrong section, etc.) | ~12 (multi-line where single-line required; body-duplicating Args; no `Examples:`) |

**Priority order for pipeline rewrites:**
1. `models.py` — most referenced; typo in `LineweauerBurkModel` class name note
2. `analyzer.py` — most complex; most inaccurate docstrings re: units and side effects
3. `calibration.py` — mostly accurate; minor omissions
4. `plotting.py` — missing class docstring; otherwise style issues only

---

## 6. Recommended Next Steps

1. **Populate `__init__.py`** — package is currently unimportable by its documented API.
2. **Fix the `LineweauerBurkModel` typo** → `LineweaverBurkModel` with a backward-compatible alias.
3. **Add Ki / Ki_std to `KineticParameters`** (as optional fields) so SI results are first-class.
4. **Guard `apply_calibration` against double application** with an idempotency flag.
5. **Run the docstring pipeline** on all four `.py` files to fix MISSING / STALE / INACCURATE items and enforce Google-style uniformly — particularly the domain-specific language from §4 above.
6. **Add unit test stubs** for: MM fitting on known analytical data, calibration round-trip (area → conc → area), and kcat calculation with and without calibration applied.
