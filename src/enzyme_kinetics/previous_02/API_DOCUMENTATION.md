# Enzyme Kinetics Analysis Package - API Documentation

## Overview

A production-quality Python package for analyzing enzyme kinetics from HPLC data. Provides modular, reusable classes for fitting Michaelis-Menten kinetics, calibrating peak area measurements, and generating publication-quality visualizations.

---

## Installation

```bash
# Add the package to your Python path
import sys
sys.path.insert(0, '/path/to/enzyme_kinetics')

# Or install in site-packages
python setup.py install
```

---

## Core Classes

### Models (`enzyme_kinetics.models`)

#### `EnzymeKineticModel` (Base Class)

Abstract base class for all kinetic models.

**Methods:**
- `__call__(substrate_conc, **params) → np.ndarray`: Evaluate model
- `fit(substrate_conc, velocity, velocity_std=None, **params) → (KineticParameters, dict)`: Fit model to data
- `calculate_r2(actual, predicted) → float`: Coefficient of determination

#### `MichaelisMentenModel(EnzymeKineticModel)`

Non-linear Michaelis-Menten kinetics: `v = (Vmax × [S]) / (Km + [S])`

**Methods:**
```python
fit(
    substrate_conc: np.ndarray,
    velocity: np.ndarray,
    velocity_std: np.ndarray | None = None,
    Vmax_init: float | None = None,
    Km_init: float | None = None,
    maxfev: int = 5000
) → (KineticParameters, dict)
```

**Returns:**
- `KineticParameters`: Fitted Km, Vmax with uncertainties and R²
- `dict`: Fit diagnostics (predictions, residuals, covariance)

**Example:**
```python
from enzyme_kinetics import MichaelisMentenModel

model = MichaelisMentenModel()
substrate = np.array([0.1, 0.2, 0.5, 1.0, 3.0])
velocity = np.array([70.97, 178.53, 297.77, 282.58, 1026.30])

params, fit_data = model.fit(substrate, velocity)
print(f"Km = {params.Km:.4f} ± {params.Km_std:.4f} µM")
print(f"Vmax = {params.Vmax:.2f} ± {params.Vmax_std:.2f}")
print(f"R² = {params.r2:.4f}")
```

#### `LineweaverBurkModel(EnzymeKineticModel)`

Lineweaver-Burk linearization: `1/v = (Km/Vmax) × (1/[S]) + 1/Vmax`

**Methods:**
```python
fit(
    substrate_conc: np.ndarray,
    velocity: np.ndarray
) → (KineticParameters, dict)
```

**Note:** More robust for visualization but can introduce bias. Use with `MichaelisMentenModel` for comparison.

#### `SubstrateInhibitionModel(EnzymeKineticModel)`

Substrate inhibition: `v = (Vmax × [S]) / (Km + [S] + [S]²/Ki)`

Use when enzyme activity decreases at high substrate concentrations.

**Methods:**
```python
fit(
    substrate_conc: np.ndarray,
    velocity: np.ndarray,
    velocity_std: np.ndarray | None = None,
    Vmax_init: float | None = None,
    Km_init: float | None = None,
    Ki_init: float | None = None,
    maxfev: int = 5000
) → (KineticParameters, dict)
```

### `KineticParameters` (Data Class)

Container for kinetic parameters and uncertainties.

**Attributes:**
- `Km: float` - Michaelis constant (µM)
- `Km_std: float` - Standard deviation
- `Vmax: float` - Maximum velocity
- `Vmax_std: float` - Standard deviation
- `r2: float` - Coefficient of determination (0–1)
- `n_points: int` - Number of data points

**Methods:**
- `__str__() → str`: Formatted string representation

---

### Calibration (`enzyme_kinetics.calibration`)

#### `CalibrationCurve`

HPLC calibration for peak area → product concentration conversion.

**Methods:**

```python
fit(
    concentrations: np.ndarray,
    peak_areas: np.ndarray,
    peak_areas_std: np.ndarray | None = None
) → CalibrationParameters
```

Fit linear calibration: `peak_area = slope × [concentration] + intercept`

```python
area_to_concentration(peak_area: float | np.ndarray) → float | np.ndarray
```

Convert peak area to concentration (µM).

```python
area_to_concentration_with_error(
    peak_area: float,
    peak_area_std: float
) → (float, float)
```

Convert with error propagation.

```python
save(filepath: str | Path) → None
```

Save calibration parameters to text file.

```python
plot(output_path: str | Path | None = None) → None
```

Generate 4-panel calibration validation plot.

**Example:**
```python
from enzyme_kinetics import CalibrationCurve

calib = CalibrationCurve()

# Your CoA standard measurements
conc = np.array([0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0])
area = np.array([0, 245.3, 612.4, 1223.1, 2451.8, 4903.5, 12258.7])
area_std = np.array([0, 12.4, 28.1, 45.2, 89.3, 156.2, 412.1])

params = calib.fit(conc, area, peak_areas_std=area_std)

# Convert HPLC areas to product concentration
product_conc = calib.area_to_concentration(peak_area=612.4)  # → 0.25 µM

# With error propagation
conc, conc_std = calib.area_to_concentration_with_error(612.4, 28.1)

# Validate
if calib.is_valid():
    print("Calibration is valid")

# Save for later use
calib.save('/path/to/calibration.txt')

# Generate plots
calib.plot(output_path='/path/to/calibration_plot.png')
```

#### `CalibrationParameters` (Data Class)

Container for calibration parameters.

**Attributes:**
- `slope: float` - area/µM
- `slope_std: float`
- `intercept: float` - area
- `intercept_std: float`
- `r2: float` - R² goodness of fit
- `rmse: float` - Root mean square error
- `n_points: int` - Number of calibration points
- `conc_range: Tuple[float, float]` - Min/max concentrations

**Methods:**
- `area_to_concentration(area) → float | np.ndarray`
- `concentration_to_area(conc) → float | np.ndarray`

---

### Analysis Engine (`enzyme_kinetics.analyzer`)

#### `EnzymeKineticsAnalyzer`

High-level analysis engine for batch processing.

**Initialization:**
```python
analyzer = EnzymeKineticsAnalyzer(
    enzyme_concentration_um=12.25,      # µM, optional
    reaction_time_seconds=10800,        # seconds, optional
    calibration=None                    # CalibrationCurve, optional
)
```

**Methods:**

```python
analyze_peak(
    peak_id: str,
    substrate_conc: np.ndarray,
    velocity: np.ndarray,
    velocity_std: np.ndarray | None = None,
    fit_lineweaver_burk: bool = True
) → EnzymeKineticsResult
```

Analyze kinetics for a single peak/enzyme.

```python
analyze_dataframe(
    df: pd.DataFrame,
    peak_id_column: str = 'peak_id',
    substrate_conc_column: str = 'substrate_conc',
    velocity_column: str = 'velocity',
    velocity_std_column: str | None = None,
    fit_lineweaver_burk: bool = True,
    min_points: int = 3
) → List[EnzymeKineticsResult]
```

Batch analyze kinetics from DataFrame. Automatically groups by peak_id.

```python
apply_calibration(calibration: CalibrationCurve) → None
```

Apply calibration curve to all results, converting peak areas to absolute concentrations.

```python
results_to_dataframe() → pd.DataFrame
```

Export all results to DataFrame.

```python
save_results(csv_path: str | Path) → None
```

Save results to CSV file.

**Example:**
```python
from enzyme_kinetics import EnzymeKineticsAnalyzer
import pandas as pd

# Load HPLC data
df = pd.read_csv('hplc_kinetics.csv')

# Initialize analyzer
analyzer = EnzymeKineticsAnalyzer(
    enzyme_concentration_um=12.25,
    reaction_time_seconds=10800  # 3 hours
)

# Batch analyze
results = analyzer.analyze_dataframe(
    df,
    peak_id_column='peak_id',
    substrate_conc_column='substrate_conc',
    velocity_column='peak_area',
    velocity_std_column='peak_area_std'
)

# Export results
df_results = analyzer.results_to_dataframe()
df_results.to_csv('kinetics_results.csv', index=False)

# Apply calibration
from enzyme_kinetics import CalibrationCurve
calib = CalibrationCurve()
calib.fit(calib_conc, calib_area)
analyzer.apply_calibration(calib)

# Re-export with absolute units
analyzer.save_results('kinetics_calibrated.csv')
```

#### `EnzymeKineticsResult` (Data Class)

Complete result for single enzyme/peak.

**Attributes:**
- `peak_id: str`
- `substrate_conc: np.ndarray`
- `velocity: np.ndarray`
- `mm_params: KineticParameters` - Michaelis-Menten fit
- `lb_params: KineticParameters | None` - Lineweaver-Burk fit
- `kcat: float | None` - Turnover number (if enzyme conc provided)
- `kcat_over_km: float | None` - Catalytic efficiency
- Plus uncertainties and diagnostic information

**Methods:**
- `to_dict() → dict`: Convert to dictionary for DataFrame

---

### Plotting (`enzyme_kinetics.plotting`)

#### `KineticsPlotter`

Static methods for generating publication-quality plots.

**Methods:**

```python
plot_michaelis_menten_curve(
    result: EnzymeKineticsResult,
    ax: plt.Axes | None = None,
    show_fit: bool = True,
    show_annotations: bool = True
) → plt.Axes
```

Plot Michaelis-Menten curve with data points, fitted curve, and Km annotation.

```python
plot_lineweaver_burk(
    result: EnzymeKineticsResult,
    ax: plt.Axes | None = None,
    show_fit: bool = True
) → plt.Axes
```

Plot Lineweaver-Burk linearization.

```python
plot_residuals(
    result: EnzymeKineticsResult,
    ax: plt.Axes | None = None
) → plt.Axes
```

Plot residuals from MM fit.

```python
plot_comparison_panel(
    results: List[EnzymeKineticsResult],
    output_path: str | Path | None = None
) → None
```

Generate multi-panel figure comparing all peaks.

```python
plot_efficiency_comparison(
    results: List[EnzymeKineticsResult],
    output_path: str | Path | None = None
) → None
```

Bar plots comparing Km, Vmax, kcat, and kcat/Km across peaks.

**Example:**
```python
from enzyme_kinetics import KineticsPlotter

# Plot individual peak
KineticsPlotter.plot_michaelis_menten_curve(
    result,
    output_path='peak_analysis.png'
)

# Multi-panel comparison
KineticsPlotter.plot_comparison_panel(
    results,
    output_path='all_peaks_comparison.png'
)

# Efficiency comparison
KineticsPlotter.plot_efficiency_comparison(
    results,
    output_path='efficiency_bars.png'
)
```

---

## Workflow Example

### Complete Analysis Pipeline

```python
import pandas as pd
import numpy as np
from enzyme_kinetics import (
    EnzymeKineticsAnalyzer,
    CalibrationCurve,
    KineticsPlotter
)

# 1. Load data
df_kinetics = pd.read_csv('enzyme_kinetics.csv')
df_calib_standards = pd.read_csv('calibration_standards.csv')

# 2. Create and fit calibration
calib = CalibrationCurve()
calib.fit(
    df_calib_standards['concentration'].values,
    df_calib_standards['peak_area'].values,
    peak_areas_std=df_calib_standards['peak_area_std'].values
)
calib.save('calibration.txt')
calib.plot('calibration_curve.png')

# 3. Analyze kinetics
analyzer = EnzymeKineticsAnalyzer(
    enzyme_concentration_um=12.25,
    reaction_time_seconds=10800
)

results = analyzer.analyze_dataframe(
    df_kinetics,
    peak_id_column='peak',
    substrate_conc_column='substrate_um',
    velocity_column='peak_area',
    velocity_std_column='peak_area_std'
)

# 4. Apply calibration for absolute units
analyzer.apply_calibration(calib)

# 5. Export and visualize
analyzer.save_results('results_calibrated.csv')
KineticsPlotter.plot_comparison_panel(results, 'curves.png')
KineticsPlotter.plot_efficiency_comparison(results, 'efficiency.png')

# 6. Summary statistics
df_results = analyzer.results_to_dataframe()
print(df_results[['peak_id', 'Km_MM', 'Vmax_MM', 'kcat_over_km']].to_string())
```

---

## Common Patterns

### Pattern 1: Analyze Single Enzyme

```python
analyzer = EnzymeKineticsAnalyzer()
result = analyzer.analyze_peak(
    'MyEnzyme',
    substrate_conc,
    velocity,
    velocity_std=velocity_std
)
print(f"Km = {result.mm_params.Km:.4f} µM")
print(f"R² = {result.mm_params.r2:.4f}")
```

### Pattern 2: Batch Analyze from CSV

```python
analyzer = EnzymeKineticsAnalyzer(enzyme_concentration_um=10)
df = pd.read_csv('data.csv')
results = analyzer.analyze_dataframe(df, 'peak_id', 'substrate_conc', 'velocity')
analyzer.save_results('results.csv')
```

### Pattern 3: Compare Multiple Models

```python
from enzyme_kinetics import MichaelisMentenModel, SubstrateInhibitionModel

mm_model = MichaelisMentenModel()
mm_params, _ = mm_model.fit(substrate_conc, velocity)

si_model = SubstrateInhibitionModel()
si_params, _ = si_model.fit(substrate_conc, velocity)

print(f"MM R² = {mm_params.r2:.4f}")
print(f"SI R² = {si_params.r2:.4f}")
print(f"Better fit: {'Substrate Inhibition' if si_params.r2 > mm_params.r2 else 'Michaelis-Menten'}")
```

---

## Error Handling

All fitting functions raise exceptions on convergence failure. Handle gracefully:

```python
try:
    result = analyzer.analyze_peak(peak_id, substrate, velocity)
except Exception as e:
    print(f"Failed to analyze {peak_id}: {e}")
    continue
```

The `analyze_dataframe()` method already handles this internally and prints status.

---

## Performance Considerations

- **Memory:** Loading large DataFrames (>100k rows) may be slow. Consider pre-filtering.
- **Fitting:** Non-linear fitting can be slow with poor initial guesses. Provide `Vmax_init` and `Km_init` if known.
- **Plotting:** Generating many plots is slow. Use `KineticsPlotter` methods for vectorized operations.

---

## References

- Michaelis, L., & Menten, M. L. (1913). Biochem. Z.
- Lineweaver, H., & Burk, D. (1934). J. Am. Chem. Soc.
- Berg, J. M., Tymoczko, J. L., & Stryer, L. (2002). Biochemistry (5th ed.)

---

## License & Attribution

This package is provided as-is for research purposes.
