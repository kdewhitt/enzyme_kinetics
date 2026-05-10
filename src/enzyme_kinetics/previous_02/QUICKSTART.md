# Enzyme Kinetics Package - Quick Start Guide

## 5-Minute Setup

### 1. Installation

```bash
# Option A: Add to your Python path
import sys
sys.path.insert(0, '/mnt/user-data/outputs')

# Option B: Copy to site-packages
cp -r enzyme_kinetics /path/to/site-packages/
```

### 2. Basic Analysis

```python
from enzyme_kinetics import EnzymeKineticsAnalyzer
import numpy as np

# Create analyzer
analyzer = EnzymeKineticsAnalyzer(
    enzyme_concentration_um=12.25,
    reaction_time_seconds=10800
)

# Define substrate concentration and velocity
substrate_conc = np.array([0.1, 0.2, 0.5, 1.0, 3.0])
velocity = np.array([71, 179, 298, 283, 1026])

# Analyze
result = analyzer.analyze_peak('MyEnzyme', substrate_conc, velocity)

# Get results
print(f"Km = {result.mm_params.Km:.4f} µM")
print(f"Vmax = {result.mm_params.Vmax:.2f}")
print(f"kcat/Km = {result.kcat_over_km:.2f}")
```

---

## Common Tasks

### Task 1: Analyze Multiple Peaks from CSV

```python
import pandas as pd
from enzyme_kinetics import EnzymeKineticsAnalyzer

# Load data
df = pd.read_csv('your_data.csv')

# Analyze
analyzer = EnzymeKineticsAnalyzer(enzyme_concentration_um=12.25)
results = analyzer.analyze_dataframe(
    df,
    peak_id_column='peak_id',
    substrate_conc_column='substrate_conc',
    velocity_column='peak_area',
    velocity_std_column='peak_area_std'
)

# Export
analyzer.save_results('results.csv')
```

### Task 2: Create Calibration Curve

```python
import numpy as np
from enzyme_kinetics import CalibrationCurve

# Your CoA standard measurements
conc = np.array([0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0])
area = np.array([0, 245, 612, 1223, 2452, 4904, 12259])

# Fit calibration
calib = CalibrationCurve()
params = calib.fit(conc, area)

# Convert peak area to concentration
product_conc = calib.area_to_concentration(peak_area=612)  # → 0.25 µM

# Save calibration
calib.save('calibration.txt')
```

### Task 3: Apply Calibration to Results

```python
from enzyme_kinetics import CalibrationCurve

# Load calibration
calib = CalibrationCurve()
calib.fit(calib_conc, calib_area)

# Apply to analyzer
analyzer.apply_calibration(calib)

# Now all results are in absolute units (µM)
analyzer.save_results('results_absolute.csv')
```

### Task 4: Generate Publication Plots

```python
from enzyme_kinetics import KineticsPlotter

# Michaelis-Menten curves for all peaks
KineticsPlotter.plot_comparison_panel(results, output_path='all_curves.png')

# Efficiency comparison bars
KineticsPlotter.plot_efficiency_comparison(results, output_path='efficiency.png')

# Single peak detailed analysis
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
KineticsPlotter.plot_michaelis_menten_curve(result, ax=axes[0])
KineticsPlotter.plot_lineweaver_burk(result, ax=axes[1])
KineticsPlotter.plot_residuals(result, ax=axes[2])
plt.tight_layout()
plt.savefig('detailed_analysis.png', dpi=300)
```

---

## Data Format Requirements

### Input: Kinetics Data

CSV with columns:
```
peak_id,substrate_conc,peak_area,peak_area_std
HTAL,0.1,70.97,4.38
HTAL,0.2,178.53,2.25
OLA,0.1,400.38,11.24
...
```

### Input: Calibration Standards

CSV with columns:
```
coa_concentration_um,peak_area_mean,peak_area_std
0.0,0.0,0.0
0.1,245.3,12.4
0.25,612.4,28.1
...
```

### Output: Results

CSV with columns:
```
peak_id,Km_MM,Km_MM_std,Vmax_MM,Vmax_MM_std,R2_MM,kcat,kcat_over_km,...
HTAL,0.1688,0.0092,316.19,6.75,0.753,0.0024,0.0142,...
OLA,0.1069,0.0050,1108.52,15.39,0.7426,0.0084,0.0784,...
```

---

## Configuration

### Enzyme Concentration

```python
# Required for calculating kcat and kcat/Km
analyzer = EnzymeKineticsAnalyzer(enzyme_concentration_um=12.25)
```

### Reaction Time

```python
# Convert peak areas to velocities per unit time
analyzer = EnzymeKineticsAnalyzer(
    reaction_time_seconds=10800  # 3 hours
)
```

### Both (Recommended)

```python
analyzer = EnzymeKineticsAnalyzer(
    enzyme_concentration_um=12.25,
    reaction_time_seconds=10800
)
```

---

## Interpreting Results

### Km (Michaelis Constant)

- **Units:** µM
- **Meaning:** Substrate concentration at 50% Vmax
- **Lower is better:** High affinity (tighter binding)
- **Typical range:** 0.01–100 µM

### Vmax (Maximum Velocity)

- **Units:** area/time (or µM/time if calibrated)
- **Meaning:** Reaction velocity at saturating substrate
- **Higher is better:** Greater catalytic power
- **Depends on:** Enzyme amount and reaction time

### kcat (Turnover Number)

- **Units:** s⁻¹ (reactions per second per enzyme)
- **Meaning:** How fast each enzyme molecule catalyzes reaction
- **Higher is better:** Faster enzyme
- **Typical range:** 1–10⁶ s⁻¹

### kcat/Km (Catalytic Efficiency)

- **Units:** M⁻¹·s⁻¹ (requires calibration)
- **Meaning:** Combined measure of affinity + turnover
- **Higher is better:** More efficient enzyme
- **Typical range:** 10³–10⁸ M⁻¹·s⁻¹
- **Best metric:** Use this to compare enzymes

### R² (Coefficient of Determination)

- **Range:** 0–1
- **R² > 0.9:** Excellent fit
- **R² > 0.7:** Good fit
- **R² < 0.5:** Poor fit, suspect data or model mismatch

---

## Troubleshooting

### "ValueError: Concentration has 0 points"

Check that `substrate_conc` and `velocity` have same length.

```python
assert len(substrate_conc) == len(velocity)
```

### "Cannot fit: All velocity values are zero"

Analyzer skips peaks with all-zero activity. Check your data.

### Fit has poor R² (< 0.5)

Possible causes:
1. **Non-Michaelis-Menten kinetics** → Try `SubstrateInhibitionModel`
2. **Data quality** → Check for outliers or measurement error
3. **Reaction not linear with time** → Run time-course experiment
4. **Substrate inhibition** → Km may be artificially high at max substrate

### "Calibration not fitted. Call fit() first."

Always call `calib.fit()` before using calibration:

```python
calib = CalibrationCurve()
calib.fit(conc, area)  # ← Must do this first!
analyzer.apply_calibration(calib)
```

### Cannot import module

Add to Python path:

```python
import sys
sys.path.insert(0, '/mnt/user-data/outputs')
from enzyme_kinetics import EnzymeKineticsAnalyzer
```

---

## Advanced: Custom Models

Add your own kinetic model by extending `EnzymeKineticModel`:

```python
from enzyme_kinetics import EnzymeKineticModel, KineticParameters
import numpy as np
from scipy.optimize import curve_fit

class MyCustomModel(EnzymeKineticModel):
    """My custom enzyme kinetics equation."""
    
    def __call__(self, substrate_conc, param1, param2):
        # Your equation here
        return param1 * substrate_conc / (param2 + substrate_conc)
    
    def fit(self, substrate_conc, velocity, velocity_std=None):
        # Fit your model
        popt, pcov = curve_fit(self, substrate_conc, velocity)
        param1, param2 = popt
        
        # Return results
        params = KineticParameters(
            Km=param2,
            Km_std=np.sqrt(pcov[1,1]),
            Vmax=param1,
            Vmax_std=np.sqrt(pcov[0,0]),
            r2=0.99,  # Calculate this
            n_points=len(substrate_conc)
        )
        
        fit_data = {'params': params, 'covariance': pcov}
        return params, fit_data
```

---

## Tips for Best Results

1. **Measure at least 5 substrate concentrations**
   - Below and above Km for good parameter definition
   - At least 3 replicates per concentration

2. **Use appropriate concentration range**
   - Min: ~0.1 × Km (if known)
   - Max: ~10 × Km
   - For Km = 0.1 µM, use 0.01–1 µM range

3. **Record velocity errors**
   - Standard deviations allow uncertainty quantification
   - Improves parameter estimates

4. **Validate with Lineweaver-Burk**
   - Compare MM and LB R² values
   - If very different, suggests kinetic complexity

5. **Always calibrate**
   - Peak area alone is not meaningful
   - Calibration enables comparison with literature

---

## Next Steps

- **Read API_DOCUMENTATION.md** for complete reference
- **Check examples/** folder for working scripts
- **See COMPLETE_SUMMARY.md** for full analysis pipeline
- **Review original ENZYME_KINETICS_REPORT.md** for scientific context
