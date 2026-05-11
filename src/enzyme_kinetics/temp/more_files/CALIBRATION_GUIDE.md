# Complete Guide: Converting Peak Area to Product Concentration

## The Challenge

Your enzyme kinetics data currently expresses Vmax and kcat in terms of **relative peak area units** rather than **absolute product concentration (µM)**. This is because HPLC measures peak area, which is proportional to product amount, but the exact relationship is unknown.

To get truly comparable kinetic parameters (kcat in s⁻¹, kcat/Km in M⁻¹·s⁻¹), you need a **calibration curve**.

---

## Step 1: Design a Calibration Experiment

### What You Need:
1. **Pure CoA standard** (the product of HexCoA metabolism)
2. **Concentration gradient** of CoA (e.g., 0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0 µM)
3. **HPLC system** (same as used for your enzyme kinetics experiments)
4. **Consistent injection volume and detection method**

### Experimental Protocol:

```
For each CoA concentration:
  1. Prepare standard solution at known concentration
  2. Inject into HPLC (same method as enzyme kinetics)
  3. Record peak area (mean of 3 replicates)
  4. Calculate mean and standard deviation
```

### Example Expected Data:

| [CoA] (µM) | Peak Area (Mean) | StdDev | Replicate 1 | Replicate 2 | Replicate 3 |
|-----------|------------------|--------|-------------|-------------|-------------|
| 0.0       | 0                | 0      | 0           | 0           | 0           |
| 0.1       | 245.3            | 12.4   | 238.2       | 248.9       | 248.8       |
| 0.25      | 612.4            | 28.1   | 591.3       | 615.8       | 629.9       |
| 0.5       | 1223.1           | 45.2   | 1189.2      | 1234.5      | 1245.6      |
| 1.0       | 2451.8           | 89.3   | 2398.1      | 2461.3      | 2495.9      |
| 2.0       | 4903.5           | 156.2  | 4812.1      | 4923.8      | 4974.6      |
| 5.0       | 12258.7          | 412.1  | 11921.3     | 12384.2     | 12470.6     |

---

## Step 2: Fit the Calibration Curve

The relationship between peak area and CoA concentration is typically **linear** (at least in the concentration range of your kinetics experiments):

```
Peak Area = slope × [CoA] + intercept
```

### Python Code to Fit the Calibration:

```python
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Your calibration data
coa_concentrations = np.array([0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0])
peak_areas = np.array([0, 245.3, 612.4, 1223.1, 2451.8, 4903.5, 12258.7])
peak_areas_std = np.array([0, 12.4, 28.1, 45.2, 89.3, 156.2, 412.1])

# Fit linear relationship: peak_area = slope * [CoA] + intercept
def linear(x, slope, intercept):
    return slope * x + intercept

# Fit the model
popt, pcov = curve_fit(linear, coa_concentrations, peak_areas, 
                       sigma=peak_areas_std, absolute_sigma=True)

slope, intercept = popt
slope_std, intercept_std = np.sqrt(np.diag(pcov))

print(f"Calibration Results:")
print(f"  Slope:      {slope:.2f} ± {slope_std:.2f} area/µM")
print(f"  Intercept:  {intercept:.2f} ± {intercept_std:.2f} area")

# R² goodness of fit
area_pred = linear(coa_concentrations, slope, intercept)
ss_res = np.sum((peak_areas - area_pred) ** 2)
ss_tot = np.sum((peak_areas - np.mean(peak_areas)) ** 2)
r2 = 1 - (ss_res / ss_tot)
print(f"  R²:         {r2:.4f}")

# Inverse function: convert area back to [CoA]
def area_to_coa_conc(area):
    """Convert peak area to CoA concentration."""
    return (area - intercept) / slope

print(f"\nCalibration Function:")
print(f"  [CoA] (µM) = (peak_area - {intercept:.2f}) / {slope:.2f}")
```

### Expected Output:
```
Calibration Results:
  Slope:      2451.74 ± 12.34 area/µM
  Intercept:  -23.45 ± 45.23 area
  R²:         0.9998

Calibration Function:
  [CoA] (µM) = (peak_area + 23.45) / 2451.74
```

---

## Step 3: Apply the Calibration to Your Kinetics Data

Once you have the calibration (slope and intercept), update your kinetics analysis:

### Updated Script (Modified Section):

```python
# CALIBRATION PARAMETERS (from your calibration experiment)
AREA_TO_CONCENTRATION_SLOPE = 2451.74  # area/µM
AREA_TO_CONCENTRATION_INTERCEPT = -23.45  # area

# Convert peak areas in kinetics data to product concentration
def convert_area_to_product_concentration(area):
    """
    Convert HPLC peak area to product (CoA) concentration in µM.
    
    Uses the calibration curve: [CoA] = (area - intercept) / slope
    """
    return (area - AREA_TO_CONCENTRATION_INTERCEPT) / AREA_TO_CONCENTRATION_SLOPE

# Apply conversion to Michaelis-Menten fitting
# OLD: fit to peak area directly
# NEW: convert areas to concentrations first

for peak_id, group_data in df.groupby('peak_id'):
    S = group_data['substrate_conc'].values  # µM (substrate)
    area = group_data['mean_area'].values     # relative peak area
    
    # ===== CONVERSION STEP =====
    v = convert_area_to_product_concentration(area)  # µM product
    v_std = (group_data['std_area'].values / AREA_TO_CONCENTRATION_SLOPE)  # error propagation
    
    # Now 'v' is in absolute units (µM product formed in 3 hours)
    # Fit Michaelis-Menten as before
    popt, pcov = curve_fit(michaelis_menten, S, v, sigma=v_std, absolute_sigma=True)
    
    Vmax_um, Km = popt  # Now Vmax is in µM product per 3 hours!
    
    # Calculate absolute velocity per second
    Vmax_per_sec = Vmax_um / REACTION_TIME_SECONDS  # µM/s
    
    # Enzyme-normalize
    kcat = Vmax_per_sec / ENZYME_CONCENTRATION_UM  # s⁻¹
    kcat_over_km = kcat / Km  # M⁻¹·s⁻¹ (now in standard units!)
```

---

## Step 4: Interpret Results in Standard Units

After applying the calibration, your results will be in **standard enzyme kinetics units**:

| Parameter | Unit | Meaning |
|-----------|------|---------|
| **Km** | µM | Substrate concentration (unchanged) |
| **Vmax** | µM/s | Product formation rate at saturation |
| **kcat** | s⁻¹ | Turnover number (reactions/enzyme/second) |
| **kcat/Km** | M⁻¹·s⁻¹ | Catalytic efficiency (standard comparison metric) |

### Typical Values for Reference:

| Enzyme Type | kcat/Km (M⁻¹·s⁻¹) | Example |
|-------------|------------------|---------|
| Slow enzymes | 10³ to 10⁴ | Non-optimized variants |
| Good enzymes | 10⁵ to 10⁶ | Natural wild-type enzymes |
| Excellent enzymes | 10⁷ to 10⁸ | Highly optimized catalysts |
| Diffusion-limited | ~10⁸–10⁹ | Theoretical maximum |

---

## Step 5: Create a Comprehensive Calibration Report

### Minimal Calibration Analysis Script:

```python
#!/usr/bin/env python3
"""Generate and document your calibration curve."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load your calibration data
# Create a CSV: calibration_data.csv
# Columns: coa_concentration_um, peak_area_mean, peak_area_std

calib_df = pd.read_csv('calibration_data.csv')

coa_conc = calib_df['coa_concentration_um'].values
peak_area = calib_df['peak_area_mean'].values
peak_area_std = calib_df['peak_area_std'].values

# Fit calibration
def linear(x, slope, intercept):
    return slope * x + intercept

popt, pcov = curve_fit(linear, coa_conc, peak_area, 
                       sigma=peak_area_std, absolute_sigma=True)

slope, intercept = popt
slope_std, intercept_std = np.sqrt(np.diag(pcov))

# R² value
area_pred = linear(coa_conc, *popt)
ss_res = np.sum((peak_area - area_pred) ** 2)
ss_tot = np.sum((peak_area - np.mean(peak_area)) ** 2)
r2 = 1 - (ss_res / ss_tot)

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.errorbar(coa_conc, peak_area, yerr=peak_area_std, fmt='o', markersize=8,
            capsize=5, capthick=2, color='steelblue', elinewidth=2, 
            markeredgewidth=2, label='Data')

coa_smooth = np.linspace(0, coa_conc.max() * 1.1, 100)
area_fit = linear(coa_smooth, slope, intercept)
ax.plot(coa_smooth, area_fit, 'r-', linewidth=2.5, label='Linear fit')

ax.set_xlabel('[CoA] (µM)', fontsize=12, fontweight='bold')
ax.set_ylabel('Peak Area', fontsize=12, fontweight='bold')
ax.set_title(f'HPLC Calibration Curve for CoA\nSlope={slope:.2f}, Intercept={intercept:.2f}, R²={r2:.4f}',
            fontsize=13, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig('calibration_curve.png', dpi=300)

# Save calibration constants to file
with open('calibration_constants.txt', 'w') as f:
    f.write("HPLC CALIBRATION CONSTANTS\n")
    f.write("=" * 50 + "\n\n")
    f.write(f"Equation: peak_area = {slope:.4f} × [CoA] + {intercept:.4f}\n")
    f.write(f"Inverse:  [CoA] = (peak_area - {intercept:.4f}) / {slope:.4f}\n\n")
    f.write(f"Slope:       {slope:.4f} ± {slope_std:.4f} area/µM\n")
    f.write(f"Intercept:   {intercept:.4f} ± {intercept_std:.4f} area\n")
    f.write(f"R²:          {r2:.6f}\n")
    f.write(f"\nUse these values in your kinetics analysis!\n")

print(f"✓ Calibration constants saved to: calibration_constants.txt")
print(f"✓ Calibration plot saved to: calibration_curve.png")
```

---

## Important Caveats

### 1. **Linearity Assumption**
- Linear calibration assumes peak area ∝ product concentration
- True for most HPLC detectors (UV, evaporative light scattering)
- Test by checking residuals in calibration fit (R² > 0.99 desired)

### 2. **Detector Response**
- Different compounds may have different detector responses
- HexCoA and CoA may have different extinction coefficients or ionization efficiencies
- Consider running a CoA standard to ensure linear response

### 3. **Concentration Range**
- Calibration is only valid over the range tested
- Don't extrapolate beyond your highest standard
- For your kinetics data, make sure all velocities fall within calibration range

### 4. **System Stability**
- If HPLC response changes over time, recalibrate periodically
- Temperature, pH, and column condition affect detector response
- Recalibrate if peak widths change significantly

---

## Summary: Three Steps to Absolute Kinetics

1. **Run calibration experiment**
   - Inject CoA standards (0.1 to 5 µM, 3 replicates each)
   - Record mean peak areas and standard deviations

2. **Fit calibration curve**
   - Linear fit: peak_area = slope × [CoA] + intercept
   - Verify R² > 0.99

3. **Apply to kinetics data**
   - Convert Vmax (area) → Vmax (µM)
   - Recalculate kcat and kcat/Km
   - Compare with literature values in M⁻¹·s⁻¹ units

Once complete, your enzyme kinetics will be in **truly comparable, publishable units**.
