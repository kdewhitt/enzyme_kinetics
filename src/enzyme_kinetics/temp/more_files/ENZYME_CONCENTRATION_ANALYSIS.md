# Why Protein Concentration Matters in Enzyme Kinetics

## The Problem with Current Analysis

The kinetic parameters I calculated earlier are **velocity measurements in terms of peak area units**, which are **enzyme-concentration dependent**. This creates several issues:

### 1. **Vmax is Not Absolute Without Enzyme Concentration**
- **Current Vmax** (2941.75 for PDAL) is meaningless without knowing the enzyme concentration
- **Correct interpretation**: This is `v_max` at 12.25 µM enzyme, not the intrinsic Vmax
- **Normalization needed**: Divide by enzyme concentration to get actual Vmax

### 2. **Km is Independent of Enzyme Concentration (GOOD)**
- Km (Michaelis constant) is an **intrinsic property** of the enzyme
- Does NOT depend on how much enzyme you used
- ✓ Current Km values are valid as-is
- ✓ PDAL's Km = 0.097 µM is the true substrate concentration at half-maximal velocity

### 3. **Cannot Compare with Literature or Other Studies WITHOUT [E]**
- Different labs use different enzyme concentrations
- Without normalization, you cannot compare your Vmax with published values
- Example: Lab A uses 1 µM enzyme, Lab B uses 10 µM → Vmax will differ 10-fold

---

## Key Metric: Turnover Number (kcat)

The **turnover number (kcat)**, also called **catalytic constant**, is the **enzyme-normalized velocity**:

```
kcat = Vmax / [E]
```

Where:
- **Vmax**: Maximum velocity (current values, in area/time units)
- **[E]**: Enzyme concentration (12.25 µM)
- **kcat**: Turnover number = catalytic events per enzyme molecule per unit time

### What kcat tells you:
- **kcat = 100 s⁻¹** means each enzyme molecule catalyzes 100 reactions per second
- This is an **intrinsic property** independent of how much enzyme you added
- Allows direct comparison across different experiments and labs

---

## Catalytic Efficiency: kcat/Km

The **most important combined metric** for enzyme performance:

```
Catalytic Efficiency = kcat / Km
```

### Why this matters:
- Reflects both **substrate affinity** (Km) and **turnover rate** (kcat)
- An enzyme with:
  - High Km (weak affinity) BUT high kcat (fast turnover)
  - Might have similar efficiency to one with low Km and slow kcat
- Units: **µM⁻¹·s⁻¹** (or M⁻¹·s⁻¹)
- Higher is better; exceptional enzymes reach 10⁵ to 10⁷ M⁻¹·s⁻¹

### Evolutionary optimization:
- Wild-type enzymes are typically optimized for:
  - **Km ≈ [S]_cellular**: Km close to actual cellular substrate concentration
  - **kcat/Km**: Often approaching diffusion-limited rates (~10⁸ M⁻¹·s⁻¹)

---

## What We Can NOW Calculate

With your enzyme concentration (12.25 µM), we can compute:

1. **kcat** (absolute turnover number per enzyme molecule)
2. **kcat/Km** (catalytic efficiency)
3. **Relative performance** (which enzyme is most efficient)
4. **Turnover frequency** (reactions/enzyme/second)

---

## Important Caveat: Peak Area ≠ Product Concentration

**Current limitation**: Peak area is still not an absolute velocity measurement.

To compute true kcat in units of **s⁻¹**, we need:
- Substrate (HexCoA) concentration → product (CoA) formation rate
- Calibration curve: peak area → µM product formed
- Reaction time for each measurement

### What we CAN do:
- Calculate **relative kcat/Km** values (enzyme comparison)
- Normalize by [E] to remove enzyme-concentration dependency
- Identify which enzyme is most efficient (dimensionless ratio)

### What we CANNOT do (without calibration):
- Absolute kcat in units of **s⁻¹**
- Comparison with literature values (which use standardized units)
- Predictive modeling of reaction rates at different [E]

---

## Recommendation Going Forward

1. **Use kcat/Km as the comparison metric** (normalizes for [E])
2. **Develop a calibration curve**:
   - Run HexCoA → CoA with known concentrations
   - Correlate HPLC peak area to product concentration
   - Convert area → product velocity (µM/min or µM/s)
   - Recalculate with absolute units

3. **Report both**:
   - Km (unchanged, intrinsic property)
   - kcat/Km (enzyme efficiency, normalized)
   - kcat (once you have product calibration)

---

## Bottom Line

✓ **Km values are already correct** — enzyme concentration doesn't affect them  
✗ **Vmax values are enzyme-concentration dependent** — need to normalize  
→ **Next step**: Calculate kcat and kcat/Km for meaningful enzyme comparison

