# HIGH priority
- [ ] Add ability to specify name of "peak_id" and related columns
- [ ] Finish comparisons between original codebase's plotting.py and refactored v1
- [ ] Add ability to remove void_peak from peak_id list in CSV file

# MEDIUM priority
- [ ] Determine units for Km and Vmax and other parameters
- [ ] Determine units for all inputs, e.g., substrate_conc (in µM)
- 
# LOW priority

- [ ] Run claude to analyze for unused parameters and attributes across codebase


# Observations/ Comments/ Notes

## Original code:
- class LineweaverBurkModel(EnzymeKineticModel).fit()
- There is no equivalent in refactored v1.

```python
# Estimate standard errors (approximate)
velocity_inv_pred = slope * substrate_conc_inv + intercept
residuals = velocity_inv - velocity_inv_pred
se = np.sqrt(np.sum(residuals ** 2) / (len(substrate_conc_inv) - 2))
```

## Original code: EnzymeKineticsAnalyzer.analyze.peak() & EnzymeKineticsAnalyzer.apply_calibration()
- Original kcat/Km calculations with error progation do not make sense to me;
- it is unclear if they are correctly implemented in refactored v1. Same comment
- also applies to apply_caligration.
- 