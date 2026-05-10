# Compare working_copy to previous_02

# Based on these comparisons, make the following changes to working_copy:

Based on your assessment, I want to migrate my codebase to using refactored v1 only after you have addressed the issues
highlighted in your functional comparison and overall assessment. for convenience, I have also listed the issues that I
would like you to address below.

1. Address assumption that Km and Vmax are uncorrelated.
2. Include intercept_se in the propagation of the delta-method formula.
3. Address fragility and opaqueness of apply_calibration and prepare_velocity only dividing by rxn_time. As you
   suggested, apply the calibration to mean_signal (raw area) directly before dividing by rxn_time.
4. Add two methods for v_sem scaling enabling the user to select between a simplified v_sem_um calculation versus the
   original's more precise approach.
4. Implement Ki as a typed field visible in KineticConstants.to_dict()
5. Add a new method to extend the plot() method that restores the LB plot, residual plot, and efficiency comparison
   panel.
6. Restore the calibration is_valid() guard.
7. Fix the semantic issue noted for the Hill FitResult regarding k_half and km.

### !!!!!!!!!! Future (for internal awareness) !!!!!!!!!!

To address the issues at a future date:

- fix original 3.6 issue