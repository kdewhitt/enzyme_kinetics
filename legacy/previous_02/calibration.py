"""
Calibration curve management for HPLC peak area to concentration conversion.

Provides classes for fitting and applying calibration curves to convert
relative peak area measurements to absolute product concentrations.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


@dataclass
class CalibrationParameters:
    """Container for calibration curve parameters and statistics."""
    
    slope: float
    slope_std: float
    intercept: float
    intercept_std: float
    r2: float
    rmse: float
    n_points: int
    conc_range: Tuple[float, float]
    
    def __str__(self) -> str:
        return (f"Slope={self.slope:.6f}±{self.slope_std:.6f} area/µM | "
                f"Intercept={self.intercept:.6f}±{self.intercept_std:.6f} | "
                f"R²={self.r2:.6f} | "
                f"RMSE={self.rmse:.4f}")
    
    def area_to_concentration(self, area: np.ndarray | float) -> np.ndarray | float:
        """Convert peak area to product concentration (µM)."""
        return (area - self.intercept) / self.slope
    
    def concentration_to_area(self, conc: np.ndarray | float) -> np.ndarray | float:
        """Convert concentration to peak area."""
        return self.slope * conc + self.intercept


class CalibrationCurve:
    """
    Manages HPLC calibration for peak area → concentration conversion.
    
    Fits a linear calibration curve using CoA standards and provides
    conversion utilities with error propagation.
    """
    
    def __init__(self):
        self.params: CalibrationParameters | None = None
        self.fit_data: dict | None = None
    
    def fit(
        self,
        concentrations: np.ndarray,
        peak_areas: np.ndarray,
        peak_areas_std: np.ndarray | None = None
    ) -> CalibrationParameters:
        """
        Fit linear calibration curve: peak_area = slope × [concentration] + intercept
        
        Parameters
        ----------
        concentrations : np.ndarray
            CoA standard concentrations (µM)
        peak_areas : np.ndarray
            Measured HPLC peak areas
        peak_areas_std : np.ndarray, optional
            Standard deviations of peak area measurements
        
        Returns
        -------
        params : CalibrationParameters
            Fitted calibration parameters and statistics
        """
        
        # Define linear model
        def linear(x, slope, intercept):
            return slope * x + intercept
        
        # Fit curve
        popt, pcov = curve_fit(
            linear,
            concentrations,
            peak_areas,
            sigma=peak_areas_std if peak_areas_std is not None and np.any(peak_areas_std > 0) else None,
            absolute_sigma=True
        )
        
        slope, intercept = popt
        slope_std, intercept_std = np.sqrt(np.diag(pcov))
        
        # Calculate R²
        area_pred = linear(concentrations, slope, intercept)
        ss_res = np.sum((peak_areas - area_pred) ** 2)
        ss_tot = np.sum((peak_areas - np.mean(peak_areas)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        # RMSE
        rmse = np.sqrt(np.mean((peak_areas - area_pred) ** 2))
        
        self.params = CalibrationParameters(
            slope=slope,
            slope_std=slope_std,
            intercept=intercept,
            intercept_std=intercept_std,
            r2=r2,
            rmse=rmse,
            n_points=len(concentrations),
            conc_range=(concentrations.min(), concentrations.max())
        )
        
        self.fit_data = {
            'concentrations': concentrations,
            'peak_areas_measured': peak_areas,
            'peak_areas_std': peak_areas_std,
            'peak_areas_predicted': area_pred,
            'residuals': peak_areas - area_pred,
            'linear_function': linear
        }
        
        return self.params
    
    def area_to_concentration(
        self,
        peak_area: np.ndarray | float
    ) -> np.ndarray | float:
        """
        Convert HPLC peak area to product concentration with error propagation.
        
        Parameters
        ----------
        peak_area : np.ndarray or float
            Peak area measurement(s)
        
        Returns
        -------
        concentration : np.ndarray or float
            Product concentration in µM
        
        Raises
        ------
        ValueError
            If calibration has not been fitted
        """
        
        if self.params is None:
            raise ValueError("Calibration curve not fitted. Call fit() first.")
        
        return self.params.area_to_concentration(peak_area)
    
    def area_to_concentration_with_error(
        self,
        peak_area: float,
        peak_area_std: float
    ) -> Tuple[float, float]:
        """
        Convert peak area to concentration with error propagation.
        
        Equation: [C] = (area - intercept) / slope
        Error: σ_C = sqrt((σ_area/slope)² + ((area - intercept)·σ_slope/slope²)²)
        
        Parameters
        ----------
        peak_area : float
            Peak area measurement
        peak_area_std : float
            Standard deviation of peak area measurement
        
        Returns
        -------
        concentration : float
            Product concentration (µM)
        concentration_std : float
            Standard deviation of concentration (µM)
        """
        
        if self.params is None:
            raise ValueError("Calibration curve not fitted. Call fit() first.")
        
        conc = self.params.area_to_concentration(peak_area)
        
        # Error propagation: [C] = (area - intercept) / slope
        numerator_error = peak_area_std
        numerator_value = peak_area - self.params.intercept
        
        # ∂[C]/∂area = 1/slope
        dC_darea = 1 / self.params.slope
        
        # ∂[C]/∂slope = -(area - intercept) / slope²
        dC_dslope = -numerator_value / (self.params.slope ** 2)
        
        # Total uncertainty
        conc_std = np.sqrt(
            (dC_darea * numerator_error) ** 2 +
            (dC_dslope * self.params.slope_std) ** 2
        )
        
        return conc, conc_std
    
    def is_valid(self) -> bool:
        """Check if calibration is fitted and valid."""
        return (self.params is not None and
                self.params.r2 > 0.95 and
                self.params.slope > 0)
    
    def save(self, filepath: str | Path) -> None:
        """
        Save calibration parameters to text file.
        
        Parameters
        ----------
        filepath : str or Path
            Output file path
        """
        
        if self.params is None:
            raise ValueError("No calibration fitted to save.")
        
        filepath = Path(filepath)
        
        with open(filepath, 'w') as f:
            f.write("="*80 + "\n")
            f.write("HPLC CALIBRATION PARAMETERS\n")
            f.write("="*80 + "\n\n")
            
            f.write("LINEAR MODEL:\n")
            f.write("  peak_area = slope × [concentration] + intercept\n\n")
            
            f.write("PARAMETERS:\n")
            f.write(f"  Slope:        {self.params.slope:.10f} ± {self.params.slope_std:.10f} area/µM\n")
            f.write(f"  Intercept:    {self.params.intercept:.10f} ± {self.params.intercept_std:.10f} area\n")
            f.write(f"  R²:           {self.params.r2:.10f}\n")
            f.write(f"  RMSE:         {self.params.rmse:.6f} area\n")
            f.write(f"  N points:     {self.params.n_points}\n")
            f.write(f"  Conc range:   {self.params.conc_range[0]:.2f}–{self.params.conc_range[1]:.2f} µM\n\n")
            
            f.write("INVERSE FUNCTION:\n")
            f.write(f"  [concentration] = (peak_area - {self.params.intercept:.10f}) / {self.params.slope:.10f}\n\n")
            
            f.write("PYTHON CODE:\n")
            f.write("```python\n")
            f.write(f"CALIB_SLOPE = {self.params.slope:.10f}\n")
            f.write(f"CALIB_INTERCEPT = {self.params.intercept:.10f}\n")
            f.write(f"CALIB_SLOPE_STD = {self.params.slope_std:.10f}\n")
            f.write(f"CALIB_INTERCEPT_STD = {self.params.intercept_std:.10f}\n\n")
            f.write("concentration = (peak_area - CALIB_INTERCEPT) / CALIB_SLOPE\n")
            f.write("```\n")
    
    def plot(self, output_path: str | Path | None = None) -> None:
        """
        Generate calibration curve plot with validation panels.
        
        Parameters
        ----------
        output_path : str or Path, optional
            Save plot to this path. If None, displays interactively.
        """
        
        if self.params is None or self.fit_data is None:
            raise ValueError("No calibration fitted to plot.")
        
        conc = self.fit_data['concentrations']
        area_measured = self.fit_data['peak_areas_measured']
        area_std = self.fit_data['peak_areas_std']
        area_pred = self.fit_data['peak_areas_predicted']
        residuals = self.fit_data['residuals']
        linear_fn = self.fit_data['linear_function']
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Plot 1: Calibration curve
        ax = axes[0, 0]
        ax.errorbar(conc, area_measured, yerr=area_std, fmt='o', markersize=10,
                    capsize=6, capthick=2.5, color='steelblue', elinewidth=2,
                    markeredgewidth=2.5, markeredgecolor='navy', label='Measurements', zorder=3)
        
        conc_smooth = np.linspace(conc.min() - conc.max()*0.1, conc.max() * 1.15, 200)
        area_fit = linear_fn(conc_smooth, self.params.slope, self.params.intercept)
        ax.plot(conc_smooth, area_fit, 'r-', linewidth=3, label='Linear fit', zorder=2)
        
        ax.set_xlabel('[CoA] (µM)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Peak Area', fontsize=12, fontweight='bold')
        ax.set_title(f'Calibration Curve (R² = {self.params.r2:.6f})', fontsize=13, fontweight='bold')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Residuals
        ax = axes[0, 1]
        ax.errorbar(conc, residuals, yerr=area_std if area_std is not None else None,
                    fmt='o', markersize=10, capsize=6, capthick=2.5, color='coral',
                    elinewidth=2, markeredgewidth=2.5, markeredgecolor='darkred', zorder=3)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=2, zorder=2)
        
        ax.set_xlabel('[CoA] (µM)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Residuals (area)', fontsize=12, fontweight='bold')
        ax.set_title('Residual Plot', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Plot 3: Normalized residuals
        ax = axes[1, 0]
        if area_std is not None and np.any(area_std > 0):
            norm_residuals = residuals / area_std
        else:
            norm_residuals = residuals / np.std(residuals)
        
        ax.scatter(conc, norm_residuals, s=150, color='green', alpha=0.6,
                   edgecolors='darkgreen', linewidth=2.5, zorder=3)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=2, zorder=2)
        ax.axhline(y=3, color='red', linestyle='--', linewidth=2, alpha=0.5, label='±3σ')
        ax.axhline(y=-3, color='red', linestyle='--', linewidth=2, alpha=0.5)
        
        ax.set_xlabel('[CoA] (µM)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Standardized Residuals', fontsize=12, fontweight='bold')
        ax.set_title('Normalized Residuals', fontsize=13, fontweight='bold')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_ylim([-4, 4])
        
        # Plot 4: Summary
        ax = axes[1, 1]
        ax.axis('off')
        
        summary = f"""
CALIBRATION SUMMARY

Linear Model:
  peak_area = slope × [CoA] + intercept

Parameters:
  Slope:       {self.params.slope:.6f} ± {self.params.slope_std:.6f} area/µM
  Intercept:   {self.params.intercept:.6f} ± {self.params.intercept_std:.6f}

Fit Quality:
  R²:          {self.params.r2:.8f}
  RMSE:        {self.params.rmse:.4f} area
  N points:    {self.params.n_points}

Concentration Range:
  Min:         {self.params.conc_range[0]:.2f} µM
  Max:         {self.params.conc_range[1]:.2f} µM
        """
        
        ax.text(0.05, 0.95, summary, transform=ax.transAxes,
                fontsize=11, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        fig.suptitle('HPLC Calibration Curve Analysis', fontsize=14, fontweight='bold')
        fig.tight_layout()
        
        if output_path is not None:
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"✓ Calibration plot saved: {output_path}")
        else:
            plt.show()
        
        plt.close(fig)
