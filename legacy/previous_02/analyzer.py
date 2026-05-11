"""
Main enzyme kinetics analysis engine.

Provides high-level interface for analyzing enzyme kinetics across
multiple samples/peaks with support for different models.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from models import MichaelisMentenModel, LineweaverBurkModel, KineticParameters
from calibration import CalibrationCurve, CalibrationParameters


# examined
@dataclass
class EnzymeKineticsResult:
    """Complete kinetics analysis result for a single enzyme/peak."""
    
    peak_id: str
    substrate_conc: np.ndarray
    velocity: np.ndarray
    velocity_std: Optional[np.ndarray]
    
    mm_params: KineticParameters
    mm_fit_data: dict
    
    lb_params: Optional[KineticParameters] = None
    lb_fit_data: Optional[dict] = None
    
    kcat: Optional[float] = None
    kcat_std: Optional[float] = None
    kcat_over_km: Optional[float] = None
    kcat_over_km_std: Optional[float] = None
    
    enzyme_concentration_um: Optional[float] = None
    reaction_time_seconds: Optional[float] = None
    
    def to_dict(self) -> dict:
        """Convert to dictionary for DataFrame serialization."""
        return {
            'peak_id': self.peak_id,
            'n_points': self.mm_params.n_points,
            'Km_MM': self.mm_params.Km,
            'Km_MM_std': self.mm_params.Km_std,
            'Vmax_MM': self.mm_params.Vmax,
            'Vmax_MM_std': self.mm_params.Vmax_std,
            'R2_MM': self.mm_params.r2,
            'Km_LB': self.lb_params.Km if self.lb_params else None,
            'Vmax_LB': self.lb_params.Vmax if self.lb_params else None,
            'R2_LB': self.lb_params.r2 if self.lb_params else None,
            'kcat': self.kcat,
            'kcat_std': self.kcat_std,
            'kcat_over_km': self.kcat_over_km,
            'kcat_over_km_std': self.kcat_over_km_std,
            'substrate_conc_min': self.substrate_conc.min(),
            'substrate_conc_max': self.substrate_conc.max(),
            'velocity_min': self.velocity.min(),
            'velocity_max': self.velocity.max(),
        }


class EnzymeKineticsAnalyzer:
    """
    Comprehensive enzyme kinetics analyzer.
    
    Handles:
    - Loading HPLC data
    - Fitting multiple kinetic models
    - Enzyme/time normalization
    - Calibration integration
    - Batch analysis across peaks
    """
    
    def __init__(
        self,
        enzyme_concentration_um: Optional[float] = None,
        reaction_time_seconds: Optional[float] = None,
        calibration: Optional[CalibrationCurve] = None
    ):
        """
        Initialize enzyme kinetics analyzer.
        
        Parameters
        ----------
        enzyme_concentration_um : float, optional
            Total enzyme concentration in the reaction (µM)
        reaction_time_seconds : float, optional
            Duration of reaction in seconds
        calibration : CalibrationCurve, optional
            Fitted calibration curve for area → concentration conversion
        """
        
        self.enzyme_concentration_um = enzyme_concentration_um
        self.reaction_time_seconds = reaction_time_seconds
        self.calibration = calibration
        
        self.mm_model = MichaelisMentenModel()
        self.lb_model = LineweaverBurkModel()
        
        self.results: List[EnzymeKineticsResult] = []
    
    def analyze_peak(
        self,
        peak_id: str,
        substrate_conc: np.ndarray,
        velocity: np.ndarray,
        velocity_std: Optional[np.ndarray] = None,
        fit_lineweaver_burk: bool = True
    ) -> EnzymeKineticsResult:
        """
        Analyze kinetics for a single peak.
        
        Parameters
        ----------
        peak_id : str
            Peak identifier
        substrate_conc : np.ndarray
            Substrate concentrations (µM)
        velocity : np.ndarray
            Measured velocities (arbitrary units or µM if calibrated)
        velocity_std : np.ndarray, optional
            Standard deviations of velocity measurements
        fit_lineweaver_burk : bool
            Whether to also fit Lineweaver-Burk model
        
        Returns
        -------
        result : EnzymeKineticsResult
            Complete analysis result
        """
        
        # Fit Michaelis-Menten
        mm_params, mm_fit_data = self.mm_model.fit(
            substrate_conc,
            velocity,
            velocity_std=velocity_std
        )
        
        # Fit Lineweaver-Burk if requested
        lb_params = None
        lb_fit_data = None
        if fit_lineweaver_burk:
            try:
                lb_params, lb_fit_data = self.lb_model.fit(substrate_conc, velocity)
            except Exception:
                pass  # LB may fail if velocity contains zeros
        
        # Calculate enzyme-normalized parameters if concentration known
        kcat = None
        kcat_std = None
        kcat_over_km = None
        kcat_over_km_std = None
        
        if self.enzyme_concentration_um is not None and self.reaction_time_seconds is not None:
            # Convert Vmax to per-second basis
            vmax_per_sec = mm_params.Vmax / self.reaction_time_seconds
            vmax_per_sec_std = mm_params.Vmax_std / self.reaction_time_seconds
            
            # kcat = Vmax / [E]
            kcat = vmax_per_sec / self.enzyme_concentration_um
            kcat_std = vmax_per_sec_std / self.enzyme_concentration_um
            
            # kcat/Km
            if mm_params.Km > 0:
                kcat_over_km = kcat / mm_params.Km
                # Error propagation
                rel_err_kcat = kcat_std / kcat if kcat > 0 else 0
                rel_err_km = mm_params.Km_std / mm_params.Km
                kcat_over_km_std = kcat_over_km * np.sqrt(rel_err_kcat**2 + rel_err_km**2)
        
        result = EnzymeKineticsResult(
            peak_id=peak_id,
            substrate_conc=substrate_conc,
            velocity=velocity,
            velocity_std=velocity_std,
            mm_params=mm_params,
            mm_fit_data=mm_fit_data,
            lb_params=lb_params,
            lb_fit_data=lb_fit_data,
            kcat=kcat,
            kcat_std=kcat_std,
            kcat_over_km=kcat_over_km,
            kcat_over_km_std=kcat_over_km_std,
            enzyme_concentration_um=self.enzyme_concentration_um,
            reaction_time_seconds=self.reaction_time_seconds
        )
        
        self.results.append(result)
        
        return result
    
    def analyze_dataframe(
        self,
        df: pd.DataFrame,
        peak_id_column: str = 'peak_id',
        substrate_conc_column: str = 'substrate_conc',
        velocity_column: str = 'velocity',
        velocity_std_column: Optional[str] = None,
        fit_lineweaver_burk: bool = True,
        min_points: int = 3
    ) -> List[EnzymeKineticsResult]:
        """
        Batch analyze kinetics from DataFrame.
        
        Parameters
        ----------
        df : pd.DataFrame
            Input data with columns for peak_id, substrate concentration, velocity
        peak_id_column : str
            Name of peak ID column
        substrate_conc_column : str
            Name of substrate concentration column
        velocity_column : str
            Name of velocity/peak area column
        velocity_std_column : str, optional
            Name of velocity standard deviation column
        fit_lineweaver_burk : bool
            Whether to fit Lineweaver-Burk model
        min_points : int
            Minimum number of data points required to fit peak
        
        Returns
        -------
        results : List[EnzymeKineticsResult]
            Results for each peak
        """
        
        results = []
        
        for peak_id, group in df.groupby(peak_id_column):
            
            group = group.sort_values(substrate_conc_column)
            
            substrate_conc = group[substrate_conc_column].values
            velocity = group[velocity_column].values
            velocity_std = (
                group[velocity_std_column].values
                if velocity_std_column and velocity_std_column in group.columns
                else None
            )
            
            # Skip if insufficient data or all zeros
            if len(substrate_conc) < min_points or np.all(velocity == 0):
                print(f"⊘ {peak_id}: Skipped (insufficient data or all zeros)")
                continue
            
            try:
                result = self.analyze_peak(
                    peak_id,
                    substrate_conc,
                    velocity,
                    velocity_std=velocity_std,
                    fit_lineweaver_burk=fit_lineweaver_burk
                )
                results.append(result)
                print(f"✓ {peak_id}: {result.mm_params}")
                
            except Exception as e:
                print(f"✗ {peak_id}: Failed - {str(e)}")
                continue
        
        return results
    
    def results_to_dataframe(self) -> pd.DataFrame:
        """Convert all results to DataFrame."""
        data = [result.to_dict() for result in self.results]
        return pd.DataFrame(data)
    
    def apply_calibration(self, calibration: CalibrationCurve) -> None:
        """
        Apply calibration curve to convert peak areas to concentrations.
        
        Parameters
        ----------
        calibration : CalibrationCurve
            Fitted calibration curve
        """
        
        self.calibration = calibration
        
        # Re-analyze all results with calibrated velocities
        for result in self.results:
            
            # Convert peak areas to concentrations
            velocity_µm = calibration.area_to_concentration(result.velocity)
            velocity_std_µm = None
            
            if result.velocity_std is not None:
                # Error propagation for each point
                velocity_std_µm = result.velocity_std / calibration.params.slope
            
            # Re-fit with calibrated values
            mm_params_calib, mm_fit_data_calib = self.mm_model.fit(
                result.substrate_conc,
                velocity_µm,
                velocity_std=velocity_std_µm
            )
            
            result.mm_params = mm_params_calib
            result.mm_fit_data = mm_fit_data_calib
            result.velocity = velocity_µm
            result.velocity_std = velocity_std_µm
            
            # Recalculate kcat/Km
            if self.enzyme_concentration_um is not None and self.reaction_time_seconds is not None:
                vmax_per_sec = mm_params_calib.Vmax / self.reaction_time_seconds
                vmax_per_sec_std = mm_params_calib.Vmax_std / self.reaction_time_seconds
                
                result.kcat = vmax_per_sec / self.enzyme_concentration_um
                result.kcat_std = vmax_per_sec_std / self.enzyme_concentration_um
                
                if mm_params_calib.Km > 0:
                    result.kcat_over_km = result.kcat / mm_params_calib.Km
                    rel_err_kcat = result.kcat_std / result.kcat if result.kcat > 0 else 0
                    rel_err_km = mm_params_calib.Km_std / mm_params_calib.Km
                    result.kcat_over_km_std = result.kcat_over_km * np.sqrt(
                        rel_err_kcat**2 + rel_err_km**2
                    )
    
    def save_results(self, csv_path: str | Path) -> None:
        """Save results to CSV file."""
        df = self.results_to_dataframe()
        df.to_csv(csv_path, index=False)
        print(f"✓ Results saved to: {csv_path}")
