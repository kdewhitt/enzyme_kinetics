"""
Core models for enzyme kinetics analysis.

Provides fundamental equations and fitting models for Michaelis-Menten
and Lineweaver-Burk kinetics.
"""

from dataclasses import dataclass
from typing import Tuple, Callable

import numpy as np
from scipy.optimize import curve_fit


@dataclass
class KineticParameters:
    """Container for kinetic parameters and their uncertainties."""
    
    Km: float
    Km_std: float
    Vmax: float
    Vmax_std: float
    r2: float
    n_points: int
    
    def __str__(self) -> str:
        return (f"Km={self.Km:.6f}±{self.Km_std:.6f} µM | "
                f"Vmax={self.Vmax:.4f}±{self.Vmax_std:.4f} | "
                f"R²={self.r2:.4f}")


class EnzymeKineticModel:
    """Base class for enzyme kinetic models."""
    
    def __call__(self, substrate_conc: np.ndarray, **params) -> np.ndarray:
        """Evaluate model at given substrate concentrations."""
        raise NotImplementedError
    
    def fit(
        self,
        substrate_conc: np.ndarray,
        velocity: np.ndarray,
        velocity_std: np.ndarray | None = None,
        **initial_params
    ) -> Tuple[KineticParameters, dict]:
        """
        Fit kinetic model to velocity data.
        
        Parameters
        ----------
        substrate_conc : np.ndarray
            Substrate concentrations (µM)
        velocity : np.ndarray
            Measured velocities (arbitrary units)
        velocity_std : np.ndarray, optional
            Standard deviations of velocity measurements
        **initial_params
            Initial parameter guesses
        
        Returns
        -------
        params : KineticParameters
            Fitted kinetic parameters
        fit_data : dict
            Additional fit information (predictions, residuals, etc.)
        """
        raise NotImplementedError
    
    def calculate_r2(
        self,
        actual: np.ndarray,
        predicted: np.ndarray
    ) -> float:
        """Calculate coefficient of determination."""
        ss_res = np.sum((actual - predicted) ** 2)
        ss_tot = np.sum((actual - np.mean(actual)) ** 2)
        return 1 - (ss_res / ss_tot) if ss_tot != 0 else 0


class MichaelisMentenModel(EnzymeKineticModel):
    """
    Michaelis-Menten enzyme kinetics model.
    
    Equation: v = (Vmax * [S]) / (Km + [S])
    """
    
    def __call__(
        self,
        substrate_conc: np.ndarray,
        Vmax: float,
        Km: float
    ) -> np.ndarray:
        """
        Michaelis-Menten equation.
        
        Parameters
        ----------
        substrate_conc : np.ndarray
            Substrate concentrations (µM)
        Vmax : float
            Maximum reaction velocity
        Km : float
            Michaelis constant
        
        Returns
        -------
        velocity : np.ndarray
            Reaction velocities
        """
        return (Vmax * substrate_conc) / (Km + substrate_conc)
    
    def fit(
        self,
        substrate_conc: np.ndarray,
        velocity: np.ndarray,
        velocity_std: np.ndarray | None = None,
        Vmax_init: float | None = None,
        Km_init: float | None = None,
        maxfev: int = 5000
    ) -> Tuple[KineticParameters, dict]:
        """
        Fit Michaelis-Menten model using Levenberg-Marquardt algorithm.
        
        Parameters
        ----------
        substrate_conc : np.ndarray
            Substrate concentrations (µM)
        velocity : np.ndarray
            Measured velocities
        velocity_std : np.ndarray, optional
            Standard deviations of measurements
        Vmax_init : float, optional
            Initial guess for Vmax (default: max velocity)
        Km_init : float, optional
            Initial guess for Km (default: median substrate concentration)
        maxfev : int
            Maximum number of function evaluations
        
        Returns
        -------
        params : KineticParameters
            Fitted parameters
        fit_data : dict
            Fit predictions, residuals, and diagnostic information
        """
        
        # Initial parameter guesses
        if Vmax_init is None:
            Vmax_init = np.max(velocity)
        if Km_init is None:
            Km_init = np.median(substrate_conc)
        
        p0 = [Vmax_init, Km_init]
        
        # Perform curve fit
        popt, pcov = curve_fit(
            self,
            substrate_conc,
            velocity,
            p0=p0,
            maxfev=maxfev,
            sigma=velocity_std if velocity_std is not None and np.any(velocity_std > 0) else None,
            absolute_sigma=True
        )
        
        Vmax_fit, Km_fit = popt
        perr = np.sqrt(np.diag(pcov))
        Vmax_std, Km_std = perr
        
        # Calculate predictions and R²
        velocity_pred = self(substrate_conc, Vmax_fit, Km_fit)
        r2 = self.calculate_r2(velocity, velocity_pred)
        
        # Calculate residuals
        residuals = velocity - velocity_pred
        
        params = KineticParameters(
            Km=Km_fit,
            Km_std=Km_std,
            Vmax=Vmax_fit,
            Vmax_std=Vmax_std,
            r2=r2,
            n_points=len(substrate_conc)
        )
        
        fit_data = {
            'predicted': velocity_pred,
            'residuals': residuals,
            'covariance': pcov,
            'substrate_conc': substrate_conc,
            'velocity_measured': velocity,
            'velocity_std': velocity_std
        }
        
        return params, fit_data


class LineweauerBurkModel(EnzymeKineticModel):
    """
    Lineweaver-Burk linearization of Michaelis-Menten kinetics.
    
    Equation: 1/v = (Km/Vmax) * (1/[S]) + 1/Vmax
    """
    
    def __call__(
        self,
        substrate_conc_inv: np.ndarray,
        Vmax: float,
        Km: float
    ) -> np.ndarray:
        """
        Lineweaver-Burk equation.
        
        Parameters
        ----------
        substrate_conc_inv : np.ndarray
            Inverse substrate concentrations (1/[S])
        Vmax : float
            Maximum reaction velocity
        Km : float
            Michaelis constant
        
        Returns
        -------
        velocity_inv : np.ndarray
            Inverse velocities (1/v)
        """
        return (Km / Vmax) * substrate_conc_inv + 1 / Vmax
    
    def fit(
        self,
        substrate_conc: np.ndarray,
        velocity: np.ndarray
    ) -> Tuple[KineticParameters, dict]:
        """
        Fit Lineweaver-Burk model using linear regression.
        
        Parameters
        ----------
        substrate_conc : np.ndarray
            Substrate concentrations (µM)
        velocity : np.ndarray
            Measured velocities
        
        Returns
        -------
        params : KineticParameters
            Fitted parameters
        fit_data : dict
            Fit information including slope, intercept, and R²
        """
        from scipy.stats import linregress
        
        # Convert to reciprocal space
        substrate_conc_inv = 1 / substrate_conc
        velocity_inv = 1 / velocity
        
        # Linear regression
        slope, intercept, r_value, _, _ = linregress(substrate_conc_inv, velocity_inv)
        r2 = r_value ** 2
        
        # Back-calculate Michaelis-Menten parameters
        if intercept != 0:
            Vmax_fit = 1 / intercept
            Km_fit = slope * Vmax_fit
        else:
            Vmax_fit = np.nan
            Km_fit = np.nan
        
        # Estimate standard errors (approximate)
        velocity_inv_pred = slope * substrate_conc_inv + intercept
        residuals = velocity_inv - velocity_inv_pred
        se = np.sqrt(np.sum(residuals ** 2) / (len(substrate_conc_inv) - 2))
        
        Vmax_std = np.nan
        Km_std = np.nan
        
        params = KineticParameters(
            Km=Km_fit,
            Km_std=Km_std,
            Vmax=Vmax_fit,
            Vmax_std=Vmax_std,
            r2=r2,
            n_points=len(substrate_conc)
        )
        
        fit_data = {
            'slope': slope,
            'intercept': intercept,
            'substrate_conc_inv': substrate_conc_inv,
            'velocity_inv_measured': velocity_inv,
            'velocity_inv_predicted': velocity_inv_pred,
            'residuals': residuals,
            'standard_error': se
        }
        
        return params, fit_data


class SubstrateInhibitionModel(EnzymeKineticModel):
    """
    Substrate inhibition kinetics model.
    
    Equation: v = (Vmax * [S]) / (Km + [S] + [S]²/Ki)
    
    Useful when enzyme activity decreases at high substrate concentrations.
    """
    
    def __call__(
        self,
        substrate_conc: np.ndarray,
        Vmax: float,
        Km: float,
        Ki: float
    ) -> np.ndarray:
        """
        Substrate inhibition equation.
        
        Parameters
        ----------
        substrate_conc : np.ndarray
            Substrate concentrations (µM)
        Vmax : float
            Maximum reaction velocity
        Km : float
            Michaelis constant
        Ki : float
            Inhibition constant
        
        Returns
        -------
        velocity : np.ndarray
            Reaction velocities
        """
        return (Vmax * substrate_conc) / (Km + substrate_conc + (substrate_conc ** 2) / Ki)
    
    def fit(
        self,
        substrate_conc: np.ndarray,
        velocity: np.ndarray,
        velocity_std: np.ndarray | None = None,
        Vmax_init: float | None = None,
        Km_init: float | None = None,
        Ki_init: float | None = None,
        maxfev: int = 5000
    ) -> Tuple[KineticParameters, dict]:
        """
        Fit substrate inhibition model.
        
        Parameters
        ----------
        substrate_conc : np.ndarray
            Substrate concentrations (µM)
        velocity : np.ndarray
            Measured velocities
        velocity_std : np.ndarray, optional
            Standard deviations
        Vmax_init : float, optional
            Initial guess for Vmax
        Km_init : float, optional
            Initial guess for Km
        Ki_init : float, optional
            Initial guess for Ki
        maxfev : int
            Maximum function evaluations
        
        Returns
        -------
        params : KineticParameters (note: Ki stored in n_points field)
        fit_data : dict
            Extended fit data including Ki parameter
        """
        
        if Vmax_init is None:
            Vmax_init = np.max(velocity)
        if Km_init is None:
            Km_init = np.median(substrate_conc)
        if Ki_init is None:
            Ki_init = np.max(substrate_conc)
        
        p0 = [Vmax_init, Km_init, Ki_init]
        
        popt, pcov = curve_fit(
            self,
            substrate_conc,
            velocity,
            p0=p0,
            maxfev=maxfev,
            sigma=velocity_std if velocity_std is not None and np.any(velocity_std > 0) else None,
            absolute_sigma=True
        )
        
        Vmax_fit, Km_fit, Ki_fit = popt
        perr = np.sqrt(np.diag(pcov))
        Vmax_std, Km_std, Ki_std = perr
        
        velocity_pred = self(substrate_conc, Vmax_fit, Km_fit, Ki_fit)
        r2 = self.calculate_r2(velocity, velocity_pred)
        residuals = velocity - velocity_pred
        
        params = KineticParameters(
            Km=Km_fit,
            Km_std=Km_std,
            Vmax=Vmax_fit,
            Vmax_std=Vmax_std,
            r2=r2,
            n_points=len(substrate_conc)
        )
        
        fit_data = {
            'predicted': velocity_pred,
            'residuals': residuals,
            'covariance': pcov,
            'Ki': Ki_fit,
            'Ki_std': Ki_std,
            'substrate_conc': substrate_conc,
            'velocity_measured': velocity,
            'velocity_std': velocity_std
        }
        
        return params, fit_data
