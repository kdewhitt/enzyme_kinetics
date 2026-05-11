"""Enzyme Kinetics Analysis Package

A comprehensive, reusable framework for analyzing enzyme kinetics from HPLC data.

Key Classes:
    - MichaelisMentenModel: Non-linear Michaelis-Menten fitting
    - LineweaverBurkModel: Lineweaver-Burk linearization
    - SubstrateInhibitionModel: Kinetics with substrate inhibition
    - CalibrationCurve: HPLC peak area to concentration calibration
    - EnzymeKineticsAnalyzer: High-level analysis interface
    - KineticsPlotter: Publication-quality visualization

Example Usage:
    >>> from enzyme_kinetics import EnzymeKineticsAnalyzer
    >>> analyzer = EnzymeKineticsAnalyzer(
    ...     enzyme_concentration_um=12.25,
    ...     reaction_time_seconds=10800  # 3 hours
    ... )
    >>> analyzer.analyze_peak('PDAL', substrate_conc, velocity)
    >>> df_results = analyzer.results_to_dataframe()
"""

from models import (
    MichaelisMentenModel,
    LineweaverBurkModel,
    SubstrateInhibitionModel,
    KineticParameters,
    EnzymeKineticModel,
)

from calibration import CalibrationCurve, CalibrationParameters

from analyzer import EnzymeKineticsAnalyzer, EnzymeKineticsResult

from plotting import KineticsPlotter

__version__ = "1.0.0"

__all__ = [
    "MichaelisMentenModel",
    "LineweaverBurkModel",
    "SubstrateInhibitionModel",
    "KineticParameters",
    "EnzymeKineticModel",
    "CalibrationCurve",
    "CalibrationParameters",
    "EnzymeKineticsAnalyzer",
    "EnzymeKineticsResult",
    "KineticsPlotter",
]
