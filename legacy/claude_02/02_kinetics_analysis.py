#!/usr/bin/env python3
"""
Example: Complete Enzyme Kinetics Analysis

Demonstrates:
1. Loading HPLC data from CSV
2. Analyzing kinetics for multiple peaks
3. Applying calibration curve
4. Generating plots
5. Exporting results
"""

from pathlib import Path
import numpy as np
import pandas as pd
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from enzyme_kinetics import (
    EnzymeKineticsAnalyzer,
    CalibrationCurve,
    KineticsPlotter
)


def example_kinetics_analysis():
    """
    Complete enzyme kinetics analysis workflow.
    """
    
    print("="*80)
    print("ENZYME KINETICS ANALYSIS EXAMPLE")
    print("="*80 + "\n")
    
    # Configuration
    ENZYME_CONCENTRATION_UM = 12.25
    REACTION_TIME_HOURS = 3
    REACTION_TIME_SECONDS = REACTION_TIME_HOURS * 3600
    
    data_path = Path('/mnt/user-data/uploads/20231225_heatmap_area_statistics_hexcoa_new_format.csv')
    output_dir = Path('/mnt/user-data/outputs')
    
    # Step 1: Initialize analyzer
    print("Step 1: Initializing analyzer...")
    
    analyzer = EnzymeKineticsAnalyzer(
        enzyme_concentration_um=ENZYME_CONCENTRATION_UM,
        reaction_time_seconds=REACTION_TIME_SECONDS
    )
    
    print(f"  Enzyme concentration: {ENZYME_CONCENTRATION_UM} µM")
    print(f"  Reaction time: {REACTION_TIME_HOURS} hours\n")
    
    # Step 2: Load HPLC data
    print("Step 2: Loading HPLC data...")
    
    df = pd.read_csv(data_path)
    df_hexcoa = df[df['substrate'] == 'HexCoA'].copy()
    
    print(f"  Loaded {len(df_hexcoa)} data points")
    print(f"  Peaks: {df_hexcoa['peak_id'].unique().tolist()}\n")
    
    # Step 3: Analyze kinetics
    print("Step 3: Analyzing kinetics for each peak...")
    
    results = analyzer.analyze_dataframe(
        df_hexcoa,
        peak_id_column='peak_id',
        substrate_conc_column='substrate_conc',
        velocity_column='mean_area',
        velocity_std_column='std_area',
        fit_lineweaver_burk=True,
        min_points=3
    )
    
    print(f"  Successfully analyzed {len(results)} peaks\n")
    
    # Step 4: Display results
    print("Step 4: Results summary...")
    
    df_results = analyzer.results_to_dataframe()
    
    # Sort by kcat/Km if available
    if 'kcat_over_km' in df_results.columns:
        df_results_sorted = df_results[df_results['kcat_over_km'].notna()].sort_values(
            'kcat_over_km',
            ascending=False
        )
    else:
        df_results_sorted = df_results.sort_values('R2_MM', ascending=False)
    
    print(df_results_sorted[['peak_id', 'Km_MM', 'Vmax_MM', 'R2_MM', 'kcat_over_km']].to_string())
    print()
    
    # Step 5: Save results
    print("Step 5: Saving results...")
    
    results_csv = output_dir / 'example_analysis_results.csv'
    analyzer.save_results(results_csv)
    
    # Step 6: Generate plots
    print("\nStep 6: Generating plots...")
    
    plot_path = output_dir / 'example_plots'
    KineticsPlotter.plot_comparison_panel(results, output_path=plot_path)
    
    efficiency_plot = output_dir / 'example_efficiency_comparison.png'
    KineticsPlotter.plot_efficiency_comparison(results, output_path=efficiency_plot)
    
    # Step 7: Apply calibration (optional)
    print("\nStep 7: Applying calibration (optional)...")
    
    # Create example calibration curve
    calib = CalibrationCurve()
    calib_conc = np.array([0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0])
    calib_area = np.array([0, 245.3, 612.4, 1223.1, 2451.8, 4903.5, 12258.7])
    calib_std = np.array([0, 12.4, 28.1, 45.2, 89.3, 156.2, 412.1])
    
    calib.fit(calib_conc, calib_area, peak_areas_std=calib_std)
    
    print(f"  Calibration fit: {calib.params}")
    
    # Apply to results
    analyzer.apply_calibration(calib)
    print("  ✓ Calibration applied\n")
    
    # Step 8: Save calibrated results
    print("Step 8: Saving calibrated results...")
    
    calibrated_results = output_dir / 'example_calibrated_results.csv'
    analyzer.save_results(calibrated_results)
    print(f"  ✓ Saved to {calibrated_results}\n")
    
    print("="*80)
    print("✓ ANALYSIS COMPLETE")
    print("="*80 + "\n")
    
    return analyzer, results


def example_single_peak_analysis():
    """
    Demonstrates analysis of a single peak with detailed diagnostic output.
    """
    
    print("\n" + "="*80)
    print("SINGLE PEAK ANALYSIS EXAMPLE")
    print("="*80 + "\n")
    
    # Create sample data
    substrate_conc = np.array([0.1, 0.2, 0.5, 1.0, 3.0])
    velocity = np.array([70.97, 178.53, 297.77, 282.58, 1026.30])
    velocity_std = np.array([4.38, 2.25, 6.69, 9.59, 40.19])
    
    # Analyze
    analyzer = EnzymeKineticsAnalyzer(
        enzyme_concentration_um=12.25,
        reaction_time_seconds=10800
    )
    
    result = analyzer.analyze_peak(
        'HTAL_example',
        substrate_conc,
        velocity,
        velocity_std=velocity_std
    )
    
    print(f"Peak ID: {result.peak_id}")
    print(f"Michaelis-Menten Fit:")
    print(f"  Km = {result.mm_params.Km:.6f} ± {result.mm_params.Km_std:.6f} µM")
    print(f"  Vmax = {result.mm_params.Vmax:.4f} ± {result.mm_params.Vmax_std:.4f}")
    print(f"  R² = {result.mm_params.r2:.4f}")
    
    if result.kcat is not None:
        print(f"\nEnzyme-Normalized Parameters:")
        print(f"  kcat = {result.kcat:.6f} ± {result.kcat_std:.6f} sec⁻¹")
        print(f"  kcat/Km = {result.kcat_over_km:.4f} ± {result.kcat_over_km_std:.4f} sec⁻¹")
    
    # Generate individual plots
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    KineticsPlotter.plot_michaelis_menten_curve(result, ax=axes[0])
    if result.lb_params:
        KineticsPlotter.plot_lineweaver_burk(result, ax=axes[1])
    KineticsPlotter.plot_residuals(result, ax=axes[2])
    
    fig.suptitle('Single Peak Detailed Analysis', fontsize=14, fontweight='bold')
    fig.tight_layout()
    
    output_path = Path('/mnt/user-data/outputs') / 'example_single_peak_analysis.png'
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Plot saved to {output_path}")
    
    plt.close(fig)


if __name__ == '__main__':
    
    # Run full analysis
    analyzer, results = example_kinetics_analysis()
    
    # Run single peak example
    example_single_peak_analysis()
