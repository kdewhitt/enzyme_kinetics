#!/usr/bin/env python3
"""
Michaelis-Menten enzyme kinetics analysis for HPLC data.

Calculates kinetic parameters (Vmax, Km, kcat) for each substrate-peak combination
using non-linear least squares fitting (Levenberg-Marquardt) and Lineweaver-Burk linearization.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import warnings

warnings.filterwarnings('ignore')


def michaelis_menten(S, Vmax, Km):
    """
    Michaelis-Menten equation: v = (Vmax * S) / (Km + S)
    
    Parameters
    ----------
    S : array-like
        Substrate concentration
    Vmax : float
        Maximum reaction velocity
    Km : float
        Michaelis constant (substrate concentration at 1/2 Vmax)
    
    Returns
    -------
    v : ndarray
        Reaction velocity
    """
    return (Vmax * S) / (Km + S)


def lineweaver_burk_line(S_inv, Vmax, Km):
    """
    Lineweaver-Burk linearization: 1/v = (Km/Vmax) * (1/S) + 1/Vmax
    
    Parameters
    ----------
    S_inv : array-like
        Inverse substrate concentration (1/S)
    Vmax : float
        Maximum reaction velocity
    Km : float
        Michaelis constant
    
    Returns
    -------
    v_inv : ndarray
        Inverse reaction velocity (1/v)
    """
    return (Km / Vmax) * S_inv + 1 / Vmax


def analyze_enzyme_kinetics(data_path, output_dir='./'):
    """
    Analyze Michaelis-Menten kinetics for each peak in the HPLC data.
    
    Parameters
    ----------
    data_path : str
        Path to the CSV file containing HPLC data
    output_dir : str
        Directory to save output plots and results
    
    Returns
    -------
    results : pd.DataFrame
        Dataframe containing kinetic parameters for each peak
    """
    
    # Load data
    df = pd.read_csv(data_path)
    
    # Filter for HexCoA substrate only
    df_hexcoa = df[df['substrate'] == 'HexCoA'].copy()
    
    results_list = []
    
    # Group by peak_id and analyze each separately
    for peak_id, group_data in df_hexcoa.groupby('peak_id'):
        
        # Sort by substrate concentration
        group_data = group_data.sort_values('substrate_conc')
        
        # Extract substrate concentrations and mean areas (velocity proxy)
        S = group_data['substrate_conc'].values
        v = group_data['mean_area'].values
        v_std = group_data['std_area'].values
        
        # Skip peaks with all zeros or insufficient data
        if np.all(v == 0) or len(S) < 3:
            print(f"Skipping {peak_id}: insufficient data or all zero values")
            continue
        
        try:
            # Non-linear fitting (Levenberg-Marquardt)
            # Initial guess: Vmax ~ max velocity, Km ~ median concentration
            p0 = [np.max(v), np.median(S)]
            popt_mm, pcov_mm = curve_fit(
                michaelis_menten, S, v, p0=p0, maxfev=5000,
                sigma=v_std if np.any(v_std > 0) else None,
                absolute_sigma=True
            )
            
            Vmax_mm, Km_mm = popt_mm
            perr_mm = np.sqrt(np.diag(pcov_mm))
            
            # Calculate R² for non-linear fit
            v_pred_mm = michaelis_menten(S, *popt_mm)
            ss_res_mm = np.sum((v - v_pred_mm) ** 2)
            ss_tot = np.sum((v - np.mean(v)) ** 2)
            r2_mm = 1 - (ss_res_mm / ss_tot) if ss_tot != 0 else 0
            
            # Lineweaver-Burk linearization (1/v vs 1/S)
            S_inv = 1 / S
            v_inv = 1 / v
            
            slope_lb, intercept_lb, r_value_lb, p_value_lb, std_err_lb = linregress(S_inv, v_inv)
            
            # Back-calculate from Lineweaver-Burk
            Vmax_lb = 1 / intercept_lb if intercept_lb != 0 else np.nan
            Km_lb = slope_lb * Vmax_lb if not np.isnan(Vmax_lb) else np.nan
            
            # Calculate kcat if enzyme concentration is known (not available here, use area units)
            # For now, we'll report Vmax in terms of peak area units
            
            results_list.append({
                'peak_id': peak_id,
                'n_points': len(S),
                'Vmax_MM': Vmax_mm,
                'Vmax_MM_std': perr_mm[0],
                'Km_MM': Km_mm,
                'Km_MM_std': perr_mm[1],
                'R2_MM': r2_mm,
                'Vmax_LB': Vmax_lb,
                'Km_LB': Km_lb,
                'R2_LB': r_value_lb ** 2,
                'LB_slope': slope_lb,
                'LB_intercept': intercept_lb,
                'min_conc': S.min(),
                'max_conc': S.max(),
                'min_velocity': v.min(),
                'max_velocity': v.max(),
            })
            
            print(f"\n✓ {peak_id}")
            print(f"  Non-linear MM: Vmax={Vmax_mm:.2f}±{perr_mm[0]:.2f}, "
                  f"Km={Km_mm:.4f}±{perr_mm[1]:.4f}, R²={r2_mm:.4f}")
            print(f"  Lineweaver-Burk: Vmax={Vmax_lb:.2f}, Km={Km_lb:.4f}, R²={r_value_lb**2:.4f}")
            
        except Exception as e:
            print(f"✗ {peak_id}: Failed to fit - {str(e)}")
            continue
    
    results_df = pd.DataFrame(results_list)
    
    # Create comprehensive plots
    create_plots(df_hexcoa, results_df, output_dir)
    
    return results_df


def create_plots(df, results_df, output_dir):
    """
    Create publication-quality plots for each peak showing:
    - Michaelis-Menten curve with data points
    - Lineweaver-Burk linearization
    - Residuals
    """
    
    n_peaks = len(results_df)
    n_cols = 2
    n_rows = (n_peaks + n_cols - 1) // n_cols
    
    # Create multi-panel figure for Michaelis-Menten curves
    fig_mm, axes_mm = plt.subplots(n_rows, n_cols, figsize=(14, 5*n_rows))
    axes_mm = np.atleast_1d(axes_mm).flatten()
    
    # Create multi-panel figure for Lineweaver-Burk plots
    fig_lb, axes_lb = plt.subplots(n_rows, n_cols, figsize=(14, 5*n_rows))
    axes_lb = np.atleast_1d(axes_lb).flatten()
    
    for idx, (_, result) in enumerate(results_df.iterrows()):
        peak_id = result['peak_id']
        peak_data = df[df['peak_id'] == peak_id].sort_values('substrate_conc')
        
        S = peak_data['substrate_conc'].values
        v = peak_data['mean_area'].values
        v_std = peak_data['std_area'].values
        
        # Michaelis-Menten plot
        ax_mm = axes_mm[idx]
        S_smooth = np.linspace(S.min() * 0.5, S.max() * 1.2, 200)
        v_fit = michaelis_menten(S_smooth, result['Vmax_MM'], result['Km_MM'])
        
        ax_mm.errorbar(S, v, yerr=v_std, fmt='o', markersize=8, capsize=5,
                       capthick=2, label='Data', color='steelblue', elinewidth=2, markeredgewidth=2)
        ax_mm.plot(S_smooth, v_fit, '-', linewidth=2.5, label='Fit', color='darkred')
        ax_mm.axhline(y=result['Vmax_MM']/2, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        ax_mm.axvline(x=result['Km_MM'], color='gray', linestyle='--', alpha=0.5, linewidth=1)
        
        ax_mm.set_xlabel('Substrate Concentration', fontsize=11, fontweight='bold')
        ax_mm.set_ylabel('Velocity (Peak Area)', fontsize=11, fontweight='bold')
        ax_mm.set_title(f'{peak_id}: Michaelis-Menten\nVmax={result["Vmax_MM"]:.2f}, '
                        f'Km={result["Km_MM"]:.4f}, R²={result["R2_MM"]:.4f}',
                        fontsize=12, fontweight='bold')
        ax_mm.legend(fontsize=10)
        ax_mm.grid(True, alpha=0.3)
        
        # Lineweaver-Burk plot
        ax_lb = axes_lb[idx]
        S_inv = 1 / S
        v_inv = 1 / v
        
        # Fit line
        S_inv_smooth = np.linspace(S_inv.min(), S_inv.max(), 200)
        v_inv_fit = lineweaver_burk_line(S_inv_smooth, result['Vmax_MM'], result['Km_MM'])
        
        ax_lb.scatter(S_inv, v_inv, s=100, alpha=0.6, color='steelblue', edgecolors='navy', linewidth=2)
        ax_lb.plot(S_inv_smooth, v_inv_fit, '-', linewidth=2.5, color='darkred', label='Fit')
        
        # Annotate x and y intercepts
        ax_lb.axhline(y=1/result['Vmax_MM'], color='gray', linestyle='--', alpha=0.5, linewidth=1)
        ax_lb.axvline(x=-1/result['Km_MM'], color='gray', linestyle='--', alpha=0.5, linewidth=1)
        
        ax_lb.set_xlabel('1/[S]', fontsize=11, fontweight='bold')
        ax_lb.set_ylabel('1/v', fontsize=11, fontweight='bold')
        ax_lb.set_title(f'{peak_id}: Lineweaver-Burk\nR²={result["R2_LB"]:.4f}',
                        fontsize=12, fontweight='bold')
        ax_lb.legend(fontsize=10)
        ax_lb.grid(True, alpha=0.3)
    
    # Remove extra subplots
    for j in range(idx + 1, len(axes_mm)):
        fig_mm.delaxes(axes_mm[j])
    for j in range(idx + 1, len(axes_lb)):
        fig_lb.delaxes(axes_lb[j])
    
    fig_mm.tight_layout()
    fig_lb.tight_layout()
    
    fig_mm.savefig(f'{output_dir}/michaelis_menten_curves.png', dpi=300, bbox_inches='tight')
    fig_lb.savefig(f'{output_dir}/lineweaver_burk_plots.png', dpi=300, bbox_inches='tight')
    
    print(f"\n✓ Plots saved:")
    print(f"  - michaelis_menten_curves.png")
    print(f"  - lineweaver_burk_plots.png")
    
    plt.close('all')


if __name__ == '__main__':
    # Run analysis
    data_path = '/mnt/user-data/uploads/20231225_heatmap_area_statistics_hexcoa_new_format.csv'
    output_dir = '/mnt/user-data/outputs'
    
    print("=" * 70)
    print("MICHAELIS-MENTEN ENZYME KINETICS ANALYSIS")
    print("=" * 70)
    
    # Analyze kinetics
    results = analyze_enzyme_kinetics(data_path, output_dir)
    
    # Save results to CSV
    results_csv = f'{output_dir}/enzyme_kinetics_results.csv'
    results.to_csv(results_csv, index=False)
    print(f"\n✓ Results saved to enzyme_kinetics_results.csv\n")
    
    # Print summary table
    print("=" * 70)
    print("SUMMARY OF KINETIC PARAMETERS")
    print("=" * 70)
    summary_cols = ['peak_id', 'Vmax_MM', 'Km_MM', 'R2_MM', 'Vmax_LB', 'Km_LB', 'R2_LB']
    print(results[summary_cols].to_string(index=False))
    print("=" * 70)
