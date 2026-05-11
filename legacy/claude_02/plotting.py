"""
Plotting and visualization utilities for enzyme kinetics analysis.

Provides reusable plotting functions for kinetic curves, residuals,
comparisons, and diagnostic plots.
"""

from pathlib import Path
from typing import List, Optional

import numpy as np
import matplotlib.pyplot as plt

from analyzer import EnzymeKineticsResult, EnzymeKineticsAnalyzer


class KineticsPlotter:
    """Generate publication-quality kinetics plots."""
    
    @staticmethod
    def plot_michaelis_menten_curve(
        result: EnzymeKineticsResult,
        ax: Optional[plt.Axes] = None,
        show_fit: bool = True,
        show_annotations: bool = True
    ) -> plt.Axes:
        """
        Plot Michaelis-Menten curve with data points.
        
        Parameters
        ----------
        result : EnzymeKineticsResult
            Kinetics analysis result
        ax : plt.Axes, optional
            Existing axes to plot on. Creates new figure if None.
        show_fit : bool
            Plot fitted curve
        show_annotations : bool
            Show Km and Vmax annotations
        
        Returns
        -------
        ax : plt.Axes
            Matplotlib axes with plot
        """
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        
        # Data points
        ax.errorbar(
            result.substrate_conc,
            result.velocity,
            yerr=result.velocity_std,
            fmt='o',
            markersize=8,
            capsize=5,
            capthick=2,
            color='steelblue',
            elinewidth=2,
            markeredgewidth=2,
            markeredgecolor='navy',
            label='Data',
            zorder=3
        )
        
        # Fitted curve
        if show_fit:
            S_smooth = np.linspace(
                result.substrate_conc.min() * 0.5,
                result.substrate_conc.max() * 1.2,
                200
            )
            v_fit = (
                result.mm_params.Vmax * S_smooth /
                (result.mm_params.Km + S_smooth)
            )
            ax.plot(S_smooth, v_fit, 'r-', linewidth=2.5, label='Fit', zorder=2)
        
        # Annotations
        if show_annotations:
            ax.axhline(
                y=result.mm_params.Vmax / 2,
                color='gray',
                linestyle='--',
                alpha=0.5,
                linewidth=1
            )
            ax.axvline(
                x=result.mm_params.Km,
                color='gray',
                linestyle='--',
                alpha=0.5,
                linewidth=1
            )
        
        ax.set_xlabel('Substrate Concentration (µM)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Velocity (arbitrary units)', fontsize=11, fontweight='bold')
        ax.set_title(
            f'{result.peak_id}: Michaelis-Menten\n'
            f'Km={result.mm_params.Km:.4f}±{result.mm_params.Km_std:.4f} µM, '
            f'R²={result.mm_params.r2:.4f}',
            fontsize=12,
            fontweight='bold'
        )
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        
        return ax
    
    @staticmethod
    def plot_lineweaver_burk(
        result: EnzymeKineticsResult,
        ax: Optional[plt.Axes] = None,
        show_fit: bool = True
    ) -> plt.Axes:
        """
        Plot Lineweaver-Burk linearization.
        
        Parameters
        ----------
        result : EnzymeKineticsResult
            Kinetics analysis result
        ax : plt.Axes, optional
            Existing axes. Creates new if None.
        show_fit : bool
            Plot fitted line
        
        Returns
        -------
        ax : plt.Axes
            Matplotlib axes
        """
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        
        if result.lb_params is None:
            ax.text(0.5, 0.5, 'Lineweaver-Burk fit not available',
                   ha='center', va='center', transform=ax.transAxes)
            return ax
        
        # Data points
        S_inv = 1 / result.substrate_conc
        v_inv = 1 / result.velocity
        
        ax.scatter(
            S_inv, v_inv,
            s=100,
            alpha=0.6,
            color='steelblue',
            edgecolors='navy',
            linewidth=2,
            zorder=3
        )
        
        # Fitted line
        if show_fit:
            S_inv_smooth = np.linspace(S_inv.min(), S_inv.max(), 200)
            v_inv_fit = (
                result.lb_params.Km / result.lb_params.Vmax * S_inv_smooth +
                1 / result.lb_params.Vmax
            )
            ax.plot(S_inv_smooth, v_inv_fit, 'r-', linewidth=2.5, label='Fit', zorder=2)
        
        ax.set_xlabel('1/[S] (µM⁻¹)', fontsize=11, fontweight='bold')
        ax.set_ylabel('1/v', fontsize=11, fontweight='bold')
        ax.set_title(
            f'{result.peak_id}: Lineweaver-Burk\n'
            f'Km={result.lb_params.Km:.4f} µM, '
            f'R²={result.lb_params.r2:.4f}',
            fontsize=12,
            fontweight='bold'
        )
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        
        return ax
    
    @staticmethod
    def plot_residuals(
        result: EnzymeKineticsResult,
        ax: Optional[plt.Axes] = None
    ) -> plt.Axes:
        """
        Plot residuals from Michaelis-Menten fit.
        
        Parameters
        ----------
        result : EnzymeKineticsResult
        ax : plt.Axes, optional
        
        Returns
        -------
        ax : plt.Axes
        """
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        
        residuals = result.mm_fit_data['residuals']
        S = result.substrate_conc
        
        ax.errorbar(
            S, residuals,
            yerr=result.velocity_std if result.velocity_std is not None else None,
            fmt='o',
            markersize=8,
            capsize=5,
            capthick=2,
            color='coral',
            elinewidth=2,
            markeredgewidth=2,
            markeredgecolor='darkred',
            zorder=3
        )
        ax.axhline(y=0, color='black', linestyle='-', linewidth=2, zorder=2)
        
        ax.set_xlabel('Substrate Concentration (µM)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Residuals', fontsize=11, fontweight='bold')
        ax.set_title(f'{result.peak_id}: Residual Plot', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        return ax
    
    @staticmethod
    def plot_comparison_panel(
        results: List[EnzymeKineticsResult],
        output_path: Optional[Path | str] = None
    ) -> None:
        """
        Create multi-panel comparison of all results.
        
        Parameters
        ----------
        results : List[EnzymeKineticsResult]
            List of analysis results
        output_path : str or Path, optional
            Save figure to this path
        """
        
        # Filter for valid results
        valid_results = [r for r in results if r.mm_params.Km > 0]
        
        if not valid_results:
            print("No valid results to plot")
            return
        
        n_peaks = len(valid_results)
        n_cols = 2
        n_rows = (n_peaks + n_cols - 1) // n_cols
        
        # Michaelis-Menten curves
        fig_mm, axes_mm = plt.subplots(n_rows, n_cols, figsize=(14, 5*n_rows))
        axes_mm = np.atleast_1d(axes_mm).flatten()
        
        for idx, result in enumerate(valid_results):
            KineticsPlotter.plot_michaelis_menten_curve(result, ax=axes_mm[idx])
        
        # Remove unused subplots
        for j in range(len(valid_results), len(axes_mm)):
            fig_mm.delaxes(axes_mm[j])
        
        fig_mm.suptitle('Michaelis-Menten Curves', fontsize=14, fontweight='bold')
        fig_mm.tight_layout()
        
        if output_path is not None:
            output_path = Path(output_path)
            mm_path = output_path.parent / f"{output_path.stem}_michaelis_menten.png"
            fig_mm.savefig(mm_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved: {mm_path}")
        
        plt.close(fig_mm)
        
        # Lineweaver-Burk plots
        fig_lb, axes_lb = plt.subplots(n_rows, n_cols, figsize=(14, 5*n_rows))
        axes_lb = np.atleast_1d(axes_lb).flatten()
        
        for idx, result in enumerate(valid_results):
            if result.lb_params is not None:
                KineticsPlotter.plot_lineweaver_burk(result, ax=axes_lb[idx])
        
        for j in range(len(valid_results), len(axes_lb)):
            fig_lb.delaxes(axes_lb[j])
        
        fig_lb.suptitle('Lineweaver-Burk Linearization', fontsize=14, fontweight='bold')
        fig_lb.tight_layout()
        
        if output_path is not None:
            lb_path = output_path.parent / f"{output_path.stem}_lineweaver_burk.png"
            fig_lb.savefig(lb_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved: {lb_path}")
        
        plt.close(fig_lb)
    
    @staticmethod
    def plot_efficiency_comparison(
        results: List[EnzymeKineticsResult],
        output_path: Optional[Path | str] = None
    ) -> None:
        """
        Create efficiency comparison plots.
        
        Parameters
        ----------
        results : List[EnzymeKineticsResult]
            List of analysis results
        output_path : str or Path, optional
            Save to this path
        """
        
        # Filter for valid results with kcat/Km
        valid_results = [
            r for r in results
            if r.mm_params.Km > 0 and r.kcat_over_km is not None
        ]
        
        if not valid_results:
            print("No results with valid kcat/Km")
            return
        
        # Sort by efficiency
        valid_results = sorted(valid_results, key=lambda r: r.kcat_over_km, reverse=True)
        
        peaks = [r.peak_id for r in valid_results]
        colors = plt.cm.Set3(np.linspace(0, 1, len(peaks)))
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Km comparison
        ax = axes[0, 0]
        ax.bar(peaks, [r.mm_params.Km for r in valid_results],
               color=colors, alpha=0.7, edgecolor='black', linewidth=2)
        ax.set_ylabel('Km (µM)', fontsize=11, fontweight='bold')
        ax.set_title('Substrate Affinity', fontsize=12, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        # Vmax comparison
        ax = axes[0, 1]
        ax.bar(peaks, [r.mm_params.Vmax for r in valid_results],
               color=colors, alpha=0.7, edgecolor='black', linewidth=2)
        ax.set_ylabel('Vmax', fontsize=11, fontweight='bold')
        ax.set_title('Maximum Velocity', fontsize=12, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        # kcat comparison
        ax = axes[1, 0]
        ax.bar(peaks, [r.kcat for r in valid_results if r.kcat is not None],
               color=colors, alpha=0.7, edgecolor='black', linewidth=2)
        ax.set_ylabel('kcat', fontsize=11, fontweight='bold')
        ax.set_title('Turnover Number', fontsize=12, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        # kcat/Km comparison
        ax = axes[1, 1]
        ax.bar(peaks, [r.kcat_over_km for r in valid_results],
               color=colors, alpha=0.7, edgecolor='black', linewidth=2)
        ax.set_ylabel('kcat/Km', fontsize=11, fontweight='bold')
        ax.set_title('Catalytic Efficiency', fontsize=12, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        fig.suptitle('Enzyme Kinetics Comparison', fontsize=14, fontweight='bold')
        fig.tight_layout()
        
        if output_path is not None:
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved: {output_path}")
        else:
            plt.show()
        
        plt.close(fig)
