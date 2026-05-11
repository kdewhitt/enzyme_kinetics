#!/usr/bin/env python3
"""Calibration Curve Generation and Validation.

Instructions:
1. Create a CSV file with your CoA standard data: calibration_standards.csv
2. Columns needed: coa_concentration_um, peak_area_mean, peak_area_std
3. Run this script to generate calibration curve and constants
4. Use the constants in your kinetics analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os


def generate_calibration_curve(csv_path, output_dir="./"):
    """Generate and validate HPLC calibration curve for CoA quantitation.

    Parameters
    ----------
    csv_path : str
        Path to CSV with columns: coa_concentration_um, peak_area_mean, peak_area_std
    output_dir : str
        Output directory for plots and calibration file

    Returns:
    -------
    calibration_dict : dict
        Dictionary with slope, intercept, and other calibration parameters
    """
    # Load calibration data
    df = pd.read_csv(csv_path)

    coa_conc = df["coa_concentration_um"].values
    peak_area = df["peak_area_mean"].values
    peak_area_std = df["peak_area_std"].values

    print(f"\n{'=' * 80}")
    print("HPLC CALIBRATION CURVE GENERATION")
    print(f"{'=' * 80}\n")

    print(f"Loaded {len(coa_conc)} calibration points:")
    for i, (conc, area, std) in enumerate(zip(coa_conc, peak_area, peak_area_std)):
        print(f"  {i + 1}. [CoA] = {conc:6.2f} µM  |  Area = {area:8.1f} ± {std:7.1f}")

    # Define linear model
    def linear(x, slope, intercept):
        return slope * x + intercept

    # Fit calibration curve
    popt, pcov = curve_fit(
        linear, coa_conc, peak_area, sigma=peak_area_std, absolute_sigma=True
    )

    slope, intercept = popt
    slope_std, intercept_std = np.sqrt(np.diag(pcov))

    # Calculate R²
    area_pred = linear(coa_conc, *popt)
    ss_res = np.sum((peak_area - area_pred) ** 2)
    ss_tot = np.sum((peak_area - np.mean(peak_area)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    # Calculate residuals
    residuals = peak_area - area_pred
    residuals_normalized = residuals / peak_area_std

    print(f"\n{'=' * 80}")
    print("CALIBRATION FIT RESULTS")
    print(f"{'=' * 80}\n")

    print("Linear Model: peak_area = slope × [CoA] + intercept\n")
    print(f"  Slope:        {slope:.6f} ± {slope_std:.6f} area/µM")
    print(f"  Intercept:    {intercept:.6f} ± {intercept_std:.6f} area")
    print(f"  R²:           {r2:.8f}")
    print(f"  σ (residuals): {np.std(residuals):.2f} area units\n")

    print(f"Inverse Function: [CoA] = (peak_area - {intercept:.6f}) / {slope:.6f}\n")

    # Validate fit quality
    print(f"{'=' * 80}")
    print("FIT QUALITY ASSESSMENT")
    print(f"{'=' * 80}\n")

    if r2 > 0.99:
        print("✓ Excellent fit (R² > 0.99)")
    elif r2 > 0.95:
        print("✓ Good fit (R² > 0.95)")
    else:
        print(f"⚠ Adequate fit (R² = {r2:.4f})")

    # Check for outliers (normalized residuals > 3σ)
    outliers = np.abs(residuals_normalized) > 3
    if np.any(outliers):
        print(
            f"⚠ Warning: {np.sum(outliers)} potential outliers detected (|residuals| > 3σ)"
        )
        for idx in np.where(outliers)[0]:
            print(
                f"    Point {idx + 1}: [CoA]={coa_conc[idx]:.2f}, residual={residuals[idx]:.2f}"
            )
    else:
        print("✓ No outliers detected (all residuals < 3σ)")

    # Create comprehensive calibration plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Calibration curve with fit
    ax = axes[0, 0]
    ax.errorbar(
        coa_conc,
        peak_area,
        yerr=peak_area_std,
        fmt="o",
        markersize=10,
        capsize=6,
        capthick=2.5,
        color="steelblue",
        elinewidth=2,
        markeredgewidth=2.5,
        markeredgecolor="navy",
        label="Measurements",
        zorder=3,
    )

    coa_smooth = np.linspace(
        coa_conc.min() - coa_conc.max() * 0.1, coa_conc.max() * 1.15, 200
    )
    area_fit = linear(coa_smooth, *popt)
    ax.plot(coa_smooth, area_fit, "r-", linewidth=3, label="Linear fit", zorder=2)

    # Confidence bands (95%)
    se = np.sqrt(np.sum(residuals**2) / (len(coa_conc) - 2))
    ci = (
        1.96
        * se
        * np.sqrt(
            1 / len(coa_conc)
            + (coa_smooth - np.mean(coa_conc)) ** 2
            / np.sum((coa_conc - np.mean(coa_conc)) ** 2)
        )
    )
    ax.fill_between(
        coa_smooth, area_fit - ci, area_fit + ci, alpha=0.2, color="red", label="95% CI"
    )

    ax.set_xlabel("[CoA] (µM)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Peak Area", fontsize=12, fontweight="bold")
    ax.set_title(f"Calibration Curve\nR² = {r2:.6f}", fontsize=13, fontweight="bold")
    ax.legend(fontsize=11, loc="upper left")
    ax.grid(True, alpha=0.3)

    # Plot 2: Residuals
    ax = axes[0, 1]
    ax.errorbar(
        coa_conc,
        residuals,
        yerr=peak_area_std,
        fmt="o",
        markersize=10,
        capsize=6,
        capthick=2.5,
        color="coral",
        elinewidth=2,
        markeredgewidth=2.5,
        markeredgecolor="darkred",
        zorder=3,
    )
    ax.axhline(y=0, color="black", linestyle="-", linewidth=2, zorder=2)
    ax.axhline(
        y=3 * np.std(residuals),
        color="red",
        linestyle="--",
        linewidth=2,
        alpha=0.5,
        label="±3σ bounds",
    )
    ax.axhline(
        y=-3 * np.std(residuals), color="red", linestyle="--", linewidth=2, alpha=0.5
    )

    ax.set_xlabel("[CoA] (µM)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Residuals (area)", fontsize=12, fontweight="bold")
    ax.set_title("Residual Plot", fontsize=13, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Plot 3: Normalized residuals
    ax = axes[1, 0]
    ax.scatter(
        coa_conc,
        residuals_normalized,
        s=150,
        color="green",
        alpha=0.6,
        edgecolors="darkgreen",
        linewidth=2.5,
        zorder=3,
    )
    ax.axhline(y=0, color="black", linestyle="-", linewidth=2, zorder=2)
    ax.axhline(
        y=3, color="red", linestyle="--", linewidth=2, alpha=0.5, label="±3σ threshold"
    )
    ax.axhline(y=-3, color="red", linestyle="--", linewidth=2, alpha=0.5)

    ax.set_xlabel("[CoA] (µM)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Normalized Residuals (σ)", fontsize=12, fontweight="bold")
    ax.set_title("Standardized Residuals", fontsize=13, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([-4, 4])

    # Plot 4: Summary statistics table
    ax = axes[1, 1]
    ax.axis("off")

    summary_text = f"""
CALIBRATION PARAMETERS

Model: peak_area = slope × [CoA] + intercept

Slope:           {slope:.6f} ± {slope_std:.6f} area/µM
Intercept:       {intercept:.6f} ± {intercept_std:.6f} area

R²:              {r2:.8f}
RMSE:            {np.sqrt(np.mean(residuals**2)):.4f} area
Number of points: {len(coa_conc)}

INVERSE FUNCTION

[CoA] = (peak_area - {intercept:.6f}) / {slope:.6f}

USE IN KINETICS ANALYSIS:

AREA_TO_CONC_SLOPE = {slope:.6f}
AREA_TO_CONC_INTERCEPT = {intercept:.6f}

v_µM = (peak_area - INTERCEPT) / SLOPE
    """

    ax.text(
        0.05,
        0.95,
        summary_text,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="top",
        family="monospace",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    fig.suptitle(
        "HPLC Calibration Curve Analysis for CoA Quantitation",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    fig.tight_layout()

    calib_plot = os.path.join(output_dir, "hplc_calibration_curve.png")
    fig.savefig(calib_plot, dpi=300, bbox_inches="tight")
    print(f"\n✓ Calibration plot saved: {calib_plot}")

    plt.close("all")

    # Save calibration constants to file
    calib_file = os.path.join(output_dir, "calibration_constants.txt")
    with open(calib_file, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("HPLC CALIBRATION CONSTANTS FOR CoA QUANTITATION\n")
        f.write("=" * 80 + "\n\n")

        f.write("LINEAR MODEL:\n")
        f.write("  peak_area = slope × [CoA] + intercept\n\n")

        f.write("PARAMETERS:\n")
        f.write(f"  Slope:        {slope:.10f} area/µM\n")
        f.write(f"  Slope_std:    {slope_std:.10f}\n")
        f.write(f"  Intercept:    {intercept:.10f} area\n")
        f.write(f"  Intercept_std: {intercept_std:.10f}\n")
        f.write(f"  R²:           {r2:.10f}\n\n")

        f.write("INVERSE FUNCTION (for converting area to concentration):\n")
        f.write(f"  [CoA] (µM) = (peak_area - {intercept:.10f}) / {slope:.10f}\n\n")

        f.write("PYTHON CODE TO USE:\n")
        f.write("```python\n")
        f.write(f"AREA_TO_CONC_SLOPE = {slope:.10f}\n")
        f.write(f"AREA_TO_CONC_INTERCEPT = {intercept:.10f}\n\n")
        f.write("def area_to_product_concentration(area):\n")
        f.write('    """Convert HPLC peak area to CoA concentration (µM)."""\n')
        f.write(f"    return (area - {intercept:.10f}) / {slope:.10f}\n")
        f.write("```\n")

    print(f"✓ Calibration constants saved: {calib_file}\n")

    return {
        "slope": slope,
        "slope_std": slope_std,
        "intercept": intercept,
        "intercept_std": intercept_std,
        "r2": r2,
        "n_points": len(coa_conc),
    }


if __name__ == "__main__":
    # TODO: Replace with your actual calibration data file
    calibration_csv = "/path/to/calibration_standards.csv"
    output_dir = "/mnt/user-data/outputs"

    # Check if file exists (provide template if not)
    if not os.path.exists(calibration_csv):
        print(f"\n⚠ Calibration data file not found: {calibration_csv}")
        print("\nCreate a CSV file with the following columns:")
        print("  coa_concentration_um, peak_area_mean, peak_area_std")
        print("\nExample:")
        print("  coa_concentration_um,peak_area_mean,peak_area_std")
        print("  0.0,0.0,0.0")
        print("  0.1,245.3,12.4")
        print("  0.25,612.4,28.1")
        print("  0.5,1223.1,45.2")
        print("  1.0,2451.8,89.3")
        print("  2.0,4903.5,156.2")
        print("  5.0,12258.7,412.1")
        print("\nThen run this script with the correct path.")
    else:
        # Run calibration
        calib_params = generate_calibration_curve(calibration_csv, output_dir)

        print(f"{'=' * 80}")
        print("✓ CALIBRATION COMPLETE")
        print(f"{'=' * 80}\n")
        print("Use these values in your kinetics analysis script:\n")
        print(f"  AREA_TO_CONC_SLOPE = {calib_params['slope']:.10f}")
        print(f"  AREA_TO_CONC_INTERCEPT = {calib_params['intercept']:.10f}\n")
