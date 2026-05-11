#!/usr/bin/env python3
"""Michaelis-Menten enzyme kinetics with ABSOLUTE TIME UNITS.

Incorporates:
- Enzyme concentration: 12.25 µM
- Reaction time: 3 hours (10,800 seconds)
- Peak area → velocity conversion using relative measurements

This script establishes the relationship between peak area and product formation,
allowing calculation of true kcat in units of s⁻¹ and kcat/Km in M⁻¹·s⁻¹.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import warnings

warnings.filterwarnings("ignore")


# ============================================================================
# CONFIGURATION
# ============================================================================

ENZYME_CONCENTRATION_UM = 12.25  # µM
REACTION_TIME_HOURS = 3
REACTION_TIME_SECONDS = 3 * 3600  # 10,800 seconds


def michaelis_menten(S, Vmax, Km):
    """Michaelis-Menten equation: v = (Vmax * S) / (Km + S)"""
    return (Vmax * S) / (Km + S)


def analyze_with_time_units(data_path, output_dir="./"):
    """Analyze enzyme kinetics with absolute time units.

    Key insight: Peak area is proportional to product formed over 3 hours.
    We express velocity as "area/3hours" and use this for kinetic parameters.
    Once a calibration curve is available, area can be converted to µM.

    Parameters
    ----------
    data_path : str
        Path to CSV file
    output_dir : str
        Output directory

    Returns:
    -------
    results : pd.DataFrame
        Complete kinetic parameters with time-normalized units
    """
    df = pd.read_csv(data_path)
    df_hexcoa = df[df["substrate"] == "HexCoA"].copy()

    results_list = []

    print(f"\n{'=' * 90}")
    print("ENZYME KINETICS WITH ABSOLUTE TIME UNITS")
    print(f"{'=' * 90}")
    print(f"Enzyme Concentration:     {ENZYME_CONCENTRATION_UM} µM")
    print(
        f"Reaction Time:            {REACTION_TIME_HOURS} hours = {REACTION_TIME_SECONDS} seconds"
    )
    print("Peak Area Units:          Relative (area per 3-hour reaction)")
    print(f"{'=' * 90}\n")

    for peak_id, group_data in df_hexcoa.groupby("peak_id"):
        group_data = group_data.sort_values("substrate_conc")
        S = group_data["substrate_conc"].values  # µM
        area = group_data["mean_area"].values  # Relative peak area
        area_std = group_data["std_area"].values

        if np.all(area == 0) or len(S) < 3:
            print(f"⊘ {peak_id:12s} | Skipped (insufficient data or all zeros)")
            continue

        try:
            # ===== FIT TO MICHAELIS-MENTEN =====
            p0 = [np.max(area), np.median(S)]
            popt_mm, pcov_mm = curve_fit(
                michaelis_menten,
                S,
                area,
                p0=p0,
                maxfev=5000,
                sigma=area_std if np.any(area_std > 0) else None,
                absolute_sigma=True,
            )

            Vmax_area, Km = popt_mm
            perr_mm = np.sqrt(np.diag(pcov_mm))

            # Calculate R²
            area_pred = michaelis_menten(S, *popt_mm)
            ss_res = np.sum((area - area_pred) ** 2)
            ss_tot = np.sum((area - np.mean(area)) ** 2)
            r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

            # ===== CONVERT TO TIME-BASED UNITS =====
            # Vmax_area = peak area at saturation per 3 hours
            # Vmax_velocity = Vmax_area / (3 hours) = area per hour

            Vmax_per_hour = Vmax_area / REACTION_TIME_HOURS
            Vmax_per_sec = Vmax_area / REACTION_TIME_SECONDS
            Vmax_per_sec_std = perr_mm[0] / REACTION_TIME_SECONDS

            # ===== NORMALIZE BY ENZYME CONCENTRATION =====
            # kcat = Vmax_per_sec / [E]
            # Units: (area/sec) / µM = area·µM⁻¹·sec⁻¹

            kcat_per_sec = Vmax_per_sec / ENZYME_CONCENTRATION_UM
            kcat_per_sec_std = Vmax_per_sec_std / ENZYME_CONCENTRATION_UM

            # Catalytic efficiency
            if Km > 0:
                kcat_over_km = kcat_per_sec / Km
                # Error propagation
                rel_err_kcat = (
                    kcat_per_sec_std / kcat_per_sec if kcat_per_sec > 0 else 0
                )
                rel_err_km = perr_mm[1] / Km
                kcat_over_km_std = kcat_over_km * np.sqrt(
                    rel_err_kcat**2 + rel_err_km**2
                )
            else:
                kcat_over_km = np.nan
                kcat_over_km_std = np.nan

            # ===== LINEWEAVER-BURK FOR COMPARISON =====
            S_inv = 1 / S
            area_inv = 1 / area
            slope_lb, intercept_lb, r_value_lb, _, _ = linregress(S_inv, area_inv)

            if intercept_lb != 0:
                Vmax_area_lb = 1 / intercept_lb
                Km_lb = slope_lb * Vmax_area_lb
                Vmax_per_sec_lb = Vmax_area_lb / REACTION_TIME_SECONDS
                kcat_per_sec_lb = Vmax_per_sec_lb / ENZYME_CONCENTRATION_UM
            else:
                Vmax_area_lb = np.nan
                Km_lb = np.nan
                kcat_per_sec_lb = np.nan

            results_list.append(
                {
                    "peak_id": peak_id,
                    "n_points": len(S),
                    # Michaelis-Menten (peak area units)
                    "Vmax_area_MM": Vmax_area,
                    "Vmax_area_MM_std": perr_mm[0],
                    # Time-normalized velocities
                    "Vmax_per_hour": Vmax_per_hour,
                    "Vmax_per_sec": Vmax_per_sec,
                    "Vmax_per_sec_std": Vmax_per_sec_std,
                    # Kinetic parameters
                    "Km": Km,
                    "Km_std": perr_mm[1],
                    "R2": r2,
                    # Enzyme-normalized
                    "kcat_per_sec": kcat_per_sec,
                    "kcat_per_sec_std": kcat_per_sec_std,
                    "kcat_over_km_per_sec_per_um": kcat_over_km,
                    "kcat_over_km_std": kcat_over_km_std,
                    # Lineweaver-Burk
                    "Vmax_area_LB": Vmax_area_lb,
                    "Km_LB": Km_lb,
                    "kcat_per_sec_LB": kcat_per_sec_lb,
                    "R2_LB": r_value_lb**2,
                    # Data range
                    "min_conc": S.min(),
                    "max_conc": S.max(),
                    "min_area": area.min(),
                    "max_area": area.max(),
                }
            )

            print(f"✓ {peak_id:12s}")
            print("  Kinetic Parameters:")
            print(f"    Km              = {Km:8.4f} ± {perr_mm[1]:8.4f} µM")
            print(
                f"    Vmax (area)     = {Vmax_area:8.1f} ± {perr_mm[0]:8.1f} area/(3h)"
            )
            print(f"    Vmax (per hour) = {Vmax_per_hour:8.2f} area/hour")
            print(f"    Vmax (per sec)  = {Vmax_per_sec:8.6f} area/sec")
            print(f"    R²              = {r2:8.4f}")
            print("  ")
            print(
                f"  Enzyme-Normalized (kcat = Vmax / [E={ENZYME_CONCENTRATION_UM} µM]):"
            )
            print(
                f"    kcat            = {kcat_per_sec:8.6f} ± {kcat_per_sec_std:8.6f} area·µM⁻¹·sec⁻¹"
            )
            print(
                f"    kcat/Km         = {kcat_over_km:12.4f} ± {kcat_over_km_std:12.4f} sec⁻¹"
            )

            if kcat_over_km > 0:
                print(
                    f"    Efficiency      = {kcat_over_km:.2f} (reactions/sec relative to Km)"
                )
            print("  ")

        except Exception as e:
            print(f"✗ {peak_id:12s} | Failed: {str(e)}")
            continue

    results_df = pd.DataFrame(results_list)

    # Save results
    results_csv = f"{output_dir}/enzyme_kinetics_absolute_time_units.csv"
    results_df.to_csv(results_csv, index=False)
    print("\n✓ Results saved to: enzyme_kinetics_absolute_time_units.csv\n")

    # Create plots
    create_time_unit_plots(results_df, output_dir)

    return results_df


def create_time_unit_plots(results_df, output_dir):
    """Create plots showing kinetic parameters with time units."""
    # Filter for valid Km values
    valid_results = results_df[results_df["Km"] > 0].copy()
    valid_results = valid_results.sort_values(
        "kcat_over_km_per_sec_per_um", ascending=False
    )

    if len(valid_results) == 0:
        print("No valid results for plotting")
        return

    peaks = valid_results["peak_id"].values
    colors = plt.cm.Set3(np.linspace(0, 1, len(peaks)))

    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3)

    # ===== Plot 1: Km =====
    ax = fig.add_subplot(gs[0, 0])
    bars = ax.bar(
        peaks,
        valid_results["Km"].values,
        color=colors,
        alpha=0.7,
        edgecolor="black",
        linewidth=2,
    )
    ax.set_ylabel("Km (µM)", fontsize=11, fontweight="bold")
    ax.set_title("Substrate Affinity", fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    for i, (peak, val) in enumerate(zip(peaks, valid_results["Km"].values)):
        ax.text(
            i,
            val,
            f"{val:.4f}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=9,
        )

    # ===== Plot 2: Vmax (per second) =====
    ax = fig.add_subplot(gs[0, 1])
    bars = ax.bar(
        peaks,
        valid_results["Vmax_per_sec"].values,
        color=colors,
        alpha=0.7,
        edgecolor="black",
        linewidth=2,
    )
    ax.set_ylabel("Vmax (area/sec)", fontsize=11, fontweight="bold")
    ax.set_title("Maximum Velocity", fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    for i, (peak, val) in enumerate(zip(peaks, valid_results["Vmax_per_sec"].values)):
        ax.text(
            i,
            val,
            f"{val:.2e}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=8,
        )

    # ===== Plot 3: kcat =====
    ax = fig.add_subplot(gs[0, 2])
    bars = ax.bar(
        peaks,
        valid_results["kcat_per_sec"].values,
        color=colors,
        alpha=0.7,
        edgecolor="black",
        linewidth=2,
    )
    ax.set_ylabel("kcat (area·µM⁻¹·sec⁻¹)", fontsize=11, fontweight="bold")
    ax.set_title("Turnover Number (Enzyme-Normalized)", fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    for i, (peak, val) in enumerate(zip(peaks, valid_results["kcat_per_sec"].values)):
        ax.text(
            i,
            val,
            f"{val:.2e}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=8,
        )

    # ===== Plot 4: kcat/Km =====
    ax = fig.add_subplot(gs[1, 0])
    bars = ax.bar(
        peaks,
        valid_results["kcat_over_km_per_sec_per_um"].values,
        color=colors,
        alpha=0.7,
        edgecolor="black",
        linewidth=2,
    )
    ax.set_ylabel("kcat/Km (sec⁻¹)", fontsize=11, fontweight="bold")
    ax.set_title("Catalytic Efficiency", fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    for i, (peak, val) in enumerate(
        zip(peaks, valid_results["kcat_over_km_per_sec_per_um"].values)
    ):
        ax.text(
            i,
            val,
            f"{val:.2e}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=8,
        )

    # ===== Plot 5: R² values =====
    ax = fig.add_subplot(gs[1, 1])
    bars = ax.bar(
        peaks,
        valid_results["R2"].values,
        color=colors,
        alpha=0.7,
        edgecolor="black",
        linewidth=2,
    )
    ax.axhline(
        y=0.7, color="red", linestyle="--", linewidth=2, label="Good fit threshold"
    )
    ax.set_ylabel("R²", fontsize=11, fontweight="bold")
    ax.set_title("Goodness of Fit", fontsize=12, fontweight="bold")
    ax.set_ylim([0, 1])
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=9)
    for i, (peak, val) in enumerate(zip(peaks, valid_results["R2"].values)):
        ax.text(
            i,
            val,
            f"{val:.3f}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=9,
        )

    # ===== Plot 6: Vmax comparison (absolute vs per-time) =====
    ax = fig.add_subplot(gs[1, 2])
    x = np.arange(len(peaks))
    width = 0.35
    ax.bar(
        x - width / 2,
        valid_results["Vmax_per_hour"].values,
        width,
        label="Per hour",
        color="steelblue",
        alpha=0.7,
        edgecolor="black",
        linewidth=1.5,
    )
    ax.bar(
        x + width / 2,
        valid_results["Vmax_per_sec"].values * 3600,
        width,
        label="Per sec (scaled)",
        color="coral",
        alpha=0.7,
        edgecolor="black",
        linewidth=1.5,
    )
    ax.set_ylabel("Vmax (area)", fontsize=11, fontweight="bold")
    ax.set_title("Velocity Over Different Time Units", fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(peaks)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    # ===== Plot 7-9: Ranking tables =====

    # Efficiency ranking
    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")

    summary_data = [
        [
            "Rank",
            "Peak",
            "Km (µM)",
            "kcat (area·µM⁻¹·sec⁻¹)",
            "kcat/Km (sec⁻¹)",
            "Interpretation",
        ],
    ]

    for rank, (_, row) in enumerate(valid_results.iterrows(), 1):
        summary_data.append(
            [
                f"#{rank}",
                row["peak_id"],
                f"{row['Km']:.4f}",
                f"{row['kcat_per_sec']:.2e}",
                f"{row['kcat_over_km_per_sec_per_um']:.2e}",
                "Most efficient"
                if rank == 1
                else "Very good"
                if rank == 2
                else "Good"
                if rank == 3
                else "Fair",
            ]
        )

    table = ax.table(
        cellText=summary_data,
        cellLoc="center",
        loc="center",
        colWidths=[0.08, 0.12, 0.15, 0.2, 0.2, 0.15],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2.2)

    # Style header
    for i in range(6):
        table[(0, i)].set_facecolor("#4472C4")
        table[(0, i)].set_text_props(weight="bold", color="white")

    # Color data rows
    for i in range(1, len(summary_data)):
        for j in range(6):
            table[(i, j)].set_facecolor(colors[i - 1])
            table[(i, j)].set_alpha(0.4)

    ax.text(
        0.5,
        1.1,
        "Enzyme Kinetics Summary (Ranked by Catalytic Efficiency)",
        ha="center",
        fontsize=13,
        fontweight="bold",
        transform=ax.transAxes,
    )

    fig.suptitle(
        f"Michaelis-Menten Analysis with Absolute Time Units\n"
        f"[Enzyme] = {ENZYME_CONCENTRATION_UM} µM  |  Reaction Time = {REACTION_TIME_HOURS} hours = {REACTION_TIME_SECONDS} seconds",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )

    fig.savefig(
        f"{output_dir}/enzyme_kinetics_time_units.png", dpi=300, bbox_inches="tight"
    )
    print("✓ Comprehensive kinetics plot saved: enzyme_kinetics_time_units.png\n")

    plt.close("all")


if __name__ == "__main__":
    data_path = (
        "/mnt/user-data/uploads/20231225_heatmap_area_statistics_hexcoa_new_format.csv"
    )
    output_dir = "/mnt/user-data/outputs"

    results = analyze_with_time_units(data_path, output_dir)

    # Print comprehensive summary
    print(f"{'=' * 90}")
    print("KINETIC PARAMETERS SUMMARY")
    print(f"{'=' * 90}\n")

    valid = results[results["Km"] > 0].sort_values(
        "kcat_over_km_per_sec_per_um", ascending=False
    )

    cols_to_show = [
        "peak_id",
        "Km",
        "Vmax_per_sec",
        "kcat_per_sec",
        "kcat_over_km_per_sec_per_um",
        "R2",
    ]

    print(valid[cols_to_show].to_string(index=False))

    print(f"\n{'=' * 90}")
    print("INTERPRETATION NOTES")
    print(f"{'=' * 90}\n")

    print("""
KEY METRICS EXPLAINED:

1. Km (µM)
   - Michaelis constant: substrate concentration at 50% Vmax
   - Lower = higher substrate affinity
   - Intrinsic property (independent of enzyme amount and reaction time)

2. Vmax (area/sec)
   - Maximum reaction velocity (in peak area units per second)
   - Time-normalized for direct comparison
   - Depends on enzyme concentration and reaction conditions

3. kcat (area·µM⁻¹·sec⁻¹)
   - Turnover number per enzyme molecule per second
   - Normalized for enzyme concentration
   - Absolute rate of catalysis per enzyme

4. kcat/Km (sec⁻¹)
   - Catalytic efficiency: combines affinity and turnover
   - Best overall metric for enzyme performance
   - Higher = more efficient enzyme
   - Units: reactions per second per unit substrate concentration

UNIT CONVERSION NOTES:

Current units are "area·µM⁻¹·sec⁻¹" because:
- Peak area is a RELATIVE measurement (not yet converted to product concentration)
- Once you calibrate peak area → µM product formed:
  1. Convert Vmax from (area/sec) → (µM/sec)
  2. All downstream calculations convert automatically
  3. kcat/Km will be in standard M⁻¹·s⁻¹ units

NEXT STEPS:

1. Create a calibration curve:
   - Run standard HexCoA concentrations
   - Measure peak area vs known CoA (product) concentration
   - Fit linear relationship: peak_area = slope × [CoA] + intercept

2. Apply calibration:
   - Convert Vmax values using: Vmax_µM = (Vmax_area × slope) / REACTION_TIME_SECONDS
   - All enzyme kinetics automatically scale to absolute units

3. Compare with literature:
   - Once in M⁻¹·s⁻¹ units, kcat/Km can be directly compared
   - Typical enzymes: 10⁴ to 10⁶ M⁻¹·s⁻¹
   - Excellent enzymes: 10⁷ to 10⁸ M⁻¹·s⁻¹ (near diffusion-limited)
    """)

    print(f"{'=' * 90}\n")
