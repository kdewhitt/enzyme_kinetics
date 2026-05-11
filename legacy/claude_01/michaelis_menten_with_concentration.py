#!/usr/bin/env python3
"""Michaelis-Menten enzyme kinetics analysis WITH PROTEIN CONCENTRATION NORMALIZATION.

Calculates kcat (turnover number) and kcat/Km (catalytic efficiency) for each enzyme.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import warnings

warnings.filterwarnings("ignore")


def michaelis_menten(S, Vmax, Km):
    """Michaelis-Menten equation: v = (Vmax * S) / (Km + S)"""
    return (Vmax * S) / (Km + S)


def analyze_enzyme_kinetics_with_concentration(
    data_path, enzyme_concentration_um=12.25, output_dir="./"
):
    """Analyze Michaelis-Menten kinetics with enzyme concentration normalization.

    Parameters
    ----------
    data_path : str
        Path to CSV file with HPLC data
    enzyme_concentration_um : float
        Total enzyme concentration in µM (default: 12.25)
    output_dir : str
        Output directory for results and plots

    Returns:
    -------
    results : pd.DataFrame
        Complete kinetic parameters including kcat and kcat/Km
    """
    # Load data
    df = pd.read_csv(data_path)
    df_hexcoa = df[df["substrate"] == "HexCoA"].copy()

    results_list = []

    print(f"\n{'=' * 80}")
    print("ENZYME KINETICS ANALYSIS WITH PROTEIN CONCENTRATION NORMALIZATION")
    print(f"{'=' * 80}")
    print(f"Enzyme Concentration: {enzyme_concentration_um} µM")
    print(f"{'=' * 80}\n")

    for peak_id, group_data in df_hexcoa.groupby("peak_id"):
        group_data = group_data.sort_values("substrate_conc")
        S = group_data["substrate_conc"].values
        v = group_data["mean_area"].values
        v_std = group_data["std_area"].values

        if np.all(v == 0) or len(S) < 3:
            print(f"⊘ {peak_id:12s} | Skipped (insufficient data)")
            continue

        try:
            # Non-linear fitting
            p0 = [np.max(v), np.median(S)]
            popt_mm, pcov_mm = curve_fit(
                michaelis_menten,
                S,
                v,
                p0=p0,
                maxfev=5000,
                sigma=v_std if np.any(v_std > 0) else None,
                absolute_sigma=True,
            )

            Vmax_mm, Km_mm = popt_mm
            perr_mm = np.sqrt(np.diag(pcov_mm))

            # Calculate R² for non-linear fit
            v_pred_mm = michaelis_menten(S, *popt_mm)
            ss_res_mm = np.sum((v - v_pred_mm) ** 2)
            ss_tot = np.sum((v - np.mean(v)) ** 2)
            r2_mm = 1 - (ss_res_mm / ss_tot) if ss_tot != 0 else 0

            # ===== NORMALIZE BY ENZYME CONCENTRATION =====
            # kcat = Vmax / [E]
            kcat = Vmax_mm / enzyme_concentration_um
            kcat_std = perr_mm[0] / enzyme_concentration_um

            # Catalytic efficiency
            if Km_mm > 0:  # Only valid for positive Km
                kcat_over_km = kcat / Km_mm
                kcat_over_km_std = kcat_over_km * np.sqrt(
                    (kcat_std / kcat) ** 2 + (perr_mm[1] / Km_mm) ** 2
                )
            else:
                kcat_over_km = np.nan
                kcat_over_km_std = np.nan

            # Lineweaver-Burk for comparison
            S_inv = 1 / S
            v_inv = 1 / v
            slope_lb, intercept_lb, r_value_lb, _, _ = linregress(S_inv, v_inv)

            if intercept_lb != 0:
                Vmax_lb = 1 / intercept_lb
                Km_lb = slope_lb * Vmax_lb
                kcat_lb = Vmax_lb / enzyme_concentration_um
            else:
                Vmax_lb = np.nan
                Km_lb = np.nan
                kcat_lb = np.nan

            results_list.append(
                {
                    "peak_id": peak_id,
                    "n_points": len(S),
                    # Michaelis-Menten parameters
                    "Vmax_MM": Vmax_mm,
                    "Vmax_MM_std": perr_mm[0],
                    "Km_MM": Km_mm,
                    "Km_MM_std": perr_mm[1],
                    "R2_MM": r2_mm,
                    # Normalized by enzyme concentration
                    "kcat_MM": kcat,
                    "kcat_MM_std": kcat_std,
                    "kcat_over_km_MM": kcat_over_km,
                    "kcat_over_km_MM_std": kcat_over_km_std,
                    # Lineweaver-Burk
                    "Vmax_LB": Vmax_lb,
                    "Km_LB": Km_lb,
                    "R2_LB": r_value_lb**2,
                    "kcat_LB": kcat_lb,
                    # Data range
                    "min_conc": S.min(),
                    "max_conc": S.max(),
                    "min_velocity": v.min(),
                    "max_velocity": v.max(),
                }
            )

            # Print results
            print(f"✓ {peak_id:12s}")
            print(f"  Vmax = {Vmax_mm:8.1f} ± {perr_mm[0]:6.1f} area/time")
            print(f"  Km   = {Km_mm:8.4f} ± {perr_mm[1]:8.4f} µM")
            print(f"  R²   = {r2_mm:8.4f}")
            print("  ")
            print(f"  kcat = {kcat:8.4f} ± {kcat_std:8.4f} area/(µM·enzyme·time)")
            print(f"  kcat/Km = {kcat_over_km:12.2f} (if Km > 0)")
            print(
                f"  Efficiency Ranking Factor: {kcat_over_km:.0f}"
                if kcat_over_km > 0
                else "  ⚠ Km < 0, see notes"
            )
            print()

        except Exception as e:
            print(f"✗ {peak_id:12s} | Failed: {str(e)}")
            continue

    results_df = pd.DataFrame(results_list)

    # Save results
    results_csv = f"{output_dir}/enzyme_kinetics_with_concentration.csv"
    results_df.to_csv(results_csv, index=False)

    # Create comparison plots
    create_comparison_plots(results_df, output_dir, enzyme_concentration_um)

    return results_df


def create_comparison_plots(results_df, output_dir, enzyme_conc):
    """Create publication-quality comparison plots."""
    # Filter for valid Km values
    valid_results = results_df[results_df["Km_MM"] > 0].copy()

    if len(valid_results) == 0:
        print("No valid results for plotting")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    peaks = valid_results["peak_id"].values
    colors = plt.cm.Set2(np.linspace(0, 1, len(peaks)))

    # ===== Plot 1: Km Comparison =====
    ax = axes[0, 0]
    bars1 = ax.bar(
        peaks,
        valid_results["Km_MM"].values,
        color=colors,
        alpha=0.7,
        edgecolor="black",
        linewidth=2,
    )
    ax.set_ylabel("Km (µM)", fontsize=12, fontweight="bold")
    ax.set_title(
        "Substrate Affinity (Lower = Higher Affinity)", fontsize=13, fontweight="bold"
    )
    ax.grid(axis="y", alpha=0.3)
    for i, (peak, val) in enumerate(zip(peaks, valid_results["Km_MM"].values)):
        ax.text(
            i,
            val,
            f"{val:.4f}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=10,
        )

    # ===== Plot 2: kcat Comparison =====
    ax = axes[0, 1]
    bars2 = ax.bar(
        peaks,
        valid_results["kcat_MM"].values,
        color=colors,
        alpha=0.7,
        edgecolor="black",
        linewidth=2,
    )
    ax.set_ylabel("kcat (area/(µM·enzyme·time))", fontsize=12, fontweight="bold")
    ax.set_title(
        "Turnover Number (Enzyme-Normalized Velocity)", fontsize=13, fontweight="bold"
    )
    ax.grid(axis="y", alpha=0.3)
    for i, (peak, val) in enumerate(zip(peaks, valid_results["kcat_MM"].values)):
        ax.text(
            i,
            val,
            f"{val:.2f}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=10,
        )

    # ===== Plot 3: kcat/Km Comparison =====
    ax = axes[1, 0]
    bars3 = ax.bar(
        peaks,
        valid_results["kcat_over_km_MM"].values,
        color=colors,
        alpha=0.7,
        edgecolor="black",
        linewidth=2,
    )
    ax.set_ylabel("kcat/Km (dimensionless ratio)", fontsize=12, fontweight="bold")
    ax.set_title(
        "Catalytic Efficiency (Combined Metric)", fontsize=13, fontweight="bold"
    )
    ax.grid(axis="y", alpha=0.3)
    for i, (peak, val) in enumerate(
        zip(peaks, valid_results["kcat_over_km_MM"].values)
    ):
        ax.text(
            i,
            val,
            f"{val:.0f}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=10,
        )

    # ===== Plot 4: Summary Table =====
    ax = axes[1, 1]
    ax.axis("off")

    summary_data = []
    summary_data.append(
        ["Peak", "Km\n(µM)", "kcat\n(norm)", "kcat/Km\n(eff)", "Ranking"]
    )

    # Sort by efficiency
    sorted_results = valid_results.sort_values("kcat_over_km_MM", ascending=False)
    for rank, (_, row) in enumerate(sorted_results.iterrows(), 1):
        summary_data.append(
            [
                row["peak_id"],
                f"{row['Km_MM']:.4f}",
                f"{row['kcat_MM']:.2f}",
                f"{row['kcat_over_km_MM']:.0f}",
                f"#{rank}",
            ]
        )

    table = ax.table(
        cellText=summary_data,
        cellLoc="center",
        loc="center",
        colWidths=[0.15, 0.15, 0.2, 0.2, 0.15],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2.5)

    # Style header row
    for i in range(5):
        table[(0, i)].set_facecolor("#4472C4")
        table[(0, i)].set_text_props(weight="bold", color="white")

    # Color data rows
    for i in range(1, len(summary_data)):
        for j in range(5):
            table[(i, j)].set_facecolor(colors[i - 1])
            table[(i, j)].set_alpha(0.5)

    ax.set_title("Efficiency Ranking", fontsize=13, fontweight="bold", pad=20)

    fig.suptitle(
        f"Enzyme Kinetics Comparison\n[Enzyme] = {enzyme_conc} µM",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )

    fig.tight_layout()
    fig.savefig(
        f"{output_dir}/enzyme_comparison_normalized.png", dpi=300, bbox_inches="tight"
    )

    print("✓ Comparison plot saved: enzyme_comparison_normalized.png")

    plt.close("all")


if __name__ == "__main__":
    data_path = (
        "/mnt/user-data/uploads/20231225_heatmap_area_statistics_hexcoa_new_format.csv"
    )
    output_dir = "/mnt/user-data/outputs"
    enzyme_conc = 12.25  # µM

    # Run analysis
    results = analyze_enzyme_kinetics_with_concentration(
        data_path, enzyme_concentration_um=enzyme_conc, output_dir=output_dir
    )

    print(f"\n{'=' * 80}")
    print("SUMMARY TABLE (Sorted by Catalytic Efficiency)")
    print(f"{'=' * 80}\n")

    # Display sorted results
    valid = results[results["Km_MM"] > 0].sort_values(
        "kcat_over_km_MM", ascending=False
    )

    display_cols = [
        "peak_id",
        "Km_MM",
        "Vmax_MM",
        "kcat_MM",
        "kcat_over_km_MM",
        "R2_MM",
    ]
    print(valid[display_cols].to_string(index=False))

    print(f"\n{'=' * 80}")
    print("INTERPRETATION GUIDE")
    print(f"{'=' * 80}")
    print("""
kcat = Vmax / [Enzyme]
    - Turnover number: reactions per enzyme molecule per unit time
    - Enzyme-normalized, comparable across experiments
    
kcat/Km = Catalytic Efficiency
    - Reflects both affinity (Km) and turnover rate (kcat)
    - Higher is better
    - Best comparative metric for enzyme performance
    
Enzyme with highest kcat/Km = Most efficient catalyst
    (Best combination of fast turnover + tight substrate binding)
    """)
    print(f"{'=' * 80}\n")
