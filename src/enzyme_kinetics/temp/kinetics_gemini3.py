import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Load data
df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")

# Constants
E_t = 12.25  # Enzyme concentration in microMolar
time_min = 180  # 3 hours


# Models
def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)


def threshold_mm(S, Vmax, Km, S0):
    return np.where(S > S0, (Vmax * (S - S0)) / (Km + (S - S0)), 0)


results = []
peaks = df['peak_id'].unique()

# Create faceted plot
num_peaks = len(peaks)
cols = 3
rows = int(np.ceil(num_peaks / cols))
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

for i, peak in enumerate(peaks):
    ax = axes[i]
    sub_df = df[df['peak_id'] == peak].dropna(subset=['substrate_conc', 'mean_area'])

    if sub_df.empty or sub_df['mean_area'].max() == 0:
        ax.axis('off')
        continue

    S = sub_df['substrate_conc'].values
    v = sub_df['mean_area'].values / time_min
    v_sem = (sub_df['std_area'].values / time_min) / np.sqrt(sub_df['counts_area'].values)

    ax.errorbar(S, v, yerr=v_sem, fmt='o', label='Data', color='blue', capsize=3)

    try:
        # Use Specialized Model for Olivetol, Standard for others
        if peak == 'OLV':
            popt, pcov = curve_fit(threshold_mm, S, v, p0=[np.max(v) * 2, 1.0, 0.5], bounds=(0, np.inf))
            vmax, km, s0 = popt
            fit_label = f"Threshold Fit\n$S_0$: {s0:.3f}"
            S_fit = np.linspace(0, 3.5, 200)
            v_fit = threshold_mm(S_fit, *popt)
        else:
            popt, pcov = curve_fit(michaelis_menten, S, v, p0=[np.max(v), 0.1], bounds=(0, np.inf))
            vmax, km = popt
            fit_label = f"MM Fit\n$K_m$: {km:.3f}"
            S_fit = np.linspace(0, 3.5, 200)
            v_fit = michaelis_menten(S_fit, *popt)

        ax.plot(S_fit, v_fit, '--r', label=fit_label)
        ax.legend()
    except:
        pass

    ax.set_title(f'Kinetics: {peak}')
    ax.set_xlabel('[HexCoA] ($\mu M$)')
    ax.set_ylabel('Velocity (Area/min)')

plt.tight_layout()
plt.savefig('final_kinetics_report.png', dpi=300)
