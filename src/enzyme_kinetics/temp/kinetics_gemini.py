import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load data
df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")


# # Define the Michaelis-Menten equation
# def michaelis_menten(S, Vmax, Km):
#     return (Vmax * S) / (Km + S)


results = []
peaks = df['peak_id'].unique()

# Initialize plot
fig, ax = plt.subplots(figsize=(10, 6))

for peak in peaks:
    # Filter data by peak and drop any missing values
    sub_df = df[df['peak_id'] == peak].dropna(subset=['substrate_conc', 'percent_activity'])

    # Skip if empty or if max velocity is zero
    if sub_df.empty or sub_df['percent_activity'].max() == 0:
        continue

    S = sub_df['substrate_conc'].values
    v = sub_df['percent_activity'].values  # Change to 'mean_area' if you prefer area as velocity

    try:
        # Initial guesses: Vmax is the maximum velocity, Km is the concentration at half Vmax
        vmax_guess = np.max(v)
        km_guess = S[np.argmin(np.abs(v - vmax_guess / 2))]

        # Fit curve using scipy
        popt, pcov = curve_fit(michaelis_menten, S, v, p0=[vmax_guess, km_guess], bounds=(0, np.inf))
        vmax, km = popt

        results.append({'peak_id': peak, 'Vmax': vmax, 'Km': km})

        # Plotting the data points and the fitted line
        color = next(ax._get_lines.prop_cycler)['color']
        ax.scatter(S, v, label=f'{peak} Data', color=color)

        # Generate points for a smooth fit line
        S_fit = np.linspace(0, np.max(S) * 1.1, 100)
        v_fit = michaelis_menten(S_fit, vmax, km)
        ax.plot(S_fit, v_fit, label=f'{peak} Fit ($K_m$={km:.2f})', linestyle='--', color=color)

    except Exception as e:
        print(f"Could not fit {peak}: {e}")

# Format the plot
ax.set_xlabel('Substrate Concentration [HexCoA]')
ax.set_ylabel('Velocity (Percent Activity)')
ax.set_title('Michaelis-Menten Enzyme Kinetics')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the outputs
plt.savefig('kinetics_plot.png', dpi=300)
results_df = pd.DataFrame(results)
results_df.to_csv('kinetics_results.csv', index=False)

print("Kinetics parameters:")
print(results_df)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Load data
df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")


# Define the Michaelis-Menten equation
def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)


results = []
peaks = df['peak_id'].unique()

# Create faceted plots (2 rows, 3 columns to fit 5 plots)
num_peaks = len(peaks)
cols = 3
rows = int(np.ceil(num_peaks / cols))
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

for i, peak in enumerate(peaks):
    ax = axes[i]

    # Filter data by peak and drop any missing values
    sub_df = df[df['peak_id'] == peak].dropna(subset=['substrate_conc', 'mean_area'])

    # Skip if empty or if max velocity is zero
    if sub_df.empty or sub_df['mean_area'].max() == 0:
        ax.set_title(f"{peak} - No Data/Zero Velocity")
        ax.axis('off')
        continue

    S = sub_df['substrate_conc'].values
    v = sub_df['mean_area'].values

    # Plot data points
    ax.scatter(S, v, label='Data', color='blue')
    ax.set_title(f'Michaelis-Menten Kinetics: {peak}')
    ax.set_xlabel('Substrate Concentration [HexCoA]')
    ax.set_ylabel('Velocity (Mean Area)')

    try:
        # Initial guesses
        vmax_guess = np.max(v)
        km_guess = S[np.argmin(np.abs(v - vmax_guess / 2))]

        # Fit curve using scipy
        popt, pcov = curve_fit(michaelis_menten, S, v, p0=[vmax_guess, km_guess], bounds=(0, np.inf), maxfev=10000)
        vmax, km = popt

        results.append({'peak_id': peak, 'Vmax': vmax, 'Km': km})

        # Generate points for a smooth fit line
        S_fit = np.linspace(0, np.max(S) * 1.1, 100)
        v_fit = michaelis_menten(S_fit, vmax, km)

        # Format labels nicely depending on magnitude
        vmax_label = f"{vmax:.2e}" if vmax > 10000 else f"{vmax:.2f}"
        km_label = f"{km:.2e}" if km > 10000 or km < 0.001 else f"{km:.4f}"

        ax.plot(S_fit, v_fit, label=f'Fit\n$V_{{max}}$={vmax_label}\n$K_m$={km_label}', linestyle='--', color='red')
        ax.legend()

    except Exception as e:
        print(f"Could not fit {peak}: {e}")
        ax.text(
            0.5,
            0.5,
            'Fit Failed',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
        )

# Hide any extra empty subplots in the grid
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.savefig('kinetics_facet_plot.png', dpi=300)

results_df = pd.DataFrame(results)
results_df.to_csv('kinetics_mean_area_results.csv', index=False)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load data
df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")

# Define total protein concentration (in micro Molar)
E_t = 12.25


# Define the Michaelis-Menten equation
def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)


results = []
peaks = df['peak_id'].unique()

# Create faceted plots (2 rows, 3 columns to fit 5 plots)
num_peaks = len(peaks)
cols = 3
rows = int(np.ceil(num_peaks / cols))
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

for i, peak in enumerate(peaks):
    ax = axes[i]

    # Filter data by peak and drop any missing values
    sub_df = df[df['peak_id'] == peak].dropna(subset=['substrate_conc', 'mean_area'])

    # Skip if empty or if max velocity is zero
    if sub_df.empty or sub_df['mean_area'].max() == 0:
        ax.set_title(f"{peak} - No Data")
        ax.axis('off')
        continue

    S = sub_df['substrate_conc'].values
    v = sub_df['mean_area'].values

    # Plot data points
    ax.scatter(S, v, label='Data', color='blue')
    ax.set_title(f'Michaelis-Menten: {peak}')
    ax.set_xlabel('Substrate Concentration [HexCoA]')
    ax.set_ylabel('Velocity (Mean Area)')

    try:
        # Initial guesses
        vmax_guess = np.max(v)
        km_guess = S[np.argmin(np.abs(v - vmax_guess / 2))]

        # Fit curve using scipy
        popt, pcov = curve_fit(michaelis_menten, S, v, p0=[vmax_guess, km_guess], bounds=(0, np.inf), maxfev=10000)
        vmax, km = popt

        # Calculate k_cat and Catalytic Efficiency
        kcat = vmax / E_t
        efficiency = kcat / km if km > 0 else np.nan

        results.append(
            {
                'peak_id': peak,
                'Vmax': vmax,
                'Km': km,
                'kcat': kcat,
                'kcat/Km': efficiency,
            },
        )

        # Generate points for a smooth fit line
        S_fit = np.linspace(0, np.max(S) * 1.1, 100)
        v_fit = michaelis_menten(S_fit, vmax, km)

        # Format labels nicely depending on magnitude
        vmax_label = f"{vmax:.2e}" if vmax > 10000 else f"{vmax:.2f}"
        km_label = f"{km:.2e}" if km > 10000 or km < 0.001 else f"{km:.4f}"
        kcat_label = f"{kcat:.2e}" if kcat > 10000 else f"{kcat:.2f}"

        # Plot the fit line
        ax.plot(
            S_fit,
            v_fit,
            label=f'Fit\n$V_{{max}}$={vmax_label}\n$K_m$={km_label}\n$k_{{cat}}$={kcat_label}',
            linestyle='--',
            color='red',
        )
        ax.legend()

    except Exception as e:
        print(f"Could not fit {peak}: {e}")
        ax.text(
            0.5,
            0.5,
            'Fit Failed',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
        )

# Hide any extra empty subplots in the grid
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.savefig('kinetics_with_protein_plot.png', dpi=300)

# Save the full results table
results_df = pd.DataFrame(results)
results_df.to_csv('kinetics_with_protein_results.csv', index=False)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load data
df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")

# Constants
E_t = 12.25  # Enzyme concentration in microMolar
time_min = 3 * 60  # 3 hours converted to 180 minutes


def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)


results = []
peaks = df['peak_id'].unique()

# Create faceted plots
num_peaks = len(peaks)
cols = 3
rows = int(np.ceil(num_peaks / cols))
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

for i, peak in enumerate(peaks):
    ax = axes[i]

    sub_df = df[df['peak_id'] == peak].dropna(subset=['substrate_conc', 'mean_area'])

    if sub_df.empty or sub_df['mean_area'].max() == 0:
        ax.set_title(f"{peak} - No Data")
        ax.axis('off')
        continue

    S = sub_df['substrate_conc'].values

    # Calculate reaction velocity: Area per minute
    v = sub_df['mean_area'].values / time_min

    ax.scatter(S, v, label='Data', color='blue')
    ax.set_title(f'Michaelis-Menten: {peak}')
    ax.set_xlabel('Substrate Concentration [HexCoA]')
    ax.set_ylabel('Velocity (Mean Area / min)')

    try:
        # Initial guesses
        vmax_guess = np.max(v)
        km_guess = S[np.argmin(np.abs(v - vmax_guess / 2))]

        # Fit curve using scipy
        popt, pcov = curve_fit(michaelis_menten, S, v, p0=[vmax_guess, km_guess], bounds=(0, np.inf), maxfev=10000)
        vmax, km = popt

        # Calculate k_cat and Catalytic Efficiency
        kcat = vmax / E_t
        efficiency = kcat / km if km > 0 else np.nan

        results.append(
            {
                'peak_id': peak,
                'Vmax_area_per_min': vmax,
                'Km': km,
                'kcat_per_min': kcat,
                'kcat/Km': efficiency,
            },
        )

        S_fit = np.linspace(0, np.max(S) * 1.1, 100)
        v_fit = michaelis_menten(S_fit, vmax, km)

        vmax_label = f"{vmax:.2e}" if vmax > 10000 else f"{vmax:.2f}"
        km_label = f"{km:.2e}" if km > 10000 or km < 0.001 else f"{km:.4f}"
        kcat_label = f"{kcat:.2e}" if kcat > 10000 else f"{kcat:.4f}"

        ax.plot(
            S_fit,
            v_fit,
            label=f'Fit\n$V_{{max}}$={vmax_label}\n$K_m$={km_label}\n$k_{{cat}}$={kcat_label}',
            linestyle='--',
            color='red',
        )
        ax.legend()

    except Exception as e:
        print(f"Could not fit {peak}: {e}")
        ax.text(
            0.5,
            0.5,
            'Fit Failed',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
        )

# Hide any extra empty subplots in the grid
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.savefig('kinetics_time_adjusted_plot.png', dpi=300)

# Save the full results table
results_df = pd.DataFrame(results)
results_df.to_csv('kinetics_time_adjusted_results.csv', index=False)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load data
df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")

# Constants
E_t = 12.25  # Enzyme concentration in microMolar
time_min = 3 * 60  # 180 minutes


def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)


results = []
peaks = df['peak_id'].unique()

# Create faceted plots
num_peaks = len(peaks)
cols = 3
rows = int(np.ceil(num_peaks / cols))
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

for i, peak in enumerate(peaks):
    ax = axes[i]
    sub_df = df[df['peak_id'] == peak].dropna(subset=['substrate_conc', 'mean_area'])

    if sub_df.empty or sub_df['mean_area'].max() == 0:
        ax.set_title(f"{peak} - No Data")
        ax.axis('off')
        continue

    S = sub_df['substrate_conc'].values
    v = sub_df['mean_area'].values / time_min  # Velocity as Area/min

    ax.scatter(S, v, label='Data', color='blue', alpha=0.7)
    ax.set_title(f'Michaelis-Menten: {peak}', fontsize=14, fontweight='bold')
    ax.set_xlabel('Substrate Concentration [HexCoA]', fontsize=12)
    ax.set_ylabel('Velocity (Mean Area / min)', fontsize=12)

    try:
        vmax_guess = np.max(v)
        km_guess = S[np.argmin(np.abs(v - vmax_guess / 2))]

        # Fit curve
        popt, pcov = curve_fit(michaelis_menten, S, v, p0=[vmax_guess, km_guess], bounds=(0, np.inf), maxfev=10000)
        vmax, km = popt

        # Calculate Standard Errors (sqrt of diagonal of covariance matrix)
        perr = np.sqrt(np.diag(pcov))
        vmax_se, km_se = perr

        kcat = vmax / E_t
        kcat_se = vmax_se / E_t  # Propagation of error

        efficiency = kcat / km if km > 0 else np.nan

        results.append(
            {
                'peak_id': peak,
                'Vmax': vmax,
                'Vmax_SE': vmax_se,
                'Km': km,
                'Km_SE': km_se,
                'kcat': kcat,
                'kcat_SE': kcat_se,
                'kcat/Km': efficiency,
            },
        )

        # Plot the fit line
        S_fit = np.linspace(0, np.max(S) * 1.1, 100)
        v_fit = michaelis_menten(S_fit, vmax, km)


        # Helper for label formatting
        def fmt(val, se):
            return f"{val:.3e} ± {se:.3e}" if val > 100 or val < 0.01 else f"{val:.3f} ± {se:.3f}"


        label_text = (f"Fit:\n"
                      f"$V_{{max}}$: {fmt(vmax, vmax_se)}\n"
                      f"$K_m$: {fmt(km, km_se)}\n"
                      f"$k_{{cat}}$: {fmt(kcat, kcat_se)}")

        ax.plot(S_fit, v_fit, label=label_text, linestyle='--', color='red', linewidth=2)
        ax.legend(fontsize=10, loc='lower right')

    except Exception as e:
        ax.text(0.5, 0.5, f'Fit Failed', ha='center', va='center', transform=ax.transAxes)

# Cleanup plot layout
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.savefig('kinetics_with_se_plot.png', dpi=300)

# Save results to CSV
pd.DataFrame(results).to_csv('kinetics_with_se_results.csv', index=False)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load data
df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")

# Constants
E_t = 12.25  # Enzyme concentration in microMolar
time_min = 180  # 3 hours


def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)


peaks = df['peak_id'].unique()
num_peaks = len(peaks)
cols = 3
rows = int(np.ceil(num_peaks / cols))
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

for i, peak in enumerate(peaks):
    ax = axes[i]
    sub_df = df[df['peak_id'] == peak].dropna(subset=['substrate_conc', 'mean_area', 'std_area'])

    if sub_df.empty or sub_df['mean_area'].max() == 0:
        ax.axis('off')
        continue

    S = sub_df['substrate_conc'].values
    v = sub_df['mean_area'].values / time_min

    # Calculate SEM: SD / sqrt(n)
    v_sem = (sub_df['std_area'].values / time_min) / np.sqrt(sub_df['counts_area'].values)

    # Plot data with SEM error bar
    ax.errorbar(S, v, yerr=v_sem, fmt='o', label='Data (Mean ± SEM)', color='blue', alpha=0.7, capsize=3)

    try:
        popt, pcov = curve_fit(michaelis_menten, S, v, p0=[np.max(v), 0.1], bounds=(0, np.inf))
        vmax, km = popt
        perr = np.sqrt(np.diag(pcov))

        # Fit Line
        S_fit = np.linspace(0, np.max(S) * 1.1, 100)
        ax.plot(S_fit, michaelis_menten(S_fit, vmax, km), linestyle='--', color='red', label='Fit')
        ax.legend()
    except:
        pass

    ax.set_title(f'Michaelis-Menten: {peak}')
    ax.set_xlabel('Substrate Concentration [HexCoA]')
    ax.set_ylabel('Velocity (Area / min)')

plt.tight_layout()
plt.savefig('kinetics_with_errorbars_plot.png', dpi=300)

# =========== start here ===========
# =========== start here ===========
# =========== start here ===========

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 1. Load data
df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")

# 2. Define constants
E_t = 12.25  # Enzyme concentration in microMolar
time_min = 180  # 3 hour reaction time


# 3. Define the standard Michaelis-Menten model
def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)


peaks = df['peak_id'].unique()
num_peaks = len(peaks)
cols = 3
rows = int(np.ceil(num_peaks / cols))

# 4. Initialize the faceted plot
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

results = []

for i, peak in enumerate(peaks):
    ax = axes[i]

# Filter data for the specific peak
sub_df = df[df['peak_id'] == peak].dropna(subset=['substrate_conc', 'mean_area', 'std_area'])

if sub_df.empty or sub_df['mean_area'].max() == 0:
    ax.axis('off')
    continue

S = sub_df['substrate_conc'].values
v = sub_df['mean_area'].values / time_min  # Velocity in Area/min

# Calculate Standard Error of the Mean (SEM) for the points
v_sem = (sub_df['std_area'].values / time_min) / np.sqrt(sub_df['counts_area'].values)

# Plot data points with error bar
ax.errorbar(S, v, yerr=v_sem, fmt='o', label='Data (Mean ± SEM)', color='blue', alpha=0.7, capsize=3)

try:
    # 5. Fit the standard model to EVERY peak
    vmax_guess = np.max(v)
    km_guess = 0.1
    popt, pcov = curve_fit(michaelis_menten, S, v, p0=[vmax_guess, km_guess], bounds=(0, np.inf))
    vmax, km = popt

    # Calculate Standard Errors for parameters
    perr = np.sqrt(np.diag(pcov))
    vmax_se, km_se = perr

# Plot the fit line
S_fit = np.linspace(0, np.max(S) * 1.1, 100)
v_fit = michaelis_menten(S_fit, vmax, km)

# Formatting for the legend
label_text = (f"Fit:\n"
              f"$V_{{max}}$: {vmax:.2f}±{vmax_se:.2f}\n"
              f"$K_m$: {km:.2f}±{km_se:.2f}")

ax.plot(S_fit, v_fit, label=label_text, linestyle='--', color='red')
ax.legend(fontsize=9, loc='lower right')

except:
    ax.text(0.5, 0.5, 'Fit Failed', ha='center', va='center', transform=ax.transAxes)

ax.set_title(f'Standard MM: {peak}')
ax.set_xlabel('Substrate Concentration [HexCoA]')
ax.set_ylabel('Velocity (Area / min)')

# Clean up empty subplots
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 1. Load data and setup constants
df = pd.read_csv("20231225_heatmap_area_statistics_hexcoa_new_format.csv")
E_t = 12.25  # Protein concentration in microMolar
time_min = 180  # 3 hours reaction time


# 2. Define the two specialized kinetic models
def hill_equation(S, Vmax, K_half, n):
    """Hill equation for sigmoidal/cooperative kinetics."""
    return (Vmax * S ** n) / (K_half ** n + S ** n)


def threshold_mm(S, Vmax, Km, S0):
    """Michaelis-Menten model with a dead-zone threshold (S0)."""
    return np.where(S > S0, (Vmax * (S - S0)) / (Km + (S - S0)), 0)


# 3. Prepare OLV data
olv_df = df[df['peak_id'] == 'OLV'].dropna(subset=['substrate_conc', 'mean_area'])
S_data = olv_df['substrate_conc'].values
v_data = olv_df['mean_area'].values / time_min
v_sem = (olv_df['std_area'].values / time_min) / np.sqrt(olv_df['counts_area'].values)

results = []

# --- Fit Hill Model ---
try:
    popt_h, pcov_h = curve_fit(hill_equation, S_data, v_data, p0=[np.max(v_data), 1.0, 2.0], bounds=(0, np.inf))
    perr_h = np.sqrt(np.diag(pcov_h))
    results.append({'Model': 'Hill', 'Vmax': popt_h[0], 'Vmax_SE': perr_h[0], 'Param': popt_h[2], 'Param_Name': 'n'})
except:
    print("Hill fit failed")

# --- Fit Threshold MM Model ---
try:
    popt_t, pcov_t = curve_fit(threshold_mm, S_data, v_data, p0=[np.max(v_data), 0.5, 0.5], bounds=(0, np.inf))
    perr_t = np.sqrt(np.diag(pcov_t))
    results.append(
        {'Model': 'Threshold MM', 'Vmax': popt_t[0], 'Vmax_SE': perr_t[0], 'Param': popt_t[2], 'Param_Name': 'S0'},
    )
except:
    print("Threshold MM fit failed")

# 4. Generate Comparison Plot
plt.figure(figsize=(10, 6))
plt.errorbar(S_data, v_data, yerr=v_sem, fmt='ok', label='OLV Data (Mean ± SEM)', capsize=5)

S_plot = np.linspace(0, 3.5, 300)
if len(results) >= 1:
    plt.plot(
        S_plot, hill_equation(S_plot, popt_h[0], popt_h[1], popt_h[2]),
        'b-', label=f'Hill Fit (n={popt_h[2]:.2f})',
    )
if len(results) >= 2:
    plt.plot(
        S_plot, threshold_mm(S_plot, popt_t[0], popt_t[1], popt_t[2]),
        'r--', label=f'Threshold MM ($S_0$={popt_t[2]:.2f})',
    )

plt.title('Specialized Modeling for Olivetol (OLV)')
plt.xlabel('Substrate Concentration [HexCoA] ($\mu M$)')
plt.ylabel('Velocity (Area / min)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)
plt.savefig('olv_specialized_fit.png', dpi=300)

# Save results to CSV
pd.DataFrame(results).to_csv('olv_specialized_kinetics.csv', index=False)

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

# 1. Your actual olivetol (OLV) velocity data (Area/min)
# Substrate [S]: 0.1, 0.2, 0.5, 1.0, 3.0
# Velocity  [v]: 0.0, 0.0, 0.0, 0.19097, 1.22213
S_data = np.array([0.1, 0.2, 0.5, 1.0, 3.0])
v_data = np.array([0, 0, 0, 34.375, 219.985]) / 180


# def threshold_mm(S, Vmax, Km, S0):
#     """The threshold model: v=0 if S <= S0, otherwise standard MM."""
#     return np.where(S > S0, (Vmax * (S - S0)) / (Km + (S - S0)), 0)


# --- CASE 1: Resulting in S0 = 0.815 (The Best Fit) ---
# Using a higher initial guess for Vmax pushes the solver to consider a
# wider curve that can reach the S=1.0 point.
p0_best = [2.4, 1.0, 0.5]
popt_best, _ = curve_fit(threshold_mm, S_data, v_data, p0=p0_best, bounds=(0, np.inf))
ssr_best = np.sum((v_data - threshold_mm(S_data, *popt_best)) ** 2)

# --- CASE 2: Resulting in S0 = 1.745 (A Local Minimum) ---
# A lower initial guess for Vmax causes the solver to shift S0 to the right
# to 'ignore' the lower data points.
p0_local = [1.2, 0.5, 0.5]
popt_local, _ = curve_fit(threshold_mm, S_data, v_data, p0=p0_local, bounds=(0, np.inf))
ssr_local = np.sum((v_data - threshold_mm(S_data, *popt_local)) ** 2)

print(f"Fit 1 (Best): S0 = {popt_best[2]:.3f} | SSR = {ssr_best:.4e}")
print(f"Fit 2 (Local): S0 = {popt_local[2]:.3f} | SSR = {ssr_local:.4f}")
