import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
file_path = "resistivity_block_iter999999.dat"
start_line = 407849  # lines before this are skipped

# --- Read file ---
# Each valid data line has 6 whitespace-separated columns
df = pd.read_csv(
    file_path,
    delim_whitespace=True,
    header=None,
    skiprows=start_line - 1,
    names=["index", "sensitivity", "col3", "col4", "col5", "col6"],
)

# --- Extract sensitivity values ---
sens = df["sensitivity"].astype(float)

# --- Compute statistics ---
median_val = sens.median()
p10 = np.percentile(sens, 10)
p90 = np.percentile(sens, 90)
mini = np.min(sens)
maxi = np.max(sens)

print(f"Read {len(sens)} sensitivity values.")
print(f"Median: {median_val:.3e}, min: {mini:.3e}, max: {maxi:.3e}, P10: {p10:.3e}, P90: {p90:.3e}")


# --- Plot histogram ---
plt.figure(figsize=(8,5))
bins = np.linspace(p10, p90, 200)
plt.hist(sens, bins=bins, color="#5DADE2", edgecolor="black", alpha=0.7)
plt.xlabel("Normalized sensitivity density", fontsize=13)
plt.ylabel("Frequency", fontsize=13)
plt.title("Distribution of normalized sensitivity density", fontsize=14)
plt.xlim(p10, p90)

# --- Highlight quantiles ---
plt.axvline(median_val, color="red", linestyle="--", label=f"Median = {median_val:.2e}")
plt.axvline(p10, color="green", linestyle="--", label=f"P10 = {p10:.2e}")
plt.axvline(p90, color="orange", linestyle="--", label=f"P90 = {p90:.2e}")
plt.legend(fontsize=11)
plt.grid(alpha=0.3, linestyle=":")

plt.tight_layout()
plt.savefig('norm_sensitivity_density_histogram.pdf', format='pdf', dpi=300, bbox_inches='tight')
#plt.show()
