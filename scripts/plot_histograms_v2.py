# -*- coding: utf-8 -*-
"""
Compare resistivity distributions of selected inversion parameters
between random_start and constant_start ensembles.

Each subplot shows:
 - constant_start ensemble (red histogram)
 - random_start ensemble (blue histogram)
 - vertical dashed lines at decile 1 and 9 (from constant_start)
 - station and geological labels from highlight_points

Author: charroyj
Date: 2025-11-07
"""

import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------
# Highlighted parameters (same as in profile plotting)
# ---------------------------------------------------------------------
highlight_points = {
    "Air4":       {12657: "A1: molasse", 43367: "A2: urgo."},
    "Ch2":        {25752: "C1: molasse", 49665: "C2: urgo."},
    "Canal2024":  {12208: "Canal1: molasse", 22011: "Canal2: urgo."},
    "Mand2_2024": {10096: "Mand2: urgo."},
    "MtAg":       {6643: "Ag1: urgo."},
    "Monod":      {21861: "Monod1: cond"}
}


cst_1_8_file='resistivity_block_iter13.dat'
rdm_81_file='../../random_start/random_start_81_1.5/resistivity_block_iter13.dat'
rdm_median_file='../../random_start/median_1.5/resistivity_block_iter7.dat'
cst_median_file='../median_1.5/resistivity_block_iter6.dat'

# Flatten all unique parameter IDs (up to 9)
selected_params = []
for sub in highlight_points.values():
    selected_params.extend(sub.keys())
selected_params = list(dict.fromkeys(selected_params))[:9]

def read_resistivity_block(block_file):
    """Read resistivity_block_iterXX.dat and return element → parameter + param table."""
    with open(block_file) as f:
        lines = f.readlines()

    n_elem, n_param = map(int, lines[0].split())

    mapping = np.loadtxt(lines[1 : 1 + n_elem], dtype=int)
    elem_to_param = dict(mapping)

    param_data = np.loadtxt(lines[1 + n_elem :], dtype=float)
    param_df = pd.DataFrame(param_data, columns=[
        "param_id", "resistivity", "lower", "upper", "nparam", "is_fixed"
    ])
    param_df["param_id"] = param_df["param_id"].astype(int)

    return elem_to_param, param_df


# ---------------------------------------------------------------------
# Load JSON distributions
# ---------------------------------------------------------------------
with open("../../random_start/cell_resistivities_log.json", "r") as f:
    cell_res_random = json.load(f)

with open("../cell_resistivities_log_filt.json", "r") as f:
    cell_res_const = json.load(f)

# ---------------------------------------------------------------------
# Determine global axis limits
# ---------------------------------------------------------------------
all_vals = []
for pid in selected_params:
    all_vals.extend(cell_res_random.get(str(pid), []))
    all_vals.extend(cell_res_const.get(str(pid), []))
all_vals = np.array(all_vals)

#x_min, x_max = np.percentile(all_vals, [1, 99])  # exclude outliers
x_min, x_max = [0.8,3.2]
y_max = 0  # computed dynamically per histogram

# ---------------------------------------------------------------------
# Prepare figure
# ---------------------------------------------------------------------
ncols, nrows = 3, 3
fig, axes = plt.subplots(nrows, ncols, figsize=(12, 10))
axes = axes.flatten()

for i, pid in enumerate(selected_params):
    ax = axes[i]

    # Get distributions
    r_vals = np.array(cell_res_random.get(str(pid), []))
    c_vals = np.array(cell_res_const.get(str(pid), []))

    if len(r_vals) == 0 or len(c_vals) == 0:
        ax.text(0.5, 0.5, f"No data for\nparam {pid}",
                ha="center", va="center", fontsize=12)
        ax.axis("off")
        continue

    # Common binning
    bins = np.linspace(x_min, x_max, 40)

    # Plot histograms
    n_const, _, _ = ax.hist(
        c_vals, bins=bins, color="red", alpha=0.6,
        label="Homogeneous start ensemble"
    )
    n_rand, _, _ = ax.hist(
        r_vals, bins=bins, color="green", alpha=0.4,
        label="Randomized start ensemble"
    )
    y_max = max(y_max, n_const.max(), n_rand.max())

    # Compute deciles (from constant_start)
    d1, d9 = np.percentile(c_vals, [10, 90])
    ax.axvline(d1, color="black", linestyle="--", lw=0.8,
               label="decile 1 (homog.)" if i == 0 else "")
    ax.axvline(d9, color="black", linestyle="--", lw=0.8,
               label="decile 9 (homog.)" if i == 0 else "")

    # extract and plot 1st reference resistivity - constant start model started with 10^1.8ohm.m and alpha=1.5
    _, df_ref1=read_resistivity_block(cst_1_8_file)
    ref_res_log1 = np.log10(float(df_ref1.loc[df_ref1["param_id"] == pid, "resistivity"].values[0]))
    ax.axvline(ref_res_log1,color="purple",linestyle="-",lw=0.8,zorder=20,
    label="m1: start = hom. 63 Ω·m" if i == 0 else "")
    
    # extract and plot 2nd reference resistivity - random start model 81 / alpha=1.5
    _, df_ref2=read_resistivity_block(rdm_81_file)
    ref_res_log2 = np.log10(float(df_ref2.loc[df_ref2["param_id"] == pid, "resistivity"].values[0]))
    ax.axvline(ref_res_log2,color="orange",linestyle="-",lw=0.8,zorder=20,
    label="m2: start = rdm. 81" if i == 0 else "")
    
    # extract and plot 2nd reference resistivity - output of inversions started with median of constant starts and alpha=1.5
    _, df_ref4=read_resistivity_block(cst_median_file)
    ref_res_log4 = np.log10(float(df_ref4.loc[df_ref4["param_id"] == pid, "resistivity"].values[0]))
    ax.axvline(ref_res_log4,color="royalblue",linestyle="-",lw=0.8,zorder=20,
    label="m3: start = median(hom. ens.)" if i == 0 else "")
        
    # extract and plot 2nd reference resistivity - output of inversions started with median of random starts and alpha=1.5
    _, df_ref3=read_resistivity_block(rdm_median_file)
    ref_res_log3 = np.log10(float(df_ref3.loc[df_ref3["param_id"] == pid, "resistivity"].values[0]))
    ax.axvline(ref_res_log3,color="red",linestyle="-",lw=0.8,zorder=20,
    label="m4: start = median(rdm. ens.)" if i == 0 else "")
    

    # ------------------------------------------------------------------
    # Add resistivity value annotations in top-left corner
    # ------------------------------------------------------------------
    # Values in linear Ω·m, no decimals
    text_lines = [
        (f"Decile 1: {10**d1:.0f} Ω·m", "black"),
        (f"m1 (Hom. 63 Ω·m start): {10**ref_res_log1:.0f} Ω·m", "purple"),
        (f"m2 (Rdm. #81 start): {10**ref_res_log2:.0f} Ω·m", "orange"),
        (f"m3 (median(hom.) start): {10**ref_res_log4:.0f} Ω·m", "royalblue"),
        (f"m4 (median(rdm.) start): {10**ref_res_log3:.0f} Ω·m", "red"),
        (f"Decile 9: {10**d9:.0f} Ω·m", "black")
    ]
    
    # Place as an anchored text block
    x_pos = 0.02  # left padding (fraction of axis width)
    y_start = 0.95  # start near top
    line_spacing = 0.06  # vertical spacing between lines
    
    for j, (txt, col) in enumerate(text_lines):
        ax.text(
            x_pos, y_start - j * line_spacing, txt,
            color=col, fontsize=10, transform=ax.transAxes,
            ha="left", va="top", fontweight="medium"
        )

    # Station + label in title
    st_name = next((st for st, pids in highlight_points.items() if pid in pids), "")
    label = highlight_points.get(st_name, {}).get(pid, "")
    ax.set_title(f"{st_name}: {label}\n(param {pid})", fontsize=14)

    ax.grid(True, linestyle="--", alpha=0.4)

# ---------------------------------------------------------------------
# Apply consistent axis limits and label visibility
# ---------------------------------------------------------------------
for j, ax in enumerate(axes):
    if not ax.has_data():
        ax.axis("off")
        continue

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, y_max * 1.1)

    row, col = divmod(j, ncols)

    # Show x labels only on bottom row
    if row == nrows - 1:
        ax.set_xlabel("log₁₀(Resistivity Ω·m)", fontsize=12)
    else:
        ax.set_xlabel("")
        ax.set_xticklabels([])

    # Show y labels only on leftmost column
    if col == 0:
        ax.set_ylabel("Count", fontsize=12)
    else:
        ax.set_ylabel("")
        ax.set_yticklabels([])

# ---------------------------------------------------------------------
# Shared legend and layout
# ---------------------------------------------------------------------
# --- Shared legend and layout ---
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(
    handles, labels,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.07),   # move legend further below the subplots
    ncol=4,
    fontsize=14,
    frameon=False
)

plt.tight_layout(rect=[0, 0, 1, 0.94])  # leave extra space at bottom for legend
plt.savefig("histograms_highlighted_params_v2.pdf", dpi=300, bbox_inches="tight")
plt.show()