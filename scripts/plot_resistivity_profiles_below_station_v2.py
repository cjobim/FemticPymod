# -*- coding: utf-8 -*-
"""
Plot vertical resistivity profiles (optionally with 1st/9th decile envelopes)
from FEMTIC parameter-averaged outputs (one line per inversion block).

Author: charroyj
Date: 2025-11-07
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------
df = pd.read_csv("profiles_resistivity.csv")
df2 = pd.read_csv("profiles_resistivity_rdm_81.csv")
df3 = pd.read_csv("profiles_resistivity_median_homogeneous.csv")
df4 = pd.read_csv("profiles_resistivity_median_rdm.csv")

stations = [
    "Ch2", "Air4", "Monod", "Canal2024",
    "Mand2_2024", "MtAg"
]

# ---------------------------------------------------------------------
# Parameters to highlight per station
# ---------------------------------------------------------------------
highlight_points = {
    "Air4":      {43367: "A2: urgo.", 12657: "A1: molasse"},
    "Ch2":        {25752: "C1: molasse", 49665: "C2: urgo."},
    "Canal2024":  {22011: "Canal2: urgo.", 12208: "Canal1: molasse"},
    "MtAg":       {6643: "Ag1: urgo."},
    "Monod":      {21861: "Monod1: cond"},
    "Mand2_2024": {10096: "Mand2: urgo."}
}

# ---------------------------------------------------------------------
# Compute uniform resistivity range
# ---------------------------------------------------------------------
sel = df[df["station"].isin(stations)]
res_min = sel["resistivity"].min()
res_max = sel["resistivity"].max()

# Round to nearest decade
res_min = 10 ** np.floor(np.log10(res_min))
res_max = 10 ** np.ceil(np.log10(res_max))

# ---------------------------------------------------------------------
# Detect spread columns (optional deciles)
# ---------------------------------------------------
has_dec9 = "resistivity_dec9" in df.columns
has_dec1 = "resistivity_dec1" in df.columns
has_spread = has_dec9 and has_dec1

if has_spread:
    print("→ 1st / 9th decile columns detected, will plot uncertainty envelopes.")

# ---------------------------------------------------------------------
# Create figure and subplots
# ---------------------------------------------------------------------
fig, axes = plt.subplots(
    1, len(stations), figsize=(3.2 * len(stations), 10), sharey=True
)

for ax, st in zip(axes, stations):
    d = df[df["station"] == st].sort_values("Z_km")
    d2 = df2[df2["station"] == st].sort_values("Z_km")
    d3 = df3[df3["station"] == st].sort_values("Z_km")
    d4 = df4[df4["station"] == st].sort_values("Z_km")
    if d.empty:
        print(f"⚠️  No data found for station {st}")
        continue

    Z_m = -d["Z_km"] * 1000  # convert km → m (depth)
    ax.semilogx(d["resistivity"], Z_m, color="purple", lw=1, label="m1: start = hom. 63 Ω·m")
    ax.semilogx(d2["resistivity"], Z_m, color="orange", lw=1, label="m2: start = rdm. #81")
    ax.semilogx(d3["resistivity"], Z_m, color="royalblue", lw=1, label="m3: start = median(hom. ens.)")
    ax.semilogx(d4["resistivity"], Z_m, color="red", lw=1, label="m4: start = median(rdm. ens.)")

    # -----------------------------------------------------------------
    # Spread envelopes (optional)
    # -----------------------------------------------------------------
    if has_spread:
        ax.fill_betweenx(
            Z_m,
            d["resistivity_dec1"],
            d["resistivity"],
            color="tomato",
            alpha=0.3,
            label="Range to decile 1" if ax == axes[0] else ""
        )
        ax.fill_betweenx(
            Z_m,
            d["resistivity"],
            d["resistivity_dec9"],
            color="royalblue",
            alpha=0.3,
            label="Range to decile 9" if ax == axes[0] else ""
        )

    # -----------------------------------------------------------------
    # Highlight parameters (red dots + labels)
    # -----------------------------------------------------------------
    if st in highlight_points:
        for pid, label in highlight_points[st].items():
            sel_param = d[d["param_id"] == pid]
            if not sel_param.empty:
                row = sel_param.iloc[0]
                x = row["resistivity"]
                y = -row["Z_km"] * 1000
                ax.plot(x, y, "o", color="red", markersize=6, zorder=10)
                ax.text(
                    x * 1.15, y, label,
                    color="red", fontsize=14,
                    va="center", ha="left"
                )
            else:
                print(f"⚠️  param_id {pid} not found for station {st}")

    # -----------------------------------------------------------------
    # Axes formatting
    # -----------------------------------------------------------------
    ax.set_xlim(res_min, res_max)
    ax.set_ylim(bottom=-3000, top=900)
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.set_title(st, fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.minorticks_on()

    if ax == axes[0]:
        ax.set_ylabel("Elevation (masl)", fontsize=14)
    ax.set_xlabel("Resistivity (Ω·m)", fontsize=14)

# ---------------------------------------------------------------------
# Common legend
# ---------------------------------------------------------------------
if has_spread:
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc="lower center",
        ncol=3,
        fontsize=14,
        frameon=False
    )
    fig.subplots_adjust(bottom=0.13)
else:
    fig.tight_layout(w_pad=0.3, h_pad=0.5)

plt.savefig("profiles_with_spreads_v2.pdf", format="pdf", dpi=300)
plt.show()