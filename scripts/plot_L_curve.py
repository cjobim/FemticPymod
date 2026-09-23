# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 12:37:37 2025

@author: charroyj
"""


#plot L curve by collating multiple femtic.cnv files contained in a directory structure with various alpha values
#run from "projects" directory
import os
import glob
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

base_dirs = {
    #"homogeneous starting model": "./constant_start",
    "randomized starting model": "./random_start"
}

all_data = []

for dist_flag, base in base_dirs.items():
    for folder in sorted(glob.glob(os.path.join(base, "inversion_*"))):
        print(f"processing folder: ${folder}")
        cnv_path = os.path.join(folder, "femtic.cnv")
        if not os.path.isfile(cnv_path):
            continue

        # Detect number of columns from first data line
        with open(cnv_path, "r") as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    n_cols = len(line.split())
                    break

        # Assign appropriate column names
        if n_cols == 8:
            colnames = ["Iter", "Retrial", "Alpha", "Damp", "Roughness", "Misfit", "RMS", "ObjFunc"]
        elif n_cols == 10:
            colnames = ["Iter", "Retrial", "Alpha", "Beta", "Damp",
                        "Roughness", "Distortion", "Misfit", "RMS", "ObjFunc"]
        else:
            print(f"[Warning] Unexpected column count ({n_cols}) in {cnv_path}, skipping.")
            continue

        # Read the file
        print("reading cnv file")
        try:
            df = pd.read_csv(cnv_path, delim_whitespace=True, comment="#", header=None, names=colnames, skiprows=1)
        except Exception as e:
            print(f"[Error] Failed to read {cnv_path}: {e}")
            continue

        # Extract alpha value from directory name
        match = re.search(r'inversion_alpha_([0-9.]+)', os.path.basename(folder))
        if not match:
            print(f"[Warning] Could not parse alpha from {folder}")
            continue
        alpha_val = float(match.group(1))

        for _, row in df.iterrows():
            all_data.append({
                "alpha_dir": os.path.basename(folder),
                "Distortion": dist_flag,
                "Iter": row["Iter"],
                "Alpha": alpha_val,
                "Roughness": row["Roughness"],
                "Misfit": row["Misfit"],
                "RMS": row["RMS"]
            })

# Convert to DataFrame
all_df = pd.DataFrame(all_data)

if all_df.empty:
    print("No data collected. Please check folder paths and file contents.")
else:
    # Keep only the final iteration for each (alpha_dir, Distortion)
    filtered_df = all_df.loc[
    all_df.groupby(["alpha_dir", "Distortion"])["RMS"].idxmin()
]
    filtered_df.to_csv("L_curve_data.csv", index=False)
    print("Saved CSV: L_curve_data.csv")
    # Plotting
    fig, ax = plt.subplots(figsize=(6, 4))
    for dist_flag, group in filtered_df.groupby("Distortion"):
        group["Alpha**2"] = group["Alpha"] ** 2 
        group_sorted = group.sort_values("Alpha**2")
    
        # ✅ Plot the actual points!
        # ax.plot(group_sorted["Roughness"], group_sorted["Misfit"],
        #         marker='o', linestyle='-', markersize=4, label=dist_flag)
        
        # First, draw the connecting line (no markers here)
        ax.plot(group_sorted["Roughness"], group_sorted["Misfit"],
                color='tab:blue' if dist_flag == "homogeneous starting model" else 'tab:orange',
                linestyle='-', linewidth=1)
        
        # Then overlay the styled markers
        ax.scatter(group_sorted["Roughness"], group_sorted["Misfit"],
           s=25,  # size of markers
           c='tab:blue' if dist_flag == "homogeneous starting model" else 'tab:orange',
           edgecolors='k', linewidths=0.5,
           label=dist_flag, zorder=3)
        
        
        # Annotate each point with alpha²
        # for _, row in group_sorted.iterrows():
        #     alpha_sq = row["Alpha"] ** 2
        #     if abs(alpha_sq - round(alpha_sq)) < 1e-6:
        #         alpha_str = f"{int(round(alpha_sq))}"
        #     else:
        #         alpha_str = f"{alpha_sq:.2f}"
    
        #     if dist_flag == "with distortion inversion":
        #         ax.annotate(alpha_str, (row["Roughness"], row["Misfit"]),
        #                     textcoords="offset points", xytext=(-3, -7),
        #                     ha='right', fontsize=8)
        #     else:
        #         ax.annotate(alpha_str, (row["Roughness"], row["Misfit"]),
        #                     textcoords="offset points", xytext=(3, 3),
        #                     ha='left', fontsize=8)
        
        # Annotate each point with alpha
        for _, row in group_sorted.iterrows():
            alpha = row["Alpha"]
            if abs(alpha - round(alpha)) < 1e-6:
                    alpha_str = f"{int(round(alpha))}"
            else:
                alpha_str = f"{alpha:.1f}"   
            if dist_flag == "with distortion inversion":
                ax.annotate(alpha_str, (row["Roughness"], row["Misfit"]),
                            textcoords="offset points", xytext=(-5, -9),
                            ha='right', fontsize=8)
            else:
                ax.annotate(alpha_str, (row["Roughness"], row["Misfit"]),
                            textcoords="offset points", xytext=(3, 3),
                            ha='left', fontsize=8)
        
    
    # ✅ Highlight alpha=1 case
    #highlight = filtered_df[((filtered_df["Alpha"] == 1.0) & (filtered_df["Distortion"] == "homogeneous starting model")) |
    #                        ((filtered_df["Alpha"] == 1.5) & (filtered_df["Distortion"] == "homogeneous starting model")) ]
    # if not highlight.empty:
    #     x = highlight["Roughness"].values[0]
    #     y = highlight["Misfit"].values[0]
    #     x2 = highlight["Roughness"].values[1]
    #     y2 = highlight["Misfit"].values[1]
        
    #     circle1 = mpatches.Circle((x, y), radius=500, edgecolor='black', facecolor='none',
    #                              linestyle=(0, (2, 1)), linewidth=1.2, zorder=10, alpha=1)
    #     ax.add_patch(circle1)
    #     circle2 = mpatches.Circle((x2, y2), radius=500, edgecolor='black', facecolor='none',
    #                              linestyle=(0, (2, 1)), linewidth=1.2, zorder=10, alpha=1)
    #     ax.add_patch(circle2)
        
    
    ax.set_xlabel("Model Roughness")
    ax.set_ylabel("Data Misfit")
    ax.set_title(r"Data Misfit vs Model Roughness ($\alpha$ annotated)")
    #ax.set_aspect('equal')
    ax.grid(True)
    ax.legend()
    ax.set_xlim(left=0, right=3200)
    ax.set_ylim(bottom=0, top=25000)
    plt.tight_layout()
    plt.savefig("L_curve_plot_synthetic.pdf", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()
