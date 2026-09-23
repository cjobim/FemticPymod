# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 11:30:01 2025
reads EQ epicenters to plot on MT sections using GMT
@author: charroyj
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

inversion_anchor = [5091.770697757265, 272.7923403189063]
file_eq = "C:/Users/charroyj/Documents/GIS/data/EQ_aftershocks_UTM.csv"
file_nodes = "./profile_nodes.csv"
proj_distance = 1.5  # in km

df_eq = pd.read_csv(file_eq)  # all in model space, in km, ie x=northing in general

# --- Create subplots ---
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

sky_blue = "#5DADE2"    # nice medium sky blue
fresh_green = "#58D68D" # nice green
median_red = "#E74C3C"  # soft red for median line

bin_width = 0.15  # desired bin width

# 1️⃣ Histogram of "erreur épicentre"
bins = np.arange(0, df_eq['erreur épicentre'].max() + bin_width, bin_width)
axs[0].hist(df_eq['erreur épicentre'].dropna(), bins=bins, color=sky_blue, edgecolor='black')

# Compute median and plot line
median_epi = df_eq['erreur épicentre'].median()
axs[0].axvline(median_epi, color=median_red, linestyle='--', linewidth=2, label=f'Median = {median_epi:.2f} km')

axs[0].set_xlabel("Horizontal error (km)", fontsize=12)
axs[0].set_ylabel("Frequency", fontsize=12)
axs[0].grid(alpha=0.3, linestyle=':')
axs[0].set_xlim(left=0, right=3)
axs[0].legend()

# 2️⃣ Histogram of "erreur profondeur"
bins = np.arange(0, df_eq["erreur profondeur"].max() + bin_width, bin_width)
axs[1].hist(df_eq['erreur profondeur'].dropna(), bins=bins, color=fresh_green, edgecolor='black')

# Compute median and plot line
median_depth = df_eq['erreur profondeur'].median()
axs[1].axvline(median_depth, color=median_red, linestyle='--', linewidth=2, label=f'Median = {median_depth:.2f} km')

axs[1].set_xlabel("Vertical error (km)", fontsize=12)
axs[1].set_ylabel("Frequency", fontsize=12)
axs[1].grid(alpha=0.3, linestyle=':')
axs[1].set_xlim(left=0, right=3)
axs[1].legend()

# --- Adjust layout and show ---
fig.tight_layout()
plt.savefig('Histogram_EQ_errors.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.show()