# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 15:46:51 2025
calculate coordinates of MT sites in modelspace
@author: charroyj
"""

import csv, pandas as pd

inversion_anchor = [272.7923403189063, 5091.770697757265]

file_coords="./coords_annecy_2025.dat"
file_coords_conv="./coords_annecy_2025_modelspace.dat"
df_coords = pd.read_csv(file_coords, delim_whitespace=True, names=['name', 'y', 'x', 'z']) #all in model space, in km, ie x=northing in general

with open(file_coords_conv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in df_coords.itertuples(index=False):
        #calculation of modelspace coordinates:
        y_ms = row.y - inversion_anchor[0]
        x_ms = row.x - inversion_anchor[1]
        z_ms = -row.z
        writer.writerow([row.name,y_ms, x_ms,z_ms])