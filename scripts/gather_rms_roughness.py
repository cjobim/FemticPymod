# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 13:55:41 2025

@author: charroyj
"""

import os
import glob
import csv

output_file = "RMS_roughness_sensitivitytest_Poisy.csv"
with open(output_file, "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["dirname", "logres", "res", "rms", "roughness"])
    for dir_name in sorted(glob.glob("mod*")):
        path = os.path.join(dir_name, "femtic.cnv")
        if os.path.isdir(dir_name) and os.path.isfile(path):
            try:
                with open(path) as f:
                    line = f.readline()  # skip header
                    line = f.readline()
                    parts = line.split()
                    logres = float(dir_name.replace("mod", ""))
                    res = 10**logres
                    rms=float(parts[8])
                    roughness=float(parts[5])
                    writer.writerow([dir_name, logres, res, rms, roughness])
                    
            except Exception as e:
                print(f"Error reading {path}: {e}")

