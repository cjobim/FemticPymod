# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 16:31:40 2025

@author: charroyj
"""

# -*- coding: utf-8 -*-
"""
Generate noisy MT datasets using FEMTICPy
Efficient: input data read once, reused each loop.
"""

import sys
import os
import copy
import shutil
import random
current_dir = os.path.dirname(os.path.abspath(__file__))
femticpy_path = os.path.abspath(os.path.join(current_dir, '..', '..'))
sys.path.append(os.path.join(femticpy_path, 'src'))
import femticPy


# ============================================================
# 1. Base inversion setup
# ============================================================

inversion = femticPy.DataGen(survey='synthetic_femtic_random_start_data', outdir='./inversion_template')

# Load real MT data
inversion.read_MTdata('./input_data/edi_files')
inversion.read_MTdata_coordinates('./input_data', 'coords_synth')

# Flags
inversion.invert_Z   = True
inversion.invert_VTF = True
inversion.invert_PT  = False

# Mesh-related files
inversion.topography = './input_data/topography.dat'
inversion.bathymetry = './input_data/bathymetry.dat'
inversion.coast_line = './input_data/coast_line.dat'

# Center and define analysis domain
inversion.center_data()
inversion.analysis_domain = [[-20, 20],
                                [-20, 20],
                                [-50, 70]] 

# Optional plots
# inversion.plot_data_loc(plot_ids=False, zoom_core=False)
# inversion.plot_coast_line()

# ============================================================
# 2. Apply error floors ONCE to create the baseline dataset
# ============================================================

inversion.error_floor_Z  = [0.15, 0.03, 0.03, 0.15]
inversion.error_floor_Tz = 0.03

# # Apply defined errors (produces proper VAR columns)
# inversion.data_Z = inversion.apply_defined_error(
#     inversion.data_Z, data_type='Z',
#     comps=['ZXX.VAR','ZYY.VAR'],
#     err_val=[0.15,0.15],
#     error_type=2
# )

# inversion.data_VTF = inversion.apply_defined_error(
#     inversion.data_VTF, data_type='VTF',
#     comps=['TXVAR.EXP','TYVAR.EXP'],
#     err_val=[0.03,0.03]
# )

# Filter unreasonable values once
inversion.filter_Z(max_val=1e5)
inversion.filter_VTF(max_val=1e5)

#============================================================
# Save CLEAN (noise-free) processed data for reuse
#============================================================

template_Z   = copy.deepcopy(inversion.data_Z)
template_VTF = copy.deepcopy(inversion.data_VTF)

#============================================================
# Frequency band for writing observe.dat
#============================================================
highfreq = 1000
lowfreq  = 0.035
subsampling = 1


# ============================================================
# 3. Loop to generate noisy datasets
# ============================================================

N = 120

for k in range(1, N+1):

    outdir = f'./inversion_{k:03d}'
    print(f"\n=== Writing noisy dataset {k}/{N} → {outdir} ===")

    os.makedirs(outdir, exist_ok=True)
    
    # ------------------------------------------------------------------
    # Create a fresh copy of the inversion object; keep geometry/settings
    # ------------------------------------------------------------------
    inv = copy.deepcopy(inversion)
    inv.outdir = outdir

    # Reset data to original clean (pre-error-floor + pre-filter) state
    inv.data_Z   = copy.deepcopy(template_Z)
    inv.data_VTF = copy.deepcopy(template_VTF)

    # ------------------------------------------------------------------
    # Apply Gaussian noise for this iteration only
    # ------------------------------------------------------------------
    inv.data_Z   = inv.add_gaussian_noise(inv.data_Z,   data_type='Z',   seed=None)
    inv.data_VTF = inv.add_gaussian_noise(inv.data_VTF, data_type='VTF', seed=None)

    # ------------------------------------------------------------------
    # Write observe.dat for this noisy realization
    # ------------------------------------------------------------------
    inv.write_observe(
        write=True,
        freq_bandwidth=[lowfreq, highfreq],
        subsampling=subsampling)
    
    #-------------------------------------------------------------------
    #Make randomized starting resistivity 
    #-------------------------------------------------------------------

    start_file = "resistivity_block.dat"
    output_file = f"{outdir}/resistivity_block_iter0.dat"
    
    with open(start_file, "r") as orig:
        start = orig.read()
    
    new_lines = []
    
    line = "\n         0        1.000000e+09   1.000000e-20   1.000000e+20   1.000000e+00         1\n"
    new_lines.append(line)
    for i in range(1, 44630):  # 54587 inclusive
        x = 10**(random.uniform(0, 4))
        line = f" {i:9d}        {x:.6e}   1.000000e-20   1.000000e+20   1.000000e+00         0\n"
        new_lines.append(line)
    
    with open(output_file, "w") as f:
        f.write(start)
        f.writelines(new_lines)
    
    shutil.copy('inversion/control.dat',f"{outdir}/control.dat")
print(f"Files in '{outdir}' successfully written.")
    

print(f"\nAll {N} noisy synthetic datasets generated successfully.\n")