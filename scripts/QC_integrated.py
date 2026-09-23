# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 10:23:16 2025
This creates the files called sites_vtf.txt and sites_imp.txt based on 
files result_VTF.txt and result_MT.txt as output from applying mergeResultOfFEMTIC 
to femtic inversion results
@author: charroyj
"""

for var in ['all_data', 'all_df', 'df', 'filtered_df', 'results', 'inversion']:
    try:
        del globals()[var]
    except:
        pass  # Variable doesn't exist, skip it

import sys, os, glob, subprocess, pandas as pd
current_dir = os.path.dirname(os.path.abspath(__file__))
femticpy_path = "c:/users/charroyj/documents/mt/annecy/femtic/software/femticpy"
sys.path.append(femticpy_path+'/src')
import femticPy

cnv_path = current_dir+'/femtic.cnv'
all_data = []
iter = None #specify a value to use that instead of the one with lowest RMS

merge_exe_path="C:\\Users\\charroyj\\Documents\\MT\\Annecy\\femtic\\Software\\mergeResultOfFEMTIC\mergeResult.exe"

# 1. Set up a data object (same as when created for the inversion)
inversion = femticPy.DataGen(survey = 'inversion_02', outdir = './')

# Loading the data and the data coordinates
inversion.read_observe('./observe.dat')

# below not needed when reading from observe.dat, rather than source EDIs
#inversion.read_MTdata_coordinates('../../annecy_2025/input_data', 'coords_annecy')

## Data to be inverted for
inversion.invert_Z   = True
inversion.invert_VTF = True
inversion.invert_PT  = False

# load the 3 files needed to create the mesh
inversion.topography = '../../annecy_2025/input_data/topography.dat'
inversion.bathymetry = '../../annecy_2025/input_data/bathymetry.dat'
inversion.coast_line = '../../annecy_2025/input_data/coast_line.dat'

# we center the data to a anchor point (center of the data set) which will also be the center of the future mesh
# inversion.center_data()
# inversion.anchor


# 2. Analyse and plot inversion results

results_directory = './'

results = femticPy.InvResults(inversion.survey,
                              inversion.mt_coords, 
                              inversion.ids_Z,
                              inversion.ids_VTF,
                              results_directory)


# Read cnv file to extract iteration with lowest RMS number. 
# this should be refactored in a femticpy method reading femtic.cnv
# merge results at iteration with lowest RMS, unless otherwise specified just below.
#
if iter is None:
    with open(cnv_path, "r") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                n_cols = len(line.split())
                break
    # Assign appropriate column names
    if n_cols == 8:
        colnames = ["Iter", "Retrial", "Alpha", "Damp", "Roughness", "Misfit", "RMS", "ObjFunc"]
        dist_flag = 0
    elif n_cols == 10:
        colnames = ["Iter", "Retrial", "Alpha", "Beta", "Damp",
                    "Roughness", "Distortion", "Misfit", "RMS", "ObjFunc"]
        dist_flag = 1
    else:
        print(f"[Warning] Unexpected column count ({n_cols}) in {cnv_path}, skipping.")
    
    # Read the file
    try:
        df = pd.read_csv(cnv_path, delim_whitespace=True, comment="#", header=None, names=colnames, skiprows=1)
    except Exception as e:
        print(f"[Error] Failed to read {cnv_path}: {e}")
    
    for _, row in df.iterrows():
                all_data.append({
                    "alpha_dir": os.path.basename(current_dir),
                    "Distortion": dist_flag,
                    "Iter": row["Iter"],
                    "Alpha": row["Alpha"],
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
    iter=str(int(filtered_df['Iter']))

# Determine number of MPI processes used based on number of femtic_*.log files found in directory
n_mpi=str(len(glob.glob(os.path.join(current_dir, "femtic_*.log"),recursive=False)))

subprocess.run([merge_exe_path, iter, n_mpi, "-csv"])
subprocess.run([merge_exe_path, iter, n_mpi, "-appphs"])

# Set the number of the last iteration
results.nIter = int(iter)
results.read_result_csv()

# print total RMS 
results.compute_rms_breakDown()
print('RMS TOTAL: ',results.rms_total)
print('RMS Z: ',results.rms_Z)
print('RMS VTF: ',results.rms_VTF)


# print RMS per site / per component
results.rms_breakDown_Z.style.background_gradient()

if inversion.invert_VTF:
    results.rms_breakDown_VTF.style.background_gradient()

# plot RMS map

results.plot_rms_map(
    Z=True,
    VTF=True,
    PT=False,  # Set to True if you want PT too
    xlim=[-5, 5],
    ylim=[-7, 7],
    vmin=0.5,       #RMS min and max values
    vmax=3,
    ms_size=40,     #marker size for plot
    save_plot=True,
    inversion=inversion,  # Required for topo and coastline background
    topo_vmin=0,
    topo_vmax=1.2
)

# plot Induction arrows
frequencies = results.resp_VTF['Freq[Hz]'][[9,10,15]].values
results.plot_induction_arrows(frequencies,
                              real = True,
                              imag = True,
                              inv_response = True,
                              scale = 30,
                              xlim = [-5, 5],
                              ylim = [-7, 7], 
                              save_plot = True)

# results.plot_Z_fit(plot_Z = False,
#                  xlim = [-3,2],
#                  ylim = [1,3],
#                  ylim_diag = [-3,2],
#                  save_plot = True,
#                  add_map_stats = True)

results.plot_Z_VTF_fit(plot_Z=False,
                    xlim_Z=[-3, 2],
                    ylim_Z=[1, 3],
                    ylim_diag=[-3, 2],
                    xlim_VTF=[-3, 2],
                    ylim_VTF=[-0.5, 0.5],
                    save_plot=True,
                    auto_ylim=True)


if inversion.invert_VTF:
    results.plot_VTF_fit(xlim = [-3,2],
                    ylim = [-0.5,0.5],
                      save_plot = True,
                      add_map_stats = True)


if dist_flag:
   results.plot_distorsion_map(xlim = [-5,5],
                            ylim = [-7,7], 
                            #vmin=0.5,
                            #vmax=3,
                            save_plot = True, 
                            inversion=inversion,
                            topo_vmin=0,
                            topo_vmax=1.2)

results.plot_cnv(save_plot=True)
