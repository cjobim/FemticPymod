import os
import numpy as np
import shutil
import json

input_file = "resistivity_block.dat"
dest_file = "resistivity_block_median.dat"
dest_file2 = "resistivity_block_dec1.dat"
dest_file4 = "resistivity_block_dec9.dat"
dest_file3 = "resistivity_block_spread.dat"
rms_file = "rms_random_start_data.dat"
json_file = "cell_resistivities_log.json"

#these 3 values determine whether to impose a threshold on RMS or convergence to accept/reject results
rms_overall_min = 1.009367
rms_threshold = 0 #1.1*rms_overall_min           #set to zero to not use it
filter_conv = 0     #filter out models that didn't converge in constant_start case

header_line = "\n        0        1.000000e+09   1.000000e-20   1.000000e+20   1.000000e+00         1\n"
count_output_file = "resistivity_count.dat"
count_header_line = "Cell Index    Count of Elements\n"

current_dir = os.path.dirname(os.path.abspath(__file__))
# List directories that start with "random_start_"
dirs = sorted([d for d in os.listdir(".") if os.path.isdir(d) and d.startswith("inversion_")])

cell_resistivities_log = {}
nb_elem = None



shutil.copy(os.path.join(current_dir,input_file),os.path.join(current_dir,dest_file))
shutil.copy(os.path.join(current_dir,input_file),os.path.join(current_dir,dest_file2))
shutil.copy(os.path.join(current_dir,input_file),os.path.join(current_dir,dest_file3))
shutil.copy(os.path.join(current_dir,input_file),os.path.join(current_dir,dest_file4))

for directory in dirs:
    cnv_path = os.path.join(directory, "femtic.cnv")
    if not os.path.exists(cnv_path):
        print(f"Missing: {cnv_path}")
        continue
    if filter_conv == 1 and directory.endswith("_nc"):
        continue        
    
    try:
        with open(cnv_path) as f:
            lines = f.readlines()[1:]
            rms_values = [(float(line.split()[-2]), int(line.split()[0])) for line in lines if len(line.split()) >= 9]

        if not rms_values:
            print(f"No RMS values in {directory}")
            continue

        # Find the minimum RMS value and the corresponding iteration
        min_rms_value, best_iter = min(rms_values, key=lambda x: x[0])
        block_path = os.path.join(directory, f"resistivity_block_iter{best_iter}.dat")
        if min_rms_value > rms_threshold:
            print(f"RMS higher than threshold: {rms_threshold}")
            continue
            
        if not os.path.exists(block_path):
            print(f"Missing: {block_path}")
            continue
        
        with open(rms_file, "a") as rms_file_out:
            rms_file_out.write(f"{directory} {min_rms_value} {best_iter}\n")

        with open(block_path) as bf:
            block_lines = bf.readlines()

        if nb_elem is None:
            try:
                nb_elem = int(block_lines[0].split()[0])
                nb_cells = int(block_lines[0].split()[1])
                print(f"nb_elem set to {nb_elem}")
            except (ValueError, IndexError) as e:
                print(f"Couldn't parse nb_elem in {directory}: {e}")
                continue

        for i, line in enumerate(block_lines[nb_elem + 2:], 1):
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    resistivity_log = np.log10(float(parts[1]))
                    cell_resistivities_log.setdefault(i, []).append(resistivity_log)
                except ValueError:
                    print(f"Couldn't convert resistivity value in {directory}: {line.strip()}")

    except Exception as e:
        print(f"Error processing {directory}: {e}")

if nb_elem is None:
    print("No valid nb_elem found. Aborting.")
else:
    with open(dest_file, "a") as out, open(dest_file2, "a") as out2, open(dest_file3, "a") as out3, open(dest_file4, "a") as out4:
        out.write(header_line)
        out2.write(header_line)
        out3.write(header_line)
        out4.write(header_line)
        for i in range(1, nb_cells + 1):
            values = cell_resistivities_log.get(i, [])
            if values:
                dec1, median, dec9 = 10**(np.percentile(values,[10, 50, 90]))
                spread = dec9/dec1
                out.write(f"{i:9d}        {median:.6e}   1.000000e-20   1.000000e+20   1.000000e+00         0\n")
                out2.write(f"{i:9d}        {dec1:.6e}   1.000000e-20   1.000000e+20   1.000000e+00         0\n")
                out3.write(f"{i:9d}        {spread:.6e}   1.000000e-20   1.000000e+20   1.000000e+00         0\n")
                out4.write(f"{i:9d}        {dec9:.6e}   1.000000e-20   1.000000e+20   1.000000e+00         0\n")

    with open(count_output_file, "w") as count_out:
        count_out.write(count_header_line)
        for i in range(1, nb_cells + 1):
            values = cell_resistivities_log.get(i, [])
            count_out.write(f"{i:9d}        {len(values)}\n")
            
    with open(json_file, "w") as json_out:
        json.dump(cell_resistivities_log, json_out, indent=4)