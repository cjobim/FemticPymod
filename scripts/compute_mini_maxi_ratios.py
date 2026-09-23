import os
import re
import json
import shutil
import numpy as np


# ======================================================================
# USER SETTINGS
# ======================================================================

input_file = "resistivity_block.dat"
ref_file = "median_1.5-m4/resistivity_block_iter6.dat"

dest_file = "resistivity_block_min.dat"
dest_file2 = "resistivity_block_max.dat"
dest_file3 = "resistivity_block_spreadbelow.dat"
dest_file4 = "resistivity_block_spreadabove.dat"

rms_file = "rms_random_start_data.dat"
ensemble_min_rms_file = "minimum_rms_ensemble.dat"

json_file = "cell_resistivities_log.json"
count_output_file = "resistivity_count.dat"


# RMS acceptance criterion:
#
# accepted if:
#
#     RMS <= rms_threshold_factor * minimum_RMS_of_ensemble
#
# Set to 0 or None to disable RMS filtering.
rms_threshold_factor = 1.10

# Filter out directories ending in "_nc"
filter_conv = 0


# Parameter 0 = fixed air resistivity
header_line = (
    "\n        0        1.000000e+09   1.000000e-20   "
    "1.000000e+20   1.000000e+00         1\n"
)

count_header_line = "Cell Index    Count of Elements\n"


# ======================================================================
# INITIALISATION
# ======================================================================

current_dir = os.path.dirname(os.path.abspath(__file__))

dirs = sorted([
    d for d in os.listdir(".")
    if os.path.isdir(d) and d.startswith("inversion_")
])

ensemble_rms = []
cell_resistivities_log = {}

nb_elem = None
nb_params = None
nb_cells = None


# ======================================================================
# HELPER FUNCTION
# ======================================================================

def get_member_number(directory):
    """
    Extract ensemble member number from names such as:
        inversion_001
        inversion_001_nc
    """
    match = re.match(r"inversion_(\d+)", directory)

    if match:
        return int(match.group(1))

    return -1


# ======================================================================
# PASS 1
#
# Find the best RMS of every ensemble member
# ======================================================================

print("\n============================================================")
print("PASS 1: Determining RMS values for ensemble")
print("============================================================\n")


for directory in dirs:

    if filter_conv == 1 and directory.endswith("_nc"):
        print(f"Skipping non-converged model: {directory}")
        continue

    cnv_path = os.path.join(directory, "femtic.cnv")

    if not os.path.exists(cnv_path):
        print(f"Missing: {cnv_path}")
        continue

    try:

        with open(cnv_path) as f:
            lines = f.readlines()[1:]

        rms_values = []

        for line in lines:

            parts = line.split()

            if len(parts) >= 9:

                try:
                    iteration = int(parts[0])
                    rms = float(parts[-2])

                    rms_values.append(
                        (rms, iteration)
                    )

                except ValueError:
                    pass

        if not rms_values:
            print(f"No RMS values found in {directory}")
            continue

        # Best RMS of this ensemble member
        min_rms_value, best_iter = min(
            rms_values,
            key=lambda x: x[0]
        )

        block_path = os.path.join(
            directory,
            f"resistivity_block_iter{best_iter}.dat"
        )

        # Only consider models whose resistivity file exists
        if not os.path.exists(block_path):

            print(
                f"{directory}: best RMS = {min_rms_value:.6f}, "
                f"iteration {best_iter}, but model file is missing"
            )

            continue

        member_number = get_member_number(directory)

        ensemble_rms.append({
            "directory": directory,
            "member": member_number,
            "rms": min_rms_value,
            "iteration": best_iter,
            "block_path": block_path,
        })

        print(
            f"{directory:<20s} "
            f"RMS = {min_rms_value:.6f} "
            f"iteration = {best_iter}"
        )

    except Exception as e:

        print(
            f"Error processing {directory}: {e}"
        )


# ======================================================================
# DETERMINE GLOBAL MINIMUM RMS
# ======================================================================

if not ensemble_rms:
    raise RuntimeError(
        "No valid ensemble RMS values were found."
    )


best_member = min(
    ensemble_rms,
    key=lambda x: x["rms"]
)

rms_overall_min = best_member["rms"]


if rms_threshold_factor is None or rms_threshold_factor <= 0:

    rms_threshold = np.inf

else:

    rms_threshold = (
        rms_threshold_factor
        * rms_overall_min
    )


print("\n============================================================")
print("ENSEMBLE RMS SUMMARY")
print("============================================================")

print(
    f"Minimum RMS       : {rms_overall_min:.8f}"
)

print(
    f"Directory         : {best_member['directory']}"
)

print(
    f"Member number     : {best_member['member']}"
)

print(
    f"Iteration         : {best_member['iteration']}"
)

if np.isfinite(rms_threshold):

    print(
        f"Threshold factor  : {rms_threshold_factor}"
    )

    print(
        f"RMS threshold     : {rms_threshold:.8f}"
    )

else:

    print("RMS filtering disabled")


# ======================================================================
# WRITE MINIMUM RMS INFORMATION
# ======================================================================

with open(ensemble_min_rms_file, "w") as f:

    f.write(
        "Directory          Member    RMS             "
        "Iteration    Threshold_factor    RMS_threshold\n"
    )

    f.write(
        f"{best_member['directory']:<18s} "
        f"{best_member['member']:6d}    "
        f"{rms_overall_min:.8e}    "
        f"{best_member['iteration']:9d}    "
    )

    if np.isfinite(rms_threshold):

        f.write(
            f"{rms_threshold_factor:.6f}          "
            f"{rms_threshold:.8e}\n"
        )

    else:

        f.write(
            "disabled          disabled\n"
        )


# ======================================================================
# DETERMINE ACCEPTED ENSEMBLE MEMBERS
# ======================================================================

accepted_models = [
    model
    for model in ensemble_rms
    if model["rms"] <= rms_threshold
]


print("\n============================================================")
print("MODEL SELECTION")
print("============================================================")

print(
    f"Valid ensemble members : {len(ensemble_rms)}"
)

print(
    f"Accepted members       : {len(accepted_models)}"
)


# ======================================================================
# WRITE ACCEPTED RMS FILE
# ======================================================================

with open(rms_file, "w") as rms_out:

    rms_out.write(
        "Directory Member RMS Iteration\n"
    )

    for model in accepted_models:

        rms_out.write(
            f"{model['directory']} "
            f"{model['member']} "
            f"{model['rms']:.8e} "
            f"{model['iteration']}\n"
        )


# ======================================================================
# PREPARE OUTPUT MODEL FILES
# ======================================================================

for output_file in [
    dest_file,
    dest_file2,
    dest_file3,
    dest_file4,
]:

    shutil.copy(
        os.path.join(current_dir, input_file),
        os.path.join(current_dir, output_file)
    )


# ======================================================================
# PASS 2
#
# Read resistivities only from accepted models
# ======================================================================

print("\n============================================================")
print("PASS 2: Reading accepted resistivity models")
print("============================================================\n")


for model in accepted_models:

    directory = model["directory"]
    block_path = model["block_path"]

    print(
        f"Reading {block_path} "
        f"(RMS = {model['rms']:.6f})"
    )

    try:

        with open(block_path) as bf:
            block_lines = bf.readlines()


        # --------------------------------------------------------------
        # Read model dimensions
        #
        # FEMTIC header contains:
        #
        #     number of elements
        #     total number of resistivity parameters
        #
        # Parameter 0 is air, therefore:
        #
        #     number of earth parameters = nb_params - 1
        # --------------------------------------------------------------

        model_nb_elem = int(
            block_lines[0].split()[0]
        )

        model_nb_params = int(
            block_lines[0].split()[1]
        )

        model_nb_cells = (
            model_nb_params - 1
        )


        # First accepted model defines expected structure
        if nb_elem is None:

            nb_elem = model_nb_elem
            nb_params = model_nb_params
            nb_cells = model_nb_cells

            print(
                f"nb_elem       = {nb_elem}"
            )

            print(
                f"nb_parameters = {nb_params} "
                f"(including air parameter 0)"
            )

            print(
                f"nb_cells      = {nb_cells} "
                f"(earth parameters only)"
            )

        else:

            if model_nb_elem != nb_elem:

                raise ValueError(
                    f"{directory}: number of elements differs "
                    f"({model_nb_elem} versus {nb_elem})"
                )

            if model_nb_params != nb_params:

                raise ValueError(
                    f"{directory}: number of parameters differs "
                    f"({model_nb_params} versus {nb_params})"
                )


        # --------------------------------------------------------------
        # Read earth-model resistivities
        #
        # Layout:
        #
        # line 0                   header
        # next nb_elem lines       element -> parameter mapping
        # next line                parameter 0 (air)
        # remaining lines          parameters 1 ... nb_params-1
        # --------------------------------------------------------------

        model_cell_lines = block_lines[
            model_nb_elem + 2:
            model_nb_elem + 2 + model_nb_cells
        ]


        if len(model_cell_lines) != model_nb_cells:

            raise ValueError(
                f"{directory}: expected {model_nb_cells} "
                f"earth parameter lines, "
                f"found {len(model_cell_lines)}"
            )


        for i, line in enumerate(
            model_cell_lines,
            start=1
        ):

            parts = line.strip().split()

            if len(parts) < 2:

                raise ValueError(
                    f"Invalid parameter line {i} "
                    f"in {directory}"
                )

            parameter_id = int(parts[0])
            resistivity = float(parts[1])


            # Optional consistency check
            if parameter_id != i:

                raise ValueError(
                    f"{directory}: expected parameter {i}, "
                    f"found parameter {parameter_id}"
                )


            if resistivity <= 0:

                raise ValueError(
                    f"Non-positive resistivity "
                    f"in {directory}, parameter {i}: "
                    f"{resistivity}"
                )


            resistivity_log = float(
                np.log10(resistivity)
            )


            cell_resistivities_log.setdefault(
                i, []
            ).append(
                resistivity_log
            )


    except Exception as e:

        print(
            f"Error reading {directory}: {e}"
        )


# ======================================================================
# READ REFERENCE MODEL
# ======================================================================

print("\n============================================================")
print("READING REFERENCE MODEL")
print("============================================================\n")


if nb_elem is None or nb_params is None or nb_cells is None:

    raise RuntimeError(
        "No accepted resistivity models were read."
    )


if not os.path.exists(ref_file):

    raise FileNotFoundError(
        f"Reference model not found: {ref_file}"
    )


with open(ref_file) as f:
    ref_lines = f.readlines()


ref_nb_elem = int(
    ref_lines[0].split()[0]
)

ref_nb_params = int(
    ref_lines[0].split()[1]
)

ref_nb_cells = (
    ref_nb_params - 1
)


# Verify identical FEMTIC structure
if ref_nb_elem != nb_elem:

    raise ValueError(
        f"Reference model has {ref_nb_elem} elements, "
        f"ensemble models have {nb_elem}"
    )


if ref_nb_params != nb_params:

    raise ValueError(
        f"Reference model has {ref_nb_params} parameters, "
        f"ensemble models have {nb_params}"
    )


# Skip mesh mapping and air parameter
ref_cell_lines = ref_lines[
    ref_nb_elem + 2:
    ref_nb_elem + 2 + ref_nb_cells
]


if len(ref_cell_lines) != ref_nb_cells:

    raise ValueError(
        f"Reference model: expected {ref_nb_cells} "
        f"earth parameter lines, "
        f"found {len(ref_cell_lines)}"
    )


ref_resistivities = {}


for i, line in enumerate(
    ref_cell_lines,
    start=1
):

    parts = line.strip().split()

    if len(parts) < 2:

        raise ValueError(
            f"Invalid reference parameter line {i}"
        )


    parameter_id = int(parts[0])
    ref = float(parts[1])


    if parameter_id != i:

        raise ValueError(
            f"Reference model: expected parameter {i}, "
            f"found parameter {parameter_id}"
        )


    if ref <= 0:

        raise ValueError(
            f"Invalid reference resistivity "
            f"for parameter {i}: {ref}"
        )


    ref_resistivities[i] = ref


print(
    f"Reference model read successfully:"
)

print(
    f"  elements              : {ref_nb_elem}"
)

print(
    f"  total parameters      : {ref_nb_params}"
)

print(
    f"  earth parameters      : {ref_nb_cells}"
)


# ======================================================================
# CHECK NUMBER OF SUCCESSFULLY READ MODELS
# ======================================================================

if not cell_resistivities_log:

    raise RuntimeError(
        "No resistivities were successfully read "
        "from the accepted ensemble."
    )


n_models_read = len(
    cell_resistivities_log.get(1, [])
)

print(
    f"\nSuccessfully read {n_models_read} "
    f"accepted resistivity models."
)


# ======================================================================
# CALCULATE MIN, MAX AND SPREAD
# ======================================================================

print("\n============================================================")
print("CALCULATING ENSEMBLE SPREAD")
print("============================================================\n")


with (
    open(dest_file, "a") as out,
    open(dest_file2, "a") as out2,
    open(dest_file3, "a") as out3,
    open(dest_file4, "a") as out4
):

    # Add parameter 0 = air
    out.write(header_line)
    out2.write(header_line)
    out3.write(header_line)
    out4.write(header_line)


    for i in range(
        1,
        nb_cells + 1
    ):

        values = cell_resistivities_log.get(
            i, []
        )


        if not values:

            raise ValueError(
                f"No ensemble resistivities "
                f"available for parameter {i}"
            )


        # --------------------------------------------------------------
        # Minimum and maximum
        #
        # values are stored in log10 resistivity
        # --------------------------------------------------------------

        mini_log, maxi_log = np.percentile(
            values,
            [0, 100]
        )

        mini = 10.0 ** mini_log
        maxi = 10.0 ** maxi_log


        # --------------------------------------------------------------
        # Reference resistivity
        # --------------------------------------------------------------

        ref = ref_resistivities[i]


        # --------------------------------------------------------------
        # Spread relative to reference model
        # --------------------------------------------------------------

        spreadbelow = (
            ref / mini
        )

        spreadabove = (
            maxi / ref
        )


        # --------------------------------------------------------------
        # Write output models
        # --------------------------------------------------------------

        out.write(
            f"{i:9d}        "
            f"{mini:.6e}   "
            f"1.000000e-20   "
            f"1.000000e+20   "
            f"1.000000e+00         0\n"
        )

        out2.write(
            f"{i:9d}        "
            f"{maxi:.6e}   "
            f"1.000000e-20   "
            f"1.000000e+20   "
            f"1.000000e+00         0\n"
        )

        out3.write(
            f"{i:9d}        "
            f"{spreadbelow:.6e}   "
            f"1.000000e-20   "
            f"1.000000e+20   "
            f"1.000000e+00         0\n"
        )

        out4.write(
            f"{i:9d}        "
            f"{spreadabove:.6e}   "
            f"1.000000e-20   "
            f"1.000000e+20   "
            f"1.000000e+00         0\n"
        )


# ======================================================================
# WRITE NUMBER OF ENSEMBLE MEMBERS CONTRIBUTING TO EACH PARAMETER
# ======================================================================

with open(
    count_output_file,
    "w"
) as count_out:

    count_out.write(
        count_header_line
    )

    for i in range(
        1,
        nb_cells + 1
    ):

        values = cell_resistivities_log.get(
            i, []
        )

        count_out.write(
            f"{i:9d}        "
            f"{len(values)}\n"
        )


# ======================================================================
# SAVE INDIVIDUAL LOG10 RESISTIVITIES
# ======================================================================

with open(
    json_file,
    "w"
) as json_out:

    json.dump(
        cell_resistivities_log,
        json_out,
        indent=4
    )


# ======================================================================
# FINAL SUMMARY
# ======================================================================

print("\n============================================================")
print("DONE")
print("============================================================")

print(
    f"Minimum ensemble RMS : "
    f"{rms_overall_min:.8f}"
)

print(
    f"Minimum RMS member   : "
    f"{best_member['directory']}"
)

print(
    f"Minimum RMS iteration: "
    f"{best_member['iteration']}"
)

if np.isfinite(rms_threshold):

    print(
        f"RMS threshold        : "
        f"{rms_threshold:.8f}"
    )

print(
    f"Models selected      : "
    f"{len(accepted_models)} / "
    f"{len(ensemble_rms)}"
)

print(
    f"Models successfully "
    f"read                : "
    f"{n_models_read}"
)

print(
    f"Total FEMTIC params  : "
    f"{nb_params}"
)

print(
    f"Earth parameters     : "
    f"{nb_cells}"
)

print(
    f"Minimum RMS info     : "
    f"{ensemble_min_rms_file}"
)
