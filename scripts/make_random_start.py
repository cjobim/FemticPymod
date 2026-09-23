# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 17:13:38 2025
This script randomizes starting resistivities of a model defined in .dat femtic format uniformly between two bounds (here 10^0 and 10^4 ohm.m). The number of resistivity blocks needs to be modified to suit, and resistivity_block.dat is a model file made by cutting out the resistivity block information (2nd part of the file) from a regular resistivity_block_iterXX.dat femtic output file. 
@author: charroyj
"""

import random
start_file = "resistivity_block.dat"
output_file = "resistivity_block_iter0.dat"

with open(start_file, "r") as orig:
    start = orig.read()

new_lines = []

line = "\n         0        1.000000e+09   1.000000e-20   1.000000e+20   1.000000e+00         1\n"
new_lines.append(line)
for i in range(1, 54588):  # 54587 inclusive
    x = 10**(random.uniform(0, 4))
    line = f" {i:9d}        {x:.6e}   1.000000e-20   1.000000e+20   1.000000e+00         0\n"
    new_lines.append(line)

with open(output_file, "w") as f:
    f.write(start)
    f.writelines(new_lines)

print(f"File '{output_file}' successfully written.")