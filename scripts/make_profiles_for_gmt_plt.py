# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 14:09:18 2025

@author: charroyj
"""
import csv

#center of Annecy inversion with 48 stations
inversion_anchor = [5091.770697757265, 272.7923403189063]

profiles_in = "profile_nodes.csv"

with open(profiles_in, "r") as infile:
    reader = csv.reader(infile, delimiter=',')
    next(reader)        #skip header
    for i in range(4):  #number of profiles
        row1=next(reader)
        row2=next(reader)
        p_name = row1[0]
        p1_x = round(float(row1[3]) / 1000,3)
        p1_y = round(float(row1[4]) / 1000,3)
        p2_x = round(float(row2[3]) / 1000,3)
        p2_y = round(float(row2[4]) / 1000,3)
        dist = round(float(row2[1]),0)
        p_angle = round(float(row2[2]),2)
        
        param_angle = -p_angle
        pivot_x = (p1_x+p2_x)/2 - inversion_anchor[1]
        pivot_y = (p1_y+p2_y)/2 - inversion_anchor[0]
        filename_out = 'param_V_' + p_name + '.dat'
        
        lines = f"0\n13\n0\n{pivot_y} {pivot_x} 0\n{param_angle}\n1\n0"
        with open(filename_out,"w") as file_out:
            file_out.write(lines)
        
        
        