# -*- coding: utf-8 -*-
"""
Created on Fri May  9 13:50:17 2025
This calculates the projection of stations from a loc.txt file with
x y z name
for all stations in model coordinates (in km)
the goal is to have for each profile in a profile_nodes.csv file with 
name,distance,angle,x,y in utm coordinates, the list of stations within a certain 
projection distance, their distance along profile and z for plotting in GMT
@author: charroyj
"""

import pandas as pd
import csv
import numpy as np

inversion_anchor = [5091.770697757265, 272.7923403189063]
file_loc="./loc.txt"
file_nodes="./profile_nodes.csv"
file_eq = "C:/Users/charroyj/Documents/GIS/data/EQ_aftershocks_UTM.csv"
proj_distance=1.5   #in km
nbprofiles=5    #number of profiles

df_eq = pd.read_csv(file_eq)  
df_eq['magnitude'] = df_eq['magnitude'].fillna(df_eq['magnitude'].min())

def distance_to_profile(x, y, x1, y1, x2, y2):
    numerator = abs((x2 - x1)*(y1 - y) - (x1 - x)*(y2 - y1))
    denominator = np.hypot(x2 - x1, y2 - y1)
    return numerator / denominator

def distance_along_profile(x, y, x1, y1, x2, y2):
    # Profile direction vector
    dx = x2 - x1
    dy = y2 - y1
    profile_length = np.hypot(dx, dy)
    # Unit vector in profile direction
    ux, uy = dx / profile_length, dy / profile_length
    # Midpoint of profile
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    # Vector from midpoint to point
    vx, vy = x - mx, y - my
    # Project onto profile unit vector to get signed distance
    distance = vx * ux + vy * uy
    return distance

with open(file_nodes, "r") as infile_nodes:
    reader = csv.reader(infile_nodes, delimiter=',')
    next(reader)
    for i in range(nbprofiles):
        row1=next(reader)
        row2=next(reader)
        p_name = row1[0]
        p1_east = round(float(row1[3]) / 1000,3)
        p1_north = round(float(row1[4]) / 1000,3)
        p2_east = round(float(row2[3]) / 1000,3)
        p2_north = round(float(row2[4]) / 1000,3)
        
        p1_x = p1_north - inversion_anchor[0]
        p1_y = p1_east - inversion_anchor[1]
        p2_x = p2_north - inversion_anchor[0]
        p2_y = p2_east - inversion_anchor[1]
        
        profile_EQ=[]
        filename_profile_EQ = f"{p_name}_EQ.csv"
        
        for row in df_eq.itertuples(index=False):
            #calculate model space coordinates of EQ
            eq_x_ms=row.Y_UTM/1000 - inversion_anchor[0]
            eq_y_ms = row.X_UTM/1000 - inversion_anchor[1]
            dist = distance_to_profile(eq_x_ms,eq_y_ms,p1_x,p1_y,p2_x,p2_y)
            if dist < proj_distance:
                x_profile = distance_along_profile(eq_x_ms,eq_y_ms,p1_x,p1_y,p2_x,p2_y)
                profile_EQ.append([row.date_heure,x_profile,row.elev_EQ,row.erreur_horiz,row.erreur_vert,row.magnitude,row.duree_j])
        
        with open(filename_profile_EQ, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # write header
            writer.writerow(['EQ_name', 'distance_along_profile', 'z', 'horiz_error', 'vert_error','magnitude','duree_jours'])       
            # Write each station's data
            for EQ in profile_EQ:
                writer.writerow(EQ)