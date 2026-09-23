# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 12:25:33 2025

@author: charroyj
Extract vertical resistivity profiles from a FEMTIC mesh.
Each row in the output corresponds to one tetrahedral element,
with centroid coordinates and element-averaged resistivity.

Coordinate system:
  X = northing (km)
  Y = easting  (km)
  Z = depth, positive downward from sea level (km)
"""

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# 1️⃣  Read mesh.dat
# ---------------------------------------------------------------------
def read_mesh_dat(mesh_file):
    """Read FEMTIC mesh.dat and return node coordinates + element connectivity."""
    with open(mesh_file) as f:
        lines = f.readlines()

    # Locate "TETRA" keyword
    i_tetra = lines.index("TETRA\n")

    # --- Node section ---
    n_nodes = int(lines[i_tetra + 1].strip())
    node_data = np.loadtxt(lines[i_tetra + 2 : i_tetra + 2 + n_nodes])
    node_ids = node_data[:, 0].astype(int)
    node_coords = node_data[:, 1:4] / 1000.0  # m → km
    node_index = {nid: i for i, nid in enumerate(node_ids)}

    # --- Element section ---
    start_elem = i_tetra + 2 + n_nodes
    n_elem = int(lines[start_elem].strip())
    elem_data = np.loadtxt(lines[start_elem + 1 : start_elem + 1 + n_elem], dtype=float).astype(int)

    elem_ids = elem_data[:, 0]
    elem_nodes = elem_data[:, 5:9]
    elem_nodes = np.vectorize(node_index.get)(elem_nodes)

    return node_coords, elem_ids.astype(int), elem_nodes


# ---------------------------------------------------------------------
# 2️⃣  Read resistivity block (element → parameter + resistivity)
# ---------------------------------------------------------------------
def read_resistivity_block(block_file):
    """Read resistivity_block_iterXX.dat and return element → parameter + param table."""
    with open(block_file) as f:
        lines = f.readlines()

    n_elem, n_param = map(int, lines[0].split())

    mapping = np.loadtxt(lines[1 : 1 + n_elem], dtype=int)
    elem_to_param = dict(mapping)

    param_data = np.loadtxt(lines[1 + n_elem :], dtype=float)
    param_df = pd.DataFrame(param_data, columns=[
        "param_id", "resistivity", "lower", "upper", "nparam", "is_fixed"
    ])
    param_df["param_id"] = param_df["param_id"].astype(int)

    return elem_to_param, param_df


# ---------------------------------------------------------------------
# 3️⃣  Compute element centroids
# ---------------------------------------------------------------------
def compute_centroids(node_coords, elem_nodes):
    """Compute centroid (mean of 4 node coordinates) for each tetrahedron."""
    return node_coords[elem_nodes].mean(axis=1)


# ---------------------------------------------------------------------
# 4️⃣  Extract one vertical profile
# ---------------------------------------------------------------------
def vertical_profile(
    x0, y0, zsurf_km,
    node_coords, elem_ids, elem_nodes,
    elem_to_param, param_df,
    zmin=5.0, r0=0.05, r5km=0.7, sky_res=1e9,
    param_df_dec1=None, param_df_dec9=None,
):
    """Select elements within a variable-radius cylinder under (x0, y0)
       and collapse them into one averaged coordinate per parameter ID.
    """
    # --- Compute centroids of each tetrahedron ---
    centroids = compute_centroids(node_coords, elem_nodes)
    Xc, Yc, Zc = centroids.T

    # --- Map each element to its parameter ID ---
    elem_params = np.array([elem_to_param[eid] for eid in elem_ids])
    resistivities = param_df.set_index("param_id").loc[elem_params, "resistivity"].values

    # --- Compute distance to station and variable search radius ---
    dist = np.sqrt((Xc - x0) ** 2 + (Yc - y0) ** 2)
    a = np.log(r5km / r0) / 5.0
    variable_radius = np.where(Zc <= 5.0, r0 * np.exp(a * Zc), r5km)

    # --- Select valid elements ---
    mask = (dist < variable_radius) & (Zc >= zsurf_km) & (Zc <= zmin) & (resistivities < sky_res)
    nsel = mask.sum()
    if nsel == 0:
        print(f"→ No cells selected for this station.")
        return pd.DataFrame()

    print(f"→ {nsel:6d} cells selected | Z range: {Zc[mask].min():.2f}–{Zc[mask].max():.2f} km")

    # --- Build dataframe for selected elements ---
    df = pd.DataFrame({
        "elem_id": elem_ids[mask],
        "param_id": elem_params[mask],
        "X_north_km": Xc[mask],
        "Y_east_km": Yc[mask],
        "Z_km": Zc[mask],
        "resistivity": resistivities[mask],
        "radius_used_km": variable_radius[mask],
    })

    # --- Collapse by parameter ID (average coords only) ---
    grouped = (
        df.groupby("param_id", as_index=False)
          .agg({
              "X_north_km": "mean",
              "Y_east_km": "mean",
              "Z_km": "mean",
              "resistivity": "first",      # same within param
              "radius_used_km": "mean"
          })
          .sort_values("Z_km")
    )

    # --- Add optional decile resistivities ---
    if param_df_dec1 is not None:
        grouped["resistivity_dec1"] = param_df_dec1.set_index("param_id").loc[
            grouped["param_id"], "resistivity"
        ].values
    if param_df_dec9 is not None:
        grouped["resistivity_dec9"] = param_df_dec9.set_index("param_id").loc[
            grouped["param_id"], "resistivity"
        ].values

    return grouped



# ---------------------------------------------------------------------
# 5️⃣  Driver for all stations
# ---------------------------------------------------------------------
def extract_all_profiles(
    mesh_file, block_file, loc_file, output_csv,
    zmin=5.0, r0=0.05, r5km=0.7, sky_res=1e9,
    **kwargs
):
    """Extract one element-based vertical profile per station."""
    node_coords, elem_ids, elem_nodes = read_mesh_dat(mesh_file)
    elem_to_param, param_df = read_resistivity_block(block_file)

    # optional deciles
    param_df_dec1 = None
    param_df_dec9 = None
    if "block_file_dec1" in kwargs and kwargs["block_file_dec1"]:
        _, param_df_dec1 = read_resistivity_block(kwargs["block_file_dec1"])
        print(f"→ Included dec1 from {kwargs['block_file_dec1']}")
    if "block_file_dec9" in kwargs and kwargs["block_file_dec9"]:
        _, param_df_dec9 = read_resistivity_block(kwargs["block_file_dec9"])
        print(f"→ Included dec9 from {kwargs['block_file_dec9']}")

    stations = pd.read_csv(
        loc_file,
        sep=r"\s+",
        names=["X", "Y", "Z", "station"]
    )

    all_profiles = []
    for _, row in stations.iterrows():
        print(f"\nProcessing {row.station:10s}  (N={row.X:.3f}, E={row.Y:.3f}, Z={row.Z:.3f} km)")
        df = vertical_profile(
            x0=row.X, y0=row.Y, zsurf_km=row.Z,
            node_coords=node_coords,
            elem_ids=elem_ids, elem_nodes=elem_nodes,
            elem_to_param=elem_to_param, param_df=param_df,
            zmin=zmin, r0=r0, r5km=r5km, sky_res=sky_res,
            param_df_dec1=param_df_dec1, param_df_dec9=param_df_dec9
        )
        df["station"] = row.station
        df["Zsurf_km"] = row.Z
        all_profiles.append(df)

    if all_profiles:
        out = pd.concat(all_profiles, ignore_index=True)
        out.to_csv(output_csv, index=False, lineterminator="\n")
        print(f"\n✅ Saved all profiles to {output_csv}")
    else:
        print("\n⚠️  No cells matched selection criteria.")
        out = pd.DataFrame()

    return out


# ---------------------------------------------------------------------
# 6️⃣  Example usage
# ---------------------------------------------------------------------
if __name__ == "__main__":
    df_all = extract_all_profiles(
        mesh_file="mesh.dat",
        block_file="../../random_start/median_1.5/resistivity_block_iter6.dat",
        loc_file="loc_profiles.txt",
        output_csv="profiles_resistivity_median_rdm.csv",
        zmin=5.0,   # bottom of profile (km)
        r0=0.1,     # radius near surface (km)
        r5km=1.0,   # radius at 5 km depth
        sky_res=1e9,
        block_file_dec1="../resistivity_block_decile1.dat",
        block_file_dec9="../resistivity_block_decile9.dat"
    )

    print("\nFirst few lines:\n", df_all.head())