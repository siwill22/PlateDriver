"""Per-step kinematic statistics, written by run_SDP.py every step.

step_stats.csv   one row per step: global RMS speed, net rotation, boundary lengths, timings
plate_stats.csv  one row per plate per step: angular velocity, RMS speed, area,
                 continental fraction and boundary lengths (for Forsyth & Uyeda style plots)

Velocities are those of the stage rotations computed at the step, i.e. the motion of the
step's plate configuration over the next delta_time.  Everything is spherical: speeds
come from the plate area tensors, lengths are great-circle lengths.
"""

import os

import numpy as np
import pandas as pd

import self_driving_plates as sdp

R = sdp.EARTH_RADIUS_KM


def _segment_lengths_km(boundaries_gdf):
    """Great-circle length (km) of every boundary row."""
    coords = boundaries_gdf.geometry.get_coordinates()
    xyz = sdp.lonlat_to_xyz(coords.x.values, coords.y.values)
    row = coords.index.values
    same = row[1:] == row[:-1]
    a, b = xyz[:-1][same], xyz[1:][same]
    arc = np.arctan2(np.linalg.norm(np.cross(a, b), axis=1), np.einsum('ij,ij->i', a, b))
    lengths = pd.Series(arc * R).groupby(row[1:][same]).sum()
    return lengths.reindex(boundaries_gdf.index, fill_value=0.0).values


def boundary_lengths_by_plate(boundaries_gdf, polarity_dict, craton_counts, params):
    """{pid: {perimeter, trench_down, trench_over, ridge, conservative}} in km."""
    out = {}
    if boundaries_gdf is None or not len(boundaries_gdf) or 'divergence' not in boundaries_gdf:
        return out
    length = _segment_lengths_km(boundaries_gdf)
    div = boundaries_gdf['divergence'].values
    p1 = np.nan_to_num(boundaries_gdf['plateid1_poly1'].values.astype(float)).astype(int)
    p2 = np.nan_to_num(boundaries_gdf['plateid1_poly2'].values.astype(float)).astype(int)
    kind = np.where(div > params['div_threshold_ridge'], 'ridge',
                    np.where(div < params['div_threshold_subduction'], 'trench', 'conservative'))
    downgoing = {}
    for a, b in set(zip(p1[kind == 'trench'], p2[kind == 'trench'])):
        downgoing[(a, b)] = sdp.resolve_downgoing_plate(a, b, polarity_dict, craton_counts)

    def add(pid, key, value):
        if pid <= 0:
            return
        d = out.setdefault(pid, dict(perimeter=0.0, trench_down=0.0, trench_over=0.0,
                                     ridge=0.0, conservative=0.0))
        d[key] += value
        if key != 'perimeter':
            d['perimeter'] += value

    for a, b, k, L in zip(p1, p2, kind, length):
        for pid in (a, b):
            if k == 'trench':
                add(pid, 'trench_down' if downgoing[(a, b)] == pid else 'trench_over', L)
            else:
                add(pid, k, L)
    return out


def compute(step, time_myr, rasters, stage_pole_dict, delta_time, boundaries_gdf,
            polarity_dict, params, grid_deg=1.0):
    """Return (step_row, plate_rows) for one step."""
    plate_id = rasters['plate_id']
    pids, J, areas = sdp.plate_area_tensors(plate_id, grid_deg=grid_deg)
    craton = rasters.get('craton')
    if craton is not None:
        cont = np.nan_to_num(craton.values) > 0
        cont_pids, _, cont_areas = sdp.plate_area_tensors(plate_id.where(cont), grid_deg=grid_deg)
        cont_area = dict(zip(cont_pids.tolist(), cont_areas))
    else:
        cont_area = {}
    craton_counts = sdp.craton_counts_from_raster(craton)
    lengths = boundary_lengths_by_plate(boundaries_gdf, polarity_dict, craton_counts, params)

    plate_rows = []
    sum_Jw = np.zeros(3)
    sum_wJw = 0.0
    for pid, Jp, Ap in zip(pids.tolist(), J, areas):
        w = sdp.stage_pole_omega(stage_pole_dict.get(pid), delta_time)
        rms = R * np.sqrt(max(w @ Jp @ w, 0.0) / Ap) if Ap > 0 else 0.0
        rate = np.degrees(np.linalg.norm(w))
        pole_lon, pole_lat = sdp.xyz_to_lonlat(w) if rate > 0 else (0.0, 90.0)
        L = lengths.get(pid, {})
        perim = L.get('perimeter', 0.0)
        plate_rows.append(dict(
            step=step, time_myr=time_myr, plate_id=pid,
            area_frac=Ap / (4 * np.pi), area_km2=Ap * R * R,
            continental_frac=cont_area.get(pid, 0.0) / Ap if Ap > 0 else 0.0,
            pole_lat=float(pole_lat), pole_lon=float(pole_lon), rate_deg_per_myr=rate,
            rms_speed_km_per_myr=rms,
            perimeter_km=perim,
            trench_down_km=L.get('trench_down', 0.0), trench_over_km=L.get('trench_over', 0.0),
            ridge_km=L.get('ridge', 0.0), conservative_km=L.get('conservative', 0.0),
            trench_frac=L.get('trench_down', 0.0) / perim if perim > 0 else 0.0,
            ridge_frac=L.get('ridge', 0.0) / perim if perim > 0 else 0.0,
        ))
        sum_Jw += Jp @ w
        sum_wJw += w @ Jp @ w

    # Net rotation of the lithosphere: ω_nr = (Σ J)⁻¹ Σ J ω  (Σ J = 8π/3 𝟙 for full cover)
    J_total = J.sum(axis=0) if len(J) else np.eye(3)
    w_nr = np.linalg.solve(J_total, sum_Jw) if len(J) else np.zeros(3)
    nr_rate = np.degrees(np.linalg.norm(w_nr))
    nr_lon, nr_lat = sdp.xyz_to_lonlat(w_nr) if nr_rate > 0 else (0.0, 90.0)
    total_area = areas.sum()
    step_row = dict(
        step=step, time_myr=time_myr, n_plates=len(pids),
        rms_speed_km_per_myr=R * np.sqrt(sum_wJw / total_area) if total_area > 0 else 0.0,
        max_plate_rms_km_per_myr=max((r['rms_speed_km_per_myr'] for r in plate_rows), default=0.0),
        net_rotation_deg_per_myr=nr_rate, net_rotation_pole_lat=float(nr_lat),
        net_rotation_pole_lon=float(nr_lon),
        trench_km=sum(r['trench_down_km'] for r in plate_rows),
        ridge_km=sum(r['ridge_km'] for r in plate_rows) / 2,
    )
    return step_row, plate_rows


def append_csv(path, rows):
    """Append rows (dicts) to a CSV, writing the header on first use."""
    if not rows:
        return
    df = pd.DataFrame(rows)
    df.to_csv(path, mode='a', header=not os.path.exists(path), index=False)
