"""
analyse_boundaries.py  —  compute boundary length and crust flux time series.

Reads boundaries_NNNN.geojson files from an output directory and produces
boundary_stats.csv with one row per saved timestep.

Divergence values in the boundary files are in km/Myr (= mm/yr), so
flux = length (km) × |rate| (km/Myr) is in km²/Myr.

Usage:
    python analyse_boundaries.py --output_dir ./output [--delta_time 1.0]

--delta_time (Myr per model step, as in config.yaml) adds a time_myr column.
"""

import argparse
import glob
import os
import re

import geopandas as gpd
import pandas as pd
import pygplates

EARTH_RADIUS_KM = pygplates.Earth.mean_radius_in_kms


def arc_length_km(geom):
    """Return spherical great-circle length of a LineString in km."""
    coords = list(geom.coords)
    if len(coords) < 2:
        return 0.0
    try:
        pts = [pygplates.PointOnSphere(c[1], c[0]) for c in coords]
        if len(pts) == 2:
            arc = pygplates.GreatCircleArc(pts[0], pts[1])
            return arc.get_arc_length() * EARTH_RADIUS_KM
        else:
            poly = pygplates.PolylineOnSphere(pts)
            return poly.get_arc_length() * EARTH_RADIUS_KM
    except Exception:
        return 0.0


def process_file(path):
    gdf = gpd.read_file(path)
    if gdf.empty or 'boundary_type' not in gdf.columns:
        return None

    gdf['length_km'] = gdf.geometry.apply(arc_length_km)

    # Divergence is already km/Myr
    div = gdf['divergence'].astype(float) if 'divergence' in gdf.columns else 0.0

    ridge = gdf['boundary_type'] == 'ridge'
    sub   = gdf['boundary_type'] == 'subduction'
    cons  = gdf['boundary_type'] == 'conservative'

    ridge_len  = gdf.loc[ridge, 'length_km'].sum()
    sub_len    = gdf.loc[sub,   'length_km'].sum()
    cons_len   = gdf.loc[cons,  'length_km'].sum()

    # Flux = sum(length_i × |rate_i|) in km² / Myr
    ridge_flux = (gdf.loc[ridge, 'length_km'] * div[ridge].abs()).sum()
    sub_flux   = (gdf.loc[sub,   'length_km'] * div[sub].abs()).sum()

    return {
        'ridge_length_km':           ridge_len,
        'subduction_length_km':      sub_len,
        'conservative_length_km':    cons_len,
        'ridge_flux_km2_per_myr':    ridge_flux,
        'subduction_flux_km2_per_myr': sub_flux,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', default='./output')
    parser.add_argument('--delta_time', type=float, default=None,
                        help='Myr per model step; adds a time_myr column')
    args = parser.parse_args()

    pattern = os.path.join(args.output_dir, 'boundaries_*.geojson')
    files = sorted(glob.glob(pattern))
    if not files:
        print(f'No boundary files found in {args.output_dir}')
        return

    rows = []
    for path in files:
        m = re.search(r'boundaries_(\d+)\.geojson$', path)
        if not m:
            continue
        timestep = int(m.group(1))
        result = process_file(path)
        if result is not None:
            result['timestep'] = timestep
            rows.append(result)
            print(f'  step {timestep:4d}: '
                  f'ridge {result["ridge_length_km"]:8.0f} km  '
                  f'sub {result["subduction_length_km"]:8.0f} km  '
                  f'gen {result["ridge_flux_km2_per_myr"]:.3e} km²/Myr  '
                  f'dest {result["subduction_flux_km2_per_myr"]:.3e} km²/Myr')

    if not rows:
        print('No valid boundary files processed.')
        return

    df = pd.DataFrame(rows).sort_values('timestep').reset_index(drop=True)
    cols = ['timestep', 'ridge_length_km', 'subduction_length_km',
            'conservative_length_km', 'ridge_flux_km2_per_myr',
            'subduction_flux_km2_per_myr']
    if args.delta_time is not None:
        df['time_myr'] = df['timestep'] * args.delta_time
        cols.insert(1, 'time_myr')
    df = df[cols]

    out_path = os.path.join(args.output_dir, 'boundary_stats.csv')
    df.to_csv(out_path, index=False, float_format='%.3f')
    print(f'\nWritten: {out_path}')


if __name__ == '__main__':
    main()
