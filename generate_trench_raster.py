#!/usr/bin/env python3
"""
generate_trench_raster.py

Generate a global trench-depth anomaly raster from a plate boundary GeoJSON
(as written by run_SDP.py at each csv_save_interval step).

Subduction segments are selected by boundary_type == 'subduction' (or by a
divergence threshold when that column is absent).  For every raster pixel the
great-circle distance to the nearest subduction-segment point is computed; the
trench profile is a Gaussian trough centred on the trench line:

    depth(d) = -peak_depth_m * exp( -(d / sigma_deg)^2 )

All values are negative (depressing the bathymetry).  Pixels further than
max_dist_deg from any subduction segment are set to exactly zero.

Subtract the output from a bathymetry grid to impose realistic ocean trenches.

Usage (command line):
    python generate_trench_raster.py output/boundaries_0020.geojson
    python generate_trench_raster.py output/boundaries_0020.geojson \\
        --output output/trench_0020.nc --spacing 0.5 \\
        --peak-depth 6000 --sigma 1.5 --max-dist 6

The function generate_trench_raster() can also be imported directly.
"""

import argparse
import json
import os
import sys

import numpy as np
import xarray as xr
from sklearn.neighbors import BallTree


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _densify_segment(lon1, lat1, lon2, lat2, spacing_deg=0.2):
    """Linearly interpolate between two points at spacing_deg intervals.

    Simple linear interpolation in lon/lat is adequate for the short arc
    segments (~0.5 deg) stored in the boundary GeoJSON.
    """
    dist = np.hypot(lon2 - lon1, lat2 - lat1)
    n = max(int(dist / spacing_deg), 1)
    lons = np.linspace(lon1, lon2, n, endpoint=False)
    lats = np.linspace(lat1, lat2, n, endpoint=False)
    return list(zip(lons, lats))


def _collect_subduction_points(features, div_threshold=-0.0001,
                                tessellation_deg=0.2):
    """Return (lon, lat) arrays for all points along subduction segments."""
    pts = []
    for feat in features:
        props = feat.get('properties', {})
        btype = props.get('boundary_type')
        div   = props.get('divergence', 0.0)

        is_sub = (btype == 'subduction') if btype is not None \
                 else (div is not None and div < div_threshold)
        if not is_sub:
            continue

        coords = feat['geometry'].get('coordinates', [])
        if len(coords) < 2:
            continue
        for i in range(len(coords) - 1):
            lon1, lat1 = coords[i]
            lon2, lat2 = coords[i + 1]
            pts.extend(_densify_segment(lon1, lat1, lon2, lat2, tessellation_deg))
        # Include the final endpoint
        pts.append(tuple(coords[-1][:2]))

    if not pts:
        return np.empty((0,)), np.empty((0,))
    arr = np.array(pts)
    return arr[:, 0], arr[:, 1]   # lons, lats


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def generate_trench_raster(boundaries_file,
                            spacing=0.5,
                            peak_depth_m=6000.0,
                            sigma_deg=1.5,
                            max_dist_deg=6.0,
                            div_threshold=-0.0001,
                            tessellation_deg=0.2):
    """Generate a trench-depth anomaly raster from a boundary GeoJSON.

    Parameters
    ----------
    boundaries_file : str
        Path to a boundaries_NNNN.geojson produced by run_SDP.py.
    spacing : float
        Output raster resolution in degrees.
    peak_depth_m : float
        Trench depth at the subduction line (positive metres; sign is negated
        in the output so the grid values are negative/deepening).
    sigma_deg : float
        Gaussian half-width in degrees (~100 km at 1 deg).
    max_dist_deg : float
        Pixels further than this from any subduction point are set to 0.
    div_threshold : float
        Fallback convergence threshold when 'boundary_type' is absent.
    tessellation_deg : float
        Spacing used when densifying boundary segment coordinates.

    Returns
    -------
    xr.DataArray  shape (nlat, nlon), dims ['lat', 'lon'], values in metres.
    """
    with open(boundaries_file) as fh:
        gj = json.load(fh)
    features = gj.get('features', [])

    lons_grid = np.arange(-180.0 + spacing / 2.0, 180.0, spacing)
    lats_grid = np.arange( -90.0 + spacing / 2.0,  90.0, spacing)
    nlon, nlat = len(lons_grid), len(lats_grid)

    # Empty DataArray returned when there are no subduction segments
    def _zero_da():
        return xr.DataArray(
            np.zeros((nlat, nlon), dtype=float),
            dims=['lat', 'lon'],
            coords={'lat': lats_grid, 'lon': lons_grid},
            attrs={'units': 'meters', 'long_name': 'synthetic trench depth anomaly'})

    sub_lons, sub_lats = _collect_subduction_points(
        features, div_threshold=div_threshold,
        tessellation_deg=tessellation_deg)

    if sub_lons.size == 0:
        print('No subduction segments found — returning zero raster.', file=sys.stderr)
        return _zero_da()

    # BallTree with haversine metric (expects lat, lon in radians)
    sub_pts_rad = np.deg2rad(np.column_stack([sub_lats, sub_lons]))
    tree = BallTree(sub_pts_rad, metric='haversine')

    lon2d, lat2d = np.meshgrid(lons_grid, lats_grid)   # (nlat, nlon)
    query_rad = np.deg2rad(np.column_stack([lat2d.ravel(), lon2d.ravel()]))

    dists_rad, _ = tree.query(query_rad, k=1)
    dists_deg = np.degrees(dists_rad.ravel())

    # Gaussian trough (negative = deeper)
    depth = np.where(
        dists_deg <= max_dist_deg,
        -peak_depth_m * np.exp(-(dists_deg / sigma_deg) ** 2),
        0.0)

    return xr.DataArray(
        depth.reshape(nlat, nlon),
        dims=['lat', 'lon'],
        coords={'lat': lats_grid, 'lon': lons_grid},
        attrs={
            'units': 'meters',
            'long_name': 'synthetic trench depth anomaly',
            'peak_depth_m': peak_depth_m,
            'sigma_deg': sigma_deg,
            'max_dist_deg': max_dist_deg,
        })


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('boundaries_file',
                        help='Path to boundaries_NNNN.geojson')
    parser.add_argument('--output', '-o', default=None,
                        help='Output NetCDF path (default: trench_NNNN.nc alongside input)')
    parser.add_argument('--spacing', type=float, default=0.5,
                        help='Raster spacing in degrees (default: 0.5)')
    parser.add_argument('--peak-depth', type=float, default=6000.0,
                        help='Trench depth at the subduction line in metres (default: 6000)')
    parser.add_argument('--sigma', type=float, default=1.5,
                        help='Gaussian half-width in degrees (default: 1.5 ≈ 165 km)')
    parser.add_argument('--max-dist', type=float, default=6.0,
                        help='Maximum distance from trench to compute, degrees (default: 6)')
    parser.add_argument('--div-threshold', type=float, default=-0.0001,
                        help='Fallback divergence threshold for subduction (default: -0.0001)')
    args = parser.parse_args()

    da = generate_trench_raster(
        args.boundaries_file,
        spacing=args.spacing,
        peak_depth_m=args.peak_depth,
        sigma_deg=args.sigma,
        max_dist_deg=args.max_dist,
        div_threshold=args.div_threshold)

    if args.output:
        outpath = args.output
    else:
        base = os.path.splitext(args.boundaries_file)[0]
        outpath = base.replace('boundaries', 'trench') + '.nc'

    da.to_netcdf(outpath)
    nonzero = int((da.values != 0).sum())
    print(f'Saved {outpath}')
    print(f'  Grid: {da.sizes["lat"]} × {da.sizes["lon"]} at {args.spacing}°')
    print(f'  Trench pixels: {nonzero} ({100 * nonzero / da.size:.1f} %)')
    print(f'  Min value: {float(da.min()):.0f} m   Max: {float(da.max()):.0f} m')


if __name__ == '__main__':
    main()
