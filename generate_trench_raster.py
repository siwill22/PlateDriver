#!/usr/bin/env python3
"""
generate_trench_raster.py

Generate a global trench-depth anomaly raster from a plate boundary GeoJSON
(as written by run_SDP.py at each csv_save_interval step).

Subduction segments are selected by boundary_type == 'subduction' (or by a
divergence threshold, in km/Myr = mm/yr, when that column is absent).  Segments
are densified along great circles; for every raster pixel the great-circle
distance to the nearest subduction-segment point is computed; the trench profile is a Gaussian trough centred on the trench line:

    depth(d) = -peak_depth_m * exp( -(d / sigma_deg)^2 )

All values are negative (depressing the bathymetry).  Pixels further than
max_dist_deg from any subduction segment are set to exactly zero.

Add the output to a bathymetry grid to impose ocean trenches (only on ocean
pixels — generate_topography.py does this, masking out land).

Usage (command line):
    python generate_trench_raster.py output/boundaries_0020.geojson
    python generate_trench_raster.py output/boundaries_0020.geojson \\
        --output output/trench_0020.nc --spacing 0.5 \\
        --peak-depth 6000 --sigma 1.5 --max-dist 6

The function generate_trench_raster() can also be imported directly (it is
defined in generate_topography.py and re-exported here).
"""

import argparse
import os

# The trench code lives in generate_topography.py (shared with the full
# topography pipeline); it is re-exported here so existing imports such as
# `from generate_trench_raster import generate_trench_raster` keep working.
from generate_topography import (   # noqa: F401
    _densify_segment,
    _collect_subduction_points,
    _check_not_input,
    generate_trench_raster,
)


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
    parser.add_argument('--spacing', type=float, default=0.1,
                        help='Raster spacing in degrees (default: 0.1)')
    parser.add_argument('--peak-depth', type=float, default=2000.0,
                        help='Trench depth at the subduction line in metres (default: 2000)')
    parser.add_argument('--sigma', type=float, default=1.5,
                        help='Gaussian half-width in degrees (default: 1.5 ≈ 165 km)')
    parser.add_argument('--max-dist', type=float, default=4.0,
                        help='Maximum distance from trench to compute, degrees (default: 4)')
    parser.add_argument('--div-threshold', type=float, default=-0.0001,
                        help='Fallback divergence threshold for subduction, km/Myr = mm/yr '
                             '(default: -0.0001)')
    args = parser.parse_args()

    if args.output:
        outpath = args.output
    else:
        # boundaries_0060.geojson -> trench_0060.nc in the same directory
        # (only the file name is changed, never the directory part)
        folder, fname = os.path.split(args.boundaries_file)
        stem = os.path.splitext(fname)[0]
        stem = stem.replace('boundaries', 'trench') if 'boundaries' in stem \
               else 'trench_' + stem
        outpath = os.path.join(folder, stem + '.nc')
    _check_not_input(outpath, args.boundaries_file)

    da = generate_trench_raster(
        args.boundaries_file,
        spacing=args.spacing,
        peak_depth_m=args.peak_depth,
        sigma_deg=args.sigma,
        max_dist_deg=args.max_dist,
        div_threshold=args.div_threshold)

    da.to_netcdf(outpath)
    nonzero = int((da.values != 0).sum())
    print(f'Saved {outpath}')
    print(f'  Grid: {da.sizes["lat"]} × {da.sizes["lon"]} at {args.spacing}°')
    print(f'  Trench pixels: {nonzero} ({100 * nonzero / da.size:.1f} %)')
    print(f'  Min value: {float(da.min()):.0f} m   Max: {float(da.max()):.0f} m')


if __name__ == '__main__':
    main()
