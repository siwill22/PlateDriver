#!/usr/bin/env python3
"""
blend_topobathy.py

Merge a synthetic topography/bathymetry grid with a reference (observed)
bathymetry, then impose synthetic trench depths on the result.

Pipeline
--------
1. Load three grids: synthetic topo, reference bathy, trench anomaly.
2. Build a "land/shallow" mask: synthetic topo pixels above ELEV_THRESHOLD_M.
3. Expand that mask outward by BUFFER_DEG (spherical great-circle distance),
   creating a transition band of NaN values between the two datasets.
4. Merge: synthetic topo where land mask is True, reference bathy elsewhere,
   NaN inside the transition band.
5. Fill the transition band with a Gaussian-weighted blend of surrounding
   known values — yields a smooth land-to-sea transition.
6. Add the trench depth anomaly (negative values, zero elsewhere).
7. Save result to NetCDF.

Usage (command line)
--------------------
    python blend_topobathy.py
    python blend_topobathy.py \\
        --topo  output/topography_0060.nc \\
        --bathy input/reference_bathymetry.nc \\
        --trench output/trench_0060.nc \\
        --output output/blended_0060.nc \\
        --threshold -1500 --buffer 2 --sigma 3

The function blend_topobathy() can also be imported directly.
"""

import argparse
import sys

import numpy as np
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import BallTree
import xarray as xr


# ---------------------------------------------------------------------------
# Default paths and parameters — edit these for your run
# ---------------------------------------------------------------------------
DEFAULT_TOPO_FILE        = './output/topo_bathy_0060.nc'
DEFAULT_BATHY_FILE       = '/Users/simon/GIT/AbFab/output/atlantis5_world_bathymetry.nc'
DEFAULT_TRENCH_FILE      = './output/trench_0060.nc'
DEFAULT_OUTPUT_FILE      = './output/blended_0060.nc'

DEFAULT_ELEV_THRESHOLD_M = -1500.0   # pixels above this → "land/shallow" mask
DEFAULT_BUFFER_DEG       =    2.0    # width of NaN transition band (degrees)
DEFAULT_FILL_SIGMA_DEG   =    3.0    # Gaussian fill kernel radius (degrees)
DEFAULT_OUTPUT_SPACING   =   None    # output grid spacing in degrees; None = use topo grid


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _norm_coords(da):
    """Rename x/y dims to lon/lat if needed."""
    renames = {}
    if 'x' in da.dims:
        renames['x'] = 'lon'
    if 'y' in da.dims:
        renames['y'] = 'lat'
    return da.rename(renames) if renames else da


def _expand_mask_spherical(land_mask, lats, lons, buffer_deg):
    """Return boolean array: True where land_mask is True OR within buffer_deg.

    Uses a haversine BallTree (same pattern as generate_trench_raster.py).
    """
    rows, cols = np.where(land_mask)
    if len(rows) == 0:
        return land_mask.copy()

    lon2d, lat2d = np.meshgrid(lons, lats)

    src_rad = np.deg2rad(np.column_stack([lat2d[rows, cols], lon2d[rows, cols]]))
    tree = BallTree(src_rad, metric='haversine')

    all_rad = np.deg2rad(np.column_stack([lat2d.ravel(), lon2d.ravel()]))
    dists_rad, _ = tree.query(all_rad, k=1)
    within = np.degrees(dists_rad.ravel()) <= buffer_deg
    return within.reshape(land_mask.shape)


def _gaussian_fill(vals, sigma_px):
    """Fill NaN pixels with a Gaussian-weighted average of known neighbours.

    Applies gaussian_filter separately to the data (NaN→0) and a binary
    weight mask, then divides to recover the weighted mean.  Only NaN pixels
    are replaced; known values are preserved exactly.
    """
    nan_mask = ~np.isfinite(vals)
    if not nan_mask.any():
        return vals.copy()

    data    = np.where(nan_mask, 0.0, vals).astype(float)
    weights = np.where(nan_mask, 0.0, 1.0).astype(float)

    smooth_data    = gaussian_filter(data,    sigma=sigma_px)
    smooth_weights = gaussian_filter(weights, sigma=sigma_px)

    filled = smooth_data / np.where(smooth_weights > 1e-6, smooth_weights, 1.0)
    result = vals.copy()
    result[nan_mask] = filled[nan_mask]
    return result


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def blend_topobathy(topo_file, bathy_file, trench_file, output_file,
                    elev_threshold_m=DEFAULT_ELEV_THRESHOLD_M,
                    buffer_deg=DEFAULT_BUFFER_DEG,
                    fill_sigma_deg=DEFAULT_FILL_SIGMA_DEG,
                    output_spacing=DEFAULT_OUTPUT_SPACING):
    """Blend synthetic topo with reference bathy and impose trench anomaly.

    Parameters
    ----------
    topo_file : str
        Synthetic topography/bathymetry NetCDF.
    bathy_file : str
        Reference (observed) bathymetry NetCDF, e.g. ETOPO or GEBCO.
    trench_file : str
        Trench depth anomaly NetCDF (negative values at trench locations,
        zero elsewhere) — as produced by generate_trench_raster.py.
    output_file : str
        Output NetCDF path.
    elev_threshold_m : float
        Elevation threshold in metres for the land/shallow mask (default -1500).
    buffer_deg : float
        Transition band width in degrees of spherical distance (default 2).
    fill_sigma_deg : float
        Gaussian smoothing radius in degrees used to fill the transition band
        (default 3).  Larger values give broader, smoother transitions.
    output_spacing : float or None
        Output grid resolution in degrees.  None (default) keeps the topo
        grid's native resolution.  All three input grids are interpolated onto
        the requested grid before any processing.
    """
    print('Loading grids ...')
    topo   = _norm_coords(xr.open_dataarray(topo_file))
    bathy  = _norm_coords(xr.open_dataarray(bathy_file))
    trench = _norm_coords(xr.open_dataarray(trench_file))

    print(f'  Topo   {dict(topo.sizes)}')
    print(f'  Bathy  {dict(bathy.sizes)}')
    print(f'  Trench {dict(trench.sizes)}')

    # Build output grid — either user-specified spacing or topo native
    if output_spacing is not None:
        spacing = float(output_spacing)
        lons = np.arange(-180.0 + spacing / 2.0, 180.0, spacing)
        lats = np.arange( -90.0 + spacing / 2.0,  90.0, spacing)
        print(f'Output grid: {len(lats)} × {len(lons)} at {spacing}°')
    else:
        lons = topo.lon.values
        lats = topo.lat.values
        spacing = float(abs(lons[1] - lons[0]))
        print(f'Output grid: {len(lats)} × {len(lons)} at {spacing}° (topo native)')

    # Regrid all inputs onto the output grid
    topo_r   = topo.interp(lat=lats, lon=lons, method='linear')
    bathy_r  = bathy.interp(lat=lats, lon=lons, method='linear')
    trench_r = trench.interp(lat=lats, lon=lons, method='linear').fillna(0.0)

    topo_v   = topo_r.values.astype(float)
    bathy_v  = bathy_r.values.astype(float)
    trench_v = trench_r.values.astype(float)

    # ------------------------------------------------------------------
    # Step 1: land / shallow mask
    # ------------------------------------------------------------------
    land_mask = topo_v > elev_threshold_m
    print(f'\nLand/shallow pixels (> {elev_threshold_m} m): '
          f'{land_mask.sum():,}  ({100 * land_mask.mean():.1f} %)')

    # ------------------------------------------------------------------
    # Step 2: expand mask outward by buffer_deg (spherical BallTree)
    # ------------------------------------------------------------------
    print(f'Expanding mask by {buffer_deg}° spherical buffer ...')
    expanded_mask = _expand_mask_spherical(land_mask, lats, lons, buffer_deg)
    buffer_only   = expanded_mask & ~land_mask
    print(f'Transition band: {buffer_only.sum():,}  ({100 * buffer_only.mean():.1f} %)')

    # ------------------------------------------------------------------
    # Step 3: merged grid
    #   land_mask  → synthetic topo
    #   ~expanded  → reference bathy
    #   buffer_only → NaN (transition band)
    # ------------------------------------------------------------------
    merged = np.where(land_mask, topo_v, bathy_v).astype(float)
    merged[buffer_only] = np.nan

    # ------------------------------------------------------------------
    # Step 4: fill transition band with Gaussian-weighted blend
    # ------------------------------------------------------------------
    sigma_px = fill_sigma_deg / spacing
    print(f'Filling transition band  (σ = {fill_sigma_deg}° = {sigma_px:.1f} px) ...')
    filled = _gaussian_fill(merged, sigma_px)

    # ------------------------------------------------------------------
    # Step 5: add trench anomaly
    # ------------------------------------------------------------------
    print('Imposing trench anomaly ...')
    final = filled + trench_v

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    result = xr.DataArray(
        final,
        dims=['lat', 'lon'],
        coords={'lat': lats, 'lon': lons},
        attrs={
            'units':             'meters',
            'long_name':         'blended topography-bathymetry with trench anomaly',
            'elev_threshold_m':  elev_threshold_m,
            'buffer_deg':        buffer_deg,
            'fill_sigma_deg':    fill_sigma_deg,
        })
    result.to_netcdf(output_file)

    trench_px = int((trench_v != 0).sum())
    print(f'\nSaved {output_file}')
    print(f'  Min: {float(result.min()):.0f} m   Max: {float(result.max()):.0f} m')
    print(f'  Trench pixels: {trench_px:,}  ({100 * trench_px / trench_v.size:.1f} %)')
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--topo',   default=DEFAULT_TOPO_FILE,
                        help='Synthetic topo NetCDF  (default: %(default)s)')
    parser.add_argument('--bathy',  default=DEFAULT_BATHY_FILE,
                        help='Reference bathy NetCDF (default: %(default)s)')
    parser.add_argument('--trench', default=DEFAULT_TRENCH_FILE,
                        help='Trench anomaly NetCDF  (default: %(default)s)')
    parser.add_argument('--output', '-o', default=DEFAULT_OUTPUT_FILE,
                        help='Output NetCDF          (default: %(default)s)')
    parser.add_argument('--threshold', type=float, default=DEFAULT_ELEV_THRESHOLD_M,
                        help='Elevation threshold for land mask in metres (default: %(default)s)')
    parser.add_argument('--buffer',    type=float, default=DEFAULT_BUFFER_DEG,
                        help='Transition band width in degrees          (default: %(default)s)')
    parser.add_argument('--sigma',     type=float, default=DEFAULT_FILL_SIGMA_DEG,
                        help='Gaussian fill radius in degrees           (default: %(default)s)')
    parser.add_argument('--spacing',   type=float, default=DEFAULT_OUTPUT_SPACING,
                        help='Output grid spacing in degrees; omit to use topo native resolution')
    args = parser.parse_args()

    blend_topobathy(
        topo_file        = args.topo,
        bathy_file       = args.bathy,
        trench_file      = args.trench,
        output_file      = args.output,
        elev_threshold_m = args.threshold,
        buffer_deg       = args.buffer,
        fill_sigma_deg   = args.sigma,
        output_spacing   = args.spacing,
    )


if __name__ == '__main__':
    main()
