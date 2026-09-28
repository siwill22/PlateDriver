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
6. Add the trench depth anomaly (negative values, zero elsewhere) on ocean
   pixels only (synthetic topo <= 0 m).
7. Save result to NetCDF.

generate_topography.py runs the same pipeline but builds the trench from the
boundary GeoJSON itself (and can add seamounts); the shared helpers live there.

Usage (command line)
--------------------
    python blend_topobathy.py \\
        --topo  output/topography_0060.nc \\
        --bathy input/reference_bathymetry.nc \\
        --trench output/trench_0060.nc \\
        --output output/blended_0060.nc \\
        --threshold -1500 --buffer 2 --sigma 3

The function blend_topobathy() can also be imported directly.
"""

import argparse

import numpy as np
import xarray as xr

# Shared helpers (defined once in generate_topography.py)
from generate_topography import (   # noqa: F401
    _norm_coords,
    _expand_mask_spherical,
    _gaussian_fill,
    _prepare_for_interp,
    _output_grid,
    _warn_missing,
    _check_not_input,
)


# ---------------------------------------------------------------------------
# Default parameters
# ---------------------------------------------------------------------------
DEFAULT_ELEV_THRESHOLD_M = -1500.0   # pixels above this → "land/shallow" mask
DEFAULT_BUFFER_DEG       =    4.0    # width of NaN transition band (degrees)
DEFAULT_FILL_SIGMA_DEG   =    3.0    # Gaussian fill kernel radius (degrees)
DEFAULT_OUTPUT_SPACING   =   None    # output grid spacing in degrees; None = use topo grid


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
        Transition band width in degrees of spherical distance (default 4).
    fill_sigma_deg : float
        Gaussian smoothing radius in degrees used to fill the transition band
        (default 3).  Larger values give broader, smoother transitions.
    output_spacing : float or None
        Output grid resolution in degrees.  None (default) keeps the topo
        grid's native resolution.  All three input grids are interpolated onto
        the requested grid before any processing.
    """
    _check_not_input(output_file, topo_file, bathy_file, trench_file)

    print('Loading grids ...')
    topo   = _norm_coords(xr.open_dataarray(topo_file))
    bathy  = _norm_coords(xr.open_dataarray(bathy_file))
    trench = _norm_coords(xr.open_dataarray(trench_file))

    print(f'  Topo   {dict(topo.sizes)}')
    print(f'  Bathy  {dict(bathy.sizes)}')
    print(f'  Trench {dict(trench.sizes)}')

    # Build output grid — either user-specified spacing or topo native
    lons, lats, spacing = _output_grid(topo, output_spacing)

    # Regrid all inputs onto the output grid.  Each input is first converted
    # to -180..180, sorted ascending and wrap-padded, so 0..360 or descending
    # grids interpolate correctly right up to the antimeridian.
    topo_r   = _prepare_for_interp(topo).interp(lat=lats, lon=lons, method='linear')
    bathy_r  = _prepare_for_interp(bathy).interp(lat=lats, lon=lons, method='linear')
    trench_r = _prepare_for_interp(trench).interp(lat=lats, lon=lons, method='linear')

    topo_v   = topo_r.values.astype(float)
    bathy_v  = bathy_r.values.astype(float)
    trench_v = trench_r.values.astype(float)
    _warn_missing(topo_v, 'Synthetic topo')
    _warn_missing(trench_v, 'Trench anomaly')
    trench_v = np.nan_to_num(trench_v, nan=0.0)

    # ------------------------------------------------------------------
    # Step 1: land / shallow mask
    # ------------------------------------------------------------------
    land_mask = topo_v > elev_threshold_m
    _warn_missing(bathy_v, 'Reference bathy', where=~land_mask)
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
    filled = _gaussian_fill(merged, sigma_px, lons=lons)

    # ------------------------------------------------------------------
    # Step 5: add trench anomaly — ocean pixels only, so overriding /
    # colliding continents are not carved below sea level
    # ------------------------------------------------------------------
    print('Imposing trench anomaly ...')
    trench_v[topo_v > 0] = 0.0   # no trench on land
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
    parser.add_argument('--topo',   required=True, metavar='NETCDF',
                        help='Synthetic topo NetCDF')
    parser.add_argument('--bathy',  required=True, metavar='NETCDF',
                        help='Reference bathy NetCDF')
    parser.add_argument('--trench', required=True, metavar='NETCDF',
                        help='Trench anomaly NetCDF (from generate_trench_raster.py)')
    parser.add_argument('--output', '-o', required=True, metavar='NETCDF',
                        help='Output NetCDF')
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
