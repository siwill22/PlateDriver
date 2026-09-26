#!/usr/bin/env python3
"""
generate_topography.py

Full topography/bathymetry generation pipeline:

  1. Generate a trench-depth anomaly raster from a plate boundary GeoJSON.
  2. Merge synthetic topo with reference (observed) bathy:
       - Build a land/shallow mask from the synthetic topo.
       - Expand the mask by a buffer to create a NaN transition band.
       - Fill the transition band with a Gaussian-weighted blend.
       - Add the trench anomaly (negative, deepens ocean; not applied on land).
  3. Optionally stamp conical seamounts onto the result (positive, raises
     ocean floor; not applied on land).

Grids may use either -180..180 or 0..360 longitudes and ascending or
descending axes; the output uses -180..180 with ascending axes.
  4. Save to NetCDF.

Usage
-----
    python generate_topography.py \\
        --boundaries output/boundaries_0060.geojson \\
        --topo       output/topo_bathy_0060.nc \\
        --bathy      /path/to/reference_bathy.nc \\
        --output     output/blended_0060.nc

    # with seamounts
    python generate_topography.py \\
        --boundaries output/boundaries_0060.geojson \\
        --topo       output/topo_bathy_0060.nc \\
        --bathy      /path/to/reference_bathy.nc \\
        --seamounts  output/seamounts.geojson \\
        --output     output/blended_0060.nc

All functions can also be imported directly.
"""

import argparse
import json
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import BallTree
import xarray as xr


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

# Blend
DEFAULT_ELEV_THRESHOLD_M = -1500.0
DEFAULT_BUFFER_DEG       =    4.0
DEFAULT_FILL_SIGMA_DEG   =    3.0
DEFAULT_OUTPUT_SPACING   =   None

# Trench
DEFAULT_TRENCH_SPACING   =    0.1   # standalone trench raster resolution
DEFAULT_PEAK_DEPTH_M     = 2000.0
DEFAULT_SIGMA_DEG        =    1.5
DEFAULT_MAX_DIST_DEG     =    4.0
DEFAULT_DIV_THRESHOLD    = -0.0001
DEFAULT_TESSELLATION_DEG =    0.2

# Seamounts
DEFAULT_SEAMOUNTS_FILE   = None
DEFAULT_SM_FRACTION      =  0.3
DEFAULT_SM_JITTER_DEG    =  0.5
DEFAULT_SM_DIAMETER_KM   = 200.0
DEFAULT_SM_HEIGHT_M      = 2000.0


# ---------------------------------------------------------------------------
# Trench helpers
# ---------------------------------------------------------------------------

def _lonlat_to_xyz(lon, lat):
    """Unit vector(s) on the sphere for lon/lat in degrees."""
    lon_r, lat_r = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(lat_r) * np.cos(lon_r),
                     np.cos(lat_r) * np.sin(lon_r),
                     np.sin(lat_r)], axis=-1)


def _densify_segment(lon1, lat1, lon2, lat2, spacing_deg=0.2):
    """Points every ~spacing_deg along the great-circle arc from point 1 to point 2.

    Interpolation is done on the sphere (spherical linear interpolation of
    the unit vectors), so the shorter arc is always followed — a segment from
    179.8°E to 179.8°W stays near the antimeridian instead of running the
    long way round the globe.  The end point itself is not included.
    """
    p1 = _lonlat_to_xyz(lon1, lat1)
    p2 = _lonlat_to_xyz(lon2, lat2)
    omega = np.arccos(np.clip(np.dot(p1, p2), -1.0, 1.0))   # arc length, radians
    if omega < 1e-12:
        return [(float(lon1), float(lat1))]
    n = max(int(np.degrees(omega) / spacing_deg), 1)

    t = np.arange(n)[:, None] / n
    pts = (np.sin((1.0 - t) * omega) * p1 + np.sin(t * omega) * p2) / np.sin(omega)
    lons = np.degrees(np.arctan2(pts[:, 1], pts[:, 0]))
    lats = np.degrees(np.arcsin(np.clip(pts[:, 2], -1.0, 1.0)))
    return list(zip(lons, lats))


def _line_parts(geometry):
    """Return a list of coordinate lists for a LineString or MultiLineString."""
    gtype  = geometry.get('type')
    coords = geometry.get('coordinates', [])
    if gtype == 'LineString':
        return [coords]
    if gtype == 'MultiLineString':
        return list(coords)
    return []


def _collect_subduction_points(features, div_threshold=-0.0001,
                                tessellation_deg=0.2):
    """Return (lon, lat) arrays for all points along subduction segments.

    A feature counts as subduction if boundary_type == 'subduction' or, when
    that property is absent, if its divergence (km/Myr = mm/yr) is below
    div_threshold.
    """
    pts = []
    for feat in features:
        props = feat.get('properties', {})
        btype = props.get('boundary_type')
        div   = props.get('divergence', 0.0)

        is_sub = (btype == 'subduction') if btype is not None \
                 else (div is not None and div < div_threshold)
        if not is_sub:
            continue

        for coords in _line_parts(feat.get('geometry') or {}):
            if len(coords) < 2:
                continue
            for i in range(len(coords) - 1):
                lon1, lat1 = coords[i][:2]          # ignore any z value
                lon2, lat2 = coords[i + 1][:2]
                pts.extend(_densify_segment(lon1, lat1, lon2, lat2, tessellation_deg))
            pts.append(tuple(coords[-1][:2]))       # include the final end point

    if not pts:
        return np.empty((0,)), np.empty((0,))
    arr = np.array(pts, dtype=float)
    return arr[:, 0], arr[:, 1]


# ---------------------------------------------------------------------------
# Grid / longitude helpers
# ---------------------------------------------------------------------------

def _norm_coords(da):
    """Rename x/y dims to lon/lat if needed."""
    renames = {}
    if 'x' in da.dims:
        renames['x'] = 'lon'
    if 'y' in da.dims:
        renames['y'] = 'lat'
    return da.rename(renames) if renames else da


def _wrap_lon_180(lons):
    """Map longitudes given as 0..360 (or -180..180) into -180..180.

    +180 is kept as +180, so a gridline-registered -180..180 axis is unchanged.
    """
    lons = np.asarray(lons, dtype=float)
    return np.where(lons > 180.0, lons - 360.0,
                    np.where(lons < -180.0, lons + 360.0, lons))


def _global_lon_layout(lons):
    """Describe an ascending, evenly spaced longitude axis.

    Returns (is_global, has_seam_duplicate):
      is_global           – the axis covers all 360° of longitude;
      has_seam_duplicate  – first and last columns are the same meridian
                            (gridline registration, e.g. both -180 and +180).
    """
    lons = np.asarray(lons, dtype=float)
    if lons.size < 2:
        return False, False
    dlon = abs(float(lons[1] - lons[0]))
    span = abs(float(lons[-1] - lons[0]))
    tol  = 1e-3 * dlon
    if abs(span - 360.0) < tol:
        return True, True
    if abs(span + dlon - 360.0) < tol:
        return True, False
    return False, False


def _prepare_for_interp(da):
    """Return a copy of a lon/lat grid that interpolates cleanly onto -180..180.

    - longitudes are converted to -180..180 (so 0..360 grids work);
    - both axes are sorted ascending (descending grids work);
    - a duplicated seam column (e.g. both 0 and 360) is removed;
    - for global grids one wrap-around column is added at each end, so points
      between the outermost columns and ±180 interpolate across the seam
      instead of becoming NaN;
    - if the grid stops less than one cell short of a pole (pixel
      registration, e.g. last row at 89.75°), a pole row is added holding the
      mean of the outermost row, so the ±90 output rows are not NaN.
    """
    da = _norm_coords(da)
    da = da.assign_coords(lon=((da.lon.values + 180.0) % 360.0) - 180.0)
    da = da.sortby('lat').sortby('lon')
    _, keep = np.unique(da.lon.values, return_index=True)
    da = da.isel(lon=keep)

    is_global, _ = _global_lon_layout(da.lon.values)
    if is_global:
        west = da.isel(lon=[-1]).assign_coords(lon=da.lon.values[-1:] - 360.0)
        east = da.isel(lon=[0]).assign_coords(lon=da.lon.values[:1] + 360.0)
        da = xr.concat([west, da, east], dim='lon')

        lat = da.lat.values
        dlat = abs(float(lat[1] - lat[0]))
        parts = [da]
        if 0 < lat[0] + 90.0 < dlat:
            parts.insert(0, (da.isel(lat=[0]) * 0 + da.isel(lat=0).mean())
                         .assign_coords(lat=[-90.0]))
        if 0 < 90.0 - lat[-1] < dlat:
            parts.append((da.isel(lat=[-1]) * 0 + da.isel(lat=-1).mean())
                         .assign_coords(lat=[90.0]))
        if len(parts) > 1:
            da = xr.concat(parts, dim='lat')
    return da


def _warn_missing(vals, name, where=None):
    """Print a warning if a regridded input has NaN (uncovered) pixels."""
    missing = ~np.isfinite(vals)
    if where is not None:
        missing &= where
    n = int(missing.sum())
    if n:
        print(f'WARNING: {name} has no data at {n:,} output pixels '
              f'({100 * n / vals.size:.1f} %) — check its lon/lat coverage.  '
              'These pixels are filled from neighbouring values, or 0 where '
              'no data are nearby.', file=sys.stderr)


def _output_grid(topo, output_spacing):
    """Return (lons, lats, spacing) for the blended output grid.

    output_spacing given → global pixel-registered grid at that spacing.
    None → the topo grid's own nodes, converted to -180..180 and sorted
    ascending (a gridline-registered -180..180 topo keeps both edge columns).
    """
    if output_spacing is not None:
        spacing = float(output_spacing)
        lons = np.arange(-180.0 + spacing / 2.0, 180.0, spacing)
        lats = np.arange( -90.0 + spacing / 2.0,  90.0, spacing)
        print(f'Output grid: {len(lats)} × {len(lons)} at {spacing}°')
    else:
        lons = np.unique(_wrap_lon_180(topo.lon.values))
        lats = np.unique(np.asarray(topo.lat.values, dtype=float))
        spacing = float(abs(lons[1] - lons[0]))
        print(f'Output grid: {len(lats)} × {len(lons)} at {spacing}° (topo native)')
    return lons, lats, spacing


def _check_not_input(output_file, *input_files):
    """Refuse to write the output over one of the input files."""
    out = os.path.realpath(output_file)
    for f in input_files:
        if f is not None and os.path.realpath(f) == out:
            raise ValueError(f'Output file {output_file!r} is the same as input '
                             f'{f!r}; refusing to overwrite it.')


# ---------------------------------------------------------------------------
# Blend helpers
# ---------------------------------------------------------------------------

def _expand_mask_spherical(land_mask, lats, lons, buffer_deg):
    """Return boolean array: True where land_mask is True OR within buffer_deg.

    Uses great-circle (haversine) distance via a BallTree, so the buffer has
    the same width in km at all latitudes and wraps across the antimeridian.
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


def _gaussian_fill(vals, sigma_px, lons=None):
    """Fill NaN pixels with a Gaussian-weighted average of known neighbours.

    Applies gaussian_filter separately to the data (NaN→0) and a binary
    weight mask, then divides to recover the weighted mean.  Only NaN pixels
    are replaced; known values are preserved exactly.

    If *lons* is a global longitude axis the filter wraps around in longitude
    (no cliff at the antimeridian); a duplicated ±180 column is dropped
    before filtering and restored afterwards.  Latitude edges use 'nearest'.
    The kernel is isotropic in *degrees*, so in km it is narrower east–west
    than north–south away from the equator.
    """
    nan_mask = ~np.isfinite(vals)
    if not nan_mask.any():
        return vals.copy()

    is_global, has_dup = _global_lon_layout(lons) if lons is not None else (False, False)
    work = vals[:, :-1] if has_dup else vals
    mode = ('nearest', 'wrap') if is_global else 'reflect'   # (lat, lon)

    work_nan = ~np.isfinite(work)
    data     = np.where(work_nan, 0.0, work).astype(float)
    weights  = np.where(work_nan, 0.0, 1.0).astype(float)

    smooth_data    = gaussian_filter(data,    sigma=sigma_px, mode=mode)
    smooth_weights = gaussian_filter(weights, sigma=sigma_px, mode=mode)

    filled = smooth_data / np.where(smooth_weights > 1e-6, smooth_weights, 1.0)
    filled = np.where(work_nan, filled, work)
    if has_dup:
        # the last column is the same meridian as the first
        filled = np.concatenate([filled, filled[:, :1]], axis=1)

    result = vals.copy()
    result[nan_mask] = filled[nan_mask]
    return result


# ---------------------------------------------------------------------------
# Seamount helpers
# ---------------------------------------------------------------------------

def _haversine_km(lon1, lat1, lon2d, lat2d):
    """Return haversine distance in km from (lon1, lat1) to each cell in 2-D arrays."""
    R = 6371.0
    dlat = np.deg2rad(lat2d - lat1)
    dlon = np.deg2rad(lon2d - lon1)
    a = np.sin(dlat / 2) ** 2 + (
        np.cos(np.deg2rad(lat1)) * np.cos(np.deg2rad(lat2d)) * np.sin(dlon / 2) ** 2)
    return 2.0 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def generate_trench_raster(boundaries_file,
                            spacing=DEFAULT_TRENCH_SPACING,
                            peak_depth_m=DEFAULT_PEAK_DEPTH_M,
                            sigma_deg=DEFAULT_SIGMA_DEG,
                            max_dist_deg=DEFAULT_MAX_DIST_DEG,
                            div_threshold=DEFAULT_DIV_THRESHOLD,
                            tessellation_deg=DEFAULT_TESSELLATION_DEG,
                            lons=None, lats=None):
    """Generate a trench-depth anomaly raster from a boundary GeoJSON.

    For every raster node the great-circle distance d to the nearest
    subduction-line point is found and the trench is a Gaussian trough:

        depth(d) = -peak_depth_m * exp( -(d / sigma_deg)^2 ),   d <= max_dist_deg

    Parameters
    ----------
    boundaries_file : str
        Path to a boundaries_NNNN.geojson produced by run_SDP.py.
    spacing : float
        Resolution in degrees of a global gridline-registered raster
        (-180..180 × -90..90).  Ignored if *lons* and *lats* are given.
    peak_depth_m : float
        Trench depth at the subduction line (positive metres; negated in output).
    sigma_deg : float
        Gaussian half-width in degrees (1.5° ≈ 165 km).
    max_dist_deg : float
        Pixels further than this from any subduction point are set to 0.
    div_threshold : float
        Fallback convergence threshold on 'divergence' (km/Myr = mm/yr) when
        'boundary_type' is absent.
    tessellation_deg : float
        Spacing used when densifying boundary segment coordinates.
    lons, lats : 1-D arrays or None
        Evaluate the trench directly on these grid axes (degrees) instead of
        building a grid from *spacing*.

    Returns
    -------
    xr.DataArray  shape (nlat, nlon), dims ['lat', 'lon'], values in metres.
    """
    with open(boundaries_file) as fh:
        gj = json.load(fh)
    features = gj.get('features', [])

    if lons is not None and lats is not None:
        lons_grid = np.asarray(lons, dtype=float)
        lats_grid = np.asarray(lats, dtype=float)
    else:
        lons_grid = np.linspace(-180.0, 180.0, int(round(360.0 / spacing)) + 1)
        lats_grid = np.linspace( -90.0,  90.0, int(round(180.0 / spacing)) + 1)
    nlon, nlat = len(lons_grid), len(lats_grid)

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

    sub_pts_rad = np.deg2rad(np.column_stack([sub_lats, sub_lons]))
    tree = BallTree(sub_pts_rad, metric='haversine')

    lon2d, lat2d = np.meshgrid(lons_grid, lats_grid)
    query_rad = np.deg2rad(np.column_stack([lat2d.ravel(), lon2d.ravel()]))

    dists_rad, _ = tree.query(query_rad, k=1)
    dists_deg = np.degrees(dists_rad.ravel())

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


def generate_seamount_raster(seamounts_file, lons, lats,
                              fraction=DEFAULT_SM_FRACTION,
                              jitter_deg=DEFAULT_SM_JITTER_DEG,
                              diameter_km=DEFAULT_SM_DIAMETER_KM,
                              height_m=DEFAULT_SM_HEIGHT_M,
                              seed=None):
    """Generate a seamount elevation raster on the provided lon/lat grid.

    A random subset of Point features from *seamounts_file* are selected,
    jittered, then stamped onto the grid as conical bumps:

        elev(d) = height_m * max(0,  1 − d / radius_km)

    where d is the haversine distance from the seamount centre.  Overlapping
    seamounts take the element-wise maximum so shapes are not double-counted.

    Parameters
    ----------
    seamounts_file : str
        Path to a GeoJSON file whose features are Point geometries.
    lons, lats : 1-D array-like
        Coordinate axes of the output grid (degrees).  Longitude wraps, so a
        seamount near ±180 is stamped on both sides of the antimeridian.
    fraction : float
        Fraction of seamounts to include, 0–1 (default 0.3).
    jitter_deg : float
        Maximum random coordinate offset in degrees (uniform ±jitter_deg).
    diameter_km : float
        Base diameter of each conical seamount in km (default 200).
    height_m : float
        Peak height of each seamount in metres (default 2000).
    seed : int or None
        Random seed for reproducibility.

    Returns
    -------
    xr.DataArray  shape (nlat, nlon), dims ['lat', 'lon'], values ≥ 0 metres.
    """
    lons = np.asarray(lons)
    lats = np.asarray(lats)
    nlat, nlon = len(lats), len(lons)

    zero_da = xr.DataArray(
        np.zeros((nlat, nlon), dtype=float),
        dims=['lat', 'lon'],
        coords={'lat': lats, 'lon': lons},
        attrs={'units': 'meters', 'long_name': 'synthetic seamount elevation'})

    with open(seamounts_file) as fh:
        gj = json.load(fh)

    all_pts = []
    for feat in gj.get('features', []):
        geom = feat.get('geometry', {})
        if geom.get('type') == 'Point':
            lon, lat = geom['coordinates'][:2]
            all_pts.append((float(lon), float(lat)))

    if not all_pts:
        print('No Point features found in seamounts file — returning zero raster.',
              file=sys.stderr)
        return zero_da

    all_pts = np.array(all_pts)
    rng = np.random.default_rng(seed)
    n_select = max(1, int(len(all_pts) * fraction))
    idx = rng.choice(len(all_pts), size=n_select, replace=False)
    selected = all_pts[idx].copy()

    selected[:, 0] += rng.uniform(-jitter_deg, jitter_deg, size=n_select)
    selected[:, 1] += rng.uniform(-jitter_deg, jitter_deg, size=n_select)
    selected[:, 0] = ((selected[:, 0] + 180.0) % 360.0) - 180.0
    selected[:, 1] = np.clip(selected[:, 1], -90.0, 90.0)

    radius_km = diameter_km / 2.0
    bbox_deg = radius_km / 111.32  # conservative upper bound at equator

    grid = np.zeros((nlat, nlon), dtype=float)

    for slon, slat in selected:
        # Rows/columns inside a lat/lon box around the seamount.  The longitude
        # difference is wrapped to -180..180, so a seamount at 179.9°E also
        # reaches the columns just east of 180°W.
        rows = np.where(np.abs(lats - slat) <= bbox_deg)[0]
        cos_lat = max(abs(np.cos(np.deg2rad(slat))), 0.01)
        dlon = ((lons - slon + 180.0) % 360.0) - 180.0
        cols = np.where(np.abs(dlon) <= bbox_deg / cos_lat)[0]

        if rows.size == 0 or cols.size == 0:
            continue

        sub_lon2d, sub_lat2d = np.meshgrid(lons[cols], lats[rows])
        dist_km = _haversine_km(slon, slat, sub_lon2d, sub_lat2d)

        elev = np.where(dist_km < radius_km,
                        height_m * (1.0 - dist_km / radius_km),
                        0.0)
        box = np.ix_(rows, cols)
        grid[box] = np.maximum(grid[box], elev)

    nonzero = int((grid > 0).sum())
    print(f'Seamounts: {n_select} placed ({fraction*100:.0f}% of {len(all_pts)} total), '
          f'{nonzero:,} non-zero pixels')

    return xr.DataArray(
        grid,
        dims=['lat', 'lon'],
        coords={'lat': lats, 'lon': lons},
        attrs={
            'units':       'meters',
            'long_name':   'synthetic seamount elevation',
            'n_seamounts': n_select,
            'diameter_km': diameter_km,
            'height_m':    height_m,
        })


def blend_topobathy(topo_file, bathy_file, boundaries_file, output_file,
                    elev_threshold_m=DEFAULT_ELEV_THRESHOLD_M,
                    buffer_deg=DEFAULT_BUFFER_DEG,
                    fill_sigma_deg=DEFAULT_FILL_SIGMA_DEG,
                    output_spacing=DEFAULT_OUTPUT_SPACING,
                    trench_spacing=None,
                    peak_depth_m=DEFAULT_PEAK_DEPTH_M,
                    sigma_deg=DEFAULT_SIGMA_DEG,
                    max_dist_deg=DEFAULT_MAX_DIST_DEG,
                    div_threshold=DEFAULT_DIV_THRESHOLD,
                    seamounts_file=DEFAULT_SEAMOUNTS_FILE,
                    seamount_fraction=DEFAULT_SM_FRACTION,
                    seamount_jitter_deg=DEFAULT_SM_JITTER_DEG,
                    seamount_diameter_km=DEFAULT_SM_DIAMETER_KM,
                    seamount_height_m=DEFAULT_SM_HEIGHT_M,
                    seamount_seed=None):
    """Run the full topography generation pipeline.

    Parameters
    ----------
    topo_file : str
        Synthetic topography/bathymetry NetCDF.
    bathy_file : str
        Reference (observed) bathymetry NetCDF, e.g. ETOPO or GEBCO.
    boundaries_file : str
        Plate boundary GeoJSON (produced by run_SDP.py).  Used to generate
        the trench anomaly raster internally.
    output_file : str
        Output NetCDF path.
    elev_threshold_m : float
        Elevation threshold in metres for the land/shallow mask (default -1500).
    buffer_deg : float
        Transition band width in degrees of spherical distance (default 4).
    fill_sigma_deg : float
        Gaussian smoothing radius in degrees for the transition band (default 3).
    output_spacing : float or None
        Output grid resolution in degrees.  None keeps the topo native resolution.
    trench_spacing : float or None
        None (default) evaluates the trench directly on the output grid.  A
        number first builds a global trench raster at that resolution and
        interpolates it onto the output grid (slower; older behaviour).
    peak_depth_m : float
        Trench peak depth in metres (default 2000).
    sigma_deg : float
        Trench Gaussian half-width in degrees (default 1.5).
    max_dist_deg : float
        Max distance from subduction line to compute trench, degrees (default 4).
    div_threshold : float
        Fallback convergence threshold on 'divergence' (km/Myr = mm/yr) when
        'boundary_type' is absent (default -0.0001).
    seamounts_file : str or None
        Path to a seamounts GeoJSON.  If None, no seamounts are added.
    seamount_fraction : float
        Fraction of seamounts to include (default 0.3).
    seamount_jitter_deg : float
        Random coordinate jitter for seamounts in degrees (default 0.5).
    seamount_diameter_km : float
        Conical seamount base diameter in km (default 200).
    seamount_height_m : float
        Seamount peak height in metres (default 2000).
    seamount_seed : int or None
        Random seed for seamount selection and jitter.
    """
    _check_not_input(output_file, topo_file, bathy_file, boundaries_file,
                     seamounts_file)

    print('Loading grids ...')
    topo  = _norm_coords(xr.open_dataarray(topo_file))
    bathy = _norm_coords(xr.open_dataarray(bathy_file))
    print(f'  Topo   {dict(topo.sizes)}')
    print(f'  Bathy  {dict(bathy.sizes)}')

    lons, lats, spacing = _output_grid(topo, output_spacing)

    # Inputs are converted to -180..180, sorted and wrap-padded before
    # interpolation, so 0..360 or descending grids are handled correctly.
    topo_r  = _prepare_for_interp(topo).interp(lat=lats, lon=lons, method='linear')
    bathy_r = _prepare_for_interp(bathy).interp(lat=lats, lon=lons, method='linear')
    topo_v  = topo_r.values.astype(float)
    bathy_v = bathy_r.values.astype(float)
    _warn_missing(topo_v, 'Synthetic topo')

    # ------------------------------------------------------------------
    # Trench — generated from boundaries
    # ------------------------------------------------------------------
    print('\nGenerating trench raster ...')
    trench_kw = dict(peak_depth_m=peak_depth_m, sigma_deg=sigma_deg,
                     max_dist_deg=max_dist_deg, div_threshold=div_threshold)
    if trench_spacing is None:
        # evaluate directly on the output grid
        trench_da = generate_trench_raster(boundaries_file, lons=lons, lats=lats,
                                           **trench_kw)
        trench_v = trench_da.values.astype(float)
    else:
        trench_da = generate_trench_raster(boundaries_file, spacing=trench_spacing,
                                           **trench_kw)
        trench_v = (_prepare_for_interp(trench_da)
                    .interp(lat=lats, lon=lons, method='linear')
                    .fillna(0.0).values.astype(float))

    # ------------------------------------------------------------------
    # Land/shallow mask and transition band
    # ------------------------------------------------------------------
    land_mask = topo_v > elev_threshold_m
    _warn_missing(bathy_v, 'Reference bathy', where=~land_mask)
    print(f'\nLand/shallow pixels (> {elev_threshold_m} m): '
          f'{land_mask.sum():,}  ({100 * land_mask.mean():.1f} %)')

    print(f'Expanding mask by {buffer_deg}° spherical buffer ...')
    expanded_mask = _expand_mask_spherical(land_mask, lats, lons, buffer_deg)
    buffer_only   = expanded_mask & ~land_mask
    print(f'Transition band: {buffer_only.sum():,}  ({100 * buffer_only.mean():.1f} %)')

    merged = np.where(land_mask, topo_v, bathy_v).astype(float)
    merged[buffer_only] = np.nan

    sigma_px = fill_sigma_deg / spacing
    print(f'Filling transition band  (σ = {fill_sigma_deg}° = {sigma_px:.1f} px) ...')
    filled = _gaussian_fill(merged, sigma_px, lons=lons)

    # ------------------------------------------------------------------
    # Trench and seamounts — only applied where the synthetic topo is ocean,
    # so overriding / colliding continents are not carved below sea level.
    # ------------------------------------------------------------------
    print('Imposing trench anomaly ...')
    trench_v[topo_v > 0] = 0.0   # no trench on land
    final = filled + trench_v

    if seamounts_file is not None:
        print('Adding seamount elevations ...')
        sm = generate_seamount_raster(
            seamounts_file, lons, lats,
            fraction=seamount_fraction,
            jitter_deg=seamount_jitter_deg,
            diameter_km=seamount_diameter_km,
            height_m=seamount_height_m,
            seed=seamount_seed)
        sm_v = sm.values.copy()
        sm_v[topo_v > 0] = 0.0   # suppress seamounts on land
        final = final + sm_v

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    result = xr.DataArray(
        final,
        dims=['lat', 'lon'],
        coords={'lat': lats, 'lon': lons},
        attrs={
            'units':            'meters',
            'long_name':        'blended topography-bathymetry',
            'elev_threshold_m': elev_threshold_m,
            'buffer_deg':       buffer_deg,
            'fill_sigma_deg':   fill_sigma_deg,
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

    # Required inputs
    parser.add_argument('--boundaries', required=True, metavar='GEOJSON',
                        help='Plate boundary GeoJSON (boundaries_NNNN.geojson)')
    parser.add_argument('--topo', required=True, metavar='NETCDF',
                        help='Synthetic topo NetCDF')
    parser.add_argument('--bathy', required=True, metavar='NETCDF',
                        help='Reference bathy NetCDF')
    parser.add_argument('--output', '-o', required=True, metavar='NETCDF',
                        help='Output NetCDF')

    # Blend parameters
    blend = parser.add_argument_group('blend parameters')
    blend.add_argument('--threshold', type=float, default=DEFAULT_ELEV_THRESHOLD_M,
                       help='Elevation threshold for land mask in metres (default: %(default)s)')
    blend.add_argument('--buffer', type=float, default=DEFAULT_BUFFER_DEG,
                       help='Transition band width in degrees (default: %(default)s)')
    blend.add_argument('--fill-sigma', type=float, default=DEFAULT_FILL_SIGMA_DEG,
                       help='Gaussian fill radius in degrees (default: %(default)s)')
    blend.add_argument('--out-spacing', type=float, default=DEFAULT_OUTPUT_SPACING,
                       help='Output grid spacing in degrees; omit to use topo native resolution')

    # Trench parameters
    trench = parser.add_argument_group('trench parameters')
    trench.add_argument('--trench-spacing', type=float, default=None,
                        help='Build the trench at this resolution (degrees) and '
                             'interpolate; omit to evaluate it on the output grid')
    trench.add_argument('--peak-depth', type=float, default=DEFAULT_PEAK_DEPTH_M,
                        help='Trench depth at subduction line in metres (default: %(default)s)')
    trench.add_argument('--sigma', type=float, default=DEFAULT_SIGMA_DEG,
                        help='Trench Gaussian half-width in degrees (default: %(default)s)')
    trench.add_argument('--max-dist', type=float, default=DEFAULT_MAX_DIST_DEG,
                        help='Max distance from trench to compute, degrees (default: %(default)s)')
    trench.add_argument('--div-threshold', type=float, default=DEFAULT_DIV_THRESHOLD,
                        help='Fallback divergence threshold for subduction, '
                             'km/Myr = mm/yr (default: %(default)s)')

    # Seamount parameters
    sm = parser.add_argument_group('seamount parameters')
    sm.add_argument('--seamounts', default=DEFAULT_SEAMOUNTS_FILE, metavar='GEOJSON',
                    help='Seamounts GeoJSON; omit to skip seamount generation')
    sm.add_argument('--sm-fraction', type=float, default=DEFAULT_SM_FRACTION,
                    help='Fraction of seamounts to include, 0–1 (default: %(default)s)')
    sm.add_argument('--sm-jitter', type=float, default=DEFAULT_SM_JITTER_DEG,
                    help='Random coordinate jitter in degrees (default: %(default)s)')
    sm.add_argument('--sm-diameter', type=float, default=DEFAULT_SM_DIAMETER_KM,
                    help='Seamount base diameter in km (default: %(default)s)')
    sm.add_argument('--sm-height', type=float, default=DEFAULT_SM_HEIGHT_M,
                    help='Seamount peak height in metres (default: %(default)s)')
    sm.add_argument('--sm-seed', type=int, default=None,
                    help='Random seed for seamount selection and jitter')

    args = parser.parse_args()

    blend_topobathy(
        topo_file            = args.topo,
        bathy_file           = args.bathy,
        boundaries_file      = args.boundaries,
        output_file          = args.output,
        elev_threshold_m     = args.threshold,
        buffer_deg           = args.buffer,
        fill_sigma_deg       = args.fill_sigma,
        output_spacing       = args.out_spacing,
        trench_spacing       = args.trench_spacing,
        peak_depth_m         = args.peak_depth,
        sigma_deg            = args.sigma,
        max_dist_deg         = args.max_dist,
        div_threshold        = args.div_threshold,
        seamounts_file       = args.seamounts,
        seamount_fraction    = args.sm_fraction,
        seamount_jitter_deg  = args.sm_jitter,
        seamount_diameter_km = args.sm_diameter,
        seamount_height_m    = args.sm_height,
        seamount_seed        = args.sm_seed,
    )


if __name__ == '__main__':
    main()
