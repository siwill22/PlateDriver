#!/usr/bin/env python3
"""
generate_ocean_currents.py

Estimate steady-state barotropic surface ocean currents from a topography/
bathymetry grid using the Stommel (1948) wind-driven circulation model.

The model solves a 2-D linear PDE for the barotropic streamfunction ψ:

    R · ∇²ψ  +  β(φ) · ∂ψ/∂x  =  curl(τ) / (ρ · H₀)

where
    R     – linear bottom-drag coefficient [s⁻¹]  (default 1 × 10⁻⁶)
    ∇²    – Laplacian on the sphere:  ∂²/∂x² + ∂²/∂y² − (tan φ / R_Earth) ∂/∂y
    β(φ)  – meridional gradient of Coriolis parameter, 2Ω cos(φ)/R_Earth [m⁻¹ s⁻¹]
    τ     – idealised zonal wind stress:  τ_x(φ) = −τ₀ cos(π φ / φ_max)
    ρ     – seawater density [kg m⁻³]
    H₀    – reference ocean depth [m]

Boundary condition:  ψ = 0 on all land pixels, poleward of --pole-mask, and
at the edges of a non-global grid.  Longitude is periodic for global grids.
Surface velocities:  u = −∂ψ/∂y,  v = ∂ψ/∂x  [m s⁻¹].

With the default φ_max = 60° the wind profile gives easterly trades from the
equator to 30° (strongest at the equator), westerlies from 30° to the pole
(strongest at 60°, fading to zero at 90°); there are no polar easterlies.
The β-term produces clockwise (NH) / anticlockwise (SH) subtropical gyres with
western-boundary-current intensification (Gulf Stream-like).

Numerics: the β·∂ψ/∂x term is discretised with a first-order upwind
(eastward one-sided) difference, because the solution is built up from the
east and the boundary layer sits on the western side.  This keeps ψ free of
grid-scale oscillations even when the Stommel boundary-layer width R/β is
smaller than the grid spacing (the layer is then about one cell wide).

Usage
-----
    python generate_ocean_currents.py --topo output/topo_bathy_0060.nc
    python generate_ocean_currents.py \\
        --topo   output/topo_bathy_0060.nc \\
        --output output/ocean_currents_0060.nc \\
        --plot

All functions can also be imported directly.
"""

import argparse
import os
import sys

import numpy as np
import xarray as xr
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

R_EARTH = 6.371e6   # m
OMEGA   = 7.292e-5  # rad s⁻¹
RHO     = 1025.0    # kg m⁻³

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_TAU0      = 0.1      # Pa   – peak wind-stress magnitude
DEFAULT_R_DRAG    = 1e-6     # s⁻¹  – linear bottom drag (Stommel 1948); boundary
                             #        layer width R/β ≈ 50 km at 35°
DEFAULT_PHI_MAX   = 60.0     # deg  – latitude of peak westerlies; trades turn to
                             #        westerlies at φ_max/2 (30°)
DEFAULT_H0        = 3500.0   # m    – reference ocean depth
DEFAULT_POLE_MASK = 85.0     # deg  – mask ocean poleward of this latitude
MAX_POLE_MASK     = 89.0     # deg  – upper limit for pole mask (cos φ → 0 guard)

# Abyssal (Stommel-Arons) model
DEFAULT_POLE_LAT  = 55.0     # deg  – deep-water formation latitude threshold
DEFAULT_W_DEEP    = 1e-7     # m/s  – deep-water formation velocity


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


def _load_ocean_grid(topo_file, pole_mask_deg):
    """Load a topo grid and build the ocean mask used by the Stommel solver.

    - Both axes are sorted ascending, so β and ∂/∂x have the right sign even
      for grids stored north-to-south or east-to-west.
    - A duplicated seam column (gridline registration: both -180 and +180) is
      dropped, so the periodic stencil does not treat one meridian as two
      separate nodes.  Use _restore_seam() to put it back on the outputs.
    - Ocean poleward of pole_mask_deg (capped at MAX_POLE_MASK) and the first
      and last latitude rows are treated as land (ψ = 0).
    - For a non-global longitude range the first and last columns are also
      treated as land, instead of wrapping one edge onto the other.

    Returns
    -------
    lats, lons : 1-D arrays (ascending; lons without the seam duplicate)
    ocean_mask : 2-D bool array (nlat, nlon)
    out_lons   : longitudes of the output grid (with the seam column, if any)
    has_seam   : True if the seam column was dropped
    """
    print(f'Loading {topo_file} ...', file=sys.stderr)
    topo = _norm_coords(xr.open_dataarray(topo_file)).sortby('lat').sortby('lon')
    lats = topo.lat.values.astype(float)
    out_lons = topo.lon.values.astype(float)
    dlat = float(lats[1] - lats[0])
    dlon = float(out_lons[1] - out_lons[0])
    print(f'  Grid  {len(lats)} × {len(out_lons)}  '
          f'(Δlat = {dlat:.3f}°, Δlon = {dlon:.3f}°)', file=sys.stderr)

    span = float(out_lons[-1] - out_lons[0])
    has_seam  = abs(span - 360.0) < 1e-3 * dlon
    is_global = has_seam or abs(span + dlon - 360.0) < 1e-3 * dlon
    values = topo.values
    lons   = out_lons
    if has_seam:
        values = values[:, :-1]
        lons   = out_lons[:-1]

    if pole_mask_deg > MAX_POLE_MASK:
        print(f'  Warning: --pole-mask {pole_mask_deg}° is too close to the pole '
              f'(cos φ → 0 makes the system ill-conditioned); using '
              f'{MAX_POLE_MASK}°.', file=sys.stderr)
        pole_mask_deg = MAX_POLE_MASK

    # Ocean mask — exclude high-latitude cells to avoid cos φ → 0 singularity
    ocean_mask = values < 0
    ocean_mask[np.abs(lats) > pole_mask_deg, :] = False
    ocean_mask[[0, -1], :] = False                 # domain edges: ψ = 0
    if not is_global:
        print('  Longitude range is not global: edge columns set to ψ = 0.',
              file=sys.stderr)
        ocean_mask[:, [0, -1]] = False

    n_ocean = int(ocean_mask.sum())
    print(f'  Ocean pixels: {n_ocean:,} / {ocean_mask.size:,}  '
          f'({100*n_ocean/ocean_mask.size:.1f}%)', file=sys.stderr)
    return lats, lons, ocean_mask, out_lons, has_seam


def _restore_seam(arr, has_seam):
    """Re-append the duplicated ±180 column (a copy of the first column)."""
    return np.concatenate([arr, arr[:, :1]], axis=1) if has_seam else arr


def _check_output_path(output_file, topo_file):
    """Refuse to overwrite the input topo grid."""
    if os.path.realpath(output_file) == os.path.realpath(topo_file):
        raise ValueError(f'Output file {output_file!r} is the same as the input '
                         'topo file; refusing to overwrite it.')


def _default_output(topo_file, suffix):
    """topo_bathy_0060.nc → topo_bathy_0060<suffix> (same directory)."""
    return os.path.splitext(topo_file)[0] + suffix


# ---------------------------------------------------------------------------
# Core physics
# ---------------------------------------------------------------------------

def _compute_wind_stress_curl(lats, tau0=DEFAULT_TAU0,
                               phi_max_deg=DEFAULT_PHI_MAX):
    """Compute the vertical component of curl(τ) [Pa m⁻¹] for a zonal-only wind.

    Wind stress:  τ_x(φ) = −τ₀ · cos(π · φ / φ_max),  τ_y = 0

    i.e. easterly (τ_x < 0) trades between the equator and φ_max/2, and
    westerlies (τ_x > 0) peaking at φ_max.  With φ_max = 60° the zero
    crossings are at ±30° and ±90°.

    Spherical curl:
        curl(τ) = −(1 / (R·cos φ)) · d(cos φ · τ_x) / dφ

    Parameters
    ----------
    lats : 1-D array  [degrees]
    tau0 : float  [Pa]
    phi_max_deg : float  [degrees]

    Returns
    -------
    curl_tau : 1-D ndarray  (nlat,)  [Pa m⁻¹]
    """
    lats_rad    = np.radians(lats)
    phi_max_rad = np.radians(phi_max_deg)

    tau_x  = -tau0 * np.cos(np.pi * lats_rad / phi_max_rad)  # < 0 = easterly
    f_phi  = np.cos(lats_rad) * tau_x                       # cos(φ)·τ_x
    df_dphi = np.gradient(f_phi, lats_rad)                  # numerical d/dφ
    cos_phi = np.clip(np.cos(lats_rad), 1e-6, None)
    return -df_dphi / (R_EARTH * cos_phi)


def _build_stommel_matrix(ocean_mask, lats, lons,
                           R_drag=DEFAULT_R_DRAG,
                           forcing_2d=None):
    """Build the sparse linear system A ψ = b for the Stommel equation.

    At each ocean pixel (i, j):
        R · ∇²ψ  +  β(φ) · ∂ψ/∂x  =  forcing_2d[i, j]

    where forcing_2d already incorporates all physical constants (ρ, H₀, etc.)
    and has units [s⁻²].  Land pixels receive identity rows (ψ = 0 Dirichlet BC).

    Discretisation (x = eastward, y = northward distance in metres):
      ∇²ψ     = ψ_xx + ψ_yy − (tan φ / R_Earth) ψ_y   (spherical Laplacian;
                centred differences, including the tan φ metric term)
      β ψ_x   ≈ β (ψ[i, j+1] − ψ[i, j]) / dx         (first-order upwind: the
                interior solution is set from the east, the boundary layer is
                on the west, so a one-sided eastward difference gives a
                non-oscillating ψ at any grid spacing)

    Longitude is treated as periodic, so lons must NOT contain a duplicated
    seam column (see _load_ocean_grid).  The first and last latitude rows
    must be masked as land (done in _load_ocean_grid).

    Parameters
    ----------
    ocean_mask : 2-D bool array  (nlat, nlon)
    lats, lons : 1-D arrays  [degrees], ascending
    R_drag : float  [s⁻¹]
    forcing_2d : 2-D array  (nlat, nlon)  [s⁻²]  or None (→ zero RHS)

    Returns
    -------
    A : scipy CSR sparse matrix  (N × N)
    b : 1-D ndarray  (N,)
    """
    nlat, nlon = ocean_mask.shape
    N = nlat * nlon

    dlat_rad = np.radians(abs(float(lats[1] - lats[0])))
    dlon_rad = np.radians(abs(float(lons[1] - lons[0])))
    lats_rad = np.radians(lats)

    dy      = R_EARTH * dlat_rad                             # m, constant
    cos_phi = np.clip(np.cos(lats_rad), 1e-6, None)         # (nlat,)
    tan_phi = np.sin(lats_rad) / cos_phi                    # (nlat,)
    dx_row  = R_EARTH * cos_phi * dlon_rad                  # (nlat,) – m
    beta    = 2.0 * OMEGA * cos_phi / R_EARTH               # (nlat,) – m⁻¹s⁻¹

    # Flat index → (lat, lon) index arrays for *ocean* pixels only
    all_flat   = np.arange(N, dtype=np.int64)
    ocean_flat = all_flat[ocean_mask.ravel()]
    land_flat  = all_flat[~ocean_mask.ravel()]

    oi = ocean_flat // nlon   # lat index per ocean pixel
    oj = ocean_flat % nlon    # lon index per ocean pixel

    dx_oi  = dx_row[oi]         # per-pixel dx (m)
    beta_oi = beta[oi]          # per-pixel β
    tan_oi = tan_phi[oi]        # per-pixel tan φ
    dy2    = dy ** 2
    dx2_oi = dx_oi ** 2

    # Metric term −(tan φ / R_Earth)·∂ψ/∂y, centred, multiplied by R
    metric_oi = R_drag * tan_oi / (2.0 * R_EARTH * dy)

    # -----------------------------------------------------------------------
    # Stencil contributions (vectorised COO format)
    # -----------------------------------------------------------------------

    # Diagonal: R · (−2/dy² − 2/dx²)  −  β/dx   (upwind β term)
    r_c = ocean_flat
    c_c = ocean_flat
    d_c = R_drag * (-2.0 / dy2 - 2.0 / dx2_oi) - beta_oi / dx_oi

    # North neighbor (i+1, j) — clamp at domain edge (edge rows are land)
    i_n = np.minimum(oi + 1, nlat - 1)
    r_n = ocean_flat;  c_n = i_n * nlon + oj
    d_n = R_drag / dy2 - metric_oi

    # South neighbor (i−1, j)
    i_s = np.maximum(oi - 1, 0)
    r_s = ocean_flat;  c_s = i_s * nlon + oj
    d_s = R_drag / dy2 + metric_oi

    # East neighbor (i, j+1) — periodic longitude; carries the upwind β term
    j_e = (oj + 1) % nlon
    r_e = ocean_flat;  c_e = oi * nlon + j_e
    d_e = R_drag / dx2_oi + beta_oi / dx_oi

    # West neighbor (i, j−1) — periodic longitude
    j_w = (oj - 1) % nlon
    r_w = ocean_flat;  c_w = oi * nlon + j_w
    d_w = R_drag / dx2_oi

    # Land identity rows
    r_land = land_flat;  c_land = land_flat
    d_land = np.ones(len(land_flat))

    # Assemble into a single COO matrix and convert to CSR
    rows = np.concatenate([r_c, r_n, r_s, r_e, r_w, r_land])
    cols = np.concatenate([c_c, c_n, c_s, c_e, c_w, c_land])
    data = np.concatenate([d_c, d_n, d_s, d_e, d_w, d_land])

    A = coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

    # RHS: forcing_2d values at ocean pixels, 0 at land
    b = np.zeros(N)
    if forcing_2d is not None:
        b[ocean_flat] = forcing_2d.ravel()[ocean_flat]

    return A, b


def _solve_streamfunction(A, b, nlat, nlon):
    """Solve A ψ = b with a sparse direct solver; return ψ as (nlat, nlon)."""
    print(f'  Solving {A.shape[0]:,}-DOF system  '
          f'({A.nnz:,} non-zeros) ...', file=sys.stderr)
    psi_flat = spsolve(A, b)
    return psi_flat.reshape(nlat, nlon)


def _recover_velocities(psi_2d, lats, lons, ocean_mask):
    """Recover surface velocities from ψ.

        u = −∂ψ/∂y   (eastward)  [m s⁻¹]
        v =  ∂ψ/∂x   (northward) [m s⁻¹]

    Returns u, v as (nlat, nlon) arrays, NaN on land.
    """
    lats_rad = np.radians(lats)
    dlon_rad = np.radians(abs(float(lons[1] - lons[0])))

    # u = −∂ψ/∂y  using physical y-coordinate (metres)
    lat_m = lats_rad * R_EARTH
    u = -np.gradient(psi_2d, lat_m, axis=0)

    # v = ∂ψ/∂x  with periodic longitude and cos-lat metric
    cos_lat = np.clip(np.cos(lats_rad), 1e-6, None)
    dx_2d   = R_EARTH * cos_lat[:, None] * dlon_rad    # (nlat, nlon)
    psi_e   = np.roll(psi_2d, -1, axis=1)              # ψ[i, j+1]
    psi_w   = np.roll(psi_2d,  1, axis=1)              # ψ[i, j-1]
    v = (psi_e - psi_w) / (2.0 * dx_2d)

    u[~ocean_mask] = np.nan
    v[~ocean_mask] = np.nan
    return u, v


def _make_diagnostic_plot(psi_2d, u, v, ocean_mask, lats, lons, out_nc_path):
    """Save a two-panel diagnostic PNG (streamfunction + speed/vectors)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    speed = np.sqrt(u**2 + v**2)
    land_rgba = np.where(ocean_mask[:, :, None],
                         np.zeros((*ocean_mask.shape, 4)),
                         np.array([[[0.78, 0.66, 0.48, 1.0]]]))

    _, axes = plt.subplots(2, 1, figsize=(16, 11))

    ax = axes[0]
    psi_ocean = np.where(ocean_mask, psi_2d, np.nan)
    vmax = np.nanpercentile(np.abs(psi_ocean), 99)
    cf = ax.contourf(lons, lats, psi_ocean, levels=30,
                     cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    plt.colorbar(cf, ax=ax, label='ψ  [m² s⁻¹]', shrink=0.8)
    ax.contour(lons, lats, psi_ocean, levels=20,
               colors='k', linewidths=0.35, alpha=0.45)
    ax.imshow(land_rgba, origin='lower', aspect='auto',
              extent=[lons[0], lons[-1], lats[0], lats[-1]])
    ax.set_title('Barotropic streamfunction  ψ  (Stommel 1948 model)',
                 fontsize=12)
    ax.set_xlabel('Longitude');  ax.set_ylabel('Latitude')

    ax2 = axes[1]
    cf2 = ax2.contourf(lons, lats, np.where(ocean_mask, speed, np.nan),
                       levels=20, cmap='viridis')
    plt.colorbar(cf2, ax=ax2, label='speed  [m s⁻¹]', shrink=0.8)
    stride = max(1, len(lats) // 36)
    ax2.quiver(lons[::stride], lats[::stride],
               u[::stride, ::stride], v[::stride, ::stride],
               scale=3.0, scale_units='width', width=0.0018, alpha=0.75)
    ax2.imshow(land_rgba, origin='lower', aspect='auto',
               extent=[lons[0], lons[-1], lats[0], lats[-1]])
    ax2.set_title('Surface current speed  +  velocity vectors')
    ax2.set_xlabel('Longitude');  ax2.set_ylabel('Latitude')

    plt.tight_layout()
    out_png = os.path.splitext(out_nc_path)[0] + '_preview.png'
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Diagnostic plot → {out_png}', file=sys.stderr)


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def generate_ocean_currents(topo_file,
                             output_file=None,
                             tau0=DEFAULT_TAU0,
                             R_drag=DEFAULT_R_DRAG,
                             phi_max_deg=DEFAULT_PHI_MAX,
                             H0=DEFAULT_H0,
                             pole_mask_deg=DEFAULT_POLE_MASK,
                             plot=False):
    """Estimate steady barotropic surface ocean currents from bathymetry.

    Reads a topo_bathy NetCDF, solves the Stommel (1948) wind-driven
    streamfunction equation, and writes u / v / psi to an output NetCDF.

    Parameters
    ----------
    topo_file : str
        Path to a topo_bathy_NNNN.nc from run_SDP.py.
        Convention: values > 0 = land (m elevation), < 0 = ocean (m depth).
    output_file : str or None
        Output NetCDF path.  Defaults to topo_file with suffix
        _surface_currents.nc.
    tau0 : float
        Peak wind-stress magnitude [Pa] (default 0.1).
    R_drag : float
        Linear bottom-drag coefficient [s⁻¹] (default 1 × 10⁻⁶).
    phi_max_deg : float
        Latitude of peak westerly wind stress [°] (default 60); the trades
        change to westerlies at φ_max/2.
    H0 : float
        Reference ocean depth used in RHS normalisation [m] (default 3500).
    pole_mask_deg : float
        Ocean pixels poleward of this latitude are masked (cos φ → 0 guard;
        capped at 89°).
    plot : bool
        Save a diagnostic PNG alongside the NetCDF.

    Returns
    -------
    xr.Dataset  with variables psi, u, v (ascending lat/lon).
    """
    if output_file is None:
        output_file = _default_output(topo_file, '_surface_currents.nc')
    _check_output_path(output_file, topo_file)

    lats, lons, ocean_mask, out_lons, has_seam = _load_ocean_grid(
        topo_file, pole_mask_deg)
    nlat, nlon = ocean_mask.shape

    # Wind-stress curl → 2-D forcing array  [s⁻²]
    print('Computing wind stress curl ...', file=sys.stderr)
    curl_1d = _compute_wind_stress_curl(lats, tau0=tau0, phi_max_deg=phi_max_deg)
    forcing_2d = np.broadcast_to(
        (curl_1d / (RHO * H0))[:, None], (nlat, nlon)
    ).copy()

    # Sparse linear system
    print('Building Stommel sparse matrix ...', file=sys.stderr)
    A, b = _build_stommel_matrix(ocean_mask, lats, lons,
                                  R_drag=R_drag, forcing_2d=forcing_2d)
    psi_2d = _solve_streamfunction(A, b, nlat, nlon)
    psi_2d[~ocean_mask] = 0.0

    # Surface velocities from streamfunction
    u, v = _recover_velocities(psi_2d, lats, lons, ocean_mask)

    psi_ocean = psi_2d[ocean_mask]
    speed_max = float(np.nanmax(np.sqrt(u**2 + v**2)))
    print(f'  ψ range (ocean): {psi_ocean.min():.2f} – {psi_ocean.max():.2f} m² s⁻¹',
          file=sys.stderr)
    print(f'  Max surface speed: {speed_max:.4f} m s⁻¹', file=sys.stderr)

    # Put back the duplicated ±180 column (if the input had one)
    psi_2d, u, v, ocean_mask = (_restore_seam(a, has_seam)
                                for a in (psi_2d, u, v, ocean_mask))
    lons = out_lons

    # Output Dataset
    coords = {'lat': lats, 'lon': lons}
    ds = xr.Dataset(
        {
            'psi': xr.DataArray(psi_2d, dims=['lat', 'lon'], coords=coords,
                                attrs={'units':     'm2 s-1',
                                       'long_name': 'barotropic streamfunction'}),
            'u':   xr.DataArray(u,      dims=['lat', 'lon'], coords=coords,
                                attrs={'units':     'm s-1',
                                       'long_name': 'eastward surface current',
                                       'note':      'NaN on land'}),
            'v':   xr.DataArray(v,      dims=['lat', 'lon'], coords=coords,
                                attrs={'units':     'm s-1',
                                       'long_name': 'northward surface current',
                                       'note':      'NaN on land'}),
        },
        attrs={
            'model':       'Stommel 1948 barotropic wind-driven circulation',
            'tau0_Pa':     tau0,
            'R_drag_s-1':  R_drag,
            'phi_max_deg': phi_max_deg,
            'H0_m':        H0,
        })

    ds.to_netcdf(output_file)
    print(f'Saved → {output_file}', file=sys.stderr)

    if plot:
        _make_diagnostic_plot(psi_2d, u, v, ocean_mask, lats, lons, output_file)

    return ds


def _compute_abyssal_forcing(ocean_mask, lats, lons,
                              pole_lat=DEFAULT_POLE_LAT,
                              w_deep=DEFAULT_W_DEEP):
    """Build the 2-D upwelling/downwelling field for the Stommel-Arons model.

    Deep water forms in polar ocean pixels (|φ| > pole_lat) and upwells
    uniformly across the remaining ocean, conserving volume.

    Parameters
    ----------
    ocean_mask : 2-D bool array  (nlat, nlon)
    lats, lons : 1-D arrays  [degrees]
    pole_lat : float  [degrees]  – formation threshold
    w_deep : float  [m s⁻¹]  – downwelling speed in formation zones

    Returns
    -------
    w_E : 2-D ndarray  (nlat, nlon)  [m s⁻¹]
        Positive = upwelling, negative = downwelling (formation).
    """
    nlat, nlon = ocean_mask.shape
    lat2d = np.broadcast_to(lats[:, None], (nlat, nlon))

    polar    = ocean_mask & (np.abs(lat2d) > pole_lat)
    interior = ocean_mask & (np.abs(lat2d) <= pole_lat)

    n_polar    = int(polar.sum())
    n_interior = int(interior.sum())

    if n_polar == 0:
        print(f'  Warning: no polar ocean found (|lat| > {pole_lat}°) — '
              'abyssal forcing is zero.', file=sys.stderr)
        return np.zeros((nlat, nlon))

    # Area-weighted conservation: ∑ w_E · cos(φ) · dA = 0
    cos_lat2d   = np.cos(np.radians(lats))[:, None] * np.ones(nlon)
    area_polar   = float((polar    * cos_lat2d).sum())
    area_interior = float((interior * cos_lat2d).sum())
    w_up = w_deep * area_polar / area_interior if area_interior > 0 else 0.0

    w_E = np.zeros((nlat, nlon))
    w_E[polar]    = -w_deep
    w_E[interior] =  w_up

    print(f'  Polar formation pixels : {n_polar:,}  (|lat| > {pole_lat}°)',
          file=sys.stderr)
    print(f'  Interior upwelling pixels: {n_interior:,}', file=sys.stderr)
    print(f'  w_deep = {w_deep:.1e} m/s  →  w_upwelling = {w_up:.2e} m/s',
          file=sys.stderr)
    return w_E


def generate_deep_currents(topo_file,
                            output_file=None,
                            R_drag=DEFAULT_R_DRAG,
                            pole_lat=DEFAULT_POLE_LAT,
                            w_deep=DEFAULT_W_DEEP,
                            H0=DEFAULT_H0,
                            pole_mask_deg=DEFAULT_POLE_MASK,
                            plot=False):
    """Estimate abyssal ocean circulation using the Stommel-Arons (1960) model.

    The Stommel-Arons model treats deep-water formation as a localised source
    in polar regions with uniform upwelling elsewhere.  The same Stommel PDE
    is solved but with thermohaline forcing on the RHS:

        R · ∇²ψ  +  β(φ) · ∂ψ/∂x  =  f(φ) · w_E(x,y) / H₀

    where w_E < 0 in formation zones (|φ| > pole_lat) and w_E > 0 (upwelling)
    elsewhere.  The result captures the broad pattern of the abyssal conveyor:
    equatorward deep western boundary currents and sluggish interior return flow.

    Parameters
    ----------
    topo_file : str
        Path to a topo_bathy_NNNN.nc from run_SDP.py.
    output_file : str or None
        Output NetCDF path.  Defaults to topo_file with suffix _deep_currents.nc.
    R_drag : float  [s⁻¹]
        Linear bottom-drag coefficient (default 1 × 10⁻⁶).
    pole_lat : float  [degrees]
        Latitude poleward of which deep water forms (default 55).
    w_deep : float  [m s⁻¹]
        Downwelling speed in formation zones (default 1 × 10⁻⁷).
    H0 : float  [m]
        Reference abyssal depth (default 3500).
    pole_mask_deg : float  [degrees]
        Ocean pixels beyond this latitude are masked (cos φ guard, default 85,
        capped at 89).
    plot : bool
        Save a diagnostic PNG alongside the NetCDF.

    Returns
    -------
    xr.Dataset  with variables psi, u, v (ascending lat/lon).
    """
    if output_file is None:
        output_file = _default_output(topo_file, '_deep_currents.nc')
    _check_output_path(output_file, topo_file)

    print('[Stommel-Arons abyssal model]', file=sys.stderr)
    lats, lons, ocean_mask, out_lons, has_seam = _load_ocean_grid(
        topo_file, pole_mask_deg)
    nlat, nlon = ocean_mask.shape

    # Abyssal forcing:  f(φ) · w_E(x,y) / H₀  [s⁻²]
    print('Computing abyssal forcing ...', file=sys.stderr)
    w_E    = _compute_abyssal_forcing(ocean_mask, lats, lons, pole_lat, w_deep)
    f_2d   = 2.0 * OMEGA * np.sin(np.radians(lats))[:, None] * np.ones((1, nlon))
    forcing_2d = f_2d * w_E / H0

    # Sparse linear system
    print('Building Stommel-Arons sparse matrix ...', file=sys.stderr)
    A, b = _build_stommel_matrix(ocean_mask, lats, lons,
                                  R_drag=R_drag, forcing_2d=forcing_2d)
    psi_2d = _solve_streamfunction(A, b, nlat, nlon)
    psi_2d[~ocean_mask] = 0.0

    u, v = _recover_velocities(psi_2d, lats, lons, ocean_mask)

    psi_ocean = psi_2d[ocean_mask]
    speed_max = float(np.nanmax(np.sqrt(u**2 + v**2)))
    print(f'  ψ range (ocean): {psi_ocean.min():.2f} – {psi_ocean.max():.2f} m² s⁻¹',
          file=sys.stderr)
    print(f'  Max abyssal speed: {speed_max:.4f} m s⁻¹', file=sys.stderr)

    # Put back the duplicated ±180 column (if the input had one)
    psi_2d, u, v, ocean_mask = (_restore_seam(a, has_seam)
                                for a in (psi_2d, u, v, ocean_mask))
    lons = out_lons

    coords = {'lat': lats, 'lon': lons}
    ds = xr.Dataset(
        {
            'psi': xr.DataArray(psi_2d, dims=['lat', 'lon'], coords=coords,
                                attrs={'units':     'm2 s-1',
                                       'long_name': 'abyssal streamfunction'}),
            'u':   xr.DataArray(u,      dims=['lat', 'lon'], coords=coords,
                                attrs={'units':     'm s-1',
                                       'long_name': 'eastward abyssal current',
                                       'note':      'NaN on land'}),
            'v':   xr.DataArray(v,      dims=['lat', 'lon'], coords=coords,
                                attrs={'units':     'm s-1',
                                       'long_name': 'northward abyssal current',
                                       'note':      'NaN on land'}),
        },
        attrs={
            'model':        'Stommel-Arons 1960 abyssal thermohaline circulation',
            'R_drag_s-1':   R_drag,
            'pole_lat_deg': pole_lat,
            'w_deep_m_s-1': w_deep,
            'H0_m':         H0,
        })

    ds.to_netcdf(output_file)
    print(f'Saved → {output_file}', file=sys.stderr)

    if plot:
        _make_diagnostic_plot(psi_2d, u, v, ocean_mask, lats, lons, output_file)

    return ds


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--topo', required=True, metavar='NC',
                        help='topo_bathy NetCDF (e.g. output/topo_bathy_0060.nc)')
    parser.add_argument('--output', '-o', default=None, metavar='NC',
                        help='Output NetCDF (default: auto-named beside --topo)')
    parser.add_argument('--R-drag', type=float, default=DEFAULT_R_DRAG,
                        dest='R_drag',
                        help='Linear drag [s⁻¹] (default: %(default)s)')
    parser.add_argument('--H0', type=float, default=DEFAULT_H0,
                        help='Reference ocean depth [m] (default: %(default)s)')
    parser.add_argument('--pole-mask', type=float, default=DEFAULT_POLE_MASK,
                        dest='pole_mask_deg',
                        help='Mask ocean poleward of this latitude [°] '
                             '(default: %(default)s)')
    parser.add_argument('--plot', action='store_true',
                        help='Save a diagnostic PNG preview')

    # Surface-model options
    surf = parser.add_argument_group('surface (wind-driven) model')
    surf.add_argument('--tau0', type=float, default=DEFAULT_TAU0,
                      help='Peak wind-stress [Pa] (default: %(default)s)')
    surf.add_argument('--phi-max', type=float, default=DEFAULT_PHI_MAX,
                      dest='phi_max_deg',
                      help='Latitude of peak westerlies [°]; trades change to '
                           'westerlies at half this latitude (default: %(default)s)')

    # Deep / abyssal model
    deep = parser.add_argument_group('deep (abyssal thermohaline) model')
    deep.add_argument('--deep', action='store_true',
                      help='Run the Stommel-Arons abyssal model instead of '
                           'the surface wind-driven model')
    deep.add_argument('--pole-lat', type=float, default=DEFAULT_POLE_LAT,
                      dest='pole_lat',
                      help='Deep-water formation latitude threshold [°] '
                           '(default: %(default)s)')
    deep.add_argument('--w-deep', type=float, default=DEFAULT_W_DEEP,
                      dest='w_deep',
                      help='Deep-water formation velocity [m s⁻¹] '
                           '(default: %(default)s)')

    args = parser.parse_args()

    if args.deep:
        generate_deep_currents(
            topo_file     = args.topo,
            output_file   = args.output,
            R_drag        = args.R_drag,
            pole_lat      = args.pole_lat,
            w_deep        = args.w_deep,
            H0            = args.H0,
            pole_mask_deg = args.pole_mask_deg,
            plot          = args.plot,
        )
    else:
        generate_ocean_currents(
            topo_file     = args.topo,
            output_file   = args.output,
            tau0          = args.tau0,
            R_drag        = args.R_drag,
            phi_max_deg   = args.phi_max_deg,
            H0            = args.H0,
            pole_mask_deg = args.pole_mask_deg,
            plot          = args.plot,
        )


if __name__ == '__main__':
    main()
