import time

import pygplates
from gprm.utils.create_gpml import gpml2gdf, gdf2gpml
from gprm.utils.proximity import rasterize_polygons
import plate_driver as pdriver
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from shapely.ops import unary_union
from shapely.geometry import box, LineString
import pygmt

import topojson as tp
try:
    from xrspatial.polygonize import polygonize                 # xarray-spatial < 0.4
except ImportError:
    from xrspatial.experimental.polygonize import polygonize    # xarray-spatial >= 0.4
from shapely.geometry import Point, Polygon
from sklearn.neighbors import BallTree


# Velocities throughout this module come from pygplates.calculate_velocities (or the
# equivalent numpy code below) and are in km/Myr, which is numerically equal to mm/yr.
EARTH_RADIUS_KM = pygplates.Earth.mean_radius_in_kms


# ---------------------------------------------------------------------------
# Small spherical helpers shared by the functions below
# ---------------------------------------------------------------------------

def lonlat_to_xyz(lon, lat):
    """Unit vectors (..., 3) for lon/lat in degrees."""
    lon_r, lat_r = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(lat_r) * np.cos(lon_r),
                     np.cos(lat_r) * np.sin(lon_r),
                     np.sin(lat_r)], axis=-1)


def xyz_to_lonlat(xyz):
    """lon, lat in degrees for unit vectors (..., 3)."""
    xyz = xyz / np.linalg.norm(xyz, axis=-1, keepdims=True)
    lat = np.degrees(np.arcsin(np.clip(xyz[..., 2], -1.0, 1.0)))
    lon = np.degrees(np.arctan2(xyz[..., 1], xyz[..., 0]))
    return lon, lat


class SphereTree:
    """Nearest-neighbour search on the sphere, a faster stand-in for
    sklearn BallTree(metric='haversine') with the same input and query() output.

    Points are given as (lat, lon) in radians, as for BallTree.  Searching uses a
    KD-tree on 3-D unit vectors: the straight-line (chord) distance c between two unit
    vectors is 2 sin(d/2) for great-circle distance d, so it orders neighbours exactly
    as d does, and query() converts c back to d.  About 7x faster than the haversine
    BallTree for tens of thousands of points.
    """

    def __init__(self, latlon_rad):
        from scipy.spatial import cKDTree
        latlon_rad = np.asarray(latlon_rad, dtype=float).reshape(-1, 2)
        self._tree = cKDTree(self._xyz(latlon_rad))

    @staticmethod
    def _xyz(latlon_rad):
        lat, lon = latlon_rad[:, 0], latlon_rad[:, 1]
        return np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])

    def query(self, latlon_rad, k=1, seed=0):
        """Return (great-circle distances in radians (n, k), indices (n, k)).

        With k = 1, exact ties (several points equally near, common on a regular grid)
        are resolved by a seeded hash of the query position and the candidate's index,
        so the answer is unbiased and does not depend on the search implementation or
        on the order of the queries.  The hash does not change between time steps, so
        a tied cell resolves the same way every step.
        """
        latlon_rad = np.asarray(latlon_rad, dtype=float).reshape(-1, 2)
        if k != 1:
            chord, idx = self._tree.query(self._xyz(latlon_rad), k=k)
            chord, idx = np.asarray(chord).reshape(len(latlon_rad), k), np.asarray(idx).reshape(len(latlon_rad), k)
            return 2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0)), idx
        dist, idx, tied = self.query_candidates(latlon_rad)
        choice = np.argmin(np.where(tied, tie_hash(latlon_rad, idx, seed), np.iinfo(np.uint64).max), axis=1)
        rows = np.arange(len(idx))
        return dist[:, None], idx[rows, choice][:, None]

    def query_candidates(self, latlon_rad, k=4, tol_rad=1e-12):
        """Nearest distance (n,), the k nearest indices (n, k) and a mask (n, k) of
        those whose distance equals the nearest (within tol_rad)."""
        latlon_rad = np.asarray(latlon_rad, dtype=float).reshape(-1, 2)
        k = min(k, self._tree.n)
        chord, idx = self._tree.query(self._xyz(latlon_rad), k=k)
        chord, idx = np.asarray(chord).reshape(len(latlon_rad), k), np.asarray(idx).reshape(len(latlon_rad), k)
        dist = 2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))
        return dist[:, 0], idx, dist <= dist[:, :1] + tol_rad


def tie_hash(latlon_rad, idx, seed=0):
    """Deterministic pseudo-random uint64 per (query position, candidate index) pair
    (splitmix64 mixing), used to choose between equally near candidates."""
    q = np.round(np.degrees(np.asarray(latlon_rad)) * 1e5).astype(np.int64)
    key = (q[:, 0].astype(np.uint64) * np.uint64(0x9E3779B97F4A7C15)
           ^ q[:, 1].astype(np.uint64) * np.uint64(0xC2B2AE3D27D4EB4F)
           ^ np.uint64(seed) * np.uint64(0x165667B19E3779F9))
    z = key[:, None] ^ (np.asarray(idx).astype(np.uint64) * np.uint64(0xD6E8FEB86659FD93))
    with np.errstate(over='ignore'):
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return z ^ (z >> np.uint64(31))


def rotation_axis_angle(rotation):
    """Return (unit axis xyz, angle in radians) of a pygplates FiniteRotation."""
    if rotation is None or rotation.represents_identity_rotation():
        return np.array([0.0, 0.0, 1.0]), 0.0
    lat, lon, angle_deg = rotation.get_lat_lon_euler_pole_and_angle_degrees()
    return lonlat_to_xyz(lon, lat), np.radians(angle_deg)


def rotate_lonlat(rotation, lon, lat):
    """Apply a FiniteRotation to arrays of lon/lat (degrees); Rodrigues' formula."""
    axis, angle = rotation_axis_angle(rotation)
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    if angle == 0.0:
        return lon.copy(), lat.copy()
    p = lonlat_to_xyz(lon, lat)
    rotated = (p * np.cos(angle)
               + np.cross(axis, p) * np.sin(angle)
               + np.outer(p @ axis, axis).reshape(p.shape) * (1.0 - np.cos(angle)))
    return xyz_to_lonlat(rotated)


def plate_velocity_ne(rotation, lon, lat, delta_time):
    """East and north velocity (km/Myr) of points moved by a stage rotation over delta_time.

    Numerically identical to pygplates.calculate_velocities + convert to north/east,
    but vectorised: v = omega x r, with omega = axis * angle / delta_time.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    axis, angle = rotation_axis_angle(rotation)
    if angle == 0.0:
        return np.zeros_like(lon), np.zeros_like(lat)
    omega = axis * angle / delta_time
    r = lonlat_to_xyz(lon, lat)
    v = np.cross(omega, r) * EARTH_RADIUS_KM
    lon_r, lat_r = np.radians(lon), np.radians(lat)
    east  = -np.sin(lon_r) * v[..., 0] + np.cos(lon_r) * v[..., 1]
    north = (-np.sin(lat_r) * np.cos(lon_r) * v[..., 0]
             - np.sin(lat_r) * np.sin(lon_r) * v[..., 1]
             + np.cos(lat_r) * v[..., 2])
    return east, north


def offset_along_normal(lon, lat, normal_xyz, dist_deg):
    """Move points a great-circle distance dist_deg in the direction of a unit vector
    normal_xyz that is perpendicular to the point (e.g. a great-circle arc normal)."""
    p = lonlat_to_xyz(lon, lat)
    n = np.asarray(normal_xyz, dtype=float)
    n = n - (p * n).sum(axis=-1, keepdims=True) * p      # make exactly tangent
    n = n / np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), 1e-12)
    d = np.radians(dist_deg)
    return xyz_to_lonlat(p * np.cos(d) + n * np.sin(d))


def raster_coords(da):
    """(lons, lats) coordinate arrays of a raster with x/y or lon/lat dims."""
    lons = (da.lon if 'lon' in da.coords else da.x).values
    lats = (da.lat if 'lat' in da.coords else da.y).values
    return lons, lats


def raster_nearest_index(lons, lats, lon, lat):
    """Row/column index of the nearest node of a regular global lon/lat grid.

    Longitude wraps around the globe, so points near ±180 find the right column
    (a duplicated ±180 seam column is handled).
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    dlon = lons[1] - lons[0]
    dlat = lats[1] - lats[0]
    nlon, nlat = len(lons), len(lats)
    ci = np.rint(((lon - lons[0]) % 360.0) / dlon).astype(int)
    if abs((lons[-1] - lons[0]) - 360.0) < 1e-6:        # duplicated seam column
        ci = ci % (nlon - 1)
    else:
        ci = ci % int(round(360.0 / dlon))
        ci = np.clip(ci, 0, nlon - 1)
    ri = np.clip(np.rint((lat - lats[0]) / dlat).astype(int), 0, nlat - 1)
    return ri, ci


def raster_lookup(da, lon, lat):
    """Nearest-node raster values at lon/lat (vectorised, longitude-wrapping)."""
    lons, lats = raster_coords(da)
    ri, ci = raster_nearest_index(lons, lats, lon, lat)
    return da.values[ri, ci]


def cell_area_weights(lats, shape):
    """Relative area of each cell of a regular lon/lat grid (∝ cos(lat))."""
    w = np.clip(np.cos(np.radians(lats)), 0.0, None)
    return np.broadcast_to(w[:, None], shape)


def plate_areas(plate_id_raster):
    """{plate_id: relative spherical area} from a plate-ID raster (cos(lat)-weighted pixel counts)."""
    _, lats = raster_coords(plate_id_raster)
    vals = plate_id_raster.values
    w = cell_area_weights(lats, vals.shape)
    ok = np.isfinite(vals) & (vals > 0)
    pids, inv = np.unique(vals[ok].astype(int), return_inverse=True)
    sums = np.bincount(inv, weights=w[ok])
    return {int(p): float(s) for p, s in zip(pids, sums)}


def grid_stride(raster, grid_deg):
    """Row/column slices that sample a regular global raster at about grid_deg spacing,
    taking the cell nearest the centre of each coarse cell.  Returns (slices, k)."""
    lons, _ = raster_coords(raster)
    k = max(1, int(round(grid_deg / abs(lons[1] - lons[0])))) if grid_deg else 1
    return (slice(k // 2, None, k), slice(k // 2, None, k)), k


def plate_area_tensors(plate_id_raster, weight=None, grid_deg=1.0):
    """Area tensor of each plate on the unit sphere: T_p = Σ w (𝟙 − r rᵀ) dA.

    For a rigid plate rotating at angular velocity ω, ωᵀ T ω = ∫ w |ω × r|² dA, so with
    w = 1 the tensor gives the plate's mean-square speed, and with w = a drag weight it
    is the basal drag tensor (torque = −c T ω).  dA is in steradians (cos(lat)
    weighted); the raster is sampled at about grid_deg spacing for speed.

    weight : optional array, same shape as the raster (e.g. a continental drag factor).
    Returns (pids, tensors (n, 3, 3), areas (n,) in steradians).
    """
    sl, _ = grid_stride(plate_id_raster, grid_deg)
    lons, lats = raster_coords(plate_id_raster)
    vals = plate_id_raster.values[sl]
    lon2d, lat2d = np.meshgrid(lons[sl[1]], lats[sl[0]])
    dlon = np.radians(abs(lons[sl[1]][1] - lons[sl[1]][0]))
    dlat = np.radians(abs(lats[sl[0]][1] - lats[sl[0]][0]))
    dA = np.cos(np.radians(lat2d)) * dlon * dlat
    ok = np.isfinite(vals) & (vals > 0)
    pids, inv = np.unique(vals[ok].astype(int), return_inverse=True)
    r = lonlat_to_xyz(lon2d[ok], lat2d[ok])
    area_w = dA[ok]
    w = area_w if weight is None else area_w * np.asarray(weight)[sl][ok]
    tensors = np.zeros((len(pids), 3, 3))
    for i in range(3):
        for j in range(i, 3):
            comp = np.bincount(inv, weights=w * ((i == j) - r[:, i] * r[:, j]),
                               minlength=len(pids))
            tensors[:, i, j] = comp
            tensors[:, j, i] = comp
    areas = np.bincount(inv, weights=area_w, minlength=len(pids))
    return pids, tensors, areas


def stage_pole_omega(rotation, delta_time):
    """Angular velocity vector (rad/Myr) of a stage rotation over delta_time."""
    axis, angle = rotation_axis_angle(rotation)
    return axis * angle / delta_time


def wrap_neighbours(mask):
    """4-connected neighbour shifts of a boolean/int array on a global grid:
    longitude wraps (np.roll), latitude does not (poles are edges)."""
    up    = np.vstack([mask[1:], mask[-1:]])
    down  = np.vstack([mask[:1], mask[:-1]])
    east  = np.roll(mask, -1, axis=1)
    west  = np.roll(mask,  1, axis=1)
    return up, down, east, west


def label_wrapped(mask, connectivity=1):
    """Connected-component labels of a mask on a global lon/lat grid, joining
    components across the ±180 seam (first and last columns are neighbours;
    if the grid duplicates the seam column they are the same meridian)."""
    from scipy.ndimage import label
    structure = None
    if connectivity == 2:
        structure = np.ones((3, 3), dtype=int)
    labels, n = label(mask, structure=structure)
    if n < 2:
        return labels, n
    # Union labels that touch across the seam
    parent = np.arange(n + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    first, last = labels[:, 0], labels[:, -1]
    for a, b in zip(first, last):
        if a > 0 and b > 0:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    roots = np.array([find(i) for i in range(n + 1)])
    uniq = np.unique(roots[1:])
    remap = np.zeros(n + 1, dtype=int)
    remap[1:] = np.searchsorted(uniq, roots[1:]) + 1
    return remap[labels], len(uniq)


def spherical_gaussian_smooth(field, lats, lons, sigma_deg):
    """Gaussian smoothing on the sphere of one 2-D field (nlat, nlon) or a stack
    (..., nlat, nlon), using great-circle distances and cell-area weights.

    Unlike scipy.ndimage.gaussian_filter on a lon/lat grid, the kernel has the same
    physical width everywhere, wraps across ±180 and passes smoothly over the poles.
    The sparse kernel is cached per grid/sigma.
    """
    if sigma_deg is None or sigma_deg <= 0:
        return field
    K = _spherical_kernel(tuple(np.round(lats, 6)), tuple(np.round(lons, 6)), float(sigma_deg))
    shape = field.shape
    flat = field.reshape(-1, shape[-2] * shape[-1])
    out = (K @ flat.T).T
    return out.reshape(shape)


_KERNEL_CACHE = {}


def _spherical_kernel(lats, lons, sigma_deg):
    key = (lats, lons, sigma_deg)
    if key in _KERNEL_CACHE:
        return _KERNEL_CACHE[key]
    from scipy.sparse import csr_matrix
    lon2d, lat2d = np.meshgrid(lons, lats)
    pts = np.radians(np.column_stack([lat2d.ravel(), lon2d.ravel()]))
    area = np.clip(np.cos(pts[:, 0]), 1e-6, None)
    tree = BallTree(pts, metric='haversine')   # radius query; built once per grid (cached)
    idx, dist = tree.query_radius(pts, r=np.radians(4 * sigma_deg), return_distance=True)
    rows = np.repeat(np.arange(len(pts)), [len(i) for i in idx])
    cols = np.concatenate(idx)
    w = np.exp(-0.5 * (np.degrees(np.concatenate(dist)) / sigma_deg) ** 2) * area[cols]
    K = csr_matrix((w, (rows, cols)), shape=(len(pts), len(pts)))
    row_sums = np.asarray(K.sum(axis=1)).ravel()
    K = csr_matrix(K.multiply(1.0 / row_sums[:, None]))
    _KERNEL_CACHE[key] = K
    return K


def nudge_topological_polygons(snapshot, stage_pole_dict=None):
    '''
    Return the resolved plate polygons of a pygplates snapshot as ordinary polygon
    features, optionally rotated by a {plate_id: FiniteRotation} dict.

    Features are cloned so the reconstruction model's topological features are
    not modified.
    '''
    if stage_pole_dict is None:
        stage_pole_dict = {}
    identity = pygplates.FiniteRotation()
    rpolygons = []
    for topology in snapshot.resolved_topologies:
        plate_id = topology.get_resolved_feature().get_reconstruction_plate_id()
        rotation = stage_pole_dict.get(plate_id) or identity
        nfeature = topology.get_feature().clone()
        nfeature.set_geometry(rotation * topology.get_resolved_geometry())
        rpolygons.append(nfeature)
    return rpolygons


def densify_plate_polygons(gdf, max_edge_deg=1.0):
    """Copy of a plate-polygon GeoDataFrame ready for pygplates.

    The polygons come from simplification in lon/lat, so their edges are straight
    lines in lon/lat.  pygplates reads every edge as a great-circle arc, which would
    make long edges bow poleward (and a pole-to-pole edge undefined).  Adding
    vertices every max_edge_deg keeps the intended shape.  (Boundaries are extracted
    from the undensified polygons; see topo_polygons_to_boundaries.)
    """
    gdf = gdf.copy()
    gdf.geometry = [_densify_lonlat(g, max_edge_deg) for g in gdf.geometry]
    return gdf


def nudge_polygons(polygon_features, stage_pole_dict):

    if isinstance(polygon_features, gpd.GeoDataFrame):
        polygon_features = gdf2gpml(densify_plate_polygons(polygon_features.set_crs(epsg=4326, allow_override=True)))

    identity = pygplates.FiniteRotation()
    rpolygons = []
    for polygon_feature in polygon_features:
        geometry = polygon_feature.get_geometry()
        if geometry is None:
            print('  [nudge] WARNING: feature with plate ID {} has no valid geometry — skipped'.format(
                polygon_feature.get_reconstruction_plate_id()))
            continue
        plate_id = polygon_feature.get_reconstruction_plate_id()
        rotation = stage_pole_dict.get(plate_id) or identity
        nfeature = polygon_feature.clone()
        nfeature.set_geometry(rotation * geometry)
        rpolygons.append(nfeature)

    return rpolygons


def tessellate_polygons(polygons, max_dist_deg=0.2):

    new_features = []
    for feature in polygons:
        new_feature = feature.clone()
        new_feature.set_geometry(feature.get_geometry().to_tessellated(np.radians(max_dist_deg)))
        new_features.append(new_feature)

    return new_features


def plate_polygons_to_wrapped_gdf(polygons):
    """Convert pygplates polygon features to a GeoDataFrame, splitting polygons at
    the dateline.  Interior rings (holes) are kept, so a plate that surrounds other
    plates does not cover them when overlaps are resolved."""
    rp_gdf = gpml2gdf(polygons)
    wrapper = pygplates.DateLineWrapper(0.0)

    rows = []
    for _, row in rp_gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        for poly in getattr(geom, 'geoms', [geom]):
            exterior = [(c[1], c[0]) for c in poly.exterior.coords[:-1]]
            interiors = [[(c[1], c[0]) for c in ring.coords[:-1]] for ring in poly.interiors]
            try:
                sphere_poly = pygplates.PolygonOnSphere(exterior, interiors)
            except Exception:
                sphere_poly = pygplates.PolygonOnSphere(exterior)
            for wrapped in wrapper.wrap(sphere_poly):
                shell = [p.to_lat_lon()[::-1] for p in wrapped.get_exterior_points()]
                holes = [[p.to_lat_lon()[::-1] for p in wrapped.get_interior_points(k)]
                         for k in range(wrapped.get_number_of_interior_rings())]
                new_row = row.copy()
                new_row['geometry'] = Polygon(shell, holes)
                rows.append(new_row)

    return gpd.GeoDataFrame(rows, crs=rp_gdf.crs).reset_index(drop=True)


def resolve_overlaps(gdf, priority_col):
    gdf = gdf.sort_values(priority_col).reset_index(drop=True)
    gdf['geometry'] = gdf.geometry.make_valid()

    result_geoms = []
    accumulated = None

    for i in range(len(gdf)):
        geom = gdf.geometry.iloc[i]
        if accumulated is not None:
            geom = geom.difference(accumulated)
        result_geoms.append(geom)
        # Build cumulative union of original geometries
        accumulated = gdf.geometry.iloc[i] if accumulated is None else accumulated.union(gdf.geometry.iloc[i])

    result = gdf.copy()
    result.geometry = result_geoms

    return result


def update_plate_id_raster(clipped_polygons, raster_spacing=1, tie_break='hash',
                           previous_plate_id=None, tie_break_seed=0):
    valid = clipped_polygons.dropna(subset=['geometry'])
    valid = valid[~valid.geometry.is_empty]
    plate_id_raster = rasterize_polygons(valid, zval_field='PLATEID1', sampling=raster_spacing)
    plate_id_raster = plate_id_raster.where(plate_id_raster>0, np.nan)

    plate_id_raster = fill_nan_spherical(plate_id_raster, tie_break=tie_break,
                                         previous=previous_plate_id, seed=tie_break_seed)

    return plate_id_raster


def update_plate_polygons(plate_id_raster, raster_spacing=1, topo_simplify=1):
    # Origin offset by half a cell so polygon edges align with raster cell boundaries
    transform = (raster_spacing, 0, -(180 + raster_spacing / 2),
                 0, raster_spacing, -(90 + raster_spacing / 2))
    new_plate_polygons = polygonize(plate_id_raster,
                                    return_type='geopandas', transform=transform)

    topo = tp.Topology(new_plate_polygons, prequantize=False)
    new_plate_polygons_smooth = topo.toposimplify(topo_simplify).to_gdf()

    # Fix any invalid geometries produced by toposimplify before clipping
    new_plate_polygons_smooth.geometry = new_plate_polygons_smooth.geometry.make_valid()
    new_plate_polygons_smooth = new_plate_polygons_smooth.clip(box(-180, -90, 180, 90))

    new_plate_polygons_smooth = new_plate_polygons_smooth.rename(columns={'DN':'PLATEID1'})

    return new_plate_polygons_smooth.set_crs(epsg=4326)


def _densify_coords(coords, max_edge_deg):
    """Add vertices along a lon/lat line so no edge is longer than max_edge_deg
    (each edge is interpolated from its lexicographically smaller end, so the same
    edge gives the same vertices whichever way it is traversed)."""
    coords = [tuple(c[:2]) for c in coords]
    out = [coords[0]]
    for a, b in zip(coords[:-1], coords[1:]):
        n = int(np.ceil(np.hypot(b[0] - a[0], b[1] - a[1]) / max_edge_deg))
        if n > 1:
            p, q = (a, b) if a <= b else (b, a)
            pts = [(p[0] + (q[0] - p[0]) * k / n, p[1] + (q[1] - p[1]) * k / n) for k in range(1, n)]
            out.extend(pts if p == a else pts[::-1])
        out.append(b)
    return out


def _densify_lonlat(geom, max_edge_deg):
    """Densify polygon edges (straight lines in lon/lat) to at most max_edge_deg, so
    that pygplates, which reads every edge as a great-circle arc, keeps their shape."""
    if geom is None or geom.is_empty:
        return geom
    if geom.geom_type == 'MultiPolygon':
        return type(geom)([_densify_lonlat(g, max_edge_deg) for g in geom.geoms])
    if geom.geom_type != 'Polygon':
        return geom
    return Polygon(_densify_coords(geom.exterior.coords, max_edge_deg),
                   [_densify_coords(r.coords, max_edge_deg) for r in geom.interiors])


def velocity_grid(plate_id_raster, plate_parameters, delta_time=1, raster_spacing=1):
    '''
    East and north velocity (km/Myr) of every raster node, from the stage rotation
    of the plate the node belongs to.
    '''
    lons, lats = raster_coords(plate_id_raster)
    pids = plate_id_raster.values     # (nlat, nlon)

    vel_east  = np.zeros(pids.shape, dtype=float)
    vel_north = np.zeros(pids.shape, dtype=float)
    lon2d, lat2d = np.meshgrid(lons, lats)

    for pid_val in np.unique(pids[np.isfinite(pids)]):
        rotation = plate_parameters.get(int(pid_val))
        if rotation is None:
            continue
        mask = pids == pid_val
        vel_east[mask], vel_north[mask] = plate_velocity_ne(
            rotation, lon2d[mask], lat2d[mask], delta_time)

    coords = {'y': lats, 'x': lons}
    da_east  = xr.DataArray(vel_east,  dims=['y', 'x'], coords=coords,
                            attrs={'units': 'km/Myr'})
    da_north = xr.DataArray(vel_north, dims=['y', 'x'], coords=coords,
                            attrs={'units': 'km/Myr'})
    return da_east, da_north


def spherical_divergence(vel_east, vel_north):
    """Horizontal divergence (1/Myr) of a velocity field (km/Myr) on a global lon/lat grid:

        div = 1/(R cos φ) · [ ∂u/∂λ + ∂(v cos φ)/∂φ ]

    Longitude is treated as periodic.  Rigid plate interiors give ~0; the signal is
    concentrated at plate boundaries (positive = spreading, negative = convergence).
    """
    lons, lats = raster_coords(vel_east)
    u = vel_east.values
    v = vel_north.values
    seam = abs((lons[-1] - lons[0]) - 360.0) < 1e-6
    if seam:
        u, v = u[:, :-1], v[:, :-1]
    lam = np.radians(lons[1] - lons[0])
    phi = np.radians(lats)
    cos_phi = np.cos(phi)[:, None]
    du_dlam = (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2.0 * lam)
    dvcos_dphi = np.gradient(v * cos_phi, phi, axis=0)
    # east/north components are undefined at the poles: limit 1/cos φ to half a cell away
    min_cos = np.cos(np.radians(90.0 - abs(lats[1] - lats[0]) / 2))
    div = (du_dlam + dvcos_dphi) / (EARTH_RADIUS_KM * np.clip(cos_phi, min_cos, None))
    if seam:
        div = np.hstack([div, div[:, :1]])
    return xr.DataArray(div, dims=vel_east.dims, coords=vel_east.coords,
                        attrs={'units': '1/Myr', 'long_name': 'horizontal divergence'})


def _velocities_by_plate(pids, lon, lat, stage_pole_dict, delta_time):
    """East/north velocity (km/Myr) at each point for the plate ID given per point."""
    east = np.zeros(len(pids))
    north = np.zeros(len(pids))
    for pid in np.unique(pids):
        rotation = stage_pole_dict.get(int(pid))
        if rotation is None:
            continue
        m = pids == pid
        east[m], north[m] = plate_velocity_ne(rotation, lon[m], lat[m], delta_time)
    return east, north


def compute_boundary_divergence(boundaries_gdf, stage_pole_dict, delta_time, plate_id_raster,
                                plate_polygons_gdf=None):
    """Compute per-arc-segment relative normal velocity.

    Explodes each boundary LineString into individual 2-vertex segments and computes
    divergence independently for each.  Returns a GeoDataFrame where each row is one
    arc segment with its own 'divergence' value in km/Myr (= mm/yr).  Positive = ridge,
    negative = subduction.  Different segments of the same plate-pair boundary can have
    different signs, correctly capturing oblique or mixed boundaries.

    Sign is determined by probing the plate-ID raster on both sides of each segment
    midpoint (a great-circle step along the arc normal) — more robust than a centroid
    hemisphere test for large or concave plates.  Segments whose side cannot be
    determined are kept with divergence 0.
    """
    OFFSET_DEG = 0.5

    # Pass 1: split the boundaries into 2-vertex segments and collect their geometry
    new_rows = []
    seg_row, seg_mid, seg_normal, seg_pids = [], [], [], []
    for _, row in boundaries_gdf.iterrows():
        pid1 = row.get('plateid1_poly1')
        pid2 = row.get('plateid1_poly2')
        coords = list(row.geometry.coords)

        if not pid1 or not pid2 or pid1 == 0 or pid2 == 0 or len(coords) < 2:
            r = row.to_dict()
            r['divergence'] = 0.0
            new_rows.append(r)
            continue

        for i in range(len(coords) - 1):
            r = row.to_dict()
            r['geometry'] = LineString([coords[i][:2], coords[i+1][:2]])
            r['divergence'] = 0.0
            new_rows.append(r)
            try:
                arc = pygplates.GreatCircleArc(pygplates.PointOnSphere(coords[i][1],   coords[i][0]),
                                               pygplates.PointOnSphere(coords[i+1][1], coords[i+1][0]))
                normal = arc.get_great_circle_normal().to_xyz()
            except pygplates.IndeterminateGreatCircleArcNormalError:
                continue
            seg_row.append(len(new_rows) - 1)
            seg_mid.append(arc.get_arc_point(0.5).to_lat_lon())
            seg_normal.append(normal)
            seg_pids.append((pid1, pid2))

    if not new_rows:
        result = boundaries_gdf.copy()
        result['divergence'] = pd.Series(dtype=float)
        return result

    # Pass 2 (vectorised): probe both sides of each midpoint, then relative velocity
    if seg_row:
        mid_lat, mid_lon = np.array(seg_mid).T
        normal = np.array(seg_normal)
        pid1, pid2 = np.array(seg_pids, dtype=float).T

        # Sign: probe the raster on both sides of the segment midpoint.  Probing both
        # directions handles triple-junction edges where one probe may land on a
        # third plate.
        pos_lon, pos_lat = offset_along_normal(mid_lon, mid_lat,  normal, OFFSET_DEG)
        neg_lon, neg_lat = offset_along_normal(mid_lon, mid_lat, -normal, OFFSET_DEG)
        pos_val = raster_lookup(plate_id_raster, pos_lon, pos_lat)
        neg_val = raster_lookup(plate_id_raster, neg_lon, neg_lat)
        votes = ((pos_val == pid1).astype(int) + (neg_val == pid2)) - ((pos_val == pid2).astype(int) + (neg_val == pid1))
        sign = np.sign(votes)    # +1: normal points toward pid1; 0: side unknown (e.g. triple junction)

        # divergence = (v_pid1 - v_pid2) · n_toward_pid1
        e1, n1 = _velocities_by_plate(pid1, mid_lon, mid_lat, stage_pole_dict, delta_time)
        e2, n2 = _velocities_by_plate(pid2, mid_lon, mid_lat, stage_pole_dict, delta_time)
        normal_e, normal_n = _tangent_ne(normal.T, mid_lon, mid_lat)
        divergence = sign * ((e1 - e2) * normal_e + (n1 - n2) * normal_n)
        for k, row_index in enumerate(seg_row):
            new_rows[row_index]['divergence'] = float(divergence[k])

    return gpd.GeoDataFrame(new_rows, geometry='geometry', crs=boundaries_gdf.crs)


def _tangent_ne(vec_xyz, lon, lat):
    """East and north components of a 3-D vector at a point on the sphere."""
    lon_r, lat_r = np.radians(lon), np.radians(lat)
    x, y, z = vec_xyz
    east  = -np.sin(lon_r) * x + np.cos(lon_r) * y
    north = -np.sin(lat_r) * np.cos(lon_r) * x - np.sin(lat_r) * np.sin(lon_r) * y + np.cos(lat_r) * z
    return east, north


def _rings_with_outward_sign(polygon):
    """Yield (segments, outward_sign) for each ring of a pygplates PolygonOnSphere.

    outward_sign converts a segment's great-circle normal into the normal pointing
    out of the plate: for a clockwise ring the great-circle normal points out of the
    area the ring encloses.  For an interior ring (hole) "out of the plate" means
    into the hole, so the sign is reversed.
    """
    def _ring_sign(points):
        ring = pygplates.PolygonOnSphere(points)
        return 1 if ring.get_orientation() == pygplates.PolygonOnSphere.Orientation.clockwise else -1

    yield polygon.get_exterior_ring_segments(), _ring_sign(polygon.get_exterior_ring_points())
    for k in range(polygon.get_number_of_interior_rings()):
        yield (polygon.get_interior_ring_segments(k),
               -_ring_sign(polygon.get_interior_ring_points(k)))


# Seafloor age for slab pull and polarity is read this far inside the plate from the
# trench: the physics age grid is 1°, and the tracers nearest the trench are removed.
TRENCH_AGE_PROBE_DEG = 1.0


def polygons_to_stage_poles(polygons_gdf, plate_id_raster, forces,
                            craton_raster=None,
                            stage_pole_dict=None,
                            delta_time=1,
                            divergence_threshold_for_subduction=-1.0,
                            divergence_threshold_for_ridge_push=1.0,
                            boundary_search_radius_deg=2,
                            polarity_dict=None,
                            collision_pairs=None,
                            collision_mode='resistance',
                            return_torques=False,
                            modes=None,
                            age_grid=None,
                            slab_tracers=None,
                            force_params=None,
                            diagnostics=None):
    """Sum boundary forces around every plate and convert the net torque to a stage rotation.

    With return_torques=True, also returns the net torque of each plate (xyz array) and
    its boundary (edge) tensor Σ L (𝟙 − m mᵀ) over its polygon segments (L = arc length in
    radians, m = segment midpoint), for solve_plate_motions.

    For each polygon segment, the neighbouring plate is found in the plate-ID raster and
    the relative normal velocity (km/Myr, from stage_pole_dict) classifies the segment:
      - convergent, continent-continent: collision resistance (or coupling)
      - convergent: slab pull on the downgoing plate, slab suction on the overriding plate
      - divergent: ridge push

    modes (force_modes names; the defaults are the simple model):
      slab_pull     hard         forces['slab_pull'] on every convergent segment
                    ramp         × clip(convergence / pull_ramp_km_per_myr, 0, 1)
                    slab_length  × √(min(age, A)/A) × min(1, z / z_ref), where age is the
                                 seafloor age 1° inside the downgoing plate (age_grid; A =
                                 slab_age_saturation_myr) and z the depth of the deepest
                                 slab tracer near the segment (slab_depth_at_segments,
                                 z_ref = slab_ref_depth_km).  Any segment that is not
                                 divergent and whose plate pair has subducted (has a
                                 polarity) is pulled if there is slab: a slab pulls by
                                 hanging there, whatever the convergence rate.
      slab_suction  ratio        forces['slab_suction'] on the overriding plate
                    depth_split  forces['slab_suction_deep'] × min(1, max(0, z − z_ref)/z_ref)
                                 on both plates, toward the trench (slab in the lower mantle)
      polarity      id           a new ocean–ocean pair: the lower plate ID subducts
                    age          the plate with the older ocean along the contact subducts
      In every mode a plate with continent (cratons) overrides one without, as before.
    diagnostics, if a dict, receives 'segments': a DataFrame of the subduction segments
    and their forces (only when a mode other than the simple ones is used).
    """
    if polarity_dict is None:
        polarity_dict = {}
    if collision_pairs is None:
        collision_pairs = set()
    if stage_pole_dict is None:
        stage_pole_dict = {}
    modes = {'slab_pull': 'hard', 'slab_suction': 'ratio', 'polarity': 'id', **(modes or {})}
    force_params = force_params or {}
    identity = pygplates.FiniteRotation()

    # Craton pixel count per plate ID — used to determine overriding plate at new boundaries
    craton_pixel_counts = craton_counts_from_raster(craton_raster)

    gpml_plates = gdf2gpml(densify_plate_polygons(polygons_gdf.explode(index_parts=False).reset_index(drop=True)))

    # Pass 1: walk every polygon ring and collect the segments (midpoint, outward
    # normal and the torque of a unit force pulling the plate toward the boundary).
    torque_dict = {}
    seg_pid, seg_mid, seg_outward, seg_torque, seg_len = [], [], [], [], []
    for plate in gpml_plates:
        pid_current = plate.get_reconstruction_plate_id()
        torque_dict.setdefault(pid_current, pygplates.Vector3D(0, 0, 0))
        for geometry in plate.get_geometries():
            for segments, outward_sign in _rings_with_outward_sign(geometry):
                for segment in segments:
                    if segment.is_zero_length():
                        continue
                    seg_pid.append(pid_current)
                    seg_mid.append(segment.get_arc_point(0.5).to_lat_lon())
                    seg_outward.append(np.array(segment.get_great_circle_normal().to_xyz()) * outward_sign)
                    # pdriver.get_torque_for_segment(segment, s) pushes the plate along
                    # s * great-circle normal, so s = outward_sign pulls the plate toward
                    # the boundary; the torque is linear in s (a push is its negative).
                    seg_torque.append(pdriver.get_torque_for_segment(segment, outward_sign))
                    seg_len.append(segment.get_arc_length())

    if seg_pid:
        seg_pid = np.array(seg_pid)
        mid_lat, mid_lon = np.array(seg_mid).T
        outward = np.array(seg_outward)

        # Neighbouring plate across each segment
        ids = two_nearest_ids_batch(mid_lat, mid_lon, plate_id_raster, boundary_search_radius_deg)
        neighbour = np.where(ids[:, 0] == seg_pid, ids[:, 1], ids[:, 0])
        valid = (neighbour > 0) & (neighbour != seg_pid)

        # Relative normal velocity: positive = diverging, negative = converging.
        # divergence = (v_neighbor - v_current) · outward_normal
        out_e, out_n = _tangent_ne(outward.T, mid_lon, mid_lat)
        e_cur, n_cur = _velocities_by_plate(seg_pid, mid_lon, mid_lat, stage_pole_dict, delta_time)
        e_nbr, n_nbr = _velocities_by_plate(neighbour, mid_lon, mid_lat, stage_pole_dict, delta_time)
        divergence = (e_nbr - e_cur) * out_e + (n_nbr - n_cur) * out_n

        convergent = valid & (divergence < divergence_threshold_for_subduction)
        divergent  = valid & (divergence > divergence_threshold_for_ridge_push)
        pair_lo = np.minimum(seg_pid, neighbour).astype(int)
        pair_hi = np.maximum(seg_pid, neighbour).astype(int)

        # Seafloor age 1° inside each segment's own plate (the downgoing side when the
        # plate subducts), used by age polarity and the age term of slab pull.  NaN
        # where the probe is on continent or there is no age grid.
        inner_age = np.full(len(seg_pid), np.nan)
        if age_grid is not None:
            probe_lon, probe_lat = offset_along_normal(mid_lon, mid_lat, -outward, TRENCH_AGE_PROBE_DEG)
            inner_age = np.asarray(raster_lookup(age_grid, probe_lon, probe_lat), dtype=float)

        # Polarity of plate pairs converging for the first time (sticky — only set once).
        # Skip when craton_raster is None so the craton rule is not contaminated by an
        # arbitrary oceanic fallback before craton info is available.
        if craton_raster is not None:
            for k in np.where(convergent)[0]:
                pair = (int(pair_lo[k]), int(pair_hi[k]))
                if pair in polarity_dict:
                    continue
                downgoing = resolve_downgoing_plate(pair[0], pair[1], {}, craton_pixel_counts)
                if (modes['polarity'] == 'age' and age_grid is not None
                        and not craton_pixel_counts.get(pair[0]) and not craton_pixel_counts.get(pair[1])):
                    # ocean–ocean: the older side along this step's contact subducts
                    on_pair = convergent & (pair_lo == pair[0]) & (pair_hi == pair[1])
                    side_age = [np.nanmean(inner_age[on_pair & (seg_pid == p)])
                                if np.isfinite(inner_age[on_pair & (seg_pid == p)]).any() else np.nan
                                for p in pair]
                    if np.all(np.isfinite(side_age)) and side_age[0] != side_age[1]:
                        downgoing = pair[0] if side_age[0] > side_age[1] else pair[1]
                    print('  [polarity] pair ({},{}): seafloor age {:.0f} / {:.0f} Myr -> {} subducts'.format(
                        pair[0], pair[1], side_age[0], side_age[1], downgoing))
                polarity_dict[pair] = downgoing

        # Segments that can carry slab forces.  Classic modes: convergent segments.
        # slab_length: any non-divergent segment of a pair that has subducted.
        slab_set = convergent.copy()
        if modes['slab_pull'] == 'slab_length':
            has_polarity = np.array([(int(a), int(b)) in polarity_dict for a, b in zip(pair_lo, pair_hi)])
            slab_set |= valid & ~divergent & has_polarity

        # Continent-continent contact: two different cratons within the search radius
        collision = np.zeros(len(seg_pid), bool)
        if craton_raster is not None and slab_set.any():
            c_ids = two_nearest_ids_batch(mid_lat[slab_set], mid_lon[slab_set], craton_raster,
                                          boundary_search_radius_deg, ignore_zero=True)
            collision[slab_set] = (c_ids[:, 0] > 0) & (c_ids[:, 1] > 0)

        # Slab forces per segment: pull on the downgoing plate, suction on the overriding
        # plate (and on the downgoing plate too with depth_split).  Both draw the plate
        # toward the trench.
        subduct = slab_set & ~collision
        is_down = np.zeros(len(seg_pid), bool)
        for k in np.where(subduct)[0]:
            is_down[k] = resolve_downgoing_plate(int(seg_pid[k]), int(neighbour[k]),
                                                 polarity_dict, craton_pixel_counts) == seg_pid[k]
        w_pull = forces.get('slab_pull', 0)
        pull_w = np.where(convergent & is_down, w_pull, 0.0)
        suction_w = np.where(convergent & ~is_down, forces.get('slab_suction', 0), 0.0)
        slab_depth = np.zeros(len(seg_pid))
        z_ref = force_params.get('slab_ref_depth_km', 660.0)
        if modes['slab_pull'] == 'ramp':
            ramp = force_params.get('pull_ramp_km_per_myr', 10.0)
            pull_w = pull_w * np.clip(-divergence / ramp, 0.0, 1.0)
        if (modes['slab_pull'] == 'slab_length' or modes['slab_suction'] == 'depth_split') \
                and slab_tracers is not None and subduct.any():
            # the slab lies on the overriding side: outward from the downgoing plate
            side = np.where(is_down[:, None], outward, -outward)
            slab_depth[subduct] = slab_depth_at_segments(
                mid_lon[subduct], mid_lat[subduct], side[subduct],
                slab_tracers['lon'].values, slab_tracers['lat'].values,
                slab_tracers['depth_km'].values,
                radius_deg=force_params.get('slab_search_radius_deg', 5.0),
                line_tol_deg=force_params.get('slab_line_tol_deg', 1.0))
        if modes['slab_pull'] == 'slab_length':
            a_sat = force_params.get('slab_age_saturation_myr', 80.0)
            f_age = np.sqrt(np.clip(np.nan_to_num(inner_age, nan=a_sat), 0.0, a_sat) / a_sat)
            pull_w = np.where(subduct & is_down, w_pull * f_age * np.minimum(1.0, slab_depth / z_ref), 0.0)
        if modes['slab_suction'] == 'depth_split':
            deep = np.clip((slab_depth - z_ref) / z_ref, 0.0, 1.0)
            suction_w = np.where(subduct, forces.get('slab_suction_deep', 0.0) * deep, 0.0)
            slab_w = pull_w + suction_w
        else:
            slab_w = np.where(is_down, pull_w, suction_w)

        if diagnostics is not None and any(modes[c] != d for c, d in
                                           (('slab_pull', 'hard'), ('slab_suction', 'ratio'), ('polarity', 'id'))):
            s = np.where(slab_set)[0]
            diagnostics['segments'] = pd.DataFrame({
                'lon': mid_lon[s], 'lat': mid_lat[s], 'plate_id': seg_pid[s].astype(int),
                'neighbour_id': neighbour[s].astype(int), 'downgoing': is_down[s],
                'collision': collision[s], 'convergence_km_per_myr': -divergence[s],
                'length_km': np.array(seg_len)[s] * EARTH_RADIUS_KM,
                'slab_depth_km': slab_depth[s], 'seafloor_age_myr': inner_age[s],
                'pull_weight': pull_w[s], 'suction_weight': np.where(is_down[s], slab_w[s] - pull_w[s], slab_w[s])})

        # Pass 2: add up the torques plate by plate, in segment order
        for k in np.where(slab_set | divergent)[0]:
            pid_current, pid_neighbor = int(seg_pid[k]), int(neighbour[k])
            if slab_set[k]:
                if collision[k]:
                    if not convergent[k]:
                        continue
                    pair = (int(pair_lo[k]), int(pair_hi[k]))
                    if collision_mode == 'coupling':
                        # Contribute zero torque; record pair for post-step velocity
                        # blending.  There is no subducting slab to pull the plates.
                        collision_pairs.add(pair)
                    elif 'continental_collision' in forces:
                        # Resistance torque pushing the plate away from the collision.
                        torque_dict[pid_current] = torque_dict[pid_current] - seg_torque[k] * forces['continental_collision']
                elif slab_w[k]:
                    torque_dict[pid_current] = torque_dict[pid_current] + seg_torque[k] * float(slab_w[k])
            elif 'ridge_push' in forces:
                torque_dict[pid_current] = torque_dict[pid_current] - seg_torque[k] * forces['ridge_push']

    new_stage_pole_dict = {pid: pdriver.torque_to_stage_rotation(torque, delta_time)
                           for pid, torque in torque_dict.items()}
    if not return_torques:
        return new_stage_pole_dict, collision_pairs

    torques = {pid: np.array(t.to_xyz()) for pid, t in torque_dict.items()}
    edge_tensors = {}
    if len(seg_pid):
        m = lonlat_to_xyz(mid_lon, mid_lat)
        L = np.array(seg_len)
        outer = L[:, None, None] * (np.eye(3)[None] - m[:, :, None] * m[:, None, :])
        for pid in np.unique(seg_pid):
            edge_tensors[int(pid)] = outer[seg_pid == pid].sum(axis=0)
    return new_stage_pole_dict, collision_pairs, torques, edge_tensors


def solve_plate_motions(torques, edge_tensors, plate_id_raster, craton_raster, params,
                        drag_mode, speed_mode, speed_state, delta_time):
    """Angular velocities from net torques with basal and boundary resistance, scaled to a
    target RMS plate speed.  Returns {pid: FiniteRotation over delta_time}.

    drag_mode
      none       ω = τ (the simple model)
      area       ω = D⁻¹ τ with D = Σ w (𝟙 − r rᵀ) dA over the plate (basal drag;
                 w = continent_drag_factor on continental cells, 1 on ocean)
      area_edge  as area, plus boundary resistance: D += ℓ Σ L (𝟙 − m mᵀ) over the plate
                 boundary, ℓ = edge_drag_width_km (boundary resistance equals basal drag
                 over a strip of that width along every boundary)
    speed_mode
      none            ω scaled by speed_scale (default 1)
      calibrate_once  scale k set at the first call so the area-weighted RMS speed equals
                      target_rms_speed, then kept (speed_scale, if set, is used instead)
      per_step        k recomputed every step
    If the RMS speed would exceed max_rms_speed_factor × target_rms_speed, that step is
    scaled down to the cap.  speed_state (a dict kept by the caller) holds k and a record
    of the last step: {'k', 'last': {'k_used', 'rms_unscaled', 'rms', 'capped'}}.
    Everything is on the unit sphere; RMS speeds are exact (from area tensors, no grid).
    """
    grid_deg = params.get('physics_grid_deg', 1.0)
    target = params.get('target_rms_speed', 50.0)
    cap = params.get('max_rms_speed_factor', 3.0) * target
    fixed_k = params.get('speed_scale')

    pids, J, areas = plate_area_tensors(plate_id_raster, grid_deg=grid_deg)
    if set(torques) - set(pids.tolist()):
        # a plate too small for the coarse grid: use the full-resolution raster
        grid_deg = None
        pids, J, areas = plate_area_tensors(plate_id_raster, grid_deg=None)
    J_of = dict(zip(pids.tolist(), J))
    area_of = dict(zip(pids.tolist(), areas))

    if drag_mode == 'none':
        omega = {pid: t for pid, t in torques.items()}
    else:
        weight = None
        if craton_raster is not None:
            cont = np.nan_to_num(craton_raster.values) > 0
            weight = np.where(cont, params.get('continent_drag_factor', 3.0), 1.0)
        _, D, _ = plate_area_tensors(plate_id_raster, weight=weight, grid_deg=grid_deg)
        D_of = dict(zip(pids.tolist(), D))
        width = params.get('edge_drag_width_km', 3000.0) / EARTH_RADIUS_KM
        omega = {}
        for pid, t in torques.items():
            Dp = D_of.get(pid, np.zeros((3, 3)))
            if drag_mode == 'area_edge':
                Dp = Dp + width * edge_tensors.get(pid, np.zeros((3, 3)))
            if not np.any(Dp):
                omega[pid] = np.zeros(3)
                continue
            omega[pid] = np.linalg.lstsq(Dp, t, rcond=1e-12)[0]

    total_area = sum(area_of.get(pid, 0.0) for pid in omega)
    msq = sum(w @ J_of[pid] @ w for pid, w in omega.items() if pid in J_of)
    rms_unscaled = EARTH_RADIUS_KM * np.sqrt(msq / total_area) if total_area > 0 else 0.0

    if fixed_k is not None and speed_mode != 'per_step':
        k = float(fixed_k)
    elif speed_mode == 'none':
        k = 1.0
    elif speed_mode == 'per_step' or speed_state.get('k') is None:
        k = target / rms_unscaled if rms_unscaled > 0 else 1.0
    else:
        k = speed_state['k']
    if speed_mode == 'calibrate_once' and speed_state.get('k') is None:
        speed_state['k'] = k
        print('  [speed] calibrated scale k = {:.4g} (RMS {:.1f} -> {:.1f} km/Myr)'.format(
            k, rms_unscaled, k * rms_unscaled))
    k_used, capped = k, False
    if k * rms_unscaled > cap > 0:
        k_used, capped = cap / rms_unscaled, True
        speed_state['capped_steps'] = speed_state.get('capped_steps', 0) + 1
        print('  [speed cap] RMS {:.1f} km/Myr > cap {:.1f}; scaled down this step'.format(
            k * rms_unscaled, cap))
    speed_state['last'] = dict(k_used=k_used, rms_unscaled=rms_unscaled,
                               rms=k_used * rms_unscaled, capped=capped)

    return {pid: pdriver.torque_to_stage_rotation(pygplates.Vector3D(*(k_used * w)), delta_time)
            for pid, w in omega.items()}


def apply_collision_coupling(stage_pole_dict, collision_pairs,
                              collision_state, plate_area_dict, params):
    """Blend stage poles of actively colliding plate pairs toward each other.

    Simulates the progressive suturing of continents (India–Asia style): the
    longer two cratons have been in contact, the more their velocities are
    coupled.  Coupling is area-weighted so the smaller plate adjusts more.

    collision_state is mutated in-place: ages for active pairs are incremented,
    pairs no longer colliding are pruned.
    """
    # Update collision ages: increment active pairs, prune inactive ones
    for pair in list(collision_state.keys()):
        if pair not in collision_pairs:
            del collision_state[pair]
    for pair in collision_pairs:
        collision_state[pair] = collision_state.get(pair, 0) + 1

    coupling_max  = params.get('collision_coupling_max', 0.0)
    timescale     = params.get('collision_coupling_timescale', 30)
    if coupling_max <= 0:
        return

    def rot_to_vec(rot):
        axis, angle = rotation_axis_angle(rot)
        return axis * angle                       # radians

    def vec_to_rot(v):
        ang = float(np.linalg.norm(v))            # radians
        if ang < 1e-12:
            return pygplates.FiniteRotation()
        lon, lat = xyz_to_lonlat(v / ang)
        return pygplates.FiniteRotation((float(lat), float(lon)), ang)

    for (pid_lo, pid_hi), age in collision_state.items():
        if pid_lo not in stage_pole_dict or pid_hi not in stage_pole_dict:
            continue
        alpha    = min(coupling_max, age / timescale * coupling_max)
        area_lo  = plate_area_dict.get(pid_lo, 1.0)
        area_hi  = plate_area_dict.get(pid_hi, 1.0)
        total    = area_lo + area_hi
        vA       = rot_to_vec(stage_pole_dict[pid_lo])
        vB       = rot_to_vec(stage_pole_dict[pid_hi])
        v_target = (area_lo * vA + area_hi * vB) / total
        stage_pole_dict[pid_lo] = vec_to_rot(vA + alpha * (v_target - vA))
        stage_pole_dict[pid_hi] = vec_to_rot(vB + alpha * (v_target - vB))
        print(f'  [collision] pair ({pid_lo},{pid_hi}) age={age} coupling={alpha:.3f}')


def _neighbour_ids(vals):
    """{plate_id: set of 4-connected neighbour plate IDs}, longitude wrapping."""
    ok = np.isfinite(vals) & (vals > 0)
    filled = np.where(ok, vals, 0)
    neighbours = {}
    for shifted in wrap_neighbours(filled):
        m = ok & (shifted > 0) & (shifted != filled)
        if not m.any():
            continue
        pairs = np.unique(np.column_stack([filled[m], shifted[m]]).astype(int), axis=0)
        for a, b in pairs:
            neighbours.setdefault(int(a), set()).add(int(b))
    return neighbours


def detect_enclosed_plates(plate_id_raster):
    """Find plates that are completely surrounded by a single, larger plate.

    Neighbours are found with longitude wrapping, so a plate straddling ±180 is
    treated like any other.  Plates that reach a pole are never considered
    enclosed (a polar cap bordered by one plate is a legitimate configuration).

    Returns a list of (enclosing_pid, enclosed_pid) pairs suitable for
    apply_plate_merges.
    """
    vals = plate_id_raster.values
    at_pole = set(int(v) for v in np.concatenate([vals[0], vals[-1]]) if np.isfinite(v) and v > 0)
    areas = plate_areas(plate_id_raster)

    merge_pairs = []
    for pid, nbrs in _neighbour_ids(vals).items():
        if len(nbrs) != 1 or pid in at_pole:
            continue
        enclosing = next(iter(nbrs))
        if areas.get(enclosing, 0) <= areas.get(pid, 0):
            continue
        print(f'  [topology] Plate {pid} is enclosed by plate {enclosing} — absorbing')
        merge_pairs.append((enclosing, pid))

    return merge_pairs


def apply_plate_merges(plate_id_raster, merge_pairs):
    """Replace pixels of smaller plate IDs with larger plate IDs.
    merge_pairs: iterable of (pid_large, pid_small).
    """
    if not merge_pairs:
        return plate_id_raster
    vals = plate_id_raster.values.copy()
    for pid_large, pid_small in merge_pairs:
        mask = vals == pid_small
        if mask.any():
            vals[mask] = pid_large
            print(f'  [merge] Plate {pid_small} absorbed into {pid_large} ({int(mask.sum())} px)')
    return xr.DataArray(vals, dims=plate_id_raster.dims, coords=plate_id_raster.coords)


def split_disconnected_plates(plate_id_raster, craton_raster=None, next_pid=None,
                              min_fragment_area_deg2=5.0):
    """Give every disconnected piece of a plate its own plate ID.

    Overlap resolution or the craton override can cut a plate in two.  The raster is
    the reference, so the pieces are found by connected-component labelling (joining
    across ±180).  For each plate ID:
      - pieces containing that plate's cratons, or else the largest piece, keep the ID;
      - pieces smaller than min_fragment_area_deg2 are absorbed into the neighbouring
        plate they share most boundary with;
      - other pieces become new plates (caller copies the parent's rotation).

    Areas are spherical (cos(lat)-weighted), in equatorial-equivalent deg².
    Returns (raster, [(new_pid, parent_pid), ...], next_pid).
    """
    lons, lats = raster_coords(plate_id_raster)
    vals = plate_id_raster.values.copy()
    if next_pid is None:
        next_pid = int(np.nanmax(vals)) + 1
    cell_area = cell_area_weights(lats, vals.shape) * abs(lons[1] - lons[0]) * abs(lats[1] - lats[0])
    craton_vals = None if craton_raster is None else np.nan_to_num(craton_raster.values)

    new_plates = []
    for pid in np.unique(vals[np.isfinite(vals) & (vals > 0)]).astype(int):
        labels, n = label_wrapped(vals == pid)
        if n < 2:
            continue
        comp_area = np.bincount(labels.ravel(), weights=cell_area.ravel(), minlength=n + 1)[1:]
        if craton_vals is not None:
            craton_hits = np.bincount(labels[craton_vals == pid], minlength=n + 1)[1:]
        else:
            craton_hits = np.zeros(n, dtype=int)
        keep = set(np.where(craton_hits > 0)[0])
        if not keep:
            keep = {int(np.argmax(comp_area))}

        for k in range(n):
            if k in keep:
                continue
            comp = labels == k + 1
            if comp_area[k] < min_fragment_area_deg2:
                ring = np.zeros_like(comp)
                for shifted in wrap_neighbours(comp):
                    ring |= shifted
                ring &= ~comp
                others = vals[ring]
                others = others[np.isfinite(others) & (others != pid) & (others > 0)]
                if len(others):
                    ids, counts = np.unique(others.astype(int), return_counts=True)
                    vals[comp] = ids[np.argmax(counts)]
                    continue
            vals[comp] = next_pid
            new_plates.append((next_pid, int(pid)))
            print(f'  [topology] Disconnected piece of plate {pid} → new plate {next_pid} '
                  f'({comp_area[k]:.0f} deg²)')
            next_pid += 1

    raster = xr.DataArray(vals, dims=plate_id_raster.dims, coords=plate_id_raster.coords)
    return raster, new_plates, next_pid


def _plate_shape_metrics(mask, lon2d, lat2d, cell_area):
    """Concavity (area / convex-hull area) and aspect ratio of a plate, from its raster
    cells and measured on the sphere.

    Aspect ratio: square root of the ratio of the principal second moments of the plate
    area in an azimuthal-equidistant frame about the plate centroid (a/b for an a x b
    rectangle or ellipse).  Convex hull: taken in a gnomonic projection about the
    centroid, where great circles are straight lines, so the planar hull of the edge
    cells is the spherical convex hull; its area comes from pygplates.

    A plate reaching more than 85 degrees from its centroid has no gnomonic hull; it is
    returned as concavity 0 so that the neck search always examines it.
    """
    w = cell_area[mask]
    xyz = lonlat_to_xyz(lon2d[mask], lat2d[mask])
    centre = (xyz * w[:, None]).sum(axis=0)
    centre /= np.linalg.norm(centre)
    # Tangent-plane basis at the centroid (east, north)
    east = np.cross([0.0, 0.0, 1.0], centre)
    if np.linalg.norm(east) < 1e-9:
        east = np.array([1.0, 0.0, 0.0])
    east /= np.linalg.norm(east)
    north = np.cross(centre, east)

    cos_d = np.clip(xyz @ centre, -1.0, 1.0)
    x, y = xyz @ east, xyz @ north
    # Azimuthal-equidistant coordinates: angular distance along the azimuth
    dist = np.arccos(cos_d)
    norm_xy = np.hypot(x, y)
    scale = np.divide(dist, norm_xy, out=np.ones_like(dist), where=norm_xy > 1e-12)
    ax, ay = x * scale, y * scale
    wsum = w.sum()
    mx, my = (w * ax).sum() / wsum, (w * ay).sum() / wsum
    cxx = (w * (ax - mx) ** 2).sum() / wsum
    cyy = (w * (ay - my) ** 2).sum() / wsum
    cxy = (w * (ax - mx) * (ay - my)).sum() / wsum
    eig = np.linalg.eigvalsh([[cxx, cxy], [cxy, cyy]])
    aspect = float(np.sqrt(eig[1] / eig[0])) if eig[0] > 0 else np.inf

    if cos_d.min() < np.cos(np.radians(85.0)):
        return 0.0, aspect

    # Spherical convex hull of the edge cells, via the gnomonic projection
    from scipy.spatial import ConvexHull, QhullError
    up, down, east_n, west_n = wrap_neighbours(mask)
    edge = (mask & ~(up & down & east_n & west_n))[mask]
    gx, gy = x[edge] / cos_d[edge], y[edge] / cos_d[edge]
    if len(gx) < 3:
        return 1.0, aspect
    try:
        hull = ConvexHull(np.column_stack([gx, gy]))
    except QhullError:
        return 1.0, aspect
    verts = centre + gx[hull.vertices, None] * east + gy[hull.vertices, None] * north
    hull_lon, hull_lat = xyz_to_lonlat(verts)
    hull_area = pygplates.PolygonOnSphere(list(zip(hull_lat, hull_lon))).get_area()  # steradians
    hull_area *= np.degrees(1.0) ** 2   # equatorial deg², the units of cell_area
    if hull_area <= 0:
        return 1.0, aspect
    return float(min(1.0, wsum / hull_area)), aspect


def _find_neck(mask, lon2d, lat2d, cell_area, min_part_area):
    """Look for a neck in a plate by spherical morphological erosion.

    The plate is eroded by increasing great-circle distances r from its edge.  The
    first r at which the eroded core falls apart into two pieces that each carry at
    least min_part_area of the plate marks a neck of half-width ~r.  Every plate
    pixel is then assigned to the nearer of those two cores, which puts the cut
    through the neck.

    Returns None or dict(part_a, part_b, neck_half_width_deg, neck_ratio), where
    neck_ratio = neck half-width / inscribed radius of the smaller part.
    """
    up, down, east, west = wrap_neighbours(mask)
    outer = ~mask & (up | down | east | west)
    if not outer.any():
        return None
    tree = SphereTree(np.radians(np.column_stack([lat2d[outer], lon2d[outer]])))
    inner_pts = np.radians(np.column_stack([lat2d[mask], lon2d[mask]]))
    dist = np.zeros(mask.shape)
    dist[mask] = np.degrees(tree.query(inner_pts, k=1)[0][:, 0])
    body_radius = dist.max()

    step = abs(lat2d[1, 0] - lat2d[0, 0])
    r = step
    while r < body_radius:
        core = mask & (dist > r)
        labels, n = label_wrapped(core, connectivity=2)
        if n == 0:
            return None
        if n >= 2:
            core_area = np.bincount(labels.ravel(), weights=cell_area.ravel(), minlength=n + 1)[1:]
            biggest = np.argsort(core_area)[::-1][:2] + 1
            seeds = [labels == b for b in biggest]
            seed_pts = np.vstack([np.column_stack([lat2d[s], lon2d[s]]) for s in seeds])
            seed_lab = np.concatenate([np.full(int(s.sum()), j) for j, s in enumerate(seeds)])
            seed_tree = SphereTree(np.radians(seed_pts))
            nearest = seed_lab[seed_tree.query(inner_pts, k=1)[1][:, 0]]
            part_a = np.zeros(mask.shape, bool)
            part_a[mask] = nearest == 0
            part_b = mask & ~part_a
            area_a, area_b = cell_area[part_a].sum(), cell_area[part_b].sum()
            if min(area_a, area_b) >= min_part_area:
                ref = min(dist[part_a].max(), dist[part_b].max())
                return dict(part_a=part_a, part_b=part_b, neck_half_width_deg=r,
                            neck_ratio=r / max(ref, 1e-9))
        r *= 1.25
    return None


def detect_and_split_plates(plate_id_raster, new_plate_polygons, params,
                             craton_raster=None, cooldown_pids=None, next_pid=None):
    """Detect geometric necks and split implausible plates in the raster.

    1. Screening: plates whose concavity is below min_concavity_ratio or whose aspect
       ratio exceeds max_aspect_ratio (both measured on the sphere from the plate's raster
       cells) are examined further.
    2. Neck search (_find_neck): spherical erosion finds the narrowest constriction that
       separates two substantial parts.  The plate is cut there if
       neck half-width / inscribed radius of the smaller part < min_neck_ratio.
    3. Cratons stay with the parent plate ID: if only the smaller part holds the
       plate's cratons the parts swap roles; if both do, the plate is not cut.

    The analysis runs on a ~1° version of the raster for speed; the final cut is made
    at full resolution.  Returns (updated_raster, split_info, next_pid) where split_info
    is a list of (new_pid, parent_pid) for rotation inheritance by the caller.
    """
    if cooldown_pids is None:
        cooldown_pids = set()

    min_concavity    = params.get('min_concavity_ratio',     0.85)
    max_aspect       = params.get('max_aspect_ratio',        4.0)
    min_neck_ratio   = params.get('min_neck_ratio',          0.35)
    min_child_area   = params.get('min_child_area_deg2',     25.0)
    min_peninsula    = params.get('min_peninsula_area_deg2', 5.0)

    lons, lats = raster_coords(plate_id_raster)
    vals = plate_id_raster.values.copy()
    if next_pid is None:
        next_pid = int(np.nanmax(vals)) + 1
    rsp = abs(lons[1] - lons[0])
    lon2d, lat2d = np.meshgrid(lons, lats)
    cell_area = cell_area_weights(lats, vals.shape) * rsp * rsp
    craton_vals = None if craton_raster is None else np.nan_to_num(craton_raster.values)

    # Coarse (~1°) copy of the grid for the neck search
    k = max(1, int(round(1.0 / rsp)))
    sl = (slice(None, None, k), slice(None, None, k))

    split_info = []
    for pid, group in new_plate_polygons.groupby('PLATEID1'):
        pid = int(pid)
        if pid in cooldown_pids:
            continue
        mask = vals == pid
        area = cell_area[mask].sum()
        if area < 2 * min_peninsula:
            continue
        concavity, aspect = _plate_shape_metrics(mask, lon2d, lat2d, cell_area)
        if concavity >= min_concavity and aspect <= max_aspect:
            continue

        neck = _find_neck(mask[sl], lon2d[sl], lat2d[sl], cell_area[sl] * k * k, min_peninsula)
        if neck is None or neck['neck_ratio'] >= min_neck_ratio:
            continue
        # Thin necks (genuine peninsulas) may sever a fragment smaller than the normal minimum
        min_area = min_peninsula if neck['neck_ratio'] < 0.2 else min_child_area

        # Full-resolution cut: assign each plate pixel to the nearer coarse part
        coarse_pts = np.column_stack([lat2d[sl][neck['part_a'] | neck['part_b']],
                                      lon2d[sl][neck['part_a'] | neck['part_b']]])
        coarse_lab = neck['part_b'][neck['part_a'] | neck['part_b']].astype(int)
        tree = SphereTree(np.radians(coarse_pts))
        nearest = coarse_lab[tree.query(np.radians(np.column_stack([lat2d[mask], lon2d[mask]])),
                                        k=1)[1][:, 0]]
        part_b = np.zeros(mask.shape, bool)
        part_b[mask] = nearest == 1
        part_a = mask & ~part_b
        if cell_area[part_a].sum() < cell_area[part_b].sum():
            part_a, part_b = part_b, part_a
        if cell_area[part_b].sum() < min_area:
            continue

        # Cratons keep their plate ID (craton features carry the parent's ID)
        child = part_b
        if craton_vals is not None:
            cr_a = ((craton_vals == pid) & part_a).any()
            cr_b = ((craton_vals == pid) & part_b).any()
            if cr_a and cr_b:
                continue
            if cr_b:
                child = part_a

        vals[child] = next_pid
        split_info.append((next_pid, pid))
        print(f'  [fragment] Plate {pid} → child {next_pid} '
              f'(concavity={concavity:.2f} aspect={aspect:.1f} '
              f'neck_ratio={neck["neck_ratio"]:.2f} neck_half_width={neck["neck_half_width_deg"]:.1f}° '
              f'child_area={cell_area[child].sum():.0f} deg²)')
        next_pid += 1

    raster = xr.DataArray(vals, dims=plate_id_raster.dims, coords=plate_id_raster.coords)
    return raster, split_info, next_pid


def compute_timestep(rpolygons, forces, stage_pole_dict, delta_time=1, cratons=None, params=None,
                     polarity_dict=None, merge_pairs=None, fragmentation_cooldown_pids=None,
                     collision_state=None, craton_id_remap=None, next_plate_id=None,
                     speed_state=None, previous_plate_id=None, force_fields=None,
                     force_diag=None):
    """Rebuild the plate configuration from moved plate polygons and compute the next
    stage rotations.

    force_fields feeds the force modes that need them (see force_fields_needed):
      'tracers'       (lon, lat, age, plate_id) seafloor tracers at this step's positions
                      (plate_id: the plate each was on before the move, or None)
      'slab_tracers'  the slab tracer DataFrame
    The seafloor-age grid is built from the tracers only if a mode reads it.
    force_diag, if a dict, receives the per-segment slab-force table ('segments').

    stage_pole_dict holds the rotations that moved the plates to this configuration;
    they are used to classify boundaries and for the velocity fields.  Returns

        (plate polygons, new stage_pole_dict, rasters Dataset, boundaries_gdf,
         new_plates, collision_state, next_plate_id)

    where new_plates lists (new_pid, parent_pid) for every plate created this step
    (fragmentation or disconnected pieces), so the caller can copy the parent's
    rotation history.  New plate IDs start at next_plate_id and are never reused.
    """
    if params is None:
        params = {}

    tessellation_deg          = params.get('tessellation_deg',           0.2)
    raster_spacing            = params.get('raster_spacing',             1)
    topo_simplify             = params.get('topo_simplify_tolerance',    1)
    div_thresh_sub            = params.get('div_threshold_subduction',  -1.0)
    div_thresh_ridge          = params.get('div_threshold_ridge',        1.0)
    boundary_search_radius    = params.get('boundary_search_radius_deg', 2)
    collision_mode            = params.get('collision_mode',             'resistance')
    min_fragment_area         = params.get('min_peninsula_area_deg2',    5.0)
    drag_mode                 = params.get('drag_mode',                  'none')
    speed_mode                = params.get('speed_normalisation',        'none')
    # The simple model (ω = τ, no scaling) keeps the original code path exactly
    solve_motions = not (drag_mode == 'none' and speed_mode == 'none'
                         and params.get('speed_scale') is None)
    if speed_state is None:
        speed_state = {}
    force_modes_used = {'slab_pull':    params.get('slab_pull_mode', 'hard'),
                        'slab_suction': params.get('slab_suction_mode', 'ratio'),
                        'polarity':     params.get('polarity_mode', 'id')}
    needs = force_fields_needed(params)
    force_fields = force_fields or {}
    absorbed_into = {}      # {absorbed pid: host pid} this step, to keep tracers on their plate

    # Local copy: children of split plates are added below
    stage_pole_dict = dict(stage_pole_dict)

    _ct = {}
    _t0 = time.perf_counter()

    new_features = tessellate_polygons(rpolygons, max_dist_deg=tessellation_deg)
    _ct['tessellate'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    rp_gdf = plate_polygons_to_wrapped_gdf(new_features)
    _ct['wrap_polygons'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    clipped_polygons = resolve_overlaps(rp_gdf, priority_col='PLATEID1')
    clipped_polygons = clipped_polygons.explode(index_parts=False).reset_index(drop=True)
    _ct['resolve_overlaps'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    plate_id_raster = update_plate_id_raster(clipped_polygons, raster_spacing=raster_spacing,
                                             tie_break=params.get('tie_break', 'previous_owner'),
                                             previous_plate_id=previous_plate_id,
                                             tie_break_seed=params.get('tie_break_seed', 0))
    if next_plate_id is None:
        next_plate_id = int(np.nanmax(plate_id_raster.values)) + 1
    _ct['plate_id_raster'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    def _record_absorbed(pairs):
        # Cratons of an absorbed plate keep their original PLATEID1; remap them to
        # the host so that (a) the craton override uses the host's ID and (b) the
        # craton moves with the host via nudge_polygons in run_SDP.py.
        for host, absorbed in pairs:
            absorbed_into[absorbed] = host
        if craton_id_remap is None:
            return
        for host, absorbed in pairs:
            craton_id_remap[absorbed] = host
            for k in list(craton_id_remap):          # redirect chains
                if craton_id_remap[k] == absorbed:
                    craton_id_remap[k] = host

    # Force values within cratons to overrule any plate polygon changes
    craton_raster = None
    if cratons is not None:
        wrapped_cratons = plate_polygons_to_wrapped_gdf(cratons)
        valid_cratons = wrapped_cratons.dropna(subset=['geometry'])
        valid_cratons = valid_cratons[~valid_cratons.geometry.is_empty]
        craton_raster = rasterize_polygons(valid_cratons, zval_field='PLATEID1', sampling=raster_spacing)

    def _remapped_cratons():
        if craton_raster is None or not craton_id_remap:
            return craton_raster
        cr_vals = craton_raster.values.copy()
        for old_id, new_id in craton_id_remap.items():
            cr_vals[craton_raster.values == old_id] = new_id
        return xr.DataArray(cr_vals, dims=craton_raster.dims, coords=craton_raster.coords)

    # Apply merges requested by caller (plates whose relative motion dropped to zero)
    if merge_pairs:
        plate_id_raster = apply_plate_merges(plate_id_raster, merge_pairs)
        _record_absorbed(merge_pairs)

    craton_raster = _remapped_cratons()
    if craton_raster is not None:
        plate_id_raster = craton_raster.where(craton_raster > 0, plate_id_raster)
    _ct['craton_raster'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    new_plates = []

    def _inherit(pairs):
        for new_pid, parent_pid in pairs:
            if parent_pid in stage_pole_dict:
                stage_pole_dict[new_pid] = stage_pole_dict[parent_pid]
        new_plates.extend(pairs)

    # Disconnected pieces of one plate ID become separate plates (or are absorbed)
    plate_id_raster, _pieces, next_plate_id = split_disconnected_plates(
        plate_id_raster, craton_raster, next_plate_id, min_fragment_area)
    _inherit(_pieces)
    _ct['split_disconnected'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    def _rebuild_polygons(raster):
        polys = update_plate_polygons(raster, raster_spacing=raster_spacing,
                                      topo_simplify=topo_simplify)
        polys = polys.explode(index_parts=False).reset_index(drop=True)
        return polys[polys.geometry.geom_type.isin(['Polygon'])]

    new_plate_polygons = _rebuild_polygons(plate_id_raster)
    _ct['update_polygons'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    # Plate fragmentation: split geometrically implausible plates at their neck
    if params.get('fragmentation_enabled', False):
        plate_id_raster, _split_info, next_plate_id = detect_and_split_plates(
            plate_id_raster, new_plate_polygons, params,
            craton_raster=craton_raster,
            cooldown_pids=fragmentation_cooldown_pids or set(),
            next_pid=next_plate_id)
        if _split_info:
            _inherit(_split_info)
            new_plate_polygons = _rebuild_polygons(plate_id_raster)
        _ct['fragment'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    # Topology repair: absorb any plate that is completely enclosed within
    # a single other plate.  This can occur when a neck closes (two sides of
    # a peninsula fuse), trapping a small plate inside a larger one.  Such
    # configurations are physically invalid and produce spurious internal
    # boundaries with large divergence.
    _enclosed = detect_enclosed_plates(plate_id_raster)
    if _enclosed:
        plate_id_raster    = apply_plate_merges(plate_id_raster, _enclosed)
        new_plate_polygons = _rebuild_polygons(plate_id_raster)
        _record_absorbed(_enclosed)
        craton_raster = _remapped_cratons()
        _ct['enclosed_repair'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    boundaries_gdf = topo_polygons_to_boundaries(new_plate_polygons,
                                                  raster_spacing=raster_spacing)
    _ct['topo_boundaries'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    boundaries_gdf = compute_boundary_divergence(boundaries_gdf, stage_pole_dict, delta_time, plate_id_raster,
                                                 plate_polygons_gdf=new_plate_polygons)
    _ct['boundary_divergence'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    velocity_grid_east, velocity_grid_north = velocity_grid(plate_id_raster, stage_pole_dict, delta_time,
                                                             raster_spacing=raster_spacing)
    div = spherical_divergence(velocity_grid_east, velocity_grid_north)
    _ct['velocity_grid'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    # Fields for the force modes, built only if a selected mode reads them
    age_grid = None
    if needs['age'] and force_fields.get('tracers') is not None:
        t_lon, t_lat, t_age, t_pid = force_fields['tracers']
        if t_pid is not None:
            # a tracer on a plate absorbed this step now belongs to the host; one on a
            # plate split off this step still belongs to the parent
            for _ in range(len(absorbed_into)):
                t_pid = map_plate_ids(t_pid, absorbed_into)
        age_grid = physics_grid_age(t_lon, t_lat, t_age, plate_id_raster, craton_raster,
                                    grid_deg=params.get('physics_grid_deg', 1.0),
                                    tracer_plate_id=t_pid,
                                    raster_id_map={new: parent for new, parent in new_plates})
        _ct['age_grid'] = time.perf_counter() - _t0; _t0 = time.perf_counter()
    slab_tracers = None
    if needs['slab'] and force_fields.get('slab_tracers') is not None:
        slab_tracers = force_fields['slab_tracers']
        slab_tracers = slab_tracers[slab_tracers['active'].astype(bool)]
    force_params = {k: params[k] for k in ('slab_ref_depth_km', 'slab_age_saturation_myr',
                                           'pull_ramp_km_per_myr', 'slab_search_radius_deg',
                                           'slab_line_tol_deg') if k in params}

    _stage_out = polygons_to_stage_poles(
                                              new_plate_polygons, plate_id_raster, forces,
                                              craton_raster,
                                              stage_pole_dict=stage_pole_dict,
                                              delta_time=delta_time,
                                              divergence_threshold_for_subduction=div_thresh_sub,
                                              divergence_threshold_for_ridge_push=div_thresh_ridge,
                                              boundary_search_radius_deg=boundary_search_radius,
                                              polarity_dict=polarity_dict,
                                              collision_mode=collision_mode,
                                              return_torques=solve_motions,
                                              modes=force_modes_used,
                                              age_grid=age_grid,
                                              slab_tracers=slab_tracers,
                                              force_params=force_params,
                                              diagnostics=force_diag)
    if solve_motions:
        _, _collision_pairs, _torques, _edge_tensors = _stage_out
        new_stage_pole_dict = solve_plate_motions(_torques, _edge_tensors, plate_id_raster,
                                                  craton_raster, params, drag_mode, speed_mode,
                                                  speed_state, delta_time)
    else:
        new_stage_pole_dict, _collision_pairs = _stage_out

    if collision_mode == 'coupling' and collision_state is not None:
        apply_collision_coupling(new_stage_pole_dict, _collision_pairs,
                                 collision_state, plate_areas(plate_id_raster), params)

    _ct['stage_poles'] = time.perf_counter() - _t0

    _ct_total = sum(_ct.values())
    print('  [compute_timestep {:.1f}s]'.format(_ct_total),
          '  '.join('{}: {:.1f}s ({:.0f}%)'.format(k, v, 100*v/_ct_total)
                    for k, v in sorted(_ct.items(), key=lambda x: -x[1])))

    rasters_dict = {
        'plate_id':        plate_id_raster.rename({'x': 'lon', 'y': 'lat'}),
        'velocity_east':   velocity_grid_east,
        'velocity_north':  velocity_grid_north,
        'divergence':      div,
    }
    if craton_raster is not None:
        rasters_dict['craton'] = craton_raster.rename({'x': 'lon', 'y': 'lat'})
    if force_diag is not None and age_grid is not None:
        force_diag['age_grid'] = age_grid

    return new_plate_polygons, \
           new_stage_pole_dict, \
           xr.Dataset(rasters_dict), \
           boundaries_gdf, \
           new_plates, \
           collision_state if collision_state is not None else {}, \
           next_plate_id


def _segment_arcs(boundaries_gdf, mask=None):
    """Yield (row_position, GreatCircleArc) for every 2-vertex piece of the boundary
    LineStrings (optionally only rows where mask is True)."""
    for pos in range(len(boundaries_gdf)):
        if mask is not None and not mask[pos]:
            continue
        coords = list(boundaries_gdf.geometry.iloc[pos].coords)
        for i in range(len(coords) - 1):
            try:
                arc = pygplates.GreatCircleArc(pygplates.PointOnSphere(coords[i][1],   coords[i][0]),
                                               pygplates.PointOnSphere(coords[i+1][1], coords[i+1][0]))
                arc.get_great_circle_normal()
            except pygplates.IndeterminateGreatCircleArcNormalError:
                continue
            yield pos, arc


def _sample_arcs(boundaries_gdf, mask, sample_deg):
    """Points every ~sample_deg along boundary rows where mask is True: each 2-vertex
    piece of a row is split into n = ceil(length / sample_deg) equal great-circle steps
    and its n + 1 points (both ends included) are returned.  Pieces of zero length (or
    antipodal ends) are skipped, as for pygplates.GreatCircleArc.
    Returns lat (N,), lon (N,), row position (N,)."""
    empty = np.array([]), np.array([]), np.array([], dtype=int)
    sel = np.where(np.asarray(mask, bool))[0] if mask is not None else np.arange(len(boundaries_gdf))
    if not len(sel):
        return empty
    coords = boundaries_gdf.geometry.iloc[sel].get_coordinates()
    row = sel[np.asarray(pd.factorize(coords.index)[0])] if len(coords) else np.array([], int)
    xyz = lonlat_to_xyz(coords.x.values, coords.y.values)
    same = row[1:] == row[:-1]
    a, b, seg_row = xyz[:-1][same], xyz[1:][same], row[1:][same]
    cross = np.linalg.norm(np.cross(a, b), axis=1)
    theta = np.arctan2(cross, np.einsum('ij,ij->i', a, b))
    ok = cross > 1e-12
    a, b, seg_row, theta = a[ok], b[ok], seg_row[ok], theta[ok]
    if not len(a):
        return empty
    # step count from acos(a·b), the way pygplates measures arc length, so that pieces
    # whose length is a whole number of steps are split exactly as before
    length_deg = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', a, b), -1.0, 1.0)))
    n = np.maximum(np.ceil(length_deg / sample_deg).astype(int), 1)
    seg = np.repeat(np.arange(len(a)), n + 1)
    f = np.concatenate([np.linspace(0.0, 1.0, k + 1) for k in n])
    th = theta[seg]
    # great-circle (spherical linear) interpolation between the segment ends
    p = (np.sin((1.0 - f) * th)[:, None] * a[seg] + np.sin(f * th)[:, None] * b[seg]) / np.sin(th)[:, None]
    lon, lat = xyz_to_lonlat(p)
    return lat, lon, seg_row[seg]


def _points_every(boundaries_gdf, mask, spacing_deg):
    """Points spaced spacing_deg apart (great-circle) along boundary rows where mask
    is True, independent of how the boundaries are split into segments: the distance
    left over at the end of one segment carries into the next.
    Returns lat, lon, row position and the arc normal (N, 3) at each point."""
    lats, lons, rows, normals = [], [], [], []
    carry = 0.5 * spacing_deg
    for pos, arc in _segment_arcs(boundaries_gdf, mask):
        length = np.degrees(arc.get_arc_length())
        normal = np.array(arc.get_great_circle_normal().to_xyz())
        d = spacing_deg - carry
        while d <= length:
            la, lo = arc.get_arc_point(d / length).to_lat_lon()
            lats.append(la); lons.append(lo); rows.append(pos); normals.append(normal)
            d += spacing_deg
        carry = length - (d - spacing_deg)
    return (np.array(lats), np.array(lons), np.array(rows, dtype=int),
            np.array(normals).reshape(-1, 3))


def seed_ridge_tracers(boundaries_gdf, spacing_deg=2.0,
                       threshold=1.0, offset_deg=0.6):
    """Seed tracer point pairs offset perpendicular to divergent boundary segments.

    Two points are seeded per position (one on each side of the ridge axis) so
    that after advection they are carried apart on their respective plates.
    Positions are spaced spacing_deg along the ridge (great-circle distance).
    offset_deg must exceed half the plate_id raster cell size so each point clearly
    falls on one plate.  boundaries_gdf must have a 'divergence' column.
    """
    empty = gpd.GeoDataFrame(columns=['geometry', 'active', 'age'], geometry='geometry').set_crs(epsg=4326)
    if len(boundaries_gdf) == 0 or 'divergence' not in boundaries_gdf.columns:
        return empty

    mask = (boundaries_gdf['divergence'] > threshold).values
    lat, lon, _, normal = _points_every(boundaries_gdf, mask, spacing_deg)
    seed_lon, seed_lat = [], []
    for sign in (+1, -1):
        lo, la = offset_along_normal(lon, lat, sign * normal, offset_deg)
        seed_lon.extend(np.atleast_1d(lo)); seed_lat.extend(np.atleast_1d(la))
    if not seed_lon:
        return empty
    return gpd.GeoDataFrame({'active': True, 'age': 0.0},
                            index=range(len(seed_lon)),
                            geometry=gpd.points_from_xy(seed_lon, seed_lat)).set_crs(epsg=4326)


def generate_seafloor_age_tracers(plate_polygons, boundaries_gdf,
                                   spacing_deg=5.0, spreading_rate_km_per_myr=50.0,
                                   div_threshold=0.00025):
    """Seed initial tracers inside each plate polygon with age estimated from ridge distance.

    For each polygon, finds its own divergent boundary segments, builds a BallTree from
    densified ridge points, then assigns age = great-circle distance (km) / spreading_rate
    to each interior grid point.  Points with no divergent boundary on their plate are skipped.
    boundaries_gdf must have a 'divergence' column (from compute_boundary_divergence).
    """
    from shapely import contains_xy

    all_lon, all_lat, all_age = [], [], []
    is_ridge = (boundaries_gdf['divergence'] > div_threshold).values if len(boundaries_gdf) else []

    for _, poly_row in plate_polygons.iterrows():
        pid = poly_row['PLATEID1']
        geom = poly_row.geometry
        if geom is None or geom.is_empty:
            continue

        on_plate = ((boundaries_gdf['plateid1_poly1'] == pid) |
                    (boundaries_gdf['plateid1_poly2'] == pid)).values & is_ridge
        if not on_plate.any():
            continue
        r_lat, r_lon, _ = _sample_arcs(boundaries_gdf, on_plate, spacing_deg / 2)
        if len(r_lat) == 0:
            continue
        tree = SphereTree(np.deg2rad(np.column_stack([r_lat, r_lon])))

        # Equal-area grid: uniform in (lon, sin(lat)) so density is ~constant on the sphere
        minx, miny, maxx, maxy = geom.bounds
        lons = np.arange(minx + spacing_deg / 2, maxx, spacing_deg)
        sin_step = np.sin(np.deg2rad(spacing_deg))
        sins = np.arange(np.sin(np.deg2rad(miny)) + sin_step / 2, np.sin(np.deg2rad(maxy)), sin_step)
        lats = np.degrees(np.arcsin(np.clip(sins, -1.0, 1.0)))
        lon2d, lat2d = np.meshgrid(lons, lats)
        inside = contains_xy(geom, lon2d.ravel(), lat2d.ravel())
        if not inside.any():
            continue
        q_lat, q_lon = lat2d.ravel()[inside], lon2d.ravel()[inside]
        dist_rad, _ = tree.query(np.deg2rad(np.column_stack([q_lat, q_lon])), k=1)
        all_lon.append(q_lon); all_lat.append(q_lat)
        all_age.append(dist_rad[:, 0] * EARTH_RADIUS_KM / spreading_rate_km_per_myr)

    if not all_lon:
        return gpd.GeoDataFrame(columns=['geometry', 'active', 'age'], geometry='geometry').set_crs(epsg=4326)
    return gpd.GeoDataFrame({'active': True, 'age': np.concatenate(all_age)},
                            geometry=gpd.points_from_xy(np.concatenate(all_lon),
                                                        np.concatenate(all_lat))).set_crs(epsg=4326)


def advect_tracers(tracers, plate_id_raster, stage_pole_dict, delta_time=1,
                   boundaries_gdf=None, div_threshold=0.0, new_plate_id_raster=None):
    """Move active tracers one step using their plate's stage rotation.

    plate_id_raster must describe the plates at the tracers' current (pre-move)
    positions; each tracer is rotated with the stage pole of the plate it sits on.

    If boundaries_gdf is provided, tracers that cross a convergent boundary during
    this step (plate ID at the new position in new_plate_id_raster differs from the
    old one, for a convergent pair) are deactivated.  This catches the case where a
    tracer slips into the overriding plate's territory before the proximity check fires.
    """
    if tracers is None or len(tracers) == 0:
        return tracers
    if new_plate_id_raster is None:
        new_plate_id_raster = plate_id_raster
    active_mask = tracers.active.astype(bool)
    active = tracers[active_mask].copy()
    inactive = tracers[~active_mask]

    convergent_pairs = set()
    if boundaries_gdf is not None and 'divergence' in boundaries_gdf.columns:
        conv = boundaries_gdf[boundaries_gdf['divergence'] < div_threshold]
        for p1, p2 in zip(conv['plateid1_poly1'], conv['plateid1_poly2']):
            if pd.notna(p1) and pd.notna(p2) and p1 != 0 and p2 != 0:
                convergent_pairs.add(frozenset([int(p1), int(p2)]))

    lons = active.geometry.x.values.copy()
    lats = active.geometry.y.values.copy()
    ages = active.age.values.copy().astype(float) + delta_time
    active_flags = active.active.values.copy().astype(bool)

    old_pids = raster_lookup(plate_id_raster, lons, lats)
    new_lons, new_lats = lons.copy(), lats.copy()
    for pid_val in np.unique(old_pids[np.isfinite(old_pids)]):
        rotation = stage_pole_dict.get(int(pid_val))
        if rotation is None:
            continue
        m = old_pids == pid_val
        new_lons[m], new_lats[m] = rotate_lonlat(rotation, lons[m], lats[m])

    if convergent_pairs:
        new_pids = raster_lookup(new_plate_id_raster, new_lons, new_lats)
        moved = (new_pids != old_pids) & (old_pids > 0) & (new_pids > 0)
        for k in np.where(moved)[0]:
            if frozenset([int(old_pids[k]), int(new_pids[k])]) in convergent_pairs:
                active_flags[k] = False

    updated = active.copy()
    updated['geometry'] = gpd.points_from_xy(new_lons, new_lats)
    updated['age'] = ages
    updated['active'] = active_flags
    return pd.concat([updated, inactive], ignore_index=True)


_CONVERGENT_CACHE = {}


def _convergent_samples(boundaries_gdf, threshold, sample_deg=0.5):
    """Points every sample_deg along convergent boundaries and a SphereTree on them.
    The last result is kept, so the seafloor-tracer and seamount calls in one step
    share it (keyed on the boundaries object itself and the threshold)."""
    c = _CONVERGENT_CACHE
    if c.get('gdf') is boundaries_gdf and c.get('key') == (threshold, sample_deg):
        return c['value']
    conv = (boundaries_gdf['divergence'] < threshold).values
    c_lat, c_lon, c_row = _sample_arcs(boundaries_gdf, conv, sample_deg)
    tree = SphereTree(np.deg2rad(np.column_stack([c_lat, c_lon]))) if len(c_lat) else None
    c.update(gdf=boundaries_gdf, key=(threshold, sample_deg), value=(c_lat, c_lon, c_row, tree))
    return c['value']


def deactivate_subducted_tracers(tracers, boundaries_gdf, threshold=-1.0, max_dist_deg=1.0,
                                 plate_id_raster=None, polarity_dict=None, craton_counts=None):
    """Mark active tracers as inactive when within max_dist_deg of a convergent boundary.

    Convergent boundaries are those with divergence < threshold in boundaries_gdf.
    Distance is great-circle distance to densely sampled boundary points.  If
    plate_id_raster is given, only tracers on the downgoing plate of the nearest
    convergent segment are removed (the overriding plate's seafloor is not subducted).
    """
    if tracers is None or len(tracers) == 0:
        return tracers
    if boundaries_gdf is None or len(boundaries_gdf) == 0:
        return tracers
    if 'divergence' not in boundaries_gdf.columns:
        return tracers

    tracers = tracers.copy()
    active_mask = tracers.active.astype(bool)
    if not active_mask.any():
        return tracers

    c_lat, c_lon, c_row, tree = _convergent_samples(boundaries_gdf, threshold)
    if len(c_lat) == 0:
        return tracers

    active_tracers = tracers[active_mask]
    t_lon = active_tracers.geometry.x.values
    t_lat = active_tracers.geometry.y.values
    dists_rad, nn = tree.query(np.deg2rad(np.column_stack([t_lat, t_lon])), k=1)
    deactivate_mask = np.degrees(dists_rad[:, 0]) < max_dist_deg

    if plate_id_raster is not None and deactivate_mask.any():
        on_pid = raster_lookup(plate_id_raster, t_lon, t_lat)
        p1 = boundaries_gdf['plateid1_poly1'].values
        p2 = boundaries_gdf['plateid1_poly2'].values
        for k in np.where(deactivate_mask)[0]:
            row = c_row[nn[k, 0]]
            dg = resolve_downgoing_plate(int(p1[row]), int(p2[row]), polarity_dict or {}, craton_counts)
            deactivate_mask[k] = on_pid[k] == dg

    tracers.loc[active_tracers.index[deactivate_mask], 'active'] = False
    return tracers


# Slablets from trenches converging at least this fast (km/Myr = mm/yr) get full weight;
# slower (near-stalled) boundaries contribute proportionally less.
SLAB_FULL_WEIGHT_CONVERGENCE_KM_PER_MYR = 10.0


def update_slab_tracers(slab_tracers, boundaries_gdf, delta_time,
                         div_threshold_subduction=-0.0001,
                         slowing_factor=4,
                         sample_deg=1.0,
                         upper_mantle_km=670.0,
                         cmb_depth_km=2890.0,
                         max_sink_rate_km_per_myr=200.0,
                         sink_rate_factor=1.0):
    """Seed new slablets at subduction zones and sink existing ones downward.

    slab_tracers is a plain DataFrame with columns:
        lon, lat, depth_km, age, active, sink_rate_km_per_myr, weight
    All sinking is vertical — lon/lat stay fixed after seeding.
    Sinking rate: sink_rate_factor × convergence rate (km/Myr) in the upper mantle,
    divided by slowing_factor in the lower mantle, capped at max_sink_rate_km_per_myr.
    Each slablet's weight is delta_time × min(1, |convergence| / 10 km/Myr): the
    mantle signal per Myr of subduction does not depend on the time step, and
    near-stalled boundaries contribute little.
    Slablets that reach the CMB are removed.
    """
    slab_tracers = slab_tracers.copy()

    if len(slab_tracers) > 0:
        active = slab_tracers['active'].astype(bool)
        upper = active & (slab_tracers['depth_km'] < upper_mantle_km)
        lower = active & (slab_tracers['depth_km'] >= upper_mantle_km)
        slab_tracers.loc[upper, 'depth_km'] += (
            slab_tracers.loc[upper, 'sink_rate_km_per_myr'] * delta_time)
        slab_tracers.loc[lower, 'depth_km'] += (
            slab_tracers.loc[lower, 'sink_rate_km_per_myr'] / slowing_factor * delta_time)
        slab_tracers.loc[active, 'age'] += delta_time
        slab_tracers = slab_tracers[slab_tracers['depth_km'] < cmb_depth_km].reset_index(drop=True)

    if boundaries_gdf is None or 'divergence' not in boundaries_gdf.columns:
        return slab_tracers

    sub = (boundaries_gdf['divergence'] < div_threshold_subduction).values
    s_lat, s_lon, s_row, _ = _points_every(boundaries_gdf, sub, sample_deg)
    if len(s_lat):
        convergence = np.abs(boundaries_gdf['divergence'].values[s_row].astype(float))
        new = pd.DataFrame({
            'lon': s_lon, 'lat': s_lat, 'depth_km': 0.0, 'age': 0.0, 'active': True,
            'sink_rate_km_per_myr': np.minimum(convergence * sink_rate_factor, max_sink_rate_km_per_myr),
            'weight': delta_time * np.minimum(1.0, convergence / SLAB_FULL_WEIGHT_CONVERGENCE_KM_PER_MYR)})
        slab_tracers = new if len(slab_tracers) == 0 else pd.concat([slab_tracers, new], ignore_index=True)
    return slab_tracers


def seed_initial_slabs(boundaries_gdf, depth_km, delta_time,
                       div_threshold_subduction=-0.0001,
                       sample_deg=1.0,
                       upper_mantle_km=670.0,
                       slowing_factor=4,
                       max_sink_rate_km_per_myr=200.0,
                       sink_rate_factor=1.0):
    """Slab tracers for the trenches of the starting configuration.

    Each convergent boundary (sampled every sample_deg, as in update_slab_tracers) gets
    the column of slablets the model itself would have built had it been subducting at
    its present rate until its slab reached depth_km: one slablet per delta_time of that
    history, sinking at the same rates as update_slab_tracers.  The sink rate is floored
    at SLAB_FULL_WEIGHT_CONVERGENCE_KM_PER_MYR so a near-stalled trench does not need
    thousands of Myr of history.  Columns are vertical (lon/lat fixed, mantle frame).
    Returns a DataFrame with the update_slab_tracers columns (empty if no trenches).
    """
    columns = ['lon', 'lat', 'depth_km', 'age', 'active', 'sink_rate_km_per_myr', 'weight']
    if (depth_km <= 0 or boundaries_gdf is None or len(boundaries_gdf) == 0
            or 'divergence' not in boundaries_gdf.columns):
        return pd.DataFrame(columns=columns)
    sub = (boundaries_gdf['divergence'] < div_threshold_subduction).values
    s_lat, s_lon, s_row, _ = _points_every(boundaries_gdf, sub, sample_deg)
    if not len(s_lat):
        return pd.DataFrame(columns=columns)
    convergence = np.abs(boundaries_gdf['divergence'].values[s_row].astype(float))
    rate = np.clip(convergence * sink_rate_factor, SLAB_FULL_WEIGHT_CONVERGENCE_KM_PER_MYR,
                   max_sink_rate_km_per_myr)
    # time for the slab to reach depth_km: upper mantle at rate, lower at rate / slowing_factor
    upper = min(depth_km, upper_mantle_km)
    t_end = upper / rate + max(0.0, depth_km - upper_mantle_km) * slowing_factor / rate
    n = np.ceil(t_end / delta_time).astype(int) + 1          # ages 0, dt, ..., >= t_end
    point = np.repeat(np.arange(len(s_lat)), n)
    age = (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)) * delta_time
    r = rate[point]
    t_um = upper_mantle_km / r
    depth = np.where(age <= t_um, r * age, upper_mantle_km + (age - t_um) * r / slowing_factor)
    return pd.DataFrame({
        'lon': s_lon[point], 'lat': s_lat[point],
        'depth_km': np.minimum(depth, depth_km), 'age': age, 'active': True,
        'sink_rate_km_per_myr': r,
        'weight': delta_time * np.minimum(1.0, convergence[point] / SLAB_FULL_WEIGHT_CONVERGENCE_KM_PER_MYR)})


def slab_depth_at_segments(mid_lon, mid_lat, side_xyz, slab_lon, slab_lat, slab_depth,
                           radius_deg=5.0, line_tol_deg=1.0):
    """Depth (km) of the deepest slab tracer near each boundary segment, on one side of it.

    side_xyz (n, 3) is the segment's great-circle normal, signed to point to the side
    where the slab lies (the overriding plate).  A tracer counts if it is within
    radius_deg (great-circle) of the segment midpoint and not more than line_tol_deg on
    the other side of the segment's great circle (slablets are seeded on the trench
    line, and the line moves).  Slablets sink vertically where they were seeded, so the
    deep part of a retreating trench's slab lies under the overriding plate, as a
    dipping slab does.  0 where there is no slab.
    """
    z = np.zeros(len(mid_lon))
    if len(slab_depth) == 0 or len(mid_lon) == 0:
        return z
    from scipy.spatial import cKDTree
    s_xyz = lonlat_to_xyz(np.asarray(slab_lon, float), np.asarray(slab_lat, float))
    m_xyz = lonlat_to_xyz(np.asarray(mid_lon, float), np.asarray(mid_lat, float))
    chord = 2.0 * np.sin(np.radians(radius_deg) / 2.0)
    pairs = cKDTree(m_xyz).sparse_distance_matrix(cKDTree(s_xyz), chord, output_type='ndarray')
    if len(pairs) == 0:
        return z
    seg, idx = pairs['i'], pairs['j']
    side = np.einsum('ij,ij->i', s_xyz[idx], np.asarray(side_xyz, float)[seg])   # sin(distance from the great circle), + on the slab side
    ok = side > -np.sin(np.radians(line_tol_deg))
    np.maximum.at(z, seg[ok], np.asarray(slab_depth, float)[idx[ok]])
    return z


def _mantle_grid(spacing_deg):
    lons_g = np.arange(-180 + spacing_deg / 2, 180, spacing_deg)
    lats_g = np.arange( -90 + spacing_deg / 2,  90, spacing_deg)
    return lons_g, lats_g


def _slab_weights(slab_tracers):
    if 'weight' in slab_tracers.columns:
        return slab_tracers['weight'].values.astype(float)
    return np.ones(len(slab_tracers))


def update_cmb_heat(cmb_heat, slab_tracers, delta_time,
                    spacing_deg=5.0,
                    cmb_depth_km=2890.0,
                    cmb_layer_km=400.0,
                    heating_rate=0.02,
                    slab_suppression=0.05,
                    surface_coupling=0.1,
                    sigma_deg=10.0,
                    max_heat=3.0):
    """Accumulate CMB heat each timestep, building LLSVP-like warm regions.

    Core heating increments the 2D field (nlat × nlon) everywhere each step;
    slablets suppress it with depth-weighted coupling.  A wide spherical
    Gaussian keeps the pattern large-scale.  Positive values = warm regions
    (slow dVs); passed to generate_mantle_temp_field as the upwelling anomaly.

    Two physical mechanisms are combined:
      - Thermal (material): slablets near the CMB suppress with full weight (1.0).
      - Dynamic pressure (incompressibility): surface slablets suppress with
        surface_coupling weight, representing instantaneous Stokes-flow coupling.
    Coupling varies linearly with depth from surface_coupling at 0 km to 1.0 at
    cmb_depth_km, so the LLSVP anti-correlation with subduction appears from the
    first timestep rather than waiting ~50 Myr for material to arrive physically.

    Parameters
    ----------
    heating_rate : float
        Heat units added per Myr everywhere (core heating).
    slab_suppression : float
        Heat units removed per unit-weight slablet per Myr (at full coupling = CMB depth).
    surface_coupling : float
        Fractional coupling at depth=0 (0 = no instantaneous signal, 1 = full).
    sigma_deg : float
        Horizontal smoothing (degrees) applied each Myr; for other time steps the
        width scales with sqrt(delta_time) so the smoothing acts like diffusion.
    max_heat : float
        Safety cap on accumulated heat units.
    """
    lons_g, lats_g = _mantle_grid(spacing_deg)
    nlat, nlon = cmb_heat.shape

    # Core heating everywhere
    cmb_heat = cmb_heat + heating_rate * delta_time

    if len(slab_tracers) > 0:
        all_active = slab_tracers[slab_tracers['active'].astype(bool)]
        if len(all_active) > 0:
            lons_f   = all_active['lon'].values.astype(float)
            lats_f   = all_active['lat'].values.astype(float)
            depths_f = all_active['depth_km'].values.astype(float)

            coupling = surface_coupling + (1.0 - surface_coupling) * (depths_f / cmb_depth_km)
            weights  = slab_suppression * delta_time * coupling * _slab_weights(all_active)

            lon_idx = np.clip(np.floor((lons_f + 180) / spacing_deg).astype(int), 0, nlon - 1)
            lat_idx = np.clip(np.floor((lats_f + 90)  / spacing_deg).astype(int), 0, nlat - 1)
            suppress = np.zeros_like(cmb_heat)
            np.add.at(suppress, (lat_idx, lon_idx), weights)
            cmb_heat = np.maximum(cmb_heat - suppress, 0.0)

    cmb_heat = spherical_gaussian_smooth(cmb_heat, lats_g, lons_g, sigma_deg * np.sqrt(delta_time))

    return np.minimum(cmb_heat, max_heat)


def _remove_area_weighted_mean(field, lats):
    w = np.cos(np.radians(lats))[:, None] * np.ones(field.shape[-1])
    means = (field * w).sum(axis=(-2, -1)) / w.sum()
    return field - means[:, None, None]


def generate_mantle_temp_field(slab_tracers,
                                spacing_deg=5.0,
                                depth_interval_km=100.0,
                                max_depth_km=2890.0,
                                temp_anomaly_K=-500.0,
                                thermal_tau_myr=200.0,
                                dVs_per_K=-0.001,
                                sigma_deg=3.0,
                                sigma_depth_km=150.0,
                                cmb_heat=None,
                                cmb_dVs_scale=-0.5,
                                cmb_decay_km=400.0,
                                tomo_smooth_horiz_deg=0.0,
                                tomo_smooth_depth_km=0.0):
    """Bin slab tracers onto a 3D grid and return a dVs/Vs anomaly Dataset.

    Each slablet contributes a cold temperature anomaly (temp_anomaly_K × its weight)
    that decays exponentially with age (thermal_tau_myr).  The summed ΔT per grid cell
    is smoothed (spherical Gaussian horizontally, Gaussian in depth) and converted to
    dVs/Vs % (dVs_per_K, default -0.001 %/K, so a fresh unit-weight slablet gives
    -0.5% per cell).

    If cmb_heat (a 2D nlat×nlon array from update_cmb_heat) is provided, a warm
    negative-dVs anomaly is added in the deepest layers, decaying upward as
    exp(-(distance from CMB / cmb_decay_km)²).

    Each depth layer has zero area-weighted mean, like tomography models.
    Returns xr.Dataset with variables 'dVs' and optionally 'cmb_heat'
    (dims: depth, lat, lon).
    """
    from scipy.ndimage import gaussian_filter1d

    lons_g, lats_g = _mantle_grid(spacing_deg)
    depths_g = np.arange(depth_interval_km / 2, max_depth_km, depth_interval_km)
    nlon, nlat, ndepth = len(lons_g), len(lats_g), len(depths_g)
    field = np.zeros((ndepth, nlat, nlon), dtype=float)

    def _smooth(f, sig_h_deg, sig_d_km):
        f = spherical_gaussian_smooth(f, lats_g, lons_g, sig_h_deg)
        if sig_d_km > 0:
            f = gaussian_filter1d(f, sig_d_km / depth_interval_km, axis=0, mode='nearest')
        return f

    if len(slab_tracers) > 0:
        active = slab_tracers[slab_tracers['active'].astype(bool)]
        if len(active) > 0:
            ages_f   = active['age'].values.astype(float)
            lons_f   = active['lon'].values.astype(float)
            lats_f   = active['lat'].values.astype(float)
            depths_f = active['depth_km'].values.astype(float)

            dT  = temp_anomaly_K * np.exp(-ages_f / thermal_tau_myr) * _slab_weights(active)
            dVs = dT * dVs_per_K

            lon_idx = np.clip(
                np.floor((lons_f + 180) / spacing_deg).astype(int), 0, nlon - 1)
            lat_idx = np.clip(
                np.floor((lats_f + 90)  / spacing_deg).astype(int), 0, nlat - 1)
            dep_idx = np.clip(
                np.floor(depths_f / depth_interval_km).astype(int), 0, ndepth - 1)
            np.add.at(field, (dep_idx, lat_idx, lon_idx), dVs)

            field = _smooth(field, sigma_deg, sigma_depth_km)

    # Inject LLSVP-like warm anomaly from CMB heat accumulation.
    # Iterate deepest→shallowest; break once decay becomes negligible.
    if cmb_heat is not None and cmb_heat.shape == (nlat, nlon):
        _n_injected = 0
        for dep_i in range(ndepth - 1, -1, -1):
            dist_from_cmb = max_depth_km - depths_g[dep_i]
            decay = np.exp(-(dist_from_cmb / cmb_decay_km) ** 2)
            if decay < 1e-4:
                break
            field[dep_i] += cmb_heat * cmb_dVs_scale * decay
            _n_injected += 1
        print(f'  [mantle] LLSVP injected into {_n_injected} depth layers; '
              f'cmb_heat range {cmb_heat.min():.3f}–{cmb_heat.max():.3f}; '
              f'dVs at CMB: {field[-1].min():.4f}–{field[-1].max():.4f} %')
    elif cmb_heat is not None:
        print(f'  [mantle] WARNING: cmb_heat shape {cmb_heat.shape} != expected ({nlat},{nlon}) — LLSVP skipped')

    field = _remove_area_weighted_mean(field, lats_g)

    if tomo_smooth_horiz_deg > 0 or tomo_smooth_depth_km > 0:
        field = _smooth(field, tomo_smooth_horiz_deg, tomo_smooth_depth_km)
        field = _remove_area_weighted_mean(field, lats_g)

    data_vars = {
        'dVs': xr.DataArray(
            field, dims=['depth', 'lat', 'lon'],
            coords={'depth': depths_g, 'lat': lats_g, 'lon': lons_g},
            attrs={'units': '%',
                   'long_name': 'dVs/Vs anomaly (positive = cold fast slab; negative = warm slow region)',
                   'temp_anomaly_K': temp_anomaly_K,
                   'thermal_tau_myr': thermal_tau_myr,
                   'dVs_per_K': dVs_per_K}),
    }
    if cmb_heat is not None and cmb_heat.shape == (nlat, nlon):
        data_vars['cmb_heat'] = xr.DataArray(
            cmb_heat, dims=['lat', 'lon'],
            coords={'lat': lats_g, 'lon': lons_g},
            attrs={'units': 'heat units',
                   'long_name': 'CMB heat accumulation (proxy for LLSVP intensity)'})
    return xr.Dataset(data_vars)


def generate_hotspot_locations(n=20, locations=None, seed=None):
    """Return a GeoDataFrame of fixed hot spot source points.

    If locations is a list of (lon, lat) tuples those are used directly.
    Otherwise n points are drawn uniformly on the sphere via the normal-vector
    method (generate (x,y,z) ~ N(0,1), normalise, convert to lon/lat).
    """
    if locations is not None:
        rows = [{'hotspot_id': i, 'geometry': Point(lon, lat)}
                for i, (lon, lat) in enumerate(locations)]
    else:
        rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
        xyz = rng.standard_normal((n, 3))
        xyz /= np.linalg.norm(xyz, axis=1, keepdims=True)
        lons, lats = xyz_to_lonlat(xyz)
        rows = [{'hotspot_id': i, 'geometry': Point(float(lons[i]), float(lats[i]))}
                for i in range(n)]
    return gpd.GeoDataFrame(rows, geometry='geometry').set_crs(epsg=4326)


def seed_hotspot_seamounts(hotspot_locs, prob=0.3, rng=None):
    """Probabilistically create new seamounts at fixed hot spot locations.

    For each hot spot, a seamount is created with probability prob.
    Returns a GeoDataFrame (possibly empty) with columns matching the tracers
    GeoDataFrame so that advect_tracers and deactivate_subducted_tracers work
    without modification.
    """
    if rng is None:
        rng = np.random.default_rng()
    rows = []
    for _, hs in hotspot_locs.iterrows():
        if rng.uniform() < prob:
            rows.append({
                'geometry':    hs.geometry,
                'hotspot_id':  int(hs['hotspot_id']),
                'age':         0.0,
                'active':      True,
            })
    if not rows:
        return gpd.GeoDataFrame(columns=['geometry', 'hotspot_id', 'age', 'active'],
                                geometry='geometry').set_crs(epsg=4326)
    return gpd.GeoDataFrame(rows, geometry='geometry').set_crs(epsg=4326)


def craton_counts_from_raster(craton_raster):
    """Return {plate_id: pixel_count} from a craton xarray DataArray.  Only counts pixels > 0."""
    if craton_raster is None:
        return {}
    vals = craton_raster.values
    vals = vals[np.isfinite(vals) & (vals > 0)].astype(int)
    ids, counts = np.unique(vals, return_counts=True)
    return {int(k): int(v) for k, v in zip(ids, counts)}


def resolve_downgoing_plate(p1, p2, polarity_dict, craton_counts=None):
    """Return the downgoing plate ID for a convergent boundary between p1 and p2.

    Priority:
      1. Craton rule: plate with cratons is overriding; the other is downgoing.
         If both have cratons, the one with more pixels is overriding.
      2. polarity_dict (sticky result from previous steps).
      3. Fallback: p1 downgoing (arbitrary, consistent).
    """
    if p1 == 0 or p2 == 0:
        return p1
    if craton_counts:
        a1 = craton_counts.get(p1, 0)
        a2 = craton_counts.get(p2, 0)
        if a1 > 0 and a2 == 0:
            return p2           # p1 continental → overriding; p2 downgoing
        if a2 > 0 and a1 == 0:
            return p1           # p2 continental → overriding; p1 downgoing
        if a1 > 0 and a2 > 0:
            return p2 if a1 >= a2 else p1   # larger craton is overriding
    pair = (min(p1, p2), max(p1, p2))
    return polarity_dict.get(pair, p1)


def _empty_earthquakes():
    return gpd.GeoDataFrame(
        columns=['geometry', 'depth_km', 'boundary_type',
                 'plateid1', 'plateid2', 'downgoing_plate', 'overriding_plate'],
        geometry='geometry').set_crs(epsg=4326)


def generate_synthetic_earthquakes(boundaries_gdf, polarity_dict=None, plate_polygons_gdf=None,
                                    craton_counts=None, plate_id_raster=None,
                                    n_candidates=100000,
                                    ridge_max_dist_deg=2.0,
                                    ridge_max_depth_km=50.0,
                                    slab_dip_deg=45.0,
                                    slab_max_depth_km=600.0,
                                    div_threshold_ridge=0.1,
                                    div_threshold_subduction=-0.1,
                                    ridge_eq_prob=1.0,
                                    subduction_eq_prob=1.0,
                                    rng=None):
    """Generate a synthetic earthquake catalog from plate boundary geometry.

    Candidate points are drawn uniformly on the sphere and assigned to the nearest
    active boundary (great-circle distance).

    Ridge earthquakes: random depth in [0, ridge_max_depth_km], symmetric about the
    ridge axis within ridge_max_dist_deg.

    Subduction earthquakes: only on the overriding plate (checked against
    plate_id_raster).  Depth increases linearly with horizontal distance from the
    boundary: depth = dist_km * tan(slab_dip_deg), capped at slab_max_depth_km.
    Maximum horizontal extent = slab_max_depth_km / tan(dip).

    Returns a GeoDataFrame with columns:
        geometry, depth_km, boundary_type, plateid1, plateid2,
        downgoing_plate, overriding_plate
    """
    if boundaries_gdf is None or 'divergence' not in boundaries_gdf.columns or len(boundaries_gdf) == 0:
        return _empty_earthquakes()
    if polarity_dict is None:
        polarity_dict = {}
    if rng is None:
        rng = np.random.default_rng()

    km_per_deg = np.radians(1.0) * EARTH_RADIUS_KM
    slab_dip_rad = np.radians(slab_dip_deg)
    slab_max_dist_deg = slab_max_depth_km / np.tan(slab_dip_rad) / km_per_deg
    max_query_deg = max(ridge_max_dist_deg, slab_max_dist_deg)

    div = boundaries_gdf['divergence'].values.astype(float)
    p1 = np.nan_to_num(boundaries_gdf['plateid1_poly1'].values.astype(float)).astype(int)
    p2 = np.nan_to_num(boundaries_gdf['plateid1_poly2'].values.astype(float)).astype(int)
    is_ridge = div > div_threshold_ridge
    is_sub   = (div < div_threshold_subduction) & (p1 != 0) & (p2 != 0)
    dg = np.array([resolve_downgoing_plate(a, b, polarity_dict, craton_counts) if s else 0
                   for a, b, s in zip(p1, p2, is_sub)], dtype=int)
    ov = np.where(dg == p1, p2, p1)

    b_lat, b_lon, b_row = _sample_arcs(boundaries_gdf, is_ridge | is_sub, 0.3)
    if len(b_lat) == 0:
        return _empty_earthquakes()

    # Uniform random points on the sphere
    lons, lats = xyz_to_lonlat(rng.standard_normal((n_candidates, 3)))

    tree = SphereTree(np.deg2rad(np.column_stack([b_lat, b_lon])))
    dists_rad, nn_idx = tree.query(np.deg2rad(np.column_stack([lats, lons])), k=1)
    dist_deg = np.degrees(dists_rad[:, 0])
    row = b_row[nn_idx[:, 0]]
    near = dist_deg < max_query_deg

    # Ridge: shallow, symmetric about axis
    ridge_ok = near & is_ridge[row] & (dist_deg <= ridge_max_dist_deg)
    ridge_ok &= rng.random(n_candidates) < ridge_eq_prob

    # Subduction: overriding-plate side only, depth from slab dip
    depth_sub = dist_deg * km_per_deg * np.tan(slab_dip_rad)
    sub_ok = near & is_sub[row] & (depth_sub <= slab_max_depth_km)
    if plate_id_raster is not None:
        sub_ok &= raster_lookup(plate_id_raster, lons, lats) == ov[row]
    sub_ok &= rng.random(n_candidates) < subduction_eq_prob

    keep = ridge_ok | sub_ok
    if not keep.any():
        return _empty_earthquakes()
    depth = np.where(ridge_ok, rng.uniform(0.0, ridge_max_depth_km, n_candidates), depth_sub)
    r = row[keep]
    return gpd.GeoDataFrame({
        'depth_km':         np.round(depth[keep], 1),
        'boundary_type':    np.where(ridge_ok[keep], 'ridge', 'subduction'),
        'plateid1':         p1[r],
        'plateid2':         p2[r],
        'downgoing_plate':  np.where(ridge_ok[keep], 0, dg[r]),
        'overriding_plate': np.where(ridge_ok[keep], 0, ov[r]),
    }, geometry=gpd.points_from_xy(lons[keep], lats[keep])).set_crs(epsg=4326)


def compute_topography_raster(plate_id_raster, craton_raster, boundaries_gdf,
                               polarity_dict=None, craton_counts=None,
                               continental_base_m=500,
                               noise_fraction=0.15,
                               spectral_beta=2.0,
                               andean_max_dist_deg=5.0,
                               andean_arc_dist_deg=1.5,
                               andean_peak_m=5000,
                               himalayan_max_dist_deg=4.0,
                               himalayan_peak_m=8000,
                               div_threshold_subduction=-0.0001,
                               coastal_margin_deg=2.0,
                               noise_seed=None):
    """Generate a synthetic elevation raster from tectonic setting.

    Continental pixels (craton_raster > 0) receive:
      - Background: continental_base_m + bandlimited spherical-harmonic noise
      - Andean zone: Gaussian elevation peak at andean_arc_dist_deg from the trench,
        only on the overriding (continental) plate side of subduction zones
      - Himalayan zone: linear ramp peaked at the boundary, for convergent
        segments with continental crust on both sides of that segment (checked
        locally, up to himalayan_max_dist_deg from the boundary)
    Distances are great-circle distances.  Ocean pixels remain at 0.
    noise_seed makes the noise pattern repeatable (e.g. identical between frames).
    Returns an xr.DataArray with lon/lat dims.
    """
    import pyshtools as pysh
    from scipy.interpolate import RegularGridInterpolator

    if polarity_dict is None:
        polarity_dict = {}
    if craton_counts is None:
        craton_counts = {}

    def _ensure_lonlat(da):
        return da.rename({k: v for k, v in {'x': 'lon', 'y': 'lat'}.items()
                          if k in da.dims})

    pid_r  = _ensure_lonlat(plate_id_raster)
    crat_r = _ensure_lonlat(craton_raster)

    lons = pid_r.lon.values
    lats = pid_r.lat.values
    raster_spacing = float(abs(lons[1] - lons[0]))

    lon2d, lat2d = np.meshgrid(lons, lats)          # (nlat, nlon)
    nlat, nlon   = lon2d.shape

    craton_vals = np.nan_to_num(crat_r.values, nan=0.0).astype(float)
    is_craton   = craton_vals > 0
    elev        = np.zeros((nlat, nlon), dtype=float)

    # --- Background continental elevation (flat; noise added after all zones) ---
    elev[is_craton] = float(continental_base_m)

    # --- Classify convergent segments: Andean vs Himalayan ---
    andean_pts   = []   # (lat, lon)
    andean_ovpid = []   # overriding plate ID for each Andean BallTree point
    himalayan_pts = []

    probe_dists = np.linspace(0.5, himalayan_max_dist_deg, 6)
    TESS_DEG = 0.5
    if boundaries_gdf is not None and 'divergence' in boundaries_gdf.columns and len(boundaries_gdf):
        conv = (boundaries_gdf['divergence'] < div_threshold_subduction).values
        for pos, arc in _segment_arcs(boundaries_gdf, conv):
            row = boundaries_gdf.iloc[pos]
            _r1, _r2 = row.get('plateid1_poly1', 0), row.get('plateid1_poly2', 0)
            p1 = int(_r1) if pd.notna(_r1) else 0
            p2 = int(_r2) if pd.notna(_r2) else 0
            if p1 == 0 or p2 == 0:
                continue

            dg = resolve_downgoing_plate(p1, p2, polarity_dict, craton_counts)
            ov = p2 if dg == p1 else p1
            is_andean = craton_counts.get(ov, 0) > 0
            if not is_andean:
                continue

            # Himalayan only if both plates have continental crust close to this segment
            mid_lat, mid_lon = arc.get_arc_point(0.5).to_lat_lon()
            normal = np.array(arc.get_great_circle_normal().to_xyz())
            probe_lon, probe_lat = [], []
            for sign in (+1, -1):
                for d in probe_dists:
                    lo, la = offset_along_normal(mid_lon, mid_lat, sign * normal, d)
                    probe_lon.append(float(lo)); probe_lat.append(float(la))
            probe_cr = raster_lookup(crat_r, probe_lon, probe_lat)
            is_himalayan = (probe_cr == dg).any() and (probe_cr == ov).any()

            n = max(int(np.ceil(np.degrees(arc.get_arc_length()) / TESS_DEG)), 1)
            for f in np.linspace(0.0, 1.0, n + 1):
                lat, lon = arc.get_arc_point(f).to_lat_lon()
                andean_pts.append([lat, lon])
                andean_ovpid.append(ov)
                if is_himalayan:
                    himalayan_pts.append([lat, lon])

    # --- Andean elevation ---
    if andean_pts:
        tree = SphereTree(np.deg2rad(andean_pts))
        ri, ci = np.where(is_craton)
        pts_rad = np.deg2rad(np.column_stack([lat2d[ri, ci], lon2d[ri, ci]]))
        dists_deg, nn_idx = tree.query(pts_rad, k=1)
        dists_deg = np.degrees(dists_deg[:, 0])
        ov_of_nearest = np.array(andean_ovpid)[nn_idx[:, 0]]
        # only pixels of the overriding plate (not the downgoing plate)
        ok = (dists_deg < andean_max_dist_deg) & (craton_vals[ri, ci].astype(int) == ov_of_nearest)
        arc_elev = andean_peak_m * np.exp(
            -((dists_deg - andean_arc_dist_deg) / andean_arc_dist_deg) ** 2)
        elev[ri[ok], ci[ok]] = np.maximum(elev[ri[ok], ci[ok]], arc_elev[ok])

    # --- Himalayan elevation ---
    if himalayan_pts:
        tree = SphereTree(np.deg2rad(himalayan_pts))
        ri, ci = np.where(is_craton)
        dists_deg, _ = tree.query(np.deg2rad(np.column_stack([lat2d[ri, ci], lon2d[ri, ci]])), k=1)
        dists_deg = np.degrees(dists_deg[:, 0])
        ok = dists_deg < himalayan_max_dist_deg
        him_elev = himalayan_peak_m * (1.0 - dists_deg / himalayan_max_dist_deg)
        elev[ri[ok], ci[ok]] = np.maximum(elev[ri[ok], ci[ok]], him_elev[ok])

    # --- Elevation-scaled spherical noise via spherical harmonic synthesis ---
    # lmax chosen so DH grid spacing ≈ raster_spacing: spacing = 90/(lmax+1)
    lmax = max(int(90.0 / raster_spacing) - 1, 10)
    l = np.arange(lmax + 1, dtype=float)
    power = np.where(l > 0, (l + 1.0) ** (-spectral_beta), 0.0)

    clm = pysh.SHCoeffs.from_random(power, seed=noise_seed)
    sh_grid = clm.expand(grid='DH')
    noise_arr = sh_grid.to_array()
    noise_arr /= noise_arr.std()

    dh_lats = sh_grid.lats()   # descending degrees
    dh_lons = sh_grid.lons()   # 0 → ~360

    interp = RegularGridInterpolator(
        (dh_lats[::-1], dh_lons), noise_arr[::-1, :],
        method='linear', bounds_error=False, fill_value=0.0)

    lon2d_360 = lon2d % 360.0
    pts = np.column_stack([lat2d.flatten(), lon2d_360.flatten()])
    noise = interp(pts).reshape(nlat, nlon)
    std = noise.std()
    if std > 0:
        noise /= std

    elev[is_craton] += noise[is_craton] * elev[is_craton] * noise_fraction
    elev[is_craton]  = np.maximum(elev[is_craton], 100.0)

    # --- Coastal taper: smooth elevation to zero within coastal_margin_deg of the
    #     craton edge (great-circle distance to the nearest non-craton pixel).
    coast = ~is_craton & np.logical_or.reduce(wrap_neighbours(is_craton))
    if coast.any() and is_craton.any():
        tree = SphereTree(np.deg2rad(np.column_stack([lat2d[coast], lon2d[coast]])))
        d, _ = tree.query(np.deg2rad(np.column_stack([lat2d[is_craton], lon2d[is_craton]])), k=1)
        t = np.clip(np.degrees(d[:, 0]) / coastal_margin_deg, 0.0, 1.0)
        elev[is_craton] *= t * t * (3.0 - 2.0 * t)   # smoothstep: 0 at coast, 1 inland

    return xr.DataArray(elev, dims=['lat', 'lon'],
                        coords={'lat': lats, 'lon': lons},
                        attrs={'units': 'meters', 'long_name': 'synthetic elevation'})


def physics_grid_age(tracer_lon, tracer_lat, tracer_age, plate_id_raster, craton_raster=None,
                     grid_deg=1.0, tracer_plate_id=None, raster_id_map=None):
    """Seafloor age (Myr) on the physics grid, from the seafloor tracers, for the force model.

    The grid is the raster sampled at about grid_deg (the same cells as plate_area_tensors).
    Each ocean cell takes the age of the nearest tracer on the same plate (great-circle
    distance), so ages do not leak across plate boundaries (at a trench the overriding
    plate would otherwise take the downgoing plate's ages).  A plate with no tracers
    takes the nearest tracer of any plate.  In a gap that has just opened at a ridge the
    nearest tracer is the youngest one beside it.  If tracer_plate_id is given, tracers
    that now sit on a different plate in plate_id_raster are left out (they crossed a
    boundary, i.e. are being subducted).  raster_id_map {raster pid: tracer pid} maps a
    plate that was split off this step back to its parent, so its tracers still count.
    Continental cells (craton raster > 0) are NaN.
    Returns a DataArray (lat, lon) on the sampled grid.
    """
    sl, _ = grid_stride(plate_id_raster, grid_deg)
    lons, lats = raster_coords(plate_id_raster)
    g_lon, g_lat = lons[sl[1]], lats[sl[0]]
    lon2d, lat2d = np.meshgrid(g_lon, g_lat)
    tracer_lon, tracer_lat = np.asarray(tracer_lon, float), np.asarray(tracer_lat, float)
    tracer_age = np.asarray(tracer_age, float)
    age = np.full(lon2d.shape, np.nan)
    if len(tracer_age):
        on_pid = raster_lookup(plate_id_raster, tracer_lon, tracer_lat).astype(float)
        cell_pid = np.asarray(plate_id_raster.values, dtype=float)[sl].ravel()
        if raster_id_map:
            on_pid = map_plate_ids(on_pid, raster_id_map)
            cell_pid = map_plate_ids(cell_pid, raster_id_map)
        if tracer_plate_id is not None:
            keep = on_pid == np.asarray(tracer_plate_id)
            tracer_lon, tracer_lat, tracer_age, on_pid = (tracer_lon[keep], tracer_lat[keep],
                                                          tracer_age[keep], on_pid[keep])
        cells = np.deg2rad(np.column_stack([lat2d.ravel(), lon2d.ravel()]))
        flat = np.full(len(cells), np.nan)
        if len(tracer_age):
            done = np.zeros(len(cells), bool)
            for pid in np.unique(cell_pid[np.isfinite(cell_pid)]):
                mine = on_pid == pid
                if not mine.any():
                    continue
                c = cell_pid == pid
                _, idx = SphereTree(np.deg2rad(np.column_stack([tracer_lat[mine], tracer_lon[mine]]))).query(cells[c], k=1)
                flat[c] = tracer_age[mine][idx[:, 0]]
                done |= c
            if not done.all():
                _, idx = SphereTree(np.deg2rad(np.column_stack([tracer_lat, tracer_lon]))).query(cells[~done], k=1)
                flat[~done] = tracer_age[idx[:, 0]]
        age = flat.reshape(lon2d.shape)
    if craton_raster is not None:
        age[np.nan_to_num(np.asarray(craton_raster.values)[sl]) > 0] = np.nan
    return xr.DataArray(age, dims=('lat', 'lon'), coords={'lat': g_lat, 'lon': g_lon},
                        attrs={'units': 'Myr', 'long_name': 'seafloor age (physics grid)'})


def map_plate_ids(ids, mapping):
    """Replace plate IDs by mapping[id] where the mapping has one (vectorised)."""
    ids = np.asarray(ids, dtype=float)
    if not mapping or not len(ids):
        return ids
    uniq, inverse = np.unique(ids, return_inverse=True)
    mapped = np.array([mapping.get(int(u), u) if np.isfinite(u) else u for u in uniq], dtype=float)
    return mapped[inverse].reshape(ids.shape)


def force_fields_needed(params):
    """Which per-step fields the selected force modes read.  The simple modes need none,
    so they cost exactly what a step cost before the force-model work."""
    pull = params.get('slab_pull_mode', 'hard')
    return {'age': (pull == 'slab_length' or params.get('polarity_mode', 'id') == 'age'
                    or params.get('ridge_push_mode', 'boundary') == 'body'),
            'slab': pull == 'slab_length' or params.get('slab_suction_mode', 'ratio') == 'depth_split'}


def predict_tracers(lon, lat, age, plate_id_raster, stage_pole_dict, delta_time):
    """Tracer positions and ages one step ahead: each tracer rotated with the stage pole of
    the plate it sits on in plate_id_raster (the raster before the move), age + delta_time.
    Returns (lon, lat, age, plate_id).  No tracers are removed or added here (that needs
    the new boundaries); pass plate_id to physics_grid_age to drop those that crossed."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    pids = raster_lookup(plate_id_raster, lon, lat)
    new_lon, new_lat = lon.copy(), lat.copy()
    for pid in np.unique(pids[np.isfinite(pids)]):
        rotation = stage_pole_dict.get(int(pid))
        if rotation is None:
            continue
        m = pids == pid
        new_lon[m], new_lat[m] = rotate_lonlat(rotation, lon[m], lat[m])
    return new_lon, new_lat, np.asarray(age, float) + delta_time, pids


def generate_seafloor_age_grid(tracer_csv, craton_raster=None, region='d', spacing='1d'):
    """Interpolate seafloor age tracers to a global grid and optionally mask continents.

    The tracer CSV is expected to have columns [age, lon, lat] (the format
    written by run_SDP.py).  blockmedian reduces duplicate points per cell
    before sphinterpolate produces the final grid.

    craton_raster : xr.DataArray with plate-ID values > 0 on continental pixels.
        Where it is positive the returned age grid is set to NaN.
    """
    # pygmt spacing must be a string like '0.5d'; convert float/int if needed
    if not isinstance(spacing, str):
        spacing = f'{spacing}d'

    res = pygmt.blockmedian(tracer_csv, region=region, incols='1,2,0',
                            spacing=spacing, output_type='numpy')
    age_grid = pygmt.sphinterpolate(res, region=region, spacing=spacing)

    if craton_raster is not None:
        # Normalise dim names so interp works regardless of x/y vs lon/lat
        rename = {k: v for k, v in {'x': 'lon', 'y': 'lat'}.items()
                  if k in craton_raster.dims}
        cr = craton_raster.rename(rename) if rename else craton_raster
        cr_on_grid = cr.interp(lat=age_grid.lat, lon=age_grid.lon, method='nearest')
        age_grid = age_grid.where(cr_on_grid <= 0)

    age_grid.attrs.update({'units': 'Myr', 'long_name': 'seafloor age'})
    return age_grid


def generate_topo_bathy_grid(seafloor_age_grid, topography_raster=None,
                              default_ocean_depth_m=-3500):
    """Combine GDH1 bathymetry with synthetic topography into a seamless grid.

    GDH1 model (Stein & Stein 1992):
      depth = -(2600 + 365 * sqrt(age))          for age < 20 Myr
      depth = -(5651 - 2473 * exp(-0.0278*age))  for age >= 20 Myr

    topography_raster : xr.DataArray, optional.  Continental elevations (> 0 m).
        Where present these overwrite the bathymetry.  If None, only bathymetry
        is returned.
    default_ocean_depth_m : fill value for ocean pixels with no age data.
        Set to None to leave those pixels as NaN.
    """
    age = seafloor_age_grid.values.copy().astype(float)

    # GDH1 — work in-place on a pre-NaN'd array
    bathy = np.full_like(age, np.nan)
    young = np.isfinite(age) & (age < 20.)
    old   = np.isfinite(age) & (age >= 20.)
    bathy[young] = -(2600. + 365. * np.sqrt(age[young]))
    bathy[old]   = -(5651. - 2473. * np.exp(-0.0278 * age[old]))

    combined = bathy.copy()

    # Fill uncovered ocean pixels with a default depth
    if default_ocean_depth_m is not None:
        combined[~np.isfinite(combined)] = float(default_ocean_depth_m)

    # Overwrite continental pixels with topography
    if topography_raster is not None:
        rename = {k: v for k, v in {'x': 'lon', 'y': 'lat'}.items()
                  if k in topography_raster.dims}
        topo = topography_raster.rename(rename) if rename else topography_raster
        topo_on_grid = topo.interp(lat=seafloor_age_grid.lat,
                                   lon=seafloor_age_grid.lon, method='nearest')
        topo_vals = topo_on_grid.values
        continental = topo_vals > 0
        combined[continental] = topo_vals[continental]

    return xr.DataArray(combined,
                        dims=seafloor_age_grid.dims,
                        coords=seafloor_age_grid.coords,
                        attrs={'units': 'meters', 'long_name': 'topography-bathymetry'})


def write_rotation_file(rotation_records, path, fixed_plate_id=0):
    """Write accumulated finite rotations to a GPlates .rot file.

    rotation_records is a list of (time_myr, plate_id, FiniteRotation) tuples
    produced by accumulating stage poles across timesteps.  A t=0 identity
    entry is written automatically for every plate that appears in the records.
    """
    from collections import defaultdict

    by_plate = defaultdict(list)
    for time_myr, pid, rot in rotation_records:
        by_plate[pid].append((time_myr, rot))

    lines = []
    for pid in sorted(by_plate.keys()):
        # t=0 identity line (required by GPlates)
        lines.append('{:d}   {:8.4f}   {:9.4f}   {:10.4f}   {:9.4f}   {:d}'.format(
            pid, 0.0, 0.0, 0.0, 0.0, fixed_plate_id))
        for time_myr, rot in sorted(by_plate[pid], key=lambda x: x[0]):
            if rot.represents_identity_rotation():
                lat, lon, angle = 0.0, 0.0, 0.0
            else:
                lat, lon, angle = rot.get_lat_lon_euler_pole_and_angle_degrees()
            lines.append('{:d}   {:8.4f}   {:9.4f}   {:10.4f}   {:9.4f}   {:d}'.format(
                pid, time_myr, lat, lon, angle, fixed_plate_id))

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def remove_boundary_segments(arc, lon_threshold=None, lat_threshold=None,
                             raster_spacing=0.25):
    """Remove arc segments that run along the dateline or poles.

    Only called for arcs shared between two different plates (outer-boundary
    arcs are excluded upstream).  Uses a narrow threshold (raster_spacing/4)
    so only vertices at exactly ±180°/±90° (from clip) are filtered, while
    real plate-boundary vertices at ±179.875° (0.25° grid cell edges) are kept.
    """
    if lon_threshold is None:
        lon_threshold = 180.0 - raster_spacing / 4   # e.g. 179.9375 for 0.25° grid — catches
    if lat_threshold is None:                        # exactly ±180° but not 179.875° cell edges
        lat_threshold = 90.0 - raster_spacing / 4
    arc = np.asarray(arc)

    # Check which points are on boundaries
    on_dateline = np.abs(np.abs(arc[:, 0]) - 180) < (180 - lon_threshold)
    on_pole = np.abs(arc[:, 1]) > lat_threshold
    on_boundary = on_dateline | on_pole
    
    # Keep segments where at least one endpoint is NOT on boundary
    keep_segments = []
    for i in range(len(arc) - 1):
        if not (on_boundary[i] and on_boundary[i+1]):
            keep_segments.append(i)
    
    if not keep_segments:
        return []
    
    # Group consecutive segments into sub-arcs
    sub_arcs = []
    current = [arc[keep_segments[0]]]
    
    for i in range(len(keep_segments)):
        current.append(arc[keep_segments[i] + 1])
        if i < len(keep_segments) - 1 and keep_segments[i+1] != keep_segments[i] + 1:
            sub_arcs.append(np.array(current))
            current = [arc[keep_segments[i+1]]]
    
    sub_arcs.append(np.array(current))
    return sub_arcs


def topo_polygons_to_boundaries(plate_polygons, raster_spacing=0.25, max_edge_deg=1.0):

    topo = tp.Topology(plate_polygons, prequantize=False, toposimplify=False)

    topo_dict = topo.to_dict()
    arcs = topo_dict['arcs']

    # Build arc-to-polygon mapping with correct decoding
    arc_to_polys = {}
    feat_to_plateid = {}  # Store PLATEID1 for each feature

    for obj_data in topo_dict['objects'].values():
        for feat_idx, feat in enumerate(obj_data['geometries']):
            # Store the PLATEID1 for this feature
            feat_to_plateid[feat_idx] = feat['properties'].get('PLATEID1')

            for ring in feat['arcs']:
                for arc_idx in ring:
                    # Decode TopoJSON arc index
                    if arc_idx < 0:
                        abs_idx = -arc_idx - 1
                    else:
                        abs_idx = arc_idx
                    arc_to_polys.setdefault(abs_idx, []).append(feat_idx)

    # Apply to all arcs with polygon IDs and PLATEID1
    boundaries = []
    for arc_idx, arc in enumerate(arcs):
        poly_list = arc_to_polys.get(arc_idx, [])
        unique_polys = list(set(poly_list))

        # Only process arcs shared between exactly two polygons with different PLATEIDs.
        # Outer-boundary arcs (world edge at ±180°/±90°) have len==1 and are never real
        # plate boundaries; same-plate arcs (duplicate ID fragments) are internal artifacts.
        if len(unique_polys) < 2:
            continue
        pid1 = feat_to_plateid.get(unique_polys[0])
        pid2 = feat_to_plateid.get(unique_polys[1])
        if pid1 == pid2:
            continue

        sub_arcs = remove_boundary_segments(arc, raster_spacing=raster_spacing)

        for sub_arc in sub_arcs:
            boundaries.append({
                # Densify: every 2-vertex piece is later treated as a great-circle arc.
                'geometry': LineString(_densify_coords(sub_arc, max_edge_deg)),
                'poly1': unique_polys[0],
                'poly2': unique_polys[1],
                'plateid1_poly1': pid1,
                'plateid1_poly2': pid2
            })

    if not boundaries:
        return gpd.GeoDataFrame(columns=['geometry', 'poly1', 'poly2',
                                         'plateid1_poly1', 'plateid1_poly2'],
                                geometry='geometry').set_crs(epsg=4326)
    return gpd.GeoDataFrame(boundaries, geometry='geometry').set_crs(epsg=4326)


def fill_nan_spherical(da, tie_break='hash', previous=None, seed=0):
    """Fill NaN cells with the value of the nearest valid cell (great-circle distance).

    Only valid cells that border a NaN cell can be the nearest valid cell, so the
    search tree is built from those alone (a few thousand instead of ~10^6 cells).

    When several valid cells are equally near (common on a regular grid):
      tie_break='hash'            seeded coin flip per cell (SphereTree.query)
      tie_break='previous_owner'  the value the cell had in `previous` (an array or
                                  DataArray on the same grid) wins if it is one of the
                                  tied values; otherwise the coin flip decides
    """
    vals = da.values
    invalid = np.isnan(vals)
    if not invalid.any() or invalid.all():
        return da.copy()
    lons, lats = raster_coords(da)
    lon, lat = np.meshgrid(lons, lats)
    candidates = ~invalid & np.logical_or.reduce(wrap_neighbours(invalid))

    # BallTree with haversine needs (lat, lon) in radians
    tree = SphereTree(np.deg2rad(np.column_stack([lat[candidates], lon[candidates]])))
    q = np.deg2rad(np.column_stack([lat[invalid], lon[invalid]]))
    _, idx = tree.query(q, k=1, seed=seed)
    chosen = vals[candidates][idx[:, 0]]

    if tie_break == 'previous_owner' and previous is not None:
        prev = np.asarray(getattr(previous, 'values', previous))[invalid]
        _, cand, tied = tree.query_candidates(q)
        cand_vals = vals[candidates][cand]
        keep_prev = np.isfinite(prev) & (tied & (cand_vals == prev[:, None])).any(axis=1)
        chosen = np.where(keep_prev, prev, chosen)

    filled = vals.copy()
    filled[invalid] = chosen

    return da.copy(data=filled)


def haversine(lon1, lat1, lon2, lat2):
    """Great-circle distance in degrees."""
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))))


def _two_nearest_ids(lat, lon, lons, lats, vals, radius_deg=2.0, ignore_zero=False):
    """The two nearest distinct raster values within radius_deg (great-circle) of a point.

    The search window wraps across ±180 and widens in longitude by 1/cos(lat), so
    neighbours are found at any longitude and latitude (the whole latitude band is
    searched when the window reaches a pole).
    """
    buf = radius_deg * 1.5
    rows = np.where((lats >= lat - buf) & (lats <= lat + buf))[0]
    if len(rows) == 0:
        return None, None
    seam = abs((lons[-1] - lons[0]) - 360.0) < 1e-6
    n_unique = len(lons) - 1 if seam else len(lons)
    edge_lat = max(abs(lat - buf), abs(lat + buf))
    if edge_lat >= 89.999:
        cols = np.arange(n_unique)
    else:
        lon_buf = buf / np.cos(np.radians(edge_lat))
        if lon_buf >= 180:
            cols = np.arange(n_unique)
        else:
            _, c0 = raster_nearest_index(lons, lats, lon, lat)
            half = int(np.ceil(lon_buf / abs(lons[1] - lons[0])))
            cols = np.unique((int(c0) + np.arange(-half, half + 1)) % n_unique)

    sub = vals[np.ix_(rows, cols)]
    lon2d, lat2d = np.meshgrid(lons[cols], lats[rows])
    dists = haversine(lon, lat, lon2d, lat2d)
    ok = (dists <= radius_deg) & np.isfinite(sub)
    if ignore_zero:
        ok &= sub > 0
    if not ok.any():
        return None, None
    ids = sub[ok].astype(int)
    d = dists[ok]
    order = np.lexsort((d, ids))                 # by id, then distance
    ids, d = ids[order], d[order]
    first = np.r_[True, ids[1:] != ids[:-1]]     # nearest cell of each id
    ids, d = ids[first], d[first]
    ranked = ids[np.argsort(d, kind='stable')]
    if len(ranked) >= 2:
        return int(ranked[0]), int(ranked[1])
    return int(ranked[0]), None


def get_two_nearest_ids(coordinate, da, radius_deg=2.0, ignore_zero=False):
    """Two nearest distinct integers of a (lat, lon) raster within radius_deg of
    coordinate = (lat, lon).  See _two_nearest_ids."""
    lons, lats = raster_coords(da)
    return _two_nearest_ids(coordinate[0], coordinate[1], lons, lats, da.values,
                            radius_deg, ignore_zero)


# ---------------------------------------------------------------------------
# Batched neighbour search (numba-compiled when numba is installed)
# ---------------------------------------------------------------------------

try:
    import numba
except ImportError:          # fall back to the per-point numpy version
    numba = None


def _two_nearest_kernel(q_lat, q_lon, q_xyz, lons, lats, gx, gy, gz, vals,
                        radius_deg, ignore_zero, n_unique):
    """Same search as _two_nearest_ids, for many points at once.  Distances are
    compared as cos(great-circle distance) = dot product of unit vectors; equal
    distances are ranked by the smaller ID, as in _two_nearest_ids."""
    n = len(q_lat)
    out = np.full((n, 2), -1, np.int64)
    dlon = lons[1] - lons[0]
    dlat = lats[1] - lats[0]
    buf = radius_deg * 1.5
    cos_radius = np.cos(np.radians(radius_deg))
    eps = 1e-13
    for k in range(n):
        la = q_lat[k]
        qx, qy, qz = q_xyz[k, 0], q_xyz[k, 1], q_xyz[k, 2]
        r0 = max(int(np.ceil((la - buf - lats[0]) / dlat)), 0)
        r1 = min(int(np.floor((la + buf - lats[0]) / dlat)), len(lats) - 1)
        edge_lat = max(abs(la - buf), abs(la + buf))
        half = n_unique // 2
        if edge_lat < 89.999:
            lon_buf = buf / np.cos(np.radians(edge_lat))
            if lon_buf < 180:
                half = int(np.ceil(lon_buf / dlon))
        c0 = int(np.rint(((q_lon[k] - lons[0]) % 360.0) / dlon))
        id1, s1, id2, s2 = -1, -2.0, -1, -2.0          # s = cos(distance): larger is nearer
        for r in range(r0, r1 + 1):
            for dc in range(-half, half + 1):
                c = (c0 + dc) % n_unique
                v = vals[r, c]
                if v != v or (ignore_zero and v <= 0):
                    continue
                s = gx[r, c] * qx + gy[r, c] * qy + gz[r, c] * qz
                if s < cos_radius:
                    continue
                iv = int(v)
                if iv == id1:
                    s1 = max(s1, s)
                elif iv == id2:
                    s2 = max(s2, s)
                    if s2 > s1 + eps or (abs(s2 - s1) <= eps and id2 < id1):
                        id1, id2, s1, s2 = id2, id1, s2, s1
                elif s > s1 + eps or (abs(s - s1) <= eps and iv < id1):
                    id2, s2, id1, s1 = id1, s1, iv, s
                elif id2 < 0 or s > s2 + eps or (abs(s - s2) <= eps and iv < id2):
                    id2, s2 = iv, s
        out[k, 0] = id1
        out[k, 1] = id2
    return out


if numba is not None:
    _two_nearest_kernel = numba.njit(cache=True)(_two_nearest_kernel)

_GRID_XYZ_CACHE = {}


def two_nearest_ids_batch(lat, lon, da, radius_deg=2.0, ignore_zero=False):
    """Two nearest distinct raster values within radius_deg of each point.

    Vectorised equivalent of calling _two_nearest_ids per point.  Returns an
    (N, 2) int array; -1 where fewer than one/two IDs are found.
    """
    lat = np.atleast_1d(np.asarray(lat, dtype=float))
    lon = np.atleast_1d(np.asarray(lon, dtype=float))
    lons, lats = raster_coords(da)
    vals = np.ascontiguousarray(da.values, dtype=float)
    if numba is None:
        out = np.full((len(lat), 2), -1, np.int64)
        for k in range(len(lat)):
            a, b = _two_nearest_ids(lat[k], lon[k], lons, lats, vals, radius_deg, ignore_zero)
            out[k] = (-1 if a is None else a, -1 if b is None else b)
        return out
    key = (vals.shape, float(lons[0]), float(lats[0]), float(lons[1] - lons[0]))
    if key not in _GRID_XYZ_CACHE:
        lon2d, lat2d = np.meshgrid(lons, lats)
        xyz = lonlat_to_xyz(lon2d, lat2d)
        _GRID_XYZ_CACHE[key] = tuple(np.ascontiguousarray(xyz[..., j]) for j in range(3))
    gx, gy, gz = _GRID_XYZ_CACHE[key]
    seam = abs((lons[-1] - lons[0]) - 360.0) < 1e-6
    n_unique = len(lons) - 1 if seam else len(lons)
    return _two_nearest_kernel(lat, lon, lonlat_to_xyz(lon, lat), lons.astype(float), lats.astype(float),
                               gx, gy, gz, vals, float(radius_deg), bool(ignore_zero), n_unique)
