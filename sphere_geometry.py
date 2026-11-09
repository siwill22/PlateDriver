"""Spherical helpers used throughout PlateDriver.

These are copies of slablet.sphere's lonlat_to_xyz, xyz_to_lonlat and spherical Gaussian kernel,
so that PlateDriver runs without slablet (a private package) whenever no mantle option is on.
tests/test_sphere_geometry.py checks that they still agree with slablet's when slablet is
installed; change both together.
"""
import numpy as np
from sklearn.neighbors import BallTree


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


_KERNEL_CACHE = {}


def spherical_kernel(lats, lons, sigma_deg):
    """Row-normalised sparse Gaussian kernel on a lon/lat grid (great-circle distances, cell-area
    weights), cached per grid and sigma.  lats and lons are tuples, so they can be the cache key."""
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
