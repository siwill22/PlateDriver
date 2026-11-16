"""The basal drag tensor of a spherical cap against its analytic value.

For a plate covering a cap of angular radius a about the north pole, on the unit sphere,
T = ∫ (𝟙 − r rᵀ) dA is diagonal with
    T_zz = ∫ sin²θ dA = 2π (2/3 − cos a + cos³a / 3)
    T_xx = T_yy = A − T_zz / 2,   A = 2π (1 − cos a).
The rest of the sphere is a second plate; the two tensors add up to (8π/3) 𝟙.
"""
import os
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from self_driving_plates import plate_area_tensors


def cap_raster(cap_deg, spacing=0.25):
    """A global plate-ID raster laid out like the model's (gridline-registered lon/lat):
    plate 1 within cap_deg of the north pole, plate 2 elsewhere."""
    lats = np.arange(-90, 90 + spacing / 2, spacing)
    lons = np.arange(-180, 180 + spacing / 2, spacing)
    pid = np.where(lats[:, None] >= 90 - cap_deg, 1.0, 2.0) * np.ones((1, len(lons)))
    return xr.DataArray(pid, coords={'lat': lats, 'lon': lons}, dims=('lat', 'lon'))


def test_cap_tensor_matches_the_analytic_value():
    a = np.radians(30.0)
    pids, tensors, areas = plate_area_tensors(cap_raster(30.0), grid_deg=1.0)
    cap = tensors[list(pids).index(1)]
    area = 2 * np.pi * (1 - np.cos(a))
    t_zz = 2 * np.pi * (2 / 3 - np.cos(a) + np.cos(a) ** 3 / 3)
    expected = np.diag([area - t_zz / 2, area - t_zz / 2, t_zz])
    assert np.allclose(areas[list(pids).index(1)], area, rtol=2e-3)
    assert np.allclose(cap, expected, rtol=2e-3, atol=1e-4)


def test_the_two_plates_add_up_to_the_whole_sphere():
    pids, tensors, _ = plate_area_tensors(cap_raster(30.0), grid_deg=1.0)
    assert np.allclose(tensors.sum(axis=0), 8 * np.pi / 3 * np.eye(3), rtol=2e-3, atol=1e-3)
