"""PlateDriver's side of the plate-history contract: self_driving_plates.plate_step.

    conda run -n pygmt17 python -m pytest tests

A two-plate world: plate 1 west of 0 deg longitude, plate 2 east of it, a trench along the meridian
from 30 S to 30 N where plate 1 goes down under plate 2, and plate 1 turning about the north pole.
"""
import os
import sys

import numpy as np
import pytest
import xarray as xr

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import geopandas as gpd                         # noqa: E402
import pygplates                                # noqa: E402
from shapely.geometry import LineString         # noqa: E402

import self_driving_plates as sdp               # noqa: E402

PH = pytest.importorskip('slablet.plate_history')
PARAMS = {'div_threshold_subduction': -1e-4, 'div_threshold_ridge': 0.2}


def _world():
    lons, lats = np.arange(-180, 180.25, 0.25), np.arange(-90, 90.25, 0.25)
    pid = np.where(lons[None, :] < 0, 1.0, 2.0) * np.ones((len(lats), 1))
    raster = xr.DataArray(pid, dims=('lat', 'lon'), coords={'lat': lats, 'lon': lons})
    poles = {1: pygplates.FiniteRotation(pygplates.PointOnSphere(0, 0, 1), np.radians(1.0)),
             2: pygplates.FiniteRotation()}
    bdf = gpd.GeoDataFrame({'divergence': [-50.0], 'plateid1_poly1': [1.0], 'plateid1_poly2': [2.0]},
                           geometry=[LineString([(0, -30), (0, 30)])], crs=4326)
    clons, clats = np.arange(-179.5, 180, 1.0), np.arange(-89.5, 90, 1.0)
    age = xr.DataArray(np.where(clons[None, :] < 0, 40.0, 10.0) * np.ones((len(clats), 1)),
                       dims=('lat', 'lon'), coords={'lat': clats, 'lon': clons})
    return raster, poles, bdf, age


def test_plate_step_meets_the_contract_and_keeps_the_trench_samples():
    raster, poles, bdf, age = _world()
    step = sdp.plate_step(5.0, 1.0, raster, poles, bdf, PARAMS, polarity_dict={(1, 2): 1}, age_grid=age)
    PH.check(step)
    plain = sdp.trench_samples(bdf, PARAMS['div_threshold_subduction'], 1.0)
    for k in ('lon', 'lat', 'convergence_km_per_myr'):
        np.testing.assert_array_equal(step.trench[k], plain[k])


def test_polarity_age_and_trench_velocity():
    raster, poles, bdf, age = _world()
    step = sdp.plate_step(5.0, 1.0, raster, poles, bdf, PARAMS, polarity_dict={(1, 2): 1}, age_grid=age)
    np.testing.assert_allclose(step.trench['downdip_azimuth_deg'], 90.0, atol=1e-6)   # towards plate 2, east
    np.testing.assert_array_equal(step.trench['seafloor_age_myr'], 40.0)                # plate 1's seafloor
    np.testing.assert_allclose(step.trench['trench_east_km_per_myr'], 0.0, atol=1e-9)   # plate 2 is still
    flipped = sdp.plate_step(5.0, 1.0, raster, poles, bdf, PARAMS, polarity_dict={(1, 2): 2}, age_grid=age)
    np.testing.assert_allclose(flipped.trench['downdip_azimuth_deg'], 270.0, atol=1e-6)


def test_velocities_match_plate_velocity_ne():
    raster, poles, bdf, age = _world()
    step = sdp.plate_step(5.0, 1.0, raster, poles, bdf, PARAMS, age_grid=age)
    lon, lat = np.array([-60.0, -20.0, 40.0]), np.array([10.0, -40.0, 0.0])
    east, north = PH.velocity_at(step, lon, lat)
    e1, n1 = sdp.plate_velocity_ne(poles[1], lon[:2], lat[:2], 1.0)
    # plate_velocity_ne uses pygplates' radius (6371.009 km), the contract slablet's (6371 km): 1.4e-6 apart
    np.testing.assert_allclose(east[:2], e1, rtol=2e-6)
    np.testing.assert_allclose(east[2], 0.0, atol=1e-12)
