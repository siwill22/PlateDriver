"""PlateDriver's copies of slablet's sphere helpers (sphere_geometry.py) must stay identical to them.

    conda run -n pygmt17 python -m pytest tests

PlateDriver keeps its own copies so it runs without slablet (a private package); when slablet is
installed, these tests check that the copies have not drifted, to the last bit.
"""
import os
import sys

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import sphere_geometry as sg                    # noqa: E402

sphere = pytest.importorskip('slablet.sphere')


def test_lonlat_xyz_round_trip_is_identical_to_slablet():
    rng = np.random.default_rng(0)
    lon, lat = rng.uniform(-180, 180, 1000), rng.uniform(-90, 90, 1000)
    np.testing.assert_array_equal(sg.lonlat_to_xyz(lon, lat), sphere.lonlat_to_xyz(lon, lat))
    xyz = rng.normal(size=(1000, 3))
    for a, b in zip(sg.xyz_to_lonlat(xyz), sphere.xyz_to_lonlat(xyz)):
        np.testing.assert_array_equal(a, b)


def test_spherical_kernel_is_identical_to_slablet():
    lats = tuple(np.round(np.arange(-88.0, 90.0, 4.0), 6))
    lons = tuple(np.round(np.arange(-178.0, 180.0, 4.0), 6))
    a = sg.spherical_kernel(lats, lons, 6.0)
    b = sphere._spherical_kernel(lats, lons, 6.0)
    assert (a != b).nnz == 0
