"""PlateDriver's side of the mantle model (slablet, ~/GIT/slablet).

    conda run -n pygmt17 python -m pytest tests

slablet's own tests (in its repo) check the model on the trench samples saved in
tests/data/golden_samples.npz.  This checks that PlateDriver's adapter,
self_driving_plates.trench_samples, still makes exactly those samples from the same
boundaries, so the two halves of the golden test still meet.
"""
import os
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import geopandas as gpd                         # noqa: E402
from shapely.geometry import LineString         # noqa: E402

import self_driving_plates as sdp               # noqa: E402

SAMPLES = os.path.join(REPO, 'tests', 'data', 'golden_samples.npz')


def boundaries(step):
    # three trenches (one crossing the dateline, one at high latitude, one slow) and a ridge
    lines = [LineString([(170, -30), (-170, -10)]), LineString([(-60, 60), (-20, 72)]),
             LineString([(30, -5), (50, 5)]), LineString([(0, -50), (10, -20)])]
    div = [-60.0 - step, -25.0, -4.0, 30.0]
    return gpd.GeoDataFrame({'divergence': div}, geometry=lines, crs=4326)


def test_adapter_makes_the_golden_samples():
    with np.load(SAMPLES) as g:
        for step in range(11):
            s = sdp.trench_samples(boundaries(step))
            for k in ('lon', 'lat', 'convergence_km_per_myr'):
                np.testing.assert_array_equal(np.asarray(s[k], dtype=float), g['s{}_{}'.format(step, k)],
                                              err_msg='step {} {}'.format(step, k))


def test_adapter_returns_no_samples_without_divergence():
    s = sdp.trench_samples(None)
    assert all(len(v) == 0 for v in s.values())
