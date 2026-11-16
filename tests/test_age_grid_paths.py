"""The seafloor age grid does not depend on the run folder's name. GMT reads '=', '+' or '?'
in a file name as format flags; an output folder named like 'drag=none' once lost every age
grid."""
import os
import shutil
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import self_driving_plates as sdp  # noqa: E402


def test_age_grid_is_the_same_whatever_the_folder_is_called(tmp_path):
    rng = np.random.default_rng(0)
    n = 4000
    lat = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))       # uniform on the sphere
    lon = rng.uniform(-180, 180, n)
    age = 100 * np.abs(np.sin(np.radians(lat))) + 10 * np.cos(np.radians(lon))
    plain = tmp_path / 'plain'
    odd = tmp_path / 'drag=none+x'
    plain.mkdir()
    odd.mkdir()
    pd.DataFrame({'age': age, 'lon': lon, 'lat': lat}).to_csv(plain / 'tracers_0000.csv', index=False)
    shutil.copyfile(plain / 'tracers_0000.csv', odd / 'tracers_0000.csv')
    a = sdp.generate_seafloor_age_grid(str(plain / 'tracers_0000.csv'), spacing=2)
    b = sdp.generate_seafloor_age_grid(str(odd / 'tracers_0000.csv'), spacing=2)
    assert np.isfinite(a.values).mean() > 0.99
    assert np.array_equal(a.values, b.values, equal_nan=True)
    assert os.listdir(odd) == ['tracers_0000.csv']
