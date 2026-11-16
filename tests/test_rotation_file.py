"""Rotation files: elapsed model time, or time counted back from the saved step (GPlates' sense)."""
import os
import sys

import numpy as np
import pygplates

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from self_driving_plates import write_rotation_file


def _records(rate_deg=1.0, steps=10):
    """Plate 101 turning about the north pole at rate_deg per Myr: (t, pid, finite rotation)."""
    return [(float(t), 101, pygplates.FiniteRotation((90.0, 0.0), np.radians(rate_deg * t)))
            for t in range(1, steps + 1)]


def _lon_after(path, age, lon0=0.0):
    rot = pygplates.RotationModel(path).get_rotation(age, 101)
    return (rot * pygplates.PointOnSphere(0.0, lon0)).to_lat_lon()[1]


def test_elapsed_time_moves_the_starting_map_forward(tmp_path):
    path = str(tmp_path / 'elapsed.rot')
    write_rotation_file(_records(), path)
    assert np.isclose(_lon_after(path, 4.0), 4.0)          # 4 Myr after the start


def test_time_counted_back_moves_the_saved_step_back(tmp_path):
    path = str(tmp_path / 'back.rot')
    write_rotation_file(_records(), path, present_time=10.0)
    # at step 10 the point is at 10E; 4 Myr before that (model time 6) it was at 6E
    assert np.isclose(_lon_after(path, 4.0, lon0=10.0), 6.0)
    assert np.isclose(_lon_after(path, 10.0, lon0=10.0), 0.0)   # the starting map
    assert np.isclose(_lon_after(path, 0.0, lon0=10.0), 10.0)   # present day


def test_plates_gone_before_the_saved_step_are_left_out(tmp_path):
    path = str(tmp_path / 'back.rot')
    write_rotation_file(_records(steps=5), path, present_time=10.0)
    assert open(path).read().strip() == ''
