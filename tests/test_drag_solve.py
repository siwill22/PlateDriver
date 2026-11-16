"""The coupled drag system is solved even where the LAPACK SVD fails to converge.

tests/data/drag_system_openblas_0334.npz is the system (M, rhs) at step 32 of a physical run
with slab_suction: ratio, on which numpy's lstsq raised "SVD did not converge" with
OpenBLAS 0.3.34 (conda-forge, November 2026); it is symmetric, with condition number 1.5e4.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from self_driving_plates import _lstsq  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def test_the_drag_system_that_broke_one_lapack_build_is_solved():
    d = np.load(os.path.join(HERE, 'data', 'drag_system_openblas_0334.npz'))
    M, rhs = d['M'], d['rhs']
    w = _lstsq(M, rhs)
    assert np.all(np.isfinite(w))
    assert np.linalg.norm(M @ w - rhs) < 1e-12 * np.linalg.norm(rhs) + 1e-15
    assert np.allclose(w, np.linalg.solve(M, rhs), rtol=1e-9, atol=1e-15)
