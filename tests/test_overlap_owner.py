"""overlap: polarity gives each overlapping cell to the overriding plate of the pair there."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from self_driving_plates import overlap_owner_by_polarity


def test_the_overriding_plate_keeps_the_overlap():
    # cells: plate 111 alone, 111 and 777 overlapping, 777 alone; 111 subducts beneath 777
    lo = np.array([111, 111, 777])
    hi = np.array([111, 777, 777])
    owner, overlap, by_polarity = overlap_owner_by_polarity(lo, hi, {(111, 777): 111})
    assert owner.tolist() == [111, 777, 777]
    assert overlap.tolist() == [False, True, False]
    assert by_polarity.tolist() == [False, True, False]


def test_without_a_polarity_the_lower_id_keeps_it():
    owner, _, by_polarity = overlap_owner_by_polarity(np.array([111]), np.array([777]), {})
    assert owner.tolist() == [111]
    assert not by_polarity.any()


def test_a_cycle_of_polarities_is_settled_pair_by_pair():
    # 1 beneath 2, 2 beneath 3, 3 beneath 1: each overlap goes to its own pair's overriding plate
    lo = np.array([1, 2, 1])
    hi = np.array([2, 3, 3])
    owner, _, _ = overlap_owner_by_polarity(lo, hi, {(1, 2): 1, (2, 3): 2, (1, 3): 3})
    assert owner.tolist() == [2, 3, 1]


def test_empty_cells_stay_empty():
    owner, overlap, _ = overlap_owner_by_polarity(np.array([0, 5]), np.array([0, 5]), {(5, 6): 5})
    assert owner.tolist() == [0, 5]
    assert not overlap.any()
