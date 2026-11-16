"""overlap: polarity draws each pair's downgoing plate before its overriding plate."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from self_driving_plates import burn_order


def test_without_polarity_the_order_is_highest_id_first():
    order, conflicts = burn_order([111, 444, 777], {})
    assert order == [777, 444, 111]
    assert conflicts == 0


def test_the_overriding_plate_is_drawn_after_the_downgoing_plate():
    # 111 subducts beneath 777: with overlap: plate_id 111 would be drawn last and keep the overlap
    order, conflicts = burn_order([111, 444, 777], {(111, 777): 111})
    assert order.index(111) < order.index(777)
    assert conflicts == 0


def test_pairs_with_a_missing_plate_are_ignored():
    order, conflicts = burn_order([111, 444], {(111, 999): 111, (111, 444): 444})
    assert order == [444, 111]
    assert conflicts == 0


def test_a_cycle_is_broken_and_counted():
    # 1 beneath 2, 2 beneath 3, 3 beneath 1: one pair must be drawn the wrong way round
    order, conflicts = burn_order([1, 2, 3], {(1, 2): 1, (2, 3): 2, (1, 3): 3})
    assert sorted(order) == [1, 2, 3]
    assert conflicts == 1
