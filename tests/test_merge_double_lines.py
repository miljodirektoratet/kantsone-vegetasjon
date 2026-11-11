"""merge_delivery.py: the double lines in the lake shore, and the lines check.

On small synthetic shapes: the script is a standalone file at the root of the
repository, and these only call its pure functions.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import shapely

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import merge_delivery as md  # noqa: E402


def piece(x0, y0, x1, y1):
    return shapely.linestrings([[x0, y0], [x1, y1]])


def ring_pieces(x0, y, n, step=10.0):
    """n pieces end to end along a straight line, as one ring's shore."""
    return [piece(x0 + k * step, y, x0 + (k + 1) * step, y) for k in range(n)]


def test_hairline_slit_is_a_double_line_and_both_copies_go():
    # One run lake, two polygons that stop 2 cm apart: the same stretch twice, the
    # second copy drawn the other way round.
    ours = ring_pieces(0, 0, 4)
    other = [piece(x1, 0.02, x0, 0.02) for x0, x1 in [(10, 20), (20, 30)]]
    lines = np.array([*ours, *other])
    run = np.array(["L"] * len(lines), dtype=object)
    flag, whole = md.double_lines(lines, run, run)
    assert flag.tolist() == [False, True, True, False, True, True]
    assert whole == []


def test_copies_from_two_run_lakes_under_one_final_lake():
    a = ring_pieces(0, 0, 3)
    b = ring_pieces(0, 0.1, 3)
    lines = np.array([*a, *b])
    run = np.array(["A"] * 3 + ["B"] * 3, dtype=object)
    final = np.array(["F"] * 6, dtype=object)
    flag, whole = md.double_lines(lines, final, run)
    assert not flag.any() and whole == ["F"]  # nothing else is left of the lake
    # With some real shore around it, the pair goes and the rest stays.
    shore = [piece(50, 0, 60, 0), piece(60, 0, 60, 10)]
    lines = np.array([*a, *b, *shore])
    final = np.array(["F"] * 8, dtype=object)
    run = np.array(["A"] * 3 + ["B"] * 3 + ["A"] * 2, dtype=object)
    flag, whole = md.double_lines(lines, final, run)
    assert flag.tolist() == [True] * 6 + [False] * 2
    assert whole == []


def test_same_line_under_two_final_lakes_is_two_lakes_shore():
    # Two lakes sharing a border keep it, each its own: the rule is within one lake.
    a = ring_pieces(0, 0, 3)
    b = ring_pieces(0, 0.1, 3)
    lines = np.array([*a, *b])
    final = np.array(["F"] * 3 + ["G"] * 3, dtype=object)
    run = final.copy()
    flag, _ = md.double_lines(lines, final, run)
    assert not flag.any()


def test_touching_at_the_midpoint_only_is_real_shore():
    # A neck: another piece of the lake crosses the middle of the first.
    lines = np.array([piece(0, 0, 10, 0), piece(5, -5, 5, 5), piece(30, 0, 40, 0)])
    run = np.array(["L"] * 3, dtype=object)
    flag, _ = md.double_lines(lines, run, run)
    assert not flag.any()


def test_ring_neighbours_are_not_a_pair_even_when_short():
    # A ring's last piece, 0.8 m long: both its ends meet a neighbour. Every
    # station is within 0.5 m of a neighbour's end, but a neighbour is no pair.
    lines = np.array([piece(0, 0, 10, 0), piece(10, 0, 10.8, 0), piece(10.8, 0, 20, 0)])
    run = np.array(["L"] * 3, dtype=object)
    flag, _ = md.double_lines(lines, run, run)
    assert not flag.any()


def test_a_hairpin_is_a_line_drawn_twice_though_its_legs_meet_end_to_end():
    # A ring folding back on itself: B starts where A ends and runs back along it.
    lines = np.array([piece(0, 0, 10, 0), piece(10, 0, 0, 0.2), piece(100, 0, 110, 0)])
    run = np.array(["L"] * 3, dtype=object)
    flag, _ = md.double_lines(lines, run, run)
    assert flag.tolist() == [True, True, False]


def test_a_copy_another_region_keeps_is_not_tested_and_not_a_partner():
    a = ring_pieces(0, 0, 2)
    b = ring_pieces(0, 0.1, 2)
    shore = [piece(50, 0, 60, 0)]
    lines = np.array([*a, *b, *shore])
    run = np.array(["L"] * 5, dtype=object)
    skip = np.array([False] * 2 + [True] * 2 + [False])
    flag, _ = md.double_lines(lines, run, run, skip)
    assert not flag.any()


def test_a_lake_of_nothing_but_a_pair_is_left_whole():
    lines = np.array([piece(0, 0, 10, 0), piece(0, 0.05, 10, 0.05)])
    run = np.array(["L", "L"], dtype=object)
    flag, whole = md.double_lines(lines, run, run)
    assert not flag.any()
    assert whole == ["L"]


def test_the_longer_copy_goes_with_the_shorter():
    # The long copy's own stations lie beyond the short one; it still goes, as the
    # pair is a line through water, because the short one stood on it at all three.
    lines = np.array([piece(0, 0, 40, 0), piece(10, 0.1, 20, 0.1), piece(20, 0.1, 30, 0.1),
                      piece(100, 0, 110, 0)])
    run = np.array(["L"] * 4, dtype=object)
    flag, _ = md.double_lines(lines, run, run)
    assert flag.tolist() == [True, True, True, False]


def test_lines_check_leaves_the_pieces_in_lakes_out():
    table = pd.DataFrame(
        {
            "line_km": [10.0, 4.0, 2.0],
            "present_km": [9.0, 1.0, 0.0],
            "missing_km": [1.0, 3.0, 2.0],  # every piece not there
            "in_lake_km": [0.4, 3.0, 2.0],  # and of that, in a lake
        }
    )
    out = md.outside_lakes(table.copy())
    assert out["missing_km"].round(6).tolist() == [0.6, 0.0, 0.0]
    assert out["in_lake_km"].tolist() == [0.4, 3.0, 2.0]
    # present_pct is of the line outside lakes; a lake's whole line has none.
    assert out["present_pct"].iloc[0] == round(9.0 / 9.6 * 100, 3)
    assert out["present_pct"].iloc[1] == 100.0
    assert np.isnan(out["present_pct"].iloc[2])
    # the totals add up to what was missing before
    assert (out["missing_km"] + out["in_lake_km"]).round(6).tolist() == [1.0, 3.0, 2.0]


def test_pieces_in_lakes_by_their_middle():
    lake = shapely.box(0, 0, 10, 10)
    a = np.array([[2.0, 5.0], [8.0, 5.0], [20.0, 20.0]])
    b = np.array([[4.0, 5.0], [14.0, 5.0], [22.0, 20.0]])
    assert md.pieces_in_lakes(a, b, [lake]).tolist() == [True, False, False]


def test_the_line_in_the_log_and_the_readme():
    c = {"line_km": 1000.0, "in_lake_km": 100.0, "present_km": 750.0, "missing_km": 150.0}
    assert md.present_line(c) == (
        "Present: 83.333 % of 900.0 km; missing outside lakes 150.00 km; "
        "inside lakes 100.00 km (delivered as lakes)"
    )
