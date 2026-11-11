"""ELVIS stretch keys: one per distinct piece of line, none lost to a repeated id.

ELVIS sometimes stores one stretch as several features under the same
elvID + strekningLnr. Keyed on that pair alone, each piece overwrote the last
and only the one read last was routed and delivered.
"""

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, MultiLineString, box

from kantsoner.geometry import elvis

A = LineString([(0, 0), (100, 0)])
B = LineString([(100, 0), (200, 0)])
C = LineString([(300, 0), (400, 0)])


def _elvis(rows):
    """An elvenett frame from (elvID, strekningLnr, objektType, geometry)."""
    return gpd.GeoDataFrame(
        [
            {"elvID": e, "strekningLnr": s, "objektType": t, "geometry": g}
            for e, s, t, g in rows
        ],
        geometry="geometry",
        crs="EPSG:25833",
    )


def _keys(rows):
    """{stretch_id: (line, objektType)}"""
    return {
        sid: (line, feat["objektType"])
        for sid, line, feat in elvis._stretches(_elvis(rows))
    }


def test_a_unique_key_stays_bare():
    assert list(_keys([("1-1", "001", "ElvBekk", A)])) == ["1-1_001"]


def test_exact_copies_collapse_to_one_even_drawn_the_other_way():
    keys = _keys(
        [
            ("1-1", "001", "ElvBekk", A),
            ("1-1", "001", "ElvBekk", A),
            ("1-1", "001", "ElvBekk", A.reverse()),
        ]
    )
    assert list(keys) == ["1-1_001"]


def test_distinct_pieces_under_one_key_are_all_kept():
    keys = _keys(
        [
            ("1-1", "001", "ElvBekk", B),
            ("1-1", "001", "InnsjøMidtlinje", A),
            ("1-1", "001", "ElvBekk", C),
        ]
    )
    assert set(keys) == {"1-1_001_p0", "1-1_001_p1", "1-1_001_p2"}
    # numbered by position, not by the order the rows were read in, and each
    # piece keeps its own feature's attributes
    assert keys["1-1_001_p0"] == (A, "InnsjøMidtlinje")
    assert keys["1-1_001_p1"] == (B, "ElvBekk")
    assert keys["1-1_001_p2"] == (C, "ElvBekk")


def test_numbering_does_not_depend_on_read_order():
    rows = [
        ("1-1", "001", "ElvBekk", C),
        ("1-1", "001", "ElvBekk", A),
        ("1-1", "001", "ElvBekk", B),
    ]
    forward = {k: v[0].wkt for k, v in _keys(rows).items()}
    backward = {k: v[0].wkt for k, v in _keys(rows[::-1]).items()}
    assert forward == backward


def test_missing_strekninglnr_is_numbered_not_lumped():
    keys = _keys(
        [
            ("1-1", None, "ElvelinjeFiktiv", A),
            ("1-1", np.nan, "InnsjøMidtlinje", C),
            ("1-2", None, "ElvBekk", B),
        ]
    )
    # two different stretches of 1-1 with no strekningLnr stay apart, and a
    # lone one is numbered too, so no key reads as an official id
    assert set(keys) == {"1-1_p0", "1-1_p1", "1-2_p0"}
    assert keys["1-1_p0"] == (A, "ElvelinjeFiktiv")


def test_multipart_parts_are_numbered():
    keys = _keys([("1-1", "001", "ElvBekk", MultiLineString([C, A]))])
    assert {k: v[0] for k, v in keys.items()} == {
        "1-1_001_p0": A,
        "1-1_001_p1": C,
    }


def test_build_delivers_every_piece_of_a_repeated_key(monkeypatch):
    # The regression itself: two consecutive pieces under one id, in one
    # river. Before, only the piece read last came out as segments.
    monkeypatch.setenv("KANTSONER_OUTER_WORKERS", "100000")  # serial zoning
    river = gpd.GeoDataFrame(
        {"id": [0]}, geometry=[box(-20, -10, 420, 10)], crs="EPSG:25833"
    )
    net = _elvis(
        [
            ("1-1", "001", "ElvBekk", LineString([(200, 0), (400, 0)])),
            ("1-1", "001", "ElvBekk", LineString([(0, 0), (200, 0)])),
        ]
    )
    segs = elvis.build(net, river, []).segments
    assert set(segs["stretch_id"]) == {"1-1_001_p0", "1-1_001_p1"}
    assert segs["segment_id"].is_unique
    assert np.isclose(segs.geometry.length.sum(), 400.0)
