"""ELVIS patched with Vann-Nett: added lines, redrawn stretches, topology kept.

Small synthetic networks, each built for one situation found in the national
data: an exact copy, a Vann-Nett-only tributary joining an ELVIS stretch
mid-way, a tributary Vann-Nett draws elsewhere joining the main stem at its
own junction, a whole junction redrawn with an ELVIS-only stretch hanging off
it, Vann-Nett with no ELVIS at all, and a Vann-Nett branch ending beside a
line instead of on its end.
"""

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import LineString, MultiLineString

from kantsoner.geometry import vannnett_patch as vp

CRS = "EPSG:25833"


def _elvis(rows):
    """ELVIS frame from (elvID, strekningLnr, objektType, Strahler, line)."""
    return gpd.GeoDataFrame(
        [
            {
                "elvID": e,
                "strekningLnr": s,
                "objektType": t,
                "elveordenStrahler": o,
                "vassdragsnummer": "001.A",
                "elvenavn": None,
                "geometry": MultiLineString([g]),
            }
            for e, s, t, o, g in rows
        ],
        geometry="geometry",
        crs=CRS,
    )


def _vn(lines, vid="001-1-R"):
    """One Vann-Nett water body made of the given lines."""
    return gpd.GeoDataFrame(
        {"vannforekomstid": [vid], "vannforekomstnavn": ["Testbekken"]},
        geometry=[MultiLineString(lines)],
        crs=CRS,
    )


def _no_elvis():
    return gpd.GeoDataFrame(
        {
            c: []
            for c in (
                "elvID",
                "strekningLnr",
                "objektType",
                "elveordenStrahler",
                "vassdragsnummer",
                "elvenavn",
            )
        },
        geometry=[],
        crs=CRS,
    )


def _rows(net):
    """{key: [(line, patch tag or ""), ...]} over the patched network."""
    out = {}
    for _, r in net.iterrows():
        key = f"{r.elvID}_{r.strekningLnr}"
        tag = r.patch if isinstance(r.patch, str) else ""
        out.setdefault(key, []).append((shapely.get_parts(r.geometry)[0], tag))
    return out


def _nodes(net):
    """End point -> number of line ends there, to the millimetre."""
    g = [shapely.get_parts(x)[0] for x in net.geometry]
    ends = [
        tuple(np.round(c, 3)) for x in g for c in (x.coords[0], x.coords[-1])
    ]
    out = {}
    for e in ends:
        out[e] = out.get(e, 0) + 1
    return out


def test_an_exact_copy_is_left_alone():
    main = LineString([(0, 0), (100, 0), (200, 0)])
    res = vp.patch(_elvis([("1-1", "001", "ElvBekk", 2, main)]), _vn([main]))
    rows = _rows(res.network)
    assert list(rows) == ["1-1_001"]
    assert rows["1-1_001"] == [(main, "")]


def test_a_vannnett_only_tributary_is_added_and_its_junction_noded():
    main = LineString([(0, 0), (200, 0)])
    trib = LineString([(100, 100), (100, 0)])
    res = vp.patch(
        _elvis([("1-1", "001", "ElvBekk", 2, main)]),
        _vn(
            [
                LineString([(0, 0), (100, 0)]),
                LineString([(100, 0), (200, 0)]),
                trib,
            ]
        ),
    )
    rows = _rows(res.network)
    # the main stem is cut where the tributary now joins, keeping its key
    assert sorted(round(g.length) for g, _ in rows["1-1_001"]) == [100, 100]
    assert all(t == vp.SPLIT for _, t in rows["1-1_001"])
    added = {k: v for k, v in rows.items() if k.startswith("VN-001-1-R_")}
    assert len(added) == 1
    [(line, tag)] = next(iter(added.values()))
    assert line.equals(trib) and tag == vp.ADDED
    assert _nodes(res.network)[(100.0, 0.0)] == 3
    new = res.network[res.network.patch == vp.ADDED].iloc[0]
    assert new.objektType == "ElvBekk" and new.elveordenStrahler == 1
    assert new.patch_vn == "001-1-R"


def test_a_redrawn_tributary_joins_at_vannnetts_junction():
    # The main stem is the same line in both; Vann-Nett draws the tributary
    # 15-20 m east, joining the main stem 20 m further down.
    m1 = LineString([(0, 0), (500, 0)])
    m2 = LineString([(500, 0), (1000, 0)])
    trib = LineString([(500, 400), (500, 0)])
    vn_trib = LineString([(515, 400), (518, 200), (520, 0)])
    res = vp.patch(
        _elvis(
            [
                ("1-1", "001", "ElvBekk", 2, m1),
                ("1-1", "002", "ElvBekk", 2, m2),
                ("1-2", "003", "ElvBekk", 1, trib),
            ]
        ),
        _vn(
            [
                LineString([(0, 0), (520, 0)]),
                LineString([(520, 0), (1000, 0)]),
                vn_trib,
            ]
        ),
    )
    rows = _rows(res.network)
    [(t, tag)] = rows["1-2_003"]
    assert vp.REDRAWN in tag
    assert t.equals(vn_trib)  # Vann-Nett's line, source to junction
    # the main stem keeps its line and node, and is cut at the new junction
    assert rows["1-1_001"] == [(m1, "")]
    assert sorted(round(g.length) for g, _ in rows["1-1_002"]) == [20, 480]
    nodes = _nodes(res.network)
    assert nodes[(520.0, 0.0)] == 3 and nodes[(500.0, 0.0)] == 2
    assert not [k for k in rows if k.startswith("VN-")]  # nothing left over


def test_a_redrawn_junction_moves_and_elvis_only_lines_follow():
    # Vann-Nett draws all three lines 10 m north of ELVIS; the stretch
    # upstream of the main stem (a lake midline, say) is ELVIS's alone.
    up = LineString([(-300, 0), (0, 0)])
    m1 = LineString([(0, 0), (500, 0)])
    m2 = LineString([(500, 0), (1000, 0)])
    trib = LineString([(500, 400), (500, 0)])
    res = vp.patch(
        _elvis(
            [
                ("1-1", "000", "InnsjøMidtlinje", 2, up),
                ("1-1", "001", "ElvBekk", 2, m1),
                ("1-1", "002", "ElvBekk", 2, m2),
                ("1-2", "003", "ElvBekk", 1, trib),
            ]
        ),
        _vn(
            [
                LineString([(0, 10), (510, 10)]),
                LineString([(510, 10), (1000, 10)]),
                LineString([(510, 400), (510, 10)]),
            ]
        ),
    )
    rows = _rows(res.network)
    for key in ("1-1_001", "1-1_002", "1-2_003"):
        [(line, tag)] = rows[key]
        assert vp.REDRAWN in tag
        assert shapely.hausdorff_distance(line, line.offset_curve(0)) == 0
    assert rows["1-1_001"][0][0].equals(LineString([(0, 10), (510, 10)]))
    # the junction moved onto Vann-Nett's, and everything meets there
    assert _nodes(res.network)[(510.0, 10.0)] == 3
    # the ELVIS-only line is still there, its end pulled to the moved node
    [(lake, tag)] = rows["1-1_000"]
    assert tag == vp.MOVED
    assert lake.coords[0] == (-300.0, 0.0) and lake.coords[-1] == (0.0, 10.0)
    # the junction, both line ends upstream and the outlet
    assert len(res.qa["moved_nodes"]) == 4


def test_vannnett_with_no_elvis_is_added_with_its_order():
    lines = [
        LineString([(0, 100), (50, 50)]),
        LineString([(100, 100), (50, 50)]),
        LineString([(50, 50), (50, 0)]),
    ]
    net = vp.patch(_no_elvis(), _vn(lines)).network
    assert len(net) == 3 and (net.patch == vp.ADDED).all()
    assert net.strekningLnr.is_unique
    order = {
        shapely.get_parts(g)[0].coords[0]: o
        for g, o in zip(net.geometry, net.elveordenStrahler, strict=True)
    }
    assert order == {(0.0, 100.0): 1, (100.0, 100.0): 1, (50.0, 50.0): 2}


def test_every_elvis_stretch_survives():
    # whatever happens to the lines, no ELVIS key is lost
    m1 = LineString([(0, 0), (500, 0)])
    trib = LineString([(250, 300), (250, 0)])
    elvis = _elvis(
        [("1-1", "001", "ElvBekk", 2, m1), ("1-2", "002", "ElvBekk", 1, trib)]
    )
    vn = _vn(
        [LineString([(0, 20), (500, 20)]), LineString([(800, 0), (900, 0)])]
    )
    keys = set(_rows(vp.patch(elvis, vn).network))
    assert {"1-1_001", "1-2_002"} <= keys


def test_a_vannnett_branch_ending_beside_a_line_is_joined_to_it():
    main = LineString([(0, 0), (100, 0)])
    branch = LineString([(50, 50), (50, 0.3)])  # ends 0.3 m short of it
    out = vp.vn_lines(_vn([main, branch]))
    ends = [
        tuple(np.round(c, 6))
        for g in out.geometry
        for c in (g.coords[0], g.coords[-1])
    ]
    assert len(out) == 3
    assert ends.count((50.0, 0.0)) == 3


def test_tributaries_vannnett_joins_above_the_main_junction_all_connect():
    # ELVIS has one four-way junction; Vann-Nett (drawn 8 m north) joins two
    # of the tributaries first and runs them 30 m on to the main junction.
    j = (0, 0)
    t1 = LineString([(-300, 100), j])
    t2 = LineString([(-300, -100), j])
    t3 = LineString([(0, 300), j])
    out = LineString([j, (300, 0)])
    a, x = (-20, 8), (10, 8)
    vn = _vn(
        [
            LineString([(-300, 108), a]),
            LineString([(-300, -92), a]),
            LineString([a, x]),
            LineString([(10, 308), x]),
            LineString([x, (300, 8)]),
        ]
    )
    elvis = _elvis(
        [
            ("1-1", "001", "ElvBekk", 1, t1),
            ("1-2", "002", "ElvBekk", 1, t2),
            ("1-3", "003", "ElvBekk", 1, t3),
            ("1-0", "004", "ElvBekk", 2, out),
        ]
    )
    res = vp.patch(elvis, vn)
    rows = _rows(res.network)
    assert not len(res.qa["redraw_failed"])
    assert all(
        vp.REDRAWN in tag for pieces in rows.values() for _, tag in pieces
    )
    # both Vann-Nett junctions are nodes of three, nothing else dangles but
    # the four outer ends
    nodes = _nodes(res.network)
    assert nodes[(-20.0, 8.0)] == 3 and nodes[(10.0, 8.0)] == 3
    assert sorted(nodes.values()).count(1) == 4
    assert not [k for k in rows if k.startswith("VN-")]
