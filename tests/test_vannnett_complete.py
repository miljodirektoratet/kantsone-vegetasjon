"""Completing a patched network with the Vann-Nett line it still lacks.

Each network is built the way the patch hands it on (one line part per row,
``patch`` and ``patch_vn`` columns), for one situation found in the national
data: an ELVIS stretch drawn straight past a Vann-Nett meander, a second
channel beside a stretch already on Vann-Nett, a tributary Vann-Nett draws
30 m off with its own junction, a Vann-Nett branch with no ELVIS line, an
ELVIS junction inside a stretch Vann-Nett draws elsewhere, ELVIS running on
past a Vann-Nett source, and the anchors that sit on a row's end.
"""

from collections import defaultdict

import geopandas as gpd
import numpy as np
import shapely
import shapely.ops
from shapely.geometry import LineString, MultiLineString

from kantsoner.geometry import vannnett_complete as vc
from kantsoner.geometry import vannnett_patch as vp

CRS = "EPSG:25833"


def _net(rows):
    """Patched network from (elvID, strekningLnr, Strahler, line)."""
    return gpd.GeoDataFrame(
        [
            {
                "OBJECTID": k + 1,
                "objektType": "ElvBekk",
                "strekningLnr": s,
                "elvID": e,
                "elveordenStrahler": o,
                "vassdragsnummer": "001.A",
                "elvenavn": None,
                "SHAPE_Length": g.length,
                "elvelengde_m": round(g.length, 2),
                vp.PATCH_VN: None,
                vp.PATCH: None,
                "geometry": MultiLineString([g]),
            }
            for k, (e, s, o, g) in enumerate(rows)
        ],
        geometry="geometry",
        crs=CRS,
    )


def _vn(lines, vid="001-1-R"):
    return gpd.GeoDataFrame(
        {"vannforekomstid": [vid], "vannforekomstnavn": ["Testbekken"]},
        geometry=[MultiLineString(lines)],
        crs=CRS,
    )


def _rows(net):
    out = {}
    for _, r in net.iterrows():
        tag = r[vp.PATCH] if isinstance(r[vp.PATCH], str) else ""
        out.setdefault(f"{r.elvID}_{r.strekningLnr}", []).append(
            (shapely.get_parts(r.geometry)[0], tag)
        )
    return out


def _covered(net, vn, tol=1.0):
    lines = shapely.get_parts(net.geometry.values)
    m = vc.missing(vp.vn_lines(vn).geometry.values, lines, tol)
    return float((m.m1 - m.m0).sum())


MEANDER = LineString(
    [(0, 0), (100, 0), (120, 25), (150, 30), (180, 25), (200, 0), (300, 0)]
)


def test_a_complete_network_is_handed_back_unchanged():
    main = LineString([(0, 0), (300, 0)])
    net = _net([("1-1", "001", 2, main)])
    res = vc.complete(net, _vn([main]))
    assert res.stats["missing_km_before"] == 0
    assert res.network.geometry.values[0] is net.geometry.values[0]


def test_a_stretch_drawn_straight_past_a_meander_takes_the_meander():
    net = _net([("1-1", "001", 2, LineString([(0, 0), (300, 0)]))])
    vn = _vn([MEANDER])
    res = vc.complete(net, vn)
    rows = _rows(res.network)
    assert list(rows) == ["1-1_001"]
    [(line, tag)] = rows["1-1_001"]
    assert tag == vc.COMPLETED
    assert line.coords[0] == (0, 0) and line.coords[-1] == (300, 0)
    assert line.hausdorff_distance(MEANDER) < 1.0
    assert res.stats["missing_km_after"] == 0
    assert res.network[vp.PATCH_VN].iloc[0] == "001-1-R"


def test_a_second_channel_beside_a_stretch_on_vannnett_is_added():
    main = LineString([(0, 0), (300, 0)])
    side = LineString([(100, 0), (120, 25), (150, 30), (180, 25), (200, 0)])
    net = _net([("1-1", "001", 2, main)])
    res = vc.complete(net, _vn([main, side]))
    rows = _rows(res.network)
    # the main stem keeps its line, cut where the channel comes within 1 m
    cuts = sorted(g.length for g, _ in rows["1-1_001"])
    assert len(cuts) == 3 and abs(sum(cuts) - 300) < 1e-6
    assert all(abs(c - 100) < 2 for c in cuts)
    assert all(t == vc.CUT for _, t in rows["1-1_001"])
    added = [v for k, v in rows.items() if k.startswith("VN-001-1-R_")]
    assert len(added) == 1
    [(line, tag)] = added[0]
    assert tag == vp.ADDED and line.hausdorff_distance(side) < 1.0
    assert res.stats["missing_km_after"] == 0


def test_a_tributary_drawn_elsewhere_moves_to_vannnetts_junction():
    m1 = LineString([(0, 0), (100, 0)])
    m2 = LineString([(100, 0), (400, 0)])
    trib = LineString([(100, 200), (100, 0)])  # ELVIS
    vtrib = LineString([(130, 200), (130, 0)])  # Vann-Nett, 30 m east
    net = _net(
        [
            ("1-1", "001", 2, m1),
            ("1-1", "002", 2, m2),
            ("1-2", "001", 1, trib),
        ]
    )
    vn = _vn(
        [
            LineString([(0, 0), (130, 0)]),
            LineString([(130, 0), (400, 0)]),
            vtrib,
        ]
    )
    res = vc.complete(net, vn)
    rows = _rows(res.network)
    [(t, tag)] = rows["1-2_001"]
    assert tag == vc.COMPLETED
    assert t.hausdorff_distance(vtrib) < 1.0
    assert t.coords[-1] == (130, 0)
    # the main stem is cut where the tributary now joins
    cuts = sorted(round(g.length) for g, _ in rows["1-1_002"])
    assert cuts == [30, 270]
    assert rows["1-1_001"] == [(m1, "")]  # untouched
    assert res.stats["missing_km_after"] == 0


def test_a_vannnett_branch_with_no_elvis_line_is_added_and_noded():
    main = LineString([(0, 0), (400, 0)])
    branch = LineString([(200, 150), (200, 0)])
    net = _net(
        [
            ("1-1", "001", 2, main),
            ("VN-001-9-R", "VN0000007", 1, LineString([(0, 50), (0, 0)])),
        ]
    )
    vn = _vn(
        [
            LineString([(0, 0), (200, 0)]),
            LineString([(200, 0), (400, 0)]),
            branch,
        ]
    )
    res = vc.complete(net, vn)
    rows = _rows(res.network)
    new = [k for k in rows if k.startswith("VN-001-1-R_")]
    assert new == ["VN-001-1-R_VN0000008"]  # numbering goes on
    [(line, tag)] = rows[new[0]]
    assert (
        tag == vp.ADDED
        and line.coords[0] == (200, 150)
        and line.coords[-1] == (200, 0)
    )
    assert sorted(round(g.length) for g, _ in rows["1-1_001"]) == [200, 200]
    added = res.network[res.network[vp.PATCH] == vp.ADDED].iloc[0]
    assert added["elveordenStrahler"] == 1
    assert res.stats["missing_km_after"] == 0


def test_rows_nothing_happens_to_are_byte_identical():
    far = LineString([(5000, 5000), (5300, 5000)])
    net = _net(
        [
            ("1-1", "001", 2, LineString([(0, 0), (300, 0)])),
            ("2-1", "001", 1, far),
        ]
    )
    res = vc.complete(net, _vn([MEANDER, far]))
    out = res.network
    k = int(np.flatnonzero(out.elvID == "2-1")[0])
    assert out.geometry.values[k] is net.geometry.values[1]
    assert out[vp.PATCH].iloc[k] is None
    assert _covered(out, _vn([MEANDER, far])) == 0


# ── topology checks ────────────────────────────────────────────────────────
def _nodes(net):
    """node (mm) -> [inflow, outflow]"""
    out = defaultdict(lambda: [0, 0])
    for g in shapely.get_parts(net.geometry.values):
        c = shapely.get_coordinates(g)
        out[tuple(np.round(c[-1], 3))][0] += 1
        out[tuple(np.round(c[0], 3))][1] += 1
    return out


def _conflicts(net):
    """Nodes of two or more rows where nothing flows in or nothing out."""
    return {
        xy for xy, (i, o) in _nodes(net).items() if i + o >= 2 and not i * o
    }


def _dangling(net):
    return {xy for xy, (i, o) in _nodes(net).items() if i + o == 1}


def _straights(net, vn, old, longer=vc.JOG_M):
    """Lengths of straight segments longer than ``longer`` drawn on neither
    Vann-Nett nor the old network (connectors and jogs)."""
    ref = shapely.union_all(
        np.r_[
            shapely.get_parts(vp.vn_lines(vn).geometry.values),
            shapely.get_parts(old.geometry.values),
        ]
    ).buffer(0.01)
    out = []
    for g in shapely.get_parts(net.geometry.values):
        c = shapely.get_coordinates(g)
        for a, b in zip(c[:-1], c[1:], strict=True):
            seg = LineString([a, b])
            if seg.length > longer + 1e-9 and seg.difference(ref).length > 1e-3:
                out.append(round(seg.length, 3))
    return sorted(out)


def test_a_path_whose_anchors_lie_on_nodes_keeps_the_row_ends_there():
    a = LineString([(0, 0), (100, 0)])
    b = LineString([(100, 0), (200, 0)])
    c = LineString([(200, 0), (300, 0)])
    net = _net(
        [("1-1", "001", 2, a), ("1-1", "002", 2, b), ("1-1", "003", 2, c)]
    )
    vn = _vn(
        [
            LineString(
                [(0, 0), (100, 0), (100, 30), (200, 30), (200, 0), (300, 0)]
            )
        ]
    )
    res = vc.complete(net, vn)
    assert _dangling(res.network) == _dangling(net)
    assert not _conflicts(res.network)
    [(line, tag)] = _rows(res.network)["1-1_002"]
    assert tag == vc.COMPLETED and line.coords[0] == (100, 0)
    assert line.coords[-1] == (200, 0)


def test_a_junction_inside_a_duplicate_moves_onto_vannnett_with_its_rows():
    """ELVIS main stem m1 + m2 meet at N, where an ELVIS-only tributary t
    joins; Vann-Nett draws the main stem 30 m south. The stem takes
    Vann-Nett's line, N goes onto it and t follows with a connector."""
    m1 = LineString([(0, 0), (150, 0)])
    m2 = LineString([(150, 0), (400, 0)])
    t = LineString([(150, 200), (150, 0)])
    net = _net(
        [("1-1", "001", 2, m1), ("1-1", "002", 2, m2), ("1-2", "001", 1, t)]
    )
    vm = LineString(
        [(0, 0), (100, 0), (110, -30), (290, -30), (300, 0), (400, 0)]
    )
    res = vc.complete(net, _vn([vm]))
    rows = _rows(res.network)
    [(l1, g1)] = rows["1-1_001"]
    [(l2, g2)] = rows["1-1_002"]
    [(lt, gt)] = rows["1-2_001"]
    assert g1 == g2 == vc.COMPLETED and gt == vc.MOVED
    assert l1.coords[-1] == l2.coords[0] == lt.coords[-1]
    assert abs(l1.coords[-1][1] + 30) < 1e-6  # on Vann-Nett's line
    assert lt.coords[-2] == (150, 0)  # t keeps its own line, plus 30 m
    assert _straights(res.network, _vn([vm]), net) == [30.0]
    assert len(res.qa["connectors"]) == 1
    assert not _conflicts(res.network)
    assert res.stats["missing_km_after"] == 0
    assert not len(res.qa.get("kept_beside", []))


def test_a_tributary_vannnett_draws_into_the_junction_takes_that_drawing():
    """As above, but Vann-Nett has the tributary too: it leaves t at
    (150, 40) and joins the main stem at J = (180, -30). t takes that line,
    N goes to J, and no connector is drawn."""
    m1 = LineString([(0, 0), (150, 0)])
    m2 = LineString([(150, 0), (400, 0)])
    t = LineString([(150, 200), (150, 0)])
    net = _net(
        [("1-1", "001", 2, m1), ("1-1", "002", 2, m2), ("1-2", "001", 1, t)]
    )
    vm1 = LineString([(0, 0), (100, 0), (110, -30), (180, -30)])
    vm2 = LineString([(180, -30), (290, -30), (300, 0), (400, 0)])
    vt = LineString([(150, 200), (150, 40), (180, -30)])
    res = vc.complete(net, _vn([vm1, vm2, vt]))
    rows = _rows(res.network)
    [(lt, gt)] = rows["1-2_001"]
    assert gt == vc.COMPLETED
    assert lt.coords[-1] == rows["1-1_001"][0][0].coords[-1]
    assert np.hypot(lt.coords[-1][0] - 180, lt.coords[-1][1] + 30) < 1e-6
    assert lt.hausdorff_distance(vt) < 1.0
    assert not _straights(res.network, _vn([vm1, vm2, vt]), net)
    assert not len(res.qa.get("connectors", []))
    assert not _conflicts(res.network)
    assert res.stats["missing_km_after"] == 0


def test_elvis_running_on_past_a_vannnett_source_keeps_its_upper_course():
    """Vann-Nett's stream starts free at F = (300, 25), 25 m beside ELVIS,
    which runs on upstream to x = 600. The row below F takes Vann-Nett's
    line, its upper course joins F by a connector: no line twice."""
    r = LineString([(600, 0), (0, 0)])
    net = _net([("1-1", "001", 2, r)])
    v = LineString([(300, 25), (200, 25), (100, 0), (0, 0)])
    res = vc.complete(net, _vn([v]))
    rows = _rows(res.network)
    [(line, tag)] = rows["1-1_001"]
    assert tag == vc.COMPLETED
    assert line.coords[0] == (600, 0) and line.coords[-1] == (0, 0)
    assert (300, 25) in line.coords
    assert _straights(res.network, _vn([v]), net) == [25.0]
    assert len(res.qa["connectors"]) == 1
    assert not [k for k in rows if k.startswith("VN-")]
    assert res.stats["missing_km_after"] == 0


def test_a_tributary_crossing_a_row_is_no_detour_of_it():
    """The Vann-Nett tributary crosses the ELVIS main stem m2 at x = 200 on
    its way to the Vann-Nett junction J 30 m south: m2 must not be redrawn
    down the tributary. The main stem takes Vann-Nett's, the tributary is
    one addition into J."""
    m1 = LineString([(0, 0), (150, 0)])
    m2 = LineString([(150, 0), (400, 0)])
    t = LineString([(150, 200), (150, 0)])
    net = _net(
        [("1-1", "001", 2, m1), ("1-1", "002", 2, m2), ("1-2", "001", 1, t)]
    )
    ms1 = LineString([(0, 0), (100, 0), (110, -30), (200, -30)])
    ms2 = LineString([(200, -30), (290, -30), (300, 0), (400, 0)])
    vt = LineString([(200, 200), (200, -30)])
    res = vc.complete(net, _vn([ms1, ms2, vt, t]))
    rows = _rows(res.network)
    for _, (g, tag) in [(0, x) for x in rows["1-1_002"]]:
        assert min(y for _, y in g.coords) < -29  # on Vann-Nett's stem
        assert vc.COMPLETED in tag
    added = [v for k, v in rows.items() if k.startswith("VN-")]
    assert len(added) == 1 and len(added[0]) == 1
    ((g, _),) = added[0]
    assert g.coords[0] == (200, 200) and g.coords[-1] == (200, -30)
    assert not _conflicts(res.network)
    assert res.stats["missing_km_after"] == 0


def test_a_loop_off_one_point_of_a_row_is_added_beside_it():
    """Vann-Nett leaves the row and comes back to the same point: the row
    stays straight, the loop is added, out and back."""
    r = LineString([(0, 0), (100, 0)])
    net = _net([("1-1", "001", 2, r)])
    v1 = LineString([(0, 0), (50, 0)])
    loop = LineString([(50, 0), (50, 20), (60, 20), (50.2, 0)])
    v2 = LineString([(50.2, 0), (100, 0)])
    res = vc.complete(net, _vn([v1, loop, v2]))
    rows = _rows(res.network)
    for g, _ in rows["1-1_001"]:
        assert g.is_simple and max(abs(y) for _, y in g.coords) < 1e-9
    assert len([k for k in rows if k.startswith("VN-")]) == 2
    assert res.stats["missing_km_after"] == 0


def test_two_detours_meeting_on_one_row_join_without_a_notch():
    """Vann-Nett touches the row only briefly between two meanders: the row
    follows Vann-Nett through, never back to its old line between them."""
    r = LineString([(0, 0), (400, 0)])
    net = _net([("1-1", "001", 2, r)])
    v = LineString(
        [
            (0, 0),
            (100, 0),
            (120, 30),
            (180, 30),
            (199.5, 0.5),
            (200.5, 0.5),
            (220, 30),
            (280, 30),
            (300, 0),
            (400, 0),
        ]
    )
    res = vc.complete(net, _vn([v]))
    [(line, tag)] = _rows(res.network)["1-1_001"]
    assert tag == vc.COMPLETED
    assert line.hausdorff_distance(v) < 1.0
    ys = [y for x, y in line.coords if 150 < x < 250]
    assert min(ys) > 0.4  # no point back on the old line at y = 0


def test_a_bump_under_the_noise_limit_is_nothing_missing():
    r = LineString([(0, 0), (100, 0)])
    net = _net([("1-1", "001", 2, r)])
    vn = _vn(
        [LineString([(0, 0), (49.9, 0), (50, 1.0003), (50.1, 0), (100, 0)])]
    )
    res = vc.complete(net, vn)
    assert res.stats["passes"] == 1 and res.stats["missing_km_before"] == 0.0


def test_passes_stop_once_a_pass_leaves_nothing_missing():
    """Pass 1 has work to do (a 1.5 m bump the row takes as a detour) and
    leaves nothing missing: it is the last."""
    r = LineString([(0, 0), (100, 0)])
    net = _net([("1-1", "001", 2, r)])
    vn = _vn([LineString([(0, 0), (45, 0), (50, 1.5), (55, 0), (100, 0)])])
    res = vc.complete(net, vn, max_passes=4)
    assert res.stats["missing_km_before"] > 0.0
    assert res.stats["passes"] == 1 and res.stats["missing_km_after"] == 0.0


def test_an_addition_with_elvis_beside_it_is_listed():
    """No rule ties the ELVIS row to the Vann-Nett line beside it (the
    Vann-Nett line starts on another row): the QA says so."""
    a = LineString([(0, 0), (300, 0)])
    q = LineString([(45, -50), (85, 10)])
    net = _net([("1-1", "001", 2, a), ("2-1", "001", 1, q)])
    v = LineString([(65, -20), (290, -20), (300, 0)])
    res = vc.complete(net, _vn([v, q, LineString([(300, 0), (400, 0)])]))
    kb = res.qa["kept_beside"]
    assert len(kb) == 1 and kb.iloc[0]["stretch"] == "1-1_001"


def test_an_earlier_addition_with_a_loose_end_goes_into_the_elvis_line_beside():
    """An earlier patch added Vann-Nett from (300, 30) on (row X, its start
    loose) beside ELVIS; the rest of Vann-Nett is missing. The ELVIS row
    beside both takes the whole Vann-Nett line and X goes."""
    e1 = LineString([(0, 0), (600, 0)])
    e2 = LineString([(600, 0), (800, 0)])
    v = LineString([(0, 0), (10, 0), (50, 30), (550, 30), (600, 0), (800, 0)])
    x = shapely.ops.substring(
        v, v.project(shapely.Point(300, 30)), v.project(shapely.Point(600, 0))
    )
    net = _net(
        [
            ("1-1", "001", 2, e1),
            ("1-1", "002", 2, e2),
            ("VN-001-1-R", "VN0000001", 2, x),
        ]
    )
    net.loc[2, vp.PATCH] = vp.ADDED
    res = vc.complete(net, _vn([v]))
    rows = _rows(res.network)
    assert "VN-001-1-R_VN0000001" not in rows
    [(line, tag)] = rows["1-1_001"]
    assert tag == vc.COMPLETED
    assert line.coords[0] == (0, 0) and line.coords[-1] == (600, 0)
    assert (
        line.hausdorff_distance(
            shapely.ops.substring(v, 0, v.project(shapely.Point(600, 0)))
        )
        < 1.0
    )
    assert len(res.qa["dropped"]) == 1
    assert res.stats["missing_km_after"] == 0
    assert not _conflicts(res.network)
