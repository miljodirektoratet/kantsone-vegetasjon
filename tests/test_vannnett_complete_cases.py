"""Round-3 counter-examples to the completion, kept as regression tests.

Each builds a small network the way the patch hands it on and asserts what
the module promises: no line drawn twice, no out-and-back, no straight line
longer than JOG_M but the listed connectors, every junction noded, no key
reused. They were written by an independent review of an earlier version.
"""

from collections import defaultdict

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import shapely
import shapely.ops
from shapely.geometry import LineString, MultiLineString, Point

from kantsoner.geometry import vannnett_complete as vc
from kantsoner.geometry import vannnett_patch as vp

CRS = "EPSG:25833"


def net_(rows):
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


def vn_(lines, vid="001-1-R"):
    return gpd.GeoDataFrame(
        {"vannforekomstid": [vid], "vannforekomstnavn": ["Testbekken"]},
        geometry=[MultiLineString(lines)],
        crs=CRS,
    )


def _lines(net):
    return shapely.get_parts(net.geometry.values)


def _keys(net):
    return (
        net["elvID"].astype(str) + "_" + net["strekningLnr"].astype(str)
    ).to_list()


def nodes(net):
    """node(mm) -> [in, out, rows]"""
    out = defaultdict(lambda: [0, 0, []])
    for k, g in zip(_keys(net), _lines(net), strict=True):
        c = shapely.get_coordinates(g)
        a = tuple(np.round(c[0], 3))
        b = tuple(np.round(c[-1], 3))
        out[a][1] += 1
        out[a][2].append(k + ":s")
        out[b][0] += 1
        out[b][2].append(k + ":e")
    return out


def dangling(net):
    return sorted(xy for xy, (i, o, _) in nodes(net).items() if i + o == 1)


def t_ends(net):
    """row ends lying on another row's interior (unnoded)"""
    g = _lines(net)
    keys = _keys(net)
    tree = shapely.STRtree(g)
    out = []
    for i, gi in enumerate(g):
        for w, p in (
            (0, shapely.get_point(gi, 0)),
            (1, shapely.get_point(gi, -1)),
        ):
            for j in tree.query(p, predicate="dwithin", distance=1e-6):
                if j == i:
                    continue
                m = g[j].project(p)
                if 1e-6 < m < g[j].length - 1e-6:
                    out.append(
                        (
                            keys[i],
                            "start" if w == 0 else "end",
                            keys[j],
                            tuple(np.round(shapely.get_coordinates(p)[0], 3)),
                        )
                    )
    return out


def long_straights(net, vn, old, jog=2.0):
    """segments longer than jog in rows whose midpoint and ends are not on VN or OLD lines"""
    V = shapely.get_parts(vp.vn_lines(vn).geometry.values)
    old_lines = _lines(old)
    ref = shapely.union_all(np.r_[V, old_lines])
    out = []
    for k, g in zip(_keys(net), _lines(net), strict=True):
        c = shapely.get_coordinates(g)
        for a, b in zip(c[:-1], c[1:], strict=True):
            L = float(np.hypot(*(b - a)))
            if L <= jog + 1e-9:
                continue
            seg = LineString([a, b])
            if seg.difference(ref.buffer(0.01)).length > 1e-3:
                out.append(
                    (
                        k,
                        round(L, 2),
                        tuple(np.round(a, 2)),
                        tuple(np.round(b, 2)),
                    )
                )
    return out


def spikes(net, min_seg=1.0, max_turn=170.0):
    out = []
    for k, g in zip(_keys(net), _lines(net), strict=True):
        c = shapely.get_coordinates(g)
        d = np.diff(c, axis=0)
        n = np.hypot(d[:, 0], d[:, 1])
        for i in range(len(d) - 1):
            if n[i] >= min_seg and n[i + 1] >= min_seg:
                cos = np.dot(d[i], d[i + 1]) / n[i] / n[i + 1]
                if np.degrees(np.arccos(np.clip(cos, -1, 1))) > max_turn:
                    out.append((k, tuple(np.round(c[i + 1], 2))))
    return out


def components(net):
    """number of connected components of the network (rows joined at shared ends, mm)"""
    parent = {}

    def find(a):
        while parent.setdefault(a, a) != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for g in _lines(net):
        c = shapely.get_coordinates(g)
        a, b = tuple(np.round(c[0], 3)), tuple(np.round(c[-1], 3))
        parent[find(a)] = find(b)
    return len({find(a) for a in list(parent)})


def overlaps(net, min_len=0.5):
    """pairs of rows (different keys or pieces) sharing a stretch of line longer than min_len"""
    g = _lines(net)
    keys = _keys(net)
    tree = shapely.STRtree(g)
    out = []
    for i, gi in enumerate(g):
        for j in tree.query(gi, predicate="intersects"):
            if j <= i:
                continue
            x = gi.intersection(g[j])
            L = sum(
                p.length
                for p in shapely.get_parts(x)
                if p.geom_type in ("LineString", "MultiLineString")
            )
            if L > min_len:
                out.append((keys[i], keys[int(j)], round(L, 3)))
    return out


def crossings(net):
    """pairs of rows whose lines cross (intersect at a point that is an end of neither)"""
    g = _lines(net)
    keys = _keys(net)
    tree = shapely.STRtree(g)
    out = []
    for i, gi in enumerate(g):
        ends_i = shapely.get_coordinates(gi)[[0, -1]]
        for j in tree.query(gi, predicate="intersects"):
            if j <= i:
                continue
            gj = g[int(j)]
            ends_j = shapely.get_coordinates(gj)[[0, -1]]
            x = gi.intersection(gj)
            for p in shapely.get_coordinates(x):
                d = min(
                    np.hypot(*(ends_i - p).T).min(),
                    np.hypot(*(ends_j - p).T).min(),
                )
                if d > 1e-6:
                    out.append((keys[i], keys[int(j)], tuple(np.round(p, 3))))
    return out


# ── the networks ───────────────────────────────────────────────────────────
def _trib_through_node(order, dx):
    """ELVIS m1, m2 (main stem) and t (tributary) meet at N=(150, 0); Vann-Nett
    draws the main stem 30 m south (junction J=(150+dx, -30)) and the
    tributary on t down to N, then on across the old main stem to J. The
    tributary's missing piece starts ~1 m south of N, nearest m1/m2 (or a
    tie), so _junction_line (e[1] != r_x) never gives it to t."""
    m1 = LineString([(0, 0), (150, 0)])
    m2 = LineString([(150, 0), (400, 0)])
    t = LineString([(150, 200), (150, 0)])
    rows = {
        "m1": ("1-1", "001", 2, m1),
        "m2": ("1-1", "002", 2, m2),
        "t": ("1-2", "001", 1, t),
    }
    net = net_([rows[k] for k in order])
    j = (150 + dx, -30)
    ms1 = LineString([(0, 0), (100, 0), (110, -30), j])
    ms2 = LineString([j, (290, -30), (300, 0), (400, 0)])
    vt = LineString([(150 + dx, 200), (150 + dx, 0), j])
    vn = vn_([ms1, ms2, vt])
    return net, vn, vc.complete(net, vn)


def two_additions_meet_on_a_redrawn_row():
    """Detour [M1, M2] through Vann-Nett junction J=(200, 30); two Vann-Nett
    tributaries (free sources, no ELVIS) end at J. _attach_all skips an
    'inner' end when two additions end there, so the row (now through J) is
    never cut at J."""
    r = LineString([(0, 0), (400, 0)])
    net = net_([("1-1", "001", 2, r)])
    m1 = LineString([(0, 0), (100, 0), (110, 30), (200, 30)])
    m2 = LineString([(200, 30), (290, 30), (300, 0), (400, 0)])
    t1 = LineString([(100, 200), (200, 30)])
    t2 = LineString([(350, 200), (200, 30)])
    vn = vn_([m1, m2, t1, t2])
    return net, vn, vc.complete(net, vn)


def two_extra_rows_at_one_node():
    """ELVIS main stem m1+m2 meets two ELVIS-only tributaries at N=(150, 0);
    Vann-Nett draws the main stem 30 m south. Both tributaries get the same
    connector N -> (150, -30): two rows drawn on one 30 m segment."""
    m1 = LineString([(0, 0), (150, 0)])
    m2 = LineString([(150, 0), (400, 0)])
    t1 = LineString([(100, 200), (150, 0)])
    t2 = LineString([(200, 200), (150, 0)])
    net = net_(
        [
            ("1-1", "001", 2, m1),
            ("1-1", "002", 2, m2),
            ("1-2", "001", 1, t1),
            ("1-3", "001", 1, t2),
        ]
    )
    vm = LineString(
        [(0, 0), (100, 0), (110, -30), (290, -30), (300, 0), (400, 0)]
    )
    vn = vn_([vm])
    return net, vn, vc.complete(net, vn)


def trib_through_node_tie():
    return _trib_through_node(["m1", "m2", "t"], 0.0)


def trib_through_node_tie_t_first():
    return _trib_through_node(["t", "m1", "m2"], 0.0)


def trib_through_node_offset():
    return _trib_through_node(["t", "m1", "m2"], 0.5)


def tributary_chain_first():
    """Main stem m0+m1+m2, tributary t joins at N=(200, 0). Vann-Nett: main
    stem L1 (leaves m0 at x=50, 30 m south) to J=(240, -30), L2 on to (400, 0);
    tributary lt leaves t at (200, 290), runs 40 m east of it, into J. One
    water body: [lt, L2] (469 m) goes before [L1, L2] (343 m). m1's end gets a
    40 m connector onto the tributary's line and L1 is added beside m0+m1."""
    m0 = LineString([(0, 0), (100, 0)])
    m1 = LineString([(100, 0), (200, 0)])
    m2 = LineString([(200, 0), (400, 0)])
    t = LineString([(200, 300), (200, 0)])
    net = net_(
        [
            ("1-1", "001", 2, m0),
            ("1-1", "002", 2, m1),
            ("1-1", "003", 2, m2),
            ("1-2", "001", 1, t),
        ]
    )
    j = (240, -30)
    l1 = LineString([(0, 0), (50, 0), (60, -30), j])
    l2 = LineString([j, (340, -30), (350, 0), (400, 0)])
    lt = LineString([(200, 300), (200, 290), (240, 250), j])
    vn = vn_([l1, l2, lt])
    return net, vn, vc.complete(net, vn)


def tributary_chain_other_water_body():
    """As tributary_chain_first, the tributary in its own water body: the
    main stem goes first (control)."""

    m0 = LineString([(0, 0), (100, 0)])
    m1 = LineString([(100, 0), (200, 0)])
    m2 = LineString([(200, 0), (400, 0)])
    t = LineString([(200, 300), (200, 0)])
    net = net_(
        [
            ("1-1", "001", 2, m0),
            ("1-1", "002", 2, m1),
            ("1-1", "003", 2, m2),
            ("1-2", "001", 1, t),
        ]
    )
    j = (240, -30)
    l1 = LineString([(0, 0), (50, 0), (60, -30), j])
    l2 = LineString([j, (340, -30), (350, 0), (400, 0)])
    lt = LineString([(200, 300), (200, 290), (240, 250), j])
    vn = pd.concat([vn_([l1, l2]), vn_([lt], vid="001-2-R")], ignore_index=True)
    return net, vn, vc.complete(net, vn)


def connector_crosses_an_addition():
    """ELVIS-only tributary t at N=(150, 0) gets a connector to (150, -40) on
    the Vann-Nett main stem; a Vann-Nett tributary (free source, no ELVIS)
    runs into the main stem at J2=(170, -40) across that connector. _crosses
    only looks at the original rows and the chain."""
    m1 = LineString([(0, 0), (150, 0)])
    m2 = LineString([(150, 0), (400, 0)])
    t = LineString([(150, 200), (150, 0)])
    net = net_(
        [("1-1", "001", 2, m1), ("1-1", "002", 2, m2), ("1-2", "001", 1, t)]
    )
    j2 = (170, -40)
    ms1 = LineString([(0, 0), (100, 0), (110, -40), j2])
    ms2 = LineString([j2, (290, -40), (300, 0), (400, 0)])
    v2 = LineString([(120, -15), j2])
    vn = vn_([ms1, ms2, v2])
    return net, vn, vc.complete(net, vn)


def merged_detours_vn_in_two_parts():
    """Two meanders of one row, Vann-Nett touching the row (0.2 m off) for
    1.2 m between them, split into two parts there: _merge_reps finds no
    single part with both ends and joins them straight (2.24 m)."""
    r = LineString([(0, 0), (400, 0)])
    net = net_([("1-1", "001", 2, r)])
    p1 = LineString(
        [(0, 0), (100, 0), (120, 30), (180, 30), (199.4, 0.2), (200, 0.2)]
    )
    p2 = LineString(
        [(200, 0.2), (200.6, 0.2), (220, 30), (280, 30), (300, 0), (400, 0)]
    )
    vn = vn_([p1, p2])
    return net, vn, vc.complete(net, vn)


def overlapping_detours_loose_anchor_unlisted():
    """(round-2 overlapping_detours_one_row) The narrower channel's anchors
    lie on the row's old line, which the wider detour moved away: _bridge
    turns both into 'free' ends and no QA layer lists them."""
    r = LineString([(0, 0), (400, 0)])
    net = net_([("1-1", "001", 2, r)])
    x = LineString([(0, 0), (100, 0), (110, 30), (290, 30), (300, 0), (400, 0)])
    y = LineString([(150, 0), (160, -20), (240, -20), (250, 0)])
    vn = vn_([x, y])
    return net, vn, vc.complete(net, vn)


def main_stem_crosses_elvis_tributary():
    """Vann-Nett draws the main stem 30 m north, on the tributary's side, and
    the tributary on ELVIS t down to (200, 30) where it joins. The main stem
    crosses t there: _crossings pairs M1 (ends on t) with M2 (starts on t),
    both steep, so both chains are dropped and M1, M2 are added beside the
    kept ELVIS main stem."""
    m1 = LineString([(0, 0), (200, 0)])
    m2 = LineString([(200, 0), (400, 0)])
    t = LineString([(200, 300), (200, 0)])
    net = net_(
        [("1-1", "001", 2, m1), ("1-1", "002", 2, m2), ("1-2", "001", 1, t)]
    )
    vm1 = LineString([(0, 0), (100, 0), (110, 30), (200, 30)])
    vm2 = LineString([(200, 30), (290, 30), (300, 0), (400, 0)])
    vt = LineString([(200, 300), (200, 30)])
    vn = vn_([vm1, vm2, vt])
    return net, vn, vc.complete(net, vn)


def main_stem_crosses_elvis_tributary_junction_beyond():
    """As above, Vann-Nett's tributary leaving t at (200, 290) and joining the
    main stem at J=(240, 30), east of where the main stem crosses t."""
    m0 = LineString([(0, 0), (100, 0)])
    m1 = LineString([(100, 0), (200, 0)])
    m2 = LineString([(200, 0), (400, 0)])
    t = LineString([(200, 300), (200, 0)])
    net = net_(
        [
            ("1-1", "001", 2, m0),
            ("1-1", "002", 2, m1),
            ("1-1", "003", 2, m2),
            ("1-2", "001", 1, t),
        ]
    )
    j = (240, 30)
    l1 = LineString([(0, 0), (50, 0), (60, 30), j])
    l2 = LineString([j, (340, 30), (350, 0), (400, 0)])
    lt = LineString([(200, 300), (200, 290), (240, 250), j])
    vn = vn_([l1, l2, lt])
    return net, vn, vc.complete(net, vn)


def loop_halves_bridged_into_one():
    """A Vann-Nett loop off one point of the row (added as two halves) whose
    anchors the main stem's detour moves away: _bridge pairs the end of one
    half with the start of the other (same Vann-Nett part, 0.5 m apart) and
    joins them along that part -- the whole loop again, reversed."""
    r = LineString([(0, 0), (400, 0)])
    net = net_([("1-1", "001", 2, r)])
    m = LineString([(0, 0), (100, 0), (110, 30), (290, 30), (300, 0), (400, 0)])
    lp = LineString([(200, 0), (200, -20), (210, -20), (200.5, 0)])
    vn = vn_([m, lp])
    return net, vn, vc.complete(net, vn)


def spur_head_out_and_back():
    """(real: line 3547 at (-41960.19, 6805688.81), 44.7 m connector)
    Vann-Nett runs along the row to x=50 and then 12 m north to its end
    F=(52, 12); ELVIS turns east there (off Vann-Nett). The 'reach' target on
    the row itself replaces 1.8 m of it (accepted through the short-line
    ratio exception) and draws the line out to F and straight back."""
    r = LineString([(0, 0), (50, 0), (100, 0)])
    net = net_([("1-1", "001", 2, r)])
    v = LineString([(0, 0), (50, 0), (52, 12)])
    vn = vn_([v])
    return net, vn, vc.complete(net, vn)


def detour_leaves_tributary_at_old_confluence():
    """(real: 123-43-663_123010516, 064-14-182_064005704) The main stem m1
    takes its Vann-Nett meander through the confluence J=(140, -30) as a
    detour; ELVIS tributary t still ends at the old confluence N=(100, 0)
    while Vann-Nett's drawing of it (leaving t at (180, -80)) is added into J:
    the tributary is drawn twice."""
    m0 = LineString([(-100, 0), (100, 0)])
    m1 = LineString([(100, 0), (400, 0)])
    t = LineString([(300, -200), (100, 0)])
    net = net_(
        [("1-1", "001", 2, m0), ("1-1", "002", 2, m1), ("1-2", "001", 1, t)]
    )
    j = (140, -30)
    vm1 = LineString([(-100, 0), (100, 0), (115, 0), (125, -15), j])
    vm2 = LineString([j, (260, -30), (300, 0), (400, 0)])
    vt = LineString([(300, -200), (180, -80), (160, -50), j])
    vn = vn_([vm1, vm2, vt])
    return net, vn, vc.complete(net, vn)


def zero_width_tiny_line():
    """(real: 087-1-737_087002795, non-simple) Vann-Nett spikes 1.3 m off the
    row and back at x=50: the missing tip is a 0.6 m line with both anchors at
    one measure, under 2 m, so no loop: step 1 inserts p -> A -> tip -> B -> p,
    out and back over the same metre."""
    r = LineString([(0, 0), (100, 0)])
    net = net_([("1-1", "001", 2, r)])
    v = LineString([(0, 0), (50, 0), (50, 1.3), (50, 0), (100, 0)])
    vn = vn_([v])
    return net, vn, vc.complete(net, vn)


def fragment_extension_both_ends_on_one_row():
    """(real, harmless only because t0 = 0: chain [2624, 2623] over
    VN-152-75-R_VN0002096 at (414146.00, 7318616.07)) An earlier addition X
    covers the middle of a Vann-Nett meander beside ELVIS row e1, loose at both
    ends; the meander's two missing ends anchor on e1 at 103 m and 697 m.
    _over_additions joins line + X + line into one chain whose two ends lie on
    e1, and _path_core walks e1 back to its start and out again: e1's first
    103 m (on Vann-Nett) are cut off and its start left loose."""
    e1 = LineString([(0, 0), (800, 0)])
    v = LineString([(0, 0), (100, 0), (110, 30), (690, 30), (700, 0), (800, 0)])
    x = shapely.ops.substring(
        v, v.project(shapely.Point(300, 30)), v.project(shapely.Point(500, 30))
    )
    net = net_([("1-1", "001", 2, e1), ("VN-001-1-R", "VN0000001", 2, x)])
    net.loc[1, vp.PATCH] = vp.ADDED
    vn = vn_([v])
    return net, vn, vc.complete(net, vn)


CASES = {
    f.__name__: f
    for f in [
        two_additions_meet_on_a_redrawn_row,
        two_extra_rows_at_one_node,
        trib_through_node_tie,
        trib_through_node_tie_t_first,
        trib_through_node_offset,
        tributary_chain_first,
        tributary_chain_other_water_body,
        connector_crosses_an_addition,
        merged_detours_vn_in_two_parts,
        overlapping_detours_loose_anchor_unlisted,
        main_stem_crosses_elvis_tributary,
        main_stem_crosses_elvis_tributary_junction_beyond,
        loop_halves_bridged_into_one,
        spur_head_out_and_back,
        detour_leaves_tributary_at_old_confluence,
        zero_width_tiny_line,
        fragment_extension_both_ends_on_one_row,
    ]
}


# ── the tests ─────────────────────────────────────────────────────────────
def _rows(res):
    out = {}
    for _, r in res.network.iterrows():
        out.setdefault(f"{r.elvID}_{r.strekningLnr}", []).append(
            (
                shapely.get_parts(r.geometry)[0],
                r["patch"] if isinstance(r["patch"], str) else "",
            )
        )
    return out


def test_two_additions_meeting_on_a_redrawn_row_are_noded():
    """_attach_all skips an inner end when two additions end there
    (vannnett_complete.py:1766), though the detoured row runs through it."""
    net, vn, res = two_additions_meet_on_a_redrawn_row()
    new_t = [t for t in t_ends(res.network) if t not in t_ends(net)]
    assert not new_t, f"addition ends on the row's interior: {new_t}"
    assert components(res.network) == 1


def test_two_rows_at_one_moved_node_get_no_shared_connector():
    """Every extra row at an inner node gets its own connector N -> x
    (vannnett_complete.py:1561-1578): two rows drawn on one segment."""
    net, vn, res = two_extra_rows_at_one_node()
    assert not overlaps(res.network), overlaps(res.network)


@pytest.mark.parametrize(
    "name", ["trib_through_node_tie", "trib_through_node_tie_t_first"]
)
def test_a_tributary_drawn_through_the_node_takes_that_drawing(name):
    """Vann-Nett draws t's stream on through the ELVIS node into J; its
    missing piece anchors on m1/m2 (tie at the node), so _junction_line
    (e[1] != r_x, :1376) misses it: t gets a 30 m connector and the same
    line is added on top of it."""
    net, vn, res = CASES[name]()
    assert not overlaps(res.network), overlaps(res.network)
    assert not len(res.qa.get("connectors", [])), (
        "connector where Vann-Nett draws the junction"
    )


def test_a_tributary_leaving_just_past_the_node_is_no_detour_of_the_main_stem():
    """The tributary's piece anchors on m2 0.5 m past the node, so the chain
    [tributary, main stem below J] is a detour of m2 (step 1, before paths):
    m2 runs down the tributary and the main stem's own line is added beside m1."""
    net, vn, res = trib_through_node_offset()
    assert not len(res.qa.get("kept_beside", [])), (
        res.qa["kept_beside"].drop(columns="geometry").to_dict("records")
    )
    rows = _rows(res)
    assert all("vn_completed" in t for _, t in rows["1-1_001"]), (
        "the main stem above the node keeps ELVIS"
    )


def test_a_connector_crosses_no_added_line():
    """_crosses (:1067) checks the original rows and the chain only; the
    Vann-Nett tributary added in step 5 crosses t's 40 m connector."""
    net, vn, res = connector_crosses_an_addition()
    new_x = [c for c in crossings(res.network) if c not in crossings(net)]
    assert not new_x, new_x


def test_vannnett_main_stem_crossing_an_elvis_tributary_replaces_the_elvis_main_stem():
    """_crossings pairs the main stem's two halves where it crosses ELVIS
    tributary t, and both chains are dropped (:587-591): the main stem is
    added beside the ELVIS one it duplicates."""
    net, vn, res = main_stem_crosses_elvis_tributary()
    rows = _rows(res)
    assert all(
        "vn_completed" in t for _, t in rows["1-1_001"] + rows["1-1_002"]
    )
    assert not len(res.qa.get("kept_beside", []))


def test_a_head_draws_no_line_out_and_back():
    """'reach' on the anchor row itself replaces 1.8 m of ELVIS (ratio
    exception for lines under near_m, :1431-1435) and returns from F by a
    12 m connector along the line it drew (real: 32 heads, 14 spikes)."""
    net, vn, res = spur_head_out_and_back()
    assert not spikes(res.network), spikes(res.network)
    assert not len(res.qa.get("connectors", []))


def test_a_detour_through_a_confluence_takes_the_tributary_drawing():
    """Step 1 has no junction rule: ELVIS tributary t keeps the old
    confluence while Vann-Nett's drawing of it, anchored on t, is added into
    J (real: 123-43-663_123010516, 064-14-182_064005704, 012-1-17326 site)."""
    net, vn, res = detour_leaves_tributary_at_old_confluence()
    assert not len(res.qa.get("kept_beside", [])), "tributary drawn twice"
    rows = _rows(res)
    assert any("vn_completed" in t for _, t in rows["1-2_001"])


def test_loop_halves_are_not_bridged_into_one():
    """_bridge pairs one half's end with the other half's start on their
    shared Vann-Nett part and joins them along it (:1644-1656): the whole loop
    again, reversed -- a self-overlapping row with a 180 deg spike."""
    net, vn, res = loop_halves_bridged_into_one()
    assert not overlaps(res.network), overlaps(res.network)
    assert not spikes(res.network), spikes(res.network)


def test_merged_detours_draw_no_straight_line_longer_than_jog_m():
    """_merge_reps joins two replacements straight when no single
    Vann-Nett part carries both ends (:976-986): 2.24 m chord."""
    net, vn, res = merged_detours_vn_in_two_parts()
    assert not long_straights(res.network, vn, net)


def test_an_anchor_end_left_loose_is_listed():
    """_bridge turns an orphaned anchor end into a 'free end' when its
    Vann-Nett part ends within BRIDGE_M (:1693, :1705-1706); loose_ends stays
    empty although two anchor ends now dangle."""
    net, vn, res = overlapping_detours_loose_anchor_unlisted()
    new_d = sorted(set(dangling(res.network)) - set(dangling(net)))
    listed = res.qa.get("loose_ends")
    pts = (
        []
        if listed is None
        else [
            tuple(np.round(shapely.get_coordinates(g)[0], 3))
            for g in listed.geometry
        ]
    )
    unlisted = [
        d for d in new_d if all(Point(d).distance(Point(p)) > 2.0 for p in pts)
    ]
    assert not unlisted, unlisted


def test_an_extension_with_both_ends_on_one_row_keeps_that_row_whole():
    """_over_additions' fragment case (:1163-1190) joins line + earlier
    addition + line into a chain anchored twice on e1; _path_core (:1205)
    then walks e1 back to its start and out again and moves e1's start to the
    anchor: e1's first 100 m (on Vann-Nett) are cut off."""
    net, vn, res = fragment_extension_both_ends_on_one_row()
    rows = _rows(res)
    [(e1, _)] = rows["1-1_001"]
    assert e1.coords[0] == (0.0, 0.0), f"e1 now starts at {e1.coords[0]}"
    assert not len(res.qa.get("loose_ends", []))


def test_an_addition_never_reuses_the_key_of_a_dropped_row():
    """_assemble numbers additions from the max VN number of the pass input
    (:1973-1981); a row dropped in pass 1 is gone from pass 2's input, so pass
    2 hands its key to a different line."""
    net, vn, res = fragment_extension_both_ends_on_one_row()
    dropped = set(res.qa["dropped"].stretch)
    keys = set((res.network.elvID + "_" + res.network.strekningLnr).to_list())
    assert not (dropped & keys), dropped & keys


# ── round 4 (an independent review of v4) and v5 ───────────────────────────
def _two(parts):
    """A Vann-Nett layer of (vannforekomstid, lines) pairs."""
    return pd.concat(
        [vn_(lines, vid) for vid, lines in parts], ignore_index=True
    )


def crossing_chain_detour_crosses_row():
    """ELVIS main stem r and a stream t ending just short of it; Vann-Nett
    draws the main stem 40 m north, across t: the chain run on across t would
    be a detour of r crossing t without a node."""
    r = LineString([(0, 0), (600, 0)])
    t = LineString([(250, 200), (250, 10)])
    net = net_([("1-1", "001", 2, r), ("1-2", "001", 1, t)])
    v = LineString([(0, 0), (100, 0), (110, 40), (390, 40), (400, 0), (600, 0)])
    vn = vn_([v])
    return net, vn, vc.complete(net, vn)


def crossing_at_a_row_end_is_a_junction():
    """(real: 123-43-524_123009394 at (354637.37, 6986465.12)) An earlier
    addition t ends on Vann-Nett's node P=(200, 40), where Vann-Nett's main
    stem, two parts meeting at P, bends: the lines on either side of P meet t
    at its end -- a junction, not a crossing."""
    r = LineString([(0, 0), (400, 0)])
    t = LineString([(200, 200), (200, 40)])
    net = net_([("1-1", "001", 2, r), ("VN-001-2-R", "VN0000005", 1, t)])
    net.loc[1, vp.PATCH] = vp.ADDED
    p1 = LineString([(0, 0), (100, 0), (110, 38), (190, 38), (200, 40)])
    p2 = LineString([(200, 40), (210, 38), (290, 38), (300, 0), (400, 0)])
    vt = LineString([(200, 200), (200, 40)])
    vn = vn_([p1, p2, vt])
    return net, vn, vc.complete(net, vn)


def any_way_detours_share_a_line():
    """Two ELVIS channels r and s; Vann-Nett's central line c joins both
    through junctions J1, J2 either way round: one row may take it, not
    both."""
    r = LineString([(0, 0), (400, 0)])
    s = LineString([(0, 100), (400, 100)])
    net = net_([("1-1", "001", 2, r), ("2-1", "001", 2, s)])
    J1, J2 = (100, 50), (300, 50)
    vn = vn_(
        [
            LineString([J1, (80, 0), (0, 0)]),
            LineString([(0, 100), (80, 100), J1]),
            LineString([J1, J2]),
            LineString([J2, (320, 0), (400, 0)]),
            LineString([(400, 100), (320, 100), J2]),
        ]
    )
    return net, vn, vc.complete(net, vn)


def stays_anchor_tie():
    """(real: 026-22-1718_026012517 at (34689.97, 6559391.21)) Vann-Nett draws
    the main stem 13 m east of node N and tributary t's stream on from N; that
    line's anchor names m_out (0.98 m), not t (1.00 m)."""
    t = LineString([(27, -90), (4, -40), (0, 0)])
    m_in = LineString([(60, -23), (25, -16), (16, -1), (8, -5), (0, 0)])
    m_out = LineString([(0, 0), (11, 7), (-40, 89)])
    net = net_(
        [
            ("1-2", "001", 1, t),
            ("1-1", "001", 2, m_in),
            ("1-1", "002", 2, m_out),
        ]
    )
    vm = LineString([(60, -23), (25, -16), (-40, 89)])
    vk = LineString([(27, -90), (4, -40), (0, 0), (-40, 89)])
    vn = _two([("001-2-R", [vm]), ("001-3-R", [vk])])
    return net, vn, vc.complete(net, vn)


def free_head_deletes_elvis_course():
    """A 5 m Vann-Nett line from free end F=(0, 5) joins row r, whose own free
    start lies 75 m upstream: r keeps that course."""
    r = LineString([(-75, 0), (0, 0), (200, 0)])
    net = net_([("1-1", "001", 2, r)])
    vn = vn_([LineString([(0, 5), (5, 1), (6, 0), (200, 0)])])
    return net, vn, vc.complete(net, vn)


def tiny_line_covered_ring():
    """(real: VN-087-223-R_VN0004908, a 3.5 mm ring) Two sub-metre spikes off
    one point of r; the first is a detour whose body carries the second."""
    r = LineString([(0, 0), (100, 0)])
    net = net_([("1-1", "001", 2, r)])
    v1 = LineString(
        [(0, 0), (50, 0), (50, 1.03), (50.02, 1.03), (50.02, 0), (100, 0)]
    )
    v2 = LineString([(50.009, 0), (50.0095, 1.0017), (50.010, 0)])
    vn = _two([("001-1-R", [v1]), ("001-2-R", [v2])])
    return net, vn, vc.complete(net, vn)


def kept_beside_reason_of_another_line():
    """r gets a detour (j) near its start; line i runs 20 m south of r and is
    added (no network path): r's stretch beside i carries i's reason."""
    u = LineString([(-30, 100), (0, 0)])
    r = LineString([(0, 0), (300, 0)])
    d = LineString([(300, 0), (400, 0)])
    q = LineString([(45, -50), (85, 10)])
    net = net_(
        [
            ("1-1", "001", 2, u),
            ("1-1", "002", 2, r),
            ("1-1", "003", 2, d),
            ("2-1", "001", 1, q),
        ]
    )
    j = LineString([(5, 0), (8, 20), (57, 20), (60, 0)])
    i = LineString([(65, -20), (290, -20), (300, 0)])
    vn = vn_([j, i, LineString([(300, 0), (400, 0)])])
    return net, vn, vc.complete(net, vn)


def kept_beside_part_of_a_long_row():
    """As above, but r runs on 700 m past the line, off Vann-Nett: the part
    beside the added line is listed although it is a fifth of r's off run."""
    u = LineString([(-30, 100), (0, 0)])
    r = LineString([(0, 0), (1000, 0)])
    q = LineString([(45, -50), (85, 10)])
    net = net_(
        [("1-1", "001", 2, u), ("1-1", "002", 2, r), ("2-1", "001", 1, q)]
    )
    i = LineString([(65, -20), (290, -20), (300, 0)])
    vn = vn_([i])
    return net, vn, vc.complete(net, vn)


def shared_connector_node_off_the_mm_grid():
    """Two extra rows at a moved node N lying off the millimetre grid share
    one connector; the row taking it is cut at N."""
    N = (150.0004, 0.0003)
    m1 = LineString([(0, 0), N])
    m2 = LineString([N, (400, 0)])
    t1 = LineString([(100, 200), N])
    t2 = LineString([(200, 200), N])
    net = net_(
        [
            ("1-1", "001", 2, m1),
            ("1-1", "002", 2, m2),
            ("1-2", "001", 1, t1),
            ("1-3", "001", 1, t2),
        ]
    )
    vm = LineString(
        [(0, 0), (100, 0), (110, -30), (290, -30), (300, 0), (400, 0)]
    )
    vn = vn_([vm])
    return net, vn, vc.complete(net, vn)


def addition_ending_on_a_row_interior():
    """(real: t13, 123-386-R beside 123-43-3683_123003609) ELVIS tributary s
    joins the main stem m at N; Vann-Nett draws the tributary 40 m east. An
    earlier addition A carries its lower part and ends on m's interior at
    (190, 0) -- a T without a node; the upper part is missing, from where
    Vann-Nett leaves ELVIS (U) to A's loose start. s is the same stream."""
    u = LineString([(110, 600), (110, 400)])
    s = LineString([(110, 400), (150, 0)])
    m1 = LineString([(0, 0), (150, 0)])
    m2 = LineString([(150, 0), (400, 0)])
    a = LineString([(190, 100), (190, 0)])
    net = net_(
        [
            ("1-2", "001", 1, u),
            ("1-2", "002", 1, s),
            ("1-1", "001", 2, m1),
            ("1-1", "002", 2, m2),
            ("VN-001-2-R", "VN0000007", 1, a),
        ]
    )
    net.loc[4, vp.PATCH] = vp.ADDED
    vt = LineString([(110, 600), (110, 400), (190, 300), (190, 0)])
    vm = LineString([(0, 0), (400, 0)])
    vn = _two([("001-2-R", [vt]), ("001-1-R", [vm])])
    return net, vn, vc.complete(net, vn)


def head_reach_runs_back():
    """ELVIS r runs 10 m beside a Vann-Nett head line and turns north at
    x=120, 30 m short of Vann-Nett's free end F=(150, -10): a connector to
    the point of r nearest F would run back along the line."""
    r = LineString([(120, 200), (120, 0), (50, 0), (0, 0)])
    net = net_([("1-1", "001", 2, r)])
    v = LineString([(150, -10), (55, -10), (50, 0), (0, 0)])
    vn = vn_([v])
    return net, vn, vc.complete(net, vn)


def small_ring_at_one_point():
    """(real: VN-064-425-R_VN0004922, a 3.1 m ring) Vann-Nett's V tip leaves r
    at one point and comes back to it, 2.4 m out and under 1 m wide: 3.2 m
    of it lies more than on_m from r (under 2 m it is a sub-metre detour)."""
    r = LineString([(0, 0), (100, 0)])
    net = net_([("1-1", "001", 2, r)])
    v = LineString(
        [(0, 0), (49.6, 0), (49.8, -2.4), (50.2, -2.4), (50.4, 0), (100, 0)]
    )
    vn = vn_([v])
    return net, vn, vc.complete(net, vn)


def test_any_way_detours_never_take_one_line_twice():
    net, vn, res = any_way_detours_share_a_line()
    assert not overlaps(res.network), overlaps(res.network)


def test_a_chain_run_across_a_row_leaves_no_unnoded_crossing():
    net, vn, res = crossing_chain_detour_crosses_row()
    new_x = [c for c in crossings(res.network) if c not in crossings(net)]
    assert not new_x, new_x


def test_a_row_ending_on_vannnett_s_junction_is_joined_there():
    net, vn, res = crossing_at_a_row_end_is_a_junction()
    assert components(res.network) == 1
    ends = {
        tuple(np.round(c, 6))
        for k, v in _rows(res).items()
        if k != "VN-001-2-R_VN0000005"
        for g, _ in v
        for c in shapely.get_coordinates(g)[[0, -1]]
    }
    assert (200.0, 40.0) in ends, "no row ends at t's end, Vann-Nett's junction"


def test_a_row_whose_stream_vannnett_carries_on_gets_no_connector():
    net, vn, res = stays_anchor_tie()
    assert not len(res.qa.get("connectors", [])), (
        res.qa["connectors"].drop(columns="geometry").to_dict("records")
    )
    [(t, tag)] = _rows(res)["1-2_001"]
    assert tag == "" and t.coords[-1] == (0.0, 0.0)


def test_a_free_elvis_end_far_beyond_a_short_line_keeps_its_course():
    net, vn, res = free_head_deletes_elvis_course()
    [(r, _)] = _rows(res)["1-1_001"]
    assert r.coords[0] == (-75.0, 0.0), f"r now starts at {r.coords[0]}"


def test_a_line_another_replacement_already_covers_is_not_added():
    net, vn, res = tiny_line_covered_ring()
    added = [
        g for k, v in _rows(res).items() if k.startswith("VN-") for g, _ in v
    ]
    assert all(g.is_simple and g.coords[0] != g.coords[-1] for g in added), [
        list(g.coords) for g in added
    ]
    assert not len(res.qa.get("flow_conflicts", []))
    mb = res.qa["missing_before"]
    assert (mb.why == "carried within on_m by a line this pass put in").any()


def test_kept_beside_gives_the_reason_of_the_line_beside_it():
    net, vn, res = kept_beside_reason_of_another_line()
    kb = res.qa["kept_beside"]
    whys = set(kb[kb.stretch == "1-1_002"].why)
    assert whys == {"no network path between the anchors"}, whys


def test_kept_beside_lists_the_part_of_a_long_row_beside_a_line():
    net, vn, res = kept_beside_part_of_a_long_row()
    kb = res.qa["kept_beside"]
    kb = kb[kb.stretch == "1-1_002"]
    assert len(kb) == 1, kb.drop(columns="geometry").to_dict("records")
    g = kb.geometry.iloc[0]
    assert 180 < g.length < 260 and g.bounds[2] < 320, (g.length, g.bounds)


def test_vn_path_joins_two_parts_at_the_nearer_node():
    j = LineString([(0, 0), (10, 0), (20, 0)])
    i = LineString([(5, 0), (5, 5), (15, 5), (15, 0)])
    vg = np.array([i, j], dtype=object)
    mid = vc._vn_path(
        vg, shapely.STRtree(vg), np.array([15.0, 0.5]), np.array([15.5, 0.0])
    )
    assert mid is not None and LineString(mid).length < 1.5


def test_rows_sharing_a_connector_stay_joined_at_a_node_off_the_mm_grid():
    net, vn, res = shared_connector_node_off_the_mm_grid()
    assert components(res.network) == 1


def test_a_tributary_drawn_through_the_node_leaves_no_duplicate():
    net, vn, res = trib_through_node_tie()
    assert not len(res.qa.get("kept_beside", [])), (
        res.qa["kept_beside"].drop(columns="geometry").to_dict("records")
    )


@pytest.mark.xfail(
    reason="the main-stem plan is kept when t's connector would cross the "
    "Vann-Nett tributary, and the head rule reroutes the main stem (round "
    "4, synthetic only: no such case nationally)",
    strict=True,
)
def test_a_connector_crossing_a_missing_line_does_not_reroute_the_main_stem():
    net, vn, res = connector_crosses_an_addition()
    conns = res.qa.get("connectors")
    assert conns is None or not (conns.kind == "head").any()
    assert not len(res.qa.get("kept_beside", []))


@pytest.mark.xfail(
    reason="the tributary chain goes first where Vann-Nett's junction lies "
    "beyond the crossing (round 3 finding 6, second variant; synthetic only)",
    strict=True,
)
def test_a_main_stem_crossing_a_tributary_before_the_junction_is_replaced():
    net, vn, res = main_stem_crosses_elvis_tributary_junction_beyond()
    assert not len(res.qa.get("kept_beside", []))


def test_an_addition_ending_on_a_row_interior_is_taken_with_the_line():
    """v4 gave up where an earlier addition ends on a row's interior (an
    unnoded T): the line was added beside s. The addition run now meets
    ELVIS there, on m2, and s takes the whole."""
    net, vn, res = addition_ending_on_a_row_interior()
    rows = _rows(res)
    assert not len(res.qa.get("kept_beside", [])), (
        res.qa["kept_beside"].drop(columns="geometry").to_dict("records")
    )
    assert "VN-001-2-R_VN0000007" in set(res.qa["dropped"].stretch)
    assert all("vn_completed" in t for _, t in rows["1-2_002"])
    assert not [t for t in t_ends(res.network) if t not in t_ends(net)]
    assert components(res.network) == 1


def test_a_head_connector_that_would_run_back_joins_elvis_further_on():
    """The connector at r's point nearest F would run back along the line;
    it joins r a little on instead, and r takes the line (no duplicate)."""
    net, vn, res = head_reach_runs_back()
    assert not len(res.qa.get("kept_beside", [])), (
        res.qa["kept_beside"].drop(columns="geometry").to_dict("records")
    )
    assert not spikes(res.network, max_turn=150.0), spikes(res.network)
    [(r, _)] = _rows(res)["1-1_001"]
    assert r.coords[0] == (120.0, 200.0)
    conns = res.qa["connectors"]
    assert len(conns) == 1 and conns.length_m.iloc[0] <= vc.CONNECT_M


def test_a_small_ring_at_one_point_is_drawn_once():
    """A Vann-Nett loop off one point that one half carries within on_m is
    one whisker, flowing back in, its tip listed -- not a ring or a 2-cycle."""
    net, vn, res = small_ring_at_one_point()
    added = [
        g for k, v in _rows(res).items() if k.startswith("VN-") for g, _ in v
    ]
    assert len(added) == 1 and added[0].coords[0] != added[0].coords[-1]
    assert not len(res.qa.get("flow_conflicts", []))
    assert "tip" in set(res.qa["loose_ends"].kind)


def test_new_sinks_and_sources_are_found_against_the_input():
    """_flow_ends: nodes where two or more rows meet with only inflow (or only
    outflow); complete() lists those the input does not have."""
    a = LineString([(0, 100), (100, 0)])
    b = LineString([(200, 100), (100, 0)])
    c = LineString([(100, 0), (100, -100)])
    before = vc._flow_ends(
        net_([("1-1", "001", 2, a), ("1-2", "001", 1, b), ("1-1", "002", 2, c)])
    )
    assert not before
    after = vc._flow_ends(net_([("1-1", "001", 2, a), ("1-2", "001", 1, b)]))
    assert [v[0] for v in after.values()] == ["sink"]


# ── round 5 (v6) ───────────────────────────────────────────────────────────
def directed_cycles(net):
    """Rows lying on a directed cycle (row start -> row end, nodes to the
    millimetre), rings included."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    g = _lines(net)
    a = shapely.get_coordinates(shapely.get_point(g, 0))
    b = shapely.get_coordinates(shapely.get_point(g, -1))
    q = np.round(np.vstack([a, b]) * 1000).astype(np.int64)
    uniq, inv = np.unique(q, axis=0, return_inverse=True)
    inv = inv.ravel()
    n = len(g)
    u, v = inv[:n], inv[n:]
    m = csr_matrix((np.ones(n), (u, v)), shape=(len(uniq), len(uniq)))
    lab = connected_components(m, directed=True, connection="strong")[1]
    return [k for k in range(n) if lab[u[k]] == lab[v[k]]]


def stays_line_covered_leaves_row_dangling():
    """(real: 064-14-690_064005751) ELVIS main stem a -> N -> b -> c, and x
    flowing west into N. Vann-Nett draws the main stem 1.6 m west of N and
    x's stream on past N into a junction of that line: the missing piece of
    it is 0.6 m. The duplicate path moves N onto the line; x keeps its end
    at N because that piece carries its stream on -- and the piece lies
    within on_m of the redrawn main stem."""
    a = LineString([(0, 40), (0, 0)])
    b = LineString([(0, 0), (0, -3.5)])
    c = LineString([(0, -3.5), (0, -40)])
    x = LineString([(40, 0), (0, 0)])
    net = net_(
        [
            ("1-1", "001", 2, a),
            ("1-1", "002", 2, b),
            ("1-1", "003", 2, c),
            ("1-2", "001", 1, x),
        ]
    )
    vm = LineString([(0, 40), (0, 6), (-1.6, 3), (-1.6, 0)])
    vm2 = LineString([(-1.6, 0), (-1.6, -3), (0, -6), (0, -40)])
    vt = LineString([(40, 0), (0, 0), (-1.6, 0)])
    vn = vn_([vm, vm2, vt])
    return net, vn, vc.complete(net, vn)


def covered_line_at_a_junction_of_additions():
    """Vann-Nett: a main line on ELVIS rows r1 -> r2 and an arch beside it
    made of a -> J, a 1.6 m link d from J to K, and b from K on, with
    tributaries e and f ending at J and K: d lies within on_m of the other
    lines' ends but no line runs along it."""
    r1 = LineString([(-100, 0), (0, 0)])
    r2 = LineString([(0, 0), (100, 0)])
    net = net_([("1-1", "001", 2, r1), ("1-1", "002", 2, r2)])
    J, K = (0.0, 1.5), (1.6, 1.5)
    vn = vn_(
        [
            LineString([(-100, 0), (100, 0)]),
            LineString([(-30, 0), (-30, 20), J]),
            LineString([J, K]),
            LineString([K, (30, 20), (30, 0)]),
            LineString([(-5, 40), J]),
            LineString([(6, 40), K]),
        ]
    )
    return net, vn, vc.complete(net, vn)


def covered_by_the_undrawn_half_of_a_whisker():
    """A narrow Vann-Nett loop off r is drawn once (its back half); a
    separate Vann-Nett line X lies within 0.9 m of the out half only."""
    r = LineString([(-50, 0), (50, 0)])
    net = net_([("1-1", "001", 2, r)])
    loop = LineString([(0, 1.0), (0, 20), (0.5, 20), (0.5, 1.0)])
    x = LineString([(-0.8, 5), (-0.8, 15)])
    vn = vn_([LineString([(-50, 0), (50, 0)]), loop, x])
    return net, vn, vc.complete(net, vn)


def follow_on_through_a_t():
    """(real: t13, 123-386-R) Main stem m runs from the north-east through A
    to N0, where tributary s ends and m2 goes on. Vann-Nett's tributary
    joins at A instead, and its main stem bends 5 m south between A and N0.
    Earlier additions carry the tributary's last 60 m (ending on m's
    interior at A -- a T without a node) and the bend (A to N0); the rest of
    the tributary is missing."""
    m = LineString([(400, 100), (140, 0), (100, 0)])
    m2 = LineString([(100, 0), (0, 0)])
    s = LineString([(100, -300), (100, 0)])
    a1 = LineString([(140, -60), (140, 0)])
    a2 = LineString([(140, 0), (120, -5), (100, 0)])
    net = net_(
        [
            ("1-1", "001", 2, m),
            ("1-1", "002", 2, m2),
            ("1-2", "001", 1, s),
            ("VN-001-2-R", "VN0000007", 1, a1),
            ("VN-001-1-R", "VN0000008", 2, a2),
        ]
    )
    net.loc[[3, 4], vp.PATCH] = vp.ADDED
    vt = LineString([(100, -300), (100, -250), (140, -200), (140, 0)])
    vm = LineString([(400, 100), (140, 0), (120, -5), (100, 0), (0, 0)])
    vn = _two([("001-2-R", [vt]), ("001-1-R", [vm])])
    return net, vn, vc.complete(net, vn)


def detour_over_a_t_another_addition_ends_at():
    """(real: VN-123-182-R_VN0000102) ELVIS r runs straight where Vann-Nett
    bends 30 m north. An earlier addition A carries the bend's last part
    and ends on r's interior at P -- a T without a node -- where another
    earlier addition B (a tributary) ends too."""
    r = LineString([(0, 0), (300, 0)])
    a = LineString([(200, 30), (250, 0)])
    b = LineString([(250, 80), (250, 0)])
    net = net_(
        [
            ("1-1", "001", 2, r),
            ("VN-001-1-R", "VN0000007", 2, a),
            ("VN-001-2-R", "VN0000008", 1, b),
        ]
    )
    net.loc[[1, 2], vp.PATCH] = vp.ADDED
    vm = LineString([(0, 0), (50, 0), (100, 30), (200, 30), (250, 0), (300, 0)])
    vt = LineString([(250, 80), (250, 0)])
    vn = _two([("001-1-R", [vm]), ("001-2-R", [vt])])
    return net, vn, vc.complete(net, vn)


def test_a_staying_row_is_joined_by_the_line_that_carries_it_on():
    """The line the stays rule relies on is added even where the redrawn
    rows carry it within on_m: x stays joined."""
    net, vn, res = stays_line_covered_leaves_row_dangling()
    assert not [d for d in dangling(res.network) if d not in dangling(net)]
    assert components(res.network) == 1


def test_a_line_other_lines_meet_at_a_junction_is_never_skipped():
    """A line others meet at a junction is added even where other lines'
    ends lie within on_m of it: the arch stays whole, nothing is drawn twice
    and the link is in the network."""
    net, vn, res = covered_line_at_a_junction_of_additions()
    assert not overlaps(res.network), overlaps(res.network)
    assert res.stats.get("covered", 0) == 0
    g = shapely.union_all(_lines(res.network))
    link = LineString([(0.0, 1.5), (1.6, 1.5)])
    assert link.difference(g.buffer(1e-6)).length < 1e-6


def test_only_the_drawn_half_of_a_loop_carries_other_lines():
    """X lies within 0.9 m of the loop's out half, which is not drawn: X is
    added in the same pass and nothing is left missing."""
    net, vn, res = covered_by_the_undrawn_half_of_a_whisker()
    assert res.stats["passes"] == 1, res.stats
    assert res.stats["missing_km_after"] == 0


def test_a_pass_whose_lines_have_no_anchor_does_not_crash():
    """An isolated Vann-Nett stream (both ends free) is the only missing
    line: it is added."""
    r = LineString([(0, 0), (100, 0)])
    net = net_([("1-1", "001", 2, r)])
    vn = vn_([LineString([(0, 0), (100, 0)]), LineString([(0, 50), (60, 50)])])
    res = vc.complete(net, vn)
    assert res.stats["missing_km_after"] == 0
    assert res.stats["added"] == 1


def test_vn_path_finds_the_way_through_a_node_either_way_round():
    """Two parts ending on a third part's interior: a -> b and b -> a both
    find Vann-Nett's way through that node."""
    i = LineString([(0, 0), (200, 0)])
    j = LineString([(100, 80), (100, 0)])
    k = LineString([(150, 60), (100, 0)])
    V = vc.vn_lines(vn_([i, j, k]))
    vg = V.geometry.values
    tree = shapely.STRtree(vg)
    a, b = np.array([60.0, 0.0]), np.array([100.0, 40.0])
    assert vc._vn_path(vg, tree, a, b) is not None
    assert vc._vn_path(vg, tree, b, a) is not None


def test_a_line_taking_additions_through_a_t_meets_that_row_at_a_node():
    """s takes the tributary, the earlier additions and the bend to N0
    (the additions go); m is cut at A, where s runs through it, so the two
    meet at a node instead of crossing."""
    net, vn, res = follow_on_through_a_t()
    new_x = [c for c in crossings(res.network) if c not in crossings(net)]
    assert not new_x, new_x
    assert not [t for t in t_ends(res.network) if t not in t_ends(net)]
    assert components(res.network) == 1
    assert res.stats["missing_km_after"] == 0


def test_a_t_where_other_additions_end_is_made_a_node():
    """r takes the bend as a detour through P and the addition A goes; B
    still ends at P, so r is cut there and B joins it at a node."""
    net, vn, res = detour_over_a_t_another_addition_ends_at()
    assert "VN-001-1-R_VN0000007" in set(res.qa["dropped"].stretch)
    assert (250.0, 0.0) not in t_ends(res.network)
    assert components(res.network) == 1


def test_crossings_near_an_end_become_nodes():
    """_node_crossings: an addition crossing a row 0.1 m before both end at
    one point drops that stretch and the row is cut where they cross; two
    additions crossing 1 m from an end are both cut there; a crossing far
    from every end is left."""
    row = LineString([(0, 0), (10, 0)])
    final = [row]
    splits = defaultdict(list)
    adds = [[np.array([(5.0, -5.0), (9.95, 0.05), (10.0, 0.0)]), 0]]
    out = vc._node_crossings(adds, final, splits)
    assert len(out) == 1 and abs(out[0][0][-1][1]) < 1e-9
    assert (
        len(splits[0]) == 1 and abs(splits[0][0][0] - out[0][0][-1][0]) < 1e-9
    )
    splits = defaultdict(list)
    adds = [
        [np.array([(0.0, 50.0), (10.0, 50.0)]), 0],
        [np.array([(9.0, 49.0), (9.0, 60.0)]), 1],
    ]
    out = vc._node_crossings(adds, final, splits)
    assert len(out) == 4
    assert (
        sum(
            np.allclose(c[0], (9, 50)) or np.allclose(c[-1], (9, 50))
            for c, _ in out
        )
        == 4
    )
    adds = [
        [np.array([(0.0, 50.0), (100.0, 50.0)]), 0],
        [np.array([(50.0, 0.0), (50.0, 100.0)]), 1],
    ]
    assert len(vc._node_crossings(adds, final, defaultdict(list))) == 2


def test_a_cycle_through_an_addition_is_broken_at_a_head():
    """An addition closing a U of rows into a directed cycle is turned so
    that the U's head feeds both (a listed source) where nothing else
    breaks the cycle; a ring of two additions off a sink is turned so that
    both flow out to their tip."""
    b = LineString([(10, 10), (10, 0)])  # from the head T
    d = LineString([(10, 0), (0, 0), (0, 10)])  # round to A
    x = LineString([(0, 10), (-50, 10)])  # out of A
    adds = [[np.array([(0.0, 10.0), (10.0, 10.0)]), 0]]  # A -> T
    left = vc._cycles(adds, [b, d, x], {})
    assert not [x for x in left if x[2] == "cycle"]
    assert np.allclose(adds[0][0][0], (10, 10))
    # the head it leaves a source is listed under the turned addition
    assert [(k, tuple(xy), kind) for k, xy, kind in left] == [
        (0, (10.0, 10.0), "source")
    ]
    sink = LineString([(-10, 0), (0, 0)])
    adds = [
        [np.array([(0.0, 0.0), (0.4, 1.0)]), 0],
        [np.array([(0.4, 1.0), (0.0, 0.0)]), 1],
    ]
    left = vc._cycles(adds, [sink], {})
    assert not [x for x in left if x[2] == "cycle"]
    assert all(np.allclose(c[0], (0, 0)) for c, _ in adds)
    assert [kind for _k, _xy, kind in left] == ["sink"]  # the ring's tip


def test_a_row_end_left_meeting_nothing_is_found():
    """_lone_ends: a row end that met another in the input and meets none in
    the output."""
    a = LineString([(0, 0), (10, 0)])
    b = LineString([(10, 0), (20, 0)])
    b2 = LineString([(11, 0), (20, 0)])
    net_in = net_([("1-1", "001", 2, a), ("1-1", "002", 2, b)])
    net_out = net_([("1-1", "001", 2, a), ("1-1", "002", 2, b2)])
    lone = vc._lone_ends(net_in, net_out)
    assert [v[0] for v in lone.values()] == ["loose"]
    assert np.allclose(list(lone.values())[0][2], (10, 0))


# ── v7: reach walks, hooks, run-on, lenses, rings, crossings, the QA ───────
def _pieces_of(res, key):
    net = res.network
    m = (net.elvID.astype(str) + "_" + net.strekningLnr.astype(str)) == key
    return list(shapely.get_parts(net.geometry.values[m.to_numpy()]))


def test_a_reach_walk_that_is_not_the_line_s_stream_keeps_the_main_stem():
    """Vann-Nett's 8.4 m head joins the ELVIS main stem at A, 28.7 m below
    the node N where the stem's upper row runs in (a through node). Walking
    from A up to N is no drawing of the head (3.4 times its length, only
    'alike' within near_m): no reach there -- the stem keeps N -> A and the
    head is added (v6 rerouted the stem N -> F -> A with a 30.6 m connector)."""
    up = LineString([(78, -124), (0, 0)])
    down = LineString([(0, 0), (-10.8, 26.6), (-300, 26.6)])
    net = net_([("1-1", "001", 3, up), ("1-1", "002", 3, down)])
    v = LineString([(-3.3, 30.4), (-7.0, 28.5), (-10.8, 26.6), (-300, 26.6)])
    res = vc.complete(net, vn_([v]))
    assert not len(res.qa.get("connectors", []))
    stem = shapely.ops.substring(down, 0.0, 28.7)
    have = shapely.union_all(_lines(res.network))
    assert stem.difference(have.buffer(0.01)).length < 1e-6
    assert res.stats["added"] == 1


def test_a_reach_into_a_row_may_take_an_elvis_glitch_beside_the_line():
    """ELVIS's row runs in from the north-east, out to a V tip and back to
    the junction J where Vann-Nett's 11 m head comes in: the walk from J up
    the row to its point nearest the head's free end is 2.7 times the line
    (the V), but it ends inside the row, at no node ELVIS runs on through:
    ELVIS drawing the same stream is enough, and the V goes."""
    r = LineString([(60, 30), (-15, 10), (0, 0), (0, -60)])
    net = net_([("1-1", "001", 2, r)])
    vn = vn_([LineString([(-4, 10), (0, 0), (0, -60)])])
    res = vc.complete(net, vn)
    assert res.stats["added"] == 0
    have = shapely.union_all(_lines(res.network))
    assert Point(-15, 10).distance(have) > 3.0  # the V tip is gone


def test_a_free_head_whose_reach_would_hook_moves_its_end_after_all():
    """ELVIS's head runs west from the junction J, Vann-Nett's 11 m head
    south-west to its free end F: ELVIS runs on 6.6 m past the point nearest
    F (more than half the line), but a connector from F to ELVIS would turn
    back off the line by more than a right angle. The free end moves to F,
    no hook drawn."""
    main = LineString([(0, 0), (300, 220)])
    h = LineString([(-17.7, 0), (0, 0)])
    net = net_([("1-1", "001", 3, main), ("1-2", "001", 1, h)])
    res = vc.complete(
        net,
        vn_(
            [
                LineString([(0, 0), (300, 220)]),
                LineString([(-8.4, -7.2), (0, 0)]),
            ]
        ),
    )
    assert not len(res.qa.get("connectors", []))
    (g,) = _pieces_of(res, "1-2_001")
    assert np.allclose(shapely.get_coordinates(g)[0], (-8.4, -7.2))


def test_a_meandering_head_ending_at_the_free_end_gets_no_connector():
    """ELVIS's head meanders beside Vann-Nett's and ends 6 m beside its free
    end F, running on past it not at all: the walk is much longer than the
    line, but that is meander -- the free end moves to F."""
    m = LineString([(0, 0), (0, 100)])
    h = LineString(
        [
            (0, 100),
            (5, 104),
            (10, 108),
            (3, 113),
            (10, 118),
            (3, 123),
            (10, 128),
        ]
        + [(3, 133), (10, 137), (6, 140)]
    )
    t = LineString([(-60, 100), (0, 100)])
    net = net_(
        [("1-1", "001", 3, m), ("1-1", "002", 2, h), ("1-2", "001", 1, t)]
    )
    vn = vn_(
        [
            LineString([(0, 0), (0, 100), (0, 140)]),
            LineString([(-60, 100), (0, 100)]),
        ]
    )
    res = vc.complete(net, vn)
    assert not len(res.qa.get("connectors", []))
    assert not long_straights(res.network, vn, net)
    (g,) = _pieces_of(res, "1-1_002")
    assert np.allclose(shapely.get_coordinates(g)[-1], (0, 140))


def test_a_head_no_longer_than_its_line_moves_its_free_end_though_elvis_runs_on():
    """ELVIS's head runs on straight 30 m past Vann-Nett's free end F, but
    the walk is no longer than Vann-Nett's meandering line: the free end
    moves to F and the run-on goes. The cap is on the walk's length -- a cap
    on the run-on alone turned clean heads into lines added beside ELVIS."""
    m = LineString([(0, 0), (0, 100)])
    h = LineString([(4, 170), (0, 100)])
    net = net_([("1-1", "001", 3, m), ("1-2", "001", 1, h)])
    vh = LineString(
        [(9, 140), (12, 132), (6, 124), (12, 116), (6, 108), (0, 100)]
    )
    res = vc.complete(net, vn_([LineString([(0, 0), (0, 100)]), vh]))
    (g,) = _pieces_of(res, "1-2_001")
    assert np.allclose(shapely.get_coordinates(g)[0], (9, 140))
    assert not len(res.qa.get("connectors", [])) and res.stats["added"] == 0


def test_a_link_vann_nett_draws_twice_is_added_once():
    """Vann-Nett runs two lines on from a sink of the rows to one point
    beyond it, 35 cm apart at most (a narrow lens): one addition draws
    both; its tip, no free end of Vann-Nett's, is listed."""
    a = LineString([(0, 100), (0, 0)])
    b = LineString([(30, 100), (0, 0)])
    net = net_([("1-1", "001", 2, a), ("1-2", "001", 2, b)])
    vn = vn_(
        [
            LineString([(0, 100), (0, 0), (-0.2, -1.5)]),
            LineString([(30, 100), (0, 0), (0.15, -1.5)]),
            LineString([(-0.2, -1.5), (0.0, -2.6)]),
            LineString([(0.15, -1.5), (0.0, -2.6)]),
        ]
    )
    res = vc.complete(net, vn)
    assert res.stats["added"] == 1 and res.stats["missing_km_after"] == 0.0
    assert list(res.qa["loose_ends"].kind) == ["tip"]
    assert not len(res.qa.get("flow_conflicts", []))


def test_a_side_channel_round_a_node_is_joined_at_both_ends():
    """A 20 m side channel leaves r1 0.75 m before the node N and rejoins r2
    0.75 m after it: both ends snap to N, so the second lands on r2 itself
    instead -- no ring, nothing loose."""
    r1 = LineString([(-100, 0), (0, 0)])
    r2 = LineString([(0, 0), (100, 0)])
    net = net_([("1-1", "001", 2, r1), ("1-1", "002", 2, r2)])
    side = LineString([(-0.75, 0), (-0.75, 10), (0.75, 10), (0.75, 0)])
    res = vc.complete(net, vn_([LineString([(-100, 0), (100, 0)]), side]))
    adds = _pieces_of(res, "VN-001-1-R_VN0000001")
    assert adds and not any(g.is_closed for g in adds)
    assert not len(res.qa.get("loose_ends", []))
    assert not sorted(set(dangling(res.network)) - set(dangling(net)))


def test_a_source_left_by_a_cycle_turn_is_listed_under_the_addition():
    """An addition closing a U of rows into a cycle is turned at the U's
    head T: T is listed as a source by the pass, under the addition."""
    b = LineString([(10, 10), (10, 0)])
    d = LineString([(10, 0), (0, 0), (0, 10)])
    x = LineString([(0, 10), (-50, 10)])
    net = net_(
        [("1-1", "001", 1, b), ("1-1", "002", 1, d), ("1-1", "003", 1, x)]
    )
    vn = vn_(
        [
            LineString([(10, 10), (10, 0), (0, 0), (0, 10), (-50, 10)]),
            LineString([(0, 10), (10, 10)]),
        ]
    )
    fc = vc.complete(net, vn).qa["flow_conflicts"]
    at_t = fc[[Point(10, 10).distance(g) < 1e-6 for g in fc.geometry]]
    assert len(at_t) == 1
    r = at_t.iloc[0]
    assert str(r.stretch).startswith("VN-") and r.kind == "source"
    assert r.vannforekomstid == "001-1-R" and r.pass_no == 1


def test_a_row_end_left_alone_at_an_old_sink_is_listed():
    """r_prev and r both flowed into n (a sink of the input); r's end moves
    away: r_prev's end left meeting nothing at n is listed (kind loose)."""
    r_prev = LineString([(0, 0), (100, 0)])
    r = LineString([(100, 60), (100, 0)])
    net = net_([("1-1", "001", 2, r_prev), ("1-2", "001", 1, r)])
    v = LineString([(70, 0), (80, 20), (95, 45), (100, 52), (100, 60)])
    res = vc.complete(net, vn_([LineString([(0, 0), (100, 0)]), v]))
    fc = res.qa["flow_conflicts"]
    loose = [
        k
        for k, g in zip(fc.kind, fc.geometry, strict=True)
        if Point(100, 0).distance(g) < 1e-6
    ]
    assert loose == ["loose"]


def test_node_crossings_cut_no_line_where_only_a_dropped_stretch_crossed():
    """g overshoots the node E where h1 ends and h2 starts, crossing h1 and
    then h2, and comes back to E: it is cut back at the h1 crossing, so h2
    (crossed only by the stretch that goes) is not cut."""
    h1 = LineString([(0, 0), (10, 0)])
    h2 = LineString([(10, 0), (10, 10)])
    g = np.array([(5.0, -3.0), (10.5, 0.5), (10.0, 0.0)])
    splits = defaultdict(list)
    out = vc._node_crossings([[g, 0]], [h1, h2], splits)
    ends = [shapely.get_coordinates(LineString(c))[[0, -1]] for c, _ in out]
    assert not splits.get(1)
    for pts in splits.values():
        for p in pts:
            assert any(np.hypot(*(e - p).T).min() <= 1e-6 for e in ends)


def test_node_crossings_drop_a_stretch_that_ends_on_the_crossed_row():
    """g crosses row h 1.1 m before the point q it was joined at (h cut
    there): that stretch goes, g ends at the crossing, and h is cut there
    instead of at q."""
    h = LineString([(0, 0), (20, 0)])
    g = np.array([(5.0, 5.0), (9.6, -0.4), (10.0, 0.0)])
    splits = defaultdict(list)
    splits[0].append(np.array([10.0, 0.0]))
    out = vc._node_crossings([[g, 0]], [h], splits)
    assert len(out) == 1
    end = shapely.get_coordinates(LineString(out[0][0]))[-1]
    assert abs(end[1]) < 1e-9 and end[0] < 10.0
    assert [tuple(np.round(p, 6)) for p in splits[0]] == [
        tuple(np.round(end, 6))
    ]


def test_node_crossings_meet_an_existing_cut_at_its_node():
    """Row h is cut at q (an addition joined there); g crosses h 5 mm from
    q: the crossing goes onto q, where g's two pieces then meet h's node."""
    h = LineString([(0, 0), (20, 0)])
    g = np.array([(10.005, 5.0), (10.005, -1.0)])
    splits = defaultdict(list)
    splits[0].append(np.array([10.0, 0.0]))
    out = vc._node_crossings([[g, 0]], [h], splits)
    nodes_h = {
        tuple(np.round(c[k], 6))
        for c in vc._cut_at(shapely.get_coordinates(h), splits[0])
        for k in (0, -1)
    }
    ends = {tuple(np.round(c[k], 6)) for c, _ in out for k in (0, -1)}
    assert not [
        e for e in ends if e not in nodes_h and h.distance(Point(e)) < 1e-6
    ]


def test_kept_beside_of_several_passes_is_listed_once_on_the_final_network():
    """A pass-1 stretch half on a row a later pass redrew keeps the half
    still there; a pass-2 stretch pass 1 already listed is not listed
    again."""
    final = net_([("1-1", "001", 2, LineString([(0, 0), (50, 0)]))])
    crs = final.crs
    p1 = gpd.GeoDataFrame(
        {
            "stretch": ["1-1_001", "1-1_001"],
            "vannforekomstid": ["001-1-R"] * 2,
            "why": ["x", "y"],
            "pass_no": [1, 1],
        },
        geometry=[
            LineString([(0, 0), (40, 0), (40, 30)]),
            LineString([(42, 0), (48, 0)]),
        ],
        crs=crs,
    )
    p2 = gpd.GeoDataFrame(
        {
            "stretch": ["1-1_001"],
            "vannforekomstid": ["001-1-R"],
            "why": ["z"],
            "pass_no": [2],
        },
        geometry=[LineString([(42, 0), (47.9, 0)])],
        crs=crs,
    )
    out = vc._kept_on([p1, p2], final)
    rows = pd.concat(out, ignore_index=True)
    assert list(rows.why) == ["x", "y"]
    assert abs(rows.geometry.length.iloc[0] - 40.0) < 0.02


# ── round 7: crossings of two additions, run-back hooks, near-through nodes,
#    ring retries, cycle labels, kept_beside over passes, exact cut joints ──
def _new_dangling(net, res):
    return sorted(set(dangling(res.network)) - set(dangling(net)))


def test_node_crossings_keep_the_junction_two_crossing_additions_share():
    """m_up, m_dn and t meet at J = (0, 0); t crosses m_dn 1.19 m from J.
    Either one's stretch to J would do as the one dropped -- only one is:
    nothing is left meeting nothing at J."""
    m_up = np.array([(-10.0, 0.0), (0.0, 0.0)])
    m_dn = np.array([(0.0, 0.0), (10.0, 0.0)])
    t = np.array([(5.0, 10.0), (1.0, -0.5), (0.0, 0.0)])
    far = [LineString([(100, 100), (200, 100)])]
    out = vc._node_crossings(
        [[m_up, 0], [m_dn, 1], [t, 2]], far, defaultdict(list)
    )
    geoms = [LineString(c) for c, _ in out]
    free = {(-10.0, 0.0), (10.0, 0.0), (5.0, 10.0)}
    lone = [
        tuple(np.round(e, 3))
        for c, _i in out
        for e in (c[0], c[-1])
        if sum(g.distance(Point(e)) <= 1e-6 for g in geoms) < 2
        and tuple(e) not in free
    ]
    assert not lone


def test_a_lens_crossing_beside_its_junction_stays_joined():
    """Vann-Nett draws one stream twice from J to K, the strands crossing
    1.6 m from J: one strand is cut back there, the other cut, and the
    network stays one piece with no new dangling end."""
    e1 = LineString([(-50, -50), (-50, 50)])
    e2 = LineString([(50, -50), (50, 50)])
    net = net_([("1-1", "001", 2, e1), ("1-2", "001", 2, e2)])
    m_up = LineString([(-50, 0), (-20, 0)])
    x = LineString(
        [(-20, 0), (-19, 0.6), (-18, 0.2), (-10, 0.4), (10, 0.4), (20, 0)]
    )
    y = LineString(
        [(-20, 0), (-19, 0.2), (-18, 0.6), (-10, 0.7), (10, 0.7), (20, 0)]
    )
    m_dn = LineString([(20, 0), (50, 0)])
    res = vc.complete(net, vn_([m_up, x, y, m_dn]))
    assert components(res.network) == 1
    assert not _new_dangling(net, res)


def test_node_crossings_undo_a_truncation_onto_a_stretch_that_goes():
    """k and j both start at J on the row; j's stretch back to J goes (it
    crosses the row), and k would be cut back onto that very stretch: k's
    truncation is undone, so k still meets j."""
    row = LineString([(-10, 0), (10, 0)])
    j_add = np.array([(0.0, 0.0), (0.5, 0.5), (1.5, -0.5), (5.0, -5.0)])
    k_add = np.array([(0.0, 0.0), (0.0, 1.0), (0.6, 0.2), (3.0, 5.0)])
    splits = defaultdict(list)
    splits[0].append(np.array([0.0, 0.0]))
    out = vc._node_crossings([[k_add, 0], [j_add, 1]], [row], splits)
    geoms = [LineString(c) for c, _ in out] + [row]
    k_pieces = [c for c, i in out if i == 0]
    assert [
        e
        for c in k_pieces
        for e in (c[0], c[-1])
        if sum(g.distance(Point(e)) <= 1e-6 for g in geoms) >= 2
    ]


def test_a_head_whose_only_reach_runs_back_moves_its_free_end():
    """The only reach connector would run back along the line (the most
    hooked of all, refused): the capped free walk, within every bound of
    the hook fallback, moves ELVIS's free end to Vann-Nett's."""
    main = LineString([(0, 0), (300, 0)])
    h = LineString(
        [(-6.5, 3), (-4, 3), (-4, 4.5), (-9.5, 4.5), (-11, 2), (-2, 2), (0, 0)]
    )
    net = net_([("1-1", "001", 3, main), ("1-2", "001", 1, h)])
    vn = vn_([LineString([(0, 0), (300, 0)]), LineString([(-16, 0), (0, 0)])])
    res = vc.complete(net, vn)
    assert res.stats["added"] == 0
    (g,) = _pieces_of(res, "1-2_001")
    assert np.allclose(shapely.get_coordinates(g)[0], (-16, 0))


def test_a_reach_two_cm_before_a_through_node_keeps_the_main_stem():
    """As test_a_reach_walk_that_is_not_the_line_s_stream_keeps_the_main_stem,
    but the walk ends 2 cm before the through node N (a vertex 1 m before N
    on the upper row): within JOG_M of N it is held to the same rule."""
    F, N = np.array([-3.3, 30.4]), np.array([0.0, 0.0])
    f = (F - N) / np.hypot(*(F - N))
    p = np.array([f[1], -f[0]])
    eps = np.arcsin(0.02 / np.hypot(*(F - N)))
    u = np.cos(eps) * p - np.sin(eps) * f
    up = LineString([(78, -124), tuple(N - u), (0, 0)])
    down = LineString([(0, 0), (-10.8, 26.6), (-300, 26.6)])
    net = net_([("1-1", "001", 3, up), ("1-1", "002", 3, down)])
    v = LineString([(-3.3, 30.4), (-7.0, 28.5), (-10.8, 26.6), (-300, 26.6)])
    res = vc.complete(net, vn_([v]))
    stem = shapely.ops.substring(down, 0.0, 28.7)
    have = shapely.union_all(_lines(res.network))
    assert not len(res.qa.get("connectors", []))
    assert stem.difference(have.buffer(0.01)).length < 1e-6


def test_a_ring_retry_cuts_no_centimetre_sliver():
    """A 1.3 m bump leaves r1 0.3 m before N and comes back 1.5 cm after it:
    landing the second end 1.5 cm from the first would cut a 1.5 cm row and
    make a ring in all but name, so it is left as a listed tip."""
    r1 = LineString([(-100, 0), (0, 0)])
    r2 = LineString([(0, 0), (100, 0)])
    net = net_([("1-1", "001", 2, r1), ("1-1", "002", 2, r2)])
    side = LineString([(-0.3, 0), (-0.3, 1.5), (0.015, 1.5), (0.015, 0)])
    res = vc.complete(net, vn_([LineString([(-100, 0), (100, 0)]), side]))
    assert not [g.length for g in _lines(res.network) if g.length < 0.1]
    adds = [
        g
        for g, t in zip(
            _lines(res.network), res.network.patch.fillna(""), strict=True
        )
        if "vn_added" in t
    ]
    assert all(
        Point(g.coords[0]).distance(Point(g.coords[-1])) > 0.1 for g in adds
    )


def test_a_cycle_turn_node_is_listed_under_the_turn_that_left_it_so():
    """Both additions into T are turned at level 2; the second turn is the
    one that leaves T a source, and T is listed under it."""
    r1 = LineString([(0, 10), (0, 0)])
    r2 = LineString([(0, 0), (10, 0)])
    r3 = LineString([(0, 0), (-10, 0)])
    adds = [
        [np.array([(10.0, 0.0), (0.0, 10.0)]), 0],
        [np.array([(-10.0, 0.0), (0.0, 10.0)]), 1],
    ]
    out = vc._cycles(adds, [r1, r2, r3], {})
    assert [(k, kind) for k, p, kind in out if np.allclose(p, (0, 10))] == [
        (1, "source")
    ]


def test_kept_beside_is_checked_against_the_final_network_after_every_pass(
    monkeypatch,
):
    """A pass-1 stretch on a row pass 2 moved goes, though pass 2 lists
    nothing itself (the passes are stubbed)."""
    r0 = LineString([(0, 0), (50, 0)])
    net0 = net_([("1-1", "001", 2, r0)])
    net2 = net_([("1-1", "001", 2, LineString([(0, 5), (50, 5)]))])
    kept1 = gpd.GeoDataFrame(
        {"stretch": ["1-1_001"], "vannforekomstid": ["001-1-R"], "why": ["x"]},
        geometry=[LineString([(10, 0), (40, 0)])],
        crs=net0.crs,
    )
    outs = [
        (net0, {"kept_beside": kept1}, 0.05),
        (net2, {}, 0.0),
    ]
    calls = []

    def stub(net, V, regine, P, next_key):
        n, q, after = outs[len(calls)]
        calls.append(1)
        st = {
            "rows_in": 1,
            "rows_out": 1,
            "missing_km_before": 0.1,
            "missing_lines": 1,
        }
        return n, q, {**st, "missing_km_after": after}, next_key

    monkeypatch.setattr(vc, "_complete_pass", stub)
    res = vc.complete(net0, vn_([LineString([(0, 0), (50, 0)])]))
    assert res.stats["passes"] == 2
    assert "kept_beside" not in res.qa or not len(res.qa["kept_beside"])


def test_kept_beside_of_a_later_pass_lists_only_what_no_earlier_one_did():
    """A pass-2 stretch overlapping a pass-1 one by half keeps the other
    half only."""
    final = net_([("1-1", "001", 2, LineString([(0, 0), (50, 0)]))])
    frames = [
        gpd.GeoDataFrame(
            {
                "stretch": ["1-1_001"],
                "vannforekomstid": ["001-1-R"],
                "why": [w],
                "pass_no": [k],
            },
            geometry=[g],
            crs=final.crs,
        )
        for k, w, g in (
            (1, "x", LineString([(10, 0), (20, 0)])),
            (2, "y", LineString([(15, 0), (25, 0)])),
        )
    ]
    rows = pd.concat(vc._kept_on(frames, final, 2), ignore_index=True)
    assert list(rows.why) == ["x", "y"]
    assert abs(rows.geometry.length.iloc[1] - 5.0) < 0.02


def test_cut_pieces_share_their_joint_exactly():
    """A cut point 4e-10 m past a vertex of the row (a T of the old network
    drawn that way): both pieces end on the cut point itself, not one on the
    vertex and one on the point."""
    v = (5.0 - 4e-10, 0.0)
    c = np.array([(0.0, 0.0), v, (10.0, 0.0)])
    a, b = vc._cut_at(c, [np.array([5.0, 0.0])])
    assert tuple(a[-1]) == tuple(b[0]) == (5.0, 0.0)


def test_kept_beside_of_a_later_pass_touching_an_earlier_one_is_kept_whole():
    """A pass-2 stretch that only touches the pass-1 one at its end lists
    nothing twice: it is kept whole, not trimmed by the 1 cm cover."""
    final = net_([("1-1", "001", 2, LineString([(0, 0), (50, 0)]))])
    frames = [
        gpd.GeoDataFrame(
            {
                "stretch": ["1-1_001"],
                "vannforekomstid": ["001-1-R"],
                "why": [w],
                "pass_no": [k],
            },
            geometry=[g],
            crs=final.crs,
        )
        for k, w, g in (
            (1, "x", LineString([(10, 0), (20, 0)])),
            (2, "y", LineString([(20, 0), (21.005, 0)])),
        )
    ]
    rows = pd.concat(vc._kept_on(frames, final), ignore_index=True)
    assert list(rows.why) == ["x", "y"]
    assert abs(rows.geometry.length.iloc[1] - 1.005) < 1e-6


def test_node_crossings_fall_back_to_a_smaller_drop_of_the_same_end():
    """k's drop to its crossing with j lies on j's dropped stretch and gives
    way; k's smaller drop to where it crosses the row (its stretch back to J
    ends on the row too) then holds: no loop is left beside the row."""
    row = LineString([(-10, 0), (10, 0)])
    k = np.array(
        [(0.0, 0.0), (-0.4, 0.4), (-0.4, -0.2), (0.7, -0.2), (3.0, -5.0)]
    )
    j = np.array([(0.0, 0.0), (0.5, -0.5), (1.5, 0.5), (5.0, 5.0)])
    splits = defaultdict(list)
    splits[0].append(np.array([0.0, 0.0]))
    out = vc._node_crossings([[k, 0], [j, 1]], [row], splits)
    assert not [
        c
        for c, _ in out
        if row.distance(Point(c[0])) <= 1e-6
        and row.distance(Point(c[-1])) <= 1e-6
    ]
    assert (-0.4, 0.0) in [tuple(np.round(c[0], 6)) for c, i in out if i == 0]
