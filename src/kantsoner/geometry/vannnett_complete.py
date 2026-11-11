"""Every Vann-Nett river line into the patched network, in its own drawing.

The patch (vannnett_patch.patch) redraws ELVIS on Vann-Nett where it can and
adds what ELVIS lacks, but it leaves Vann-Nett line out: next to an ELVIS
line that kept its own drawing (a failed redraw, the head or tail of a
redrawn stretch), in off runs shorter than min_run, in slivers under
sliver_m, and where an end pull or a route's fade moved Vann-Nett geometry.
On the 2026-10-27 network that was 147.7 km more than 1 m from any line,
104 km of it in Nidelva. The pipeline couples a segment to its water body
only within 1 m of the Vann-Nett line, so that line went undelivered.

``complete`` closes the gap on an existing patched network and changes
nothing else. Missing line -- Vann-Nett more than on_m from every network
line -- is cut out exactly and chained through Vann-Nett nodes into
components; the pieces merge into lines, each end of which is an anchor (the
Vann-Nett line is back on a network row there), a free Vann-Nett end (a
source or an outlet) or a junction inside the component. Lines that follow
one another downstream through such junctions make chains; so do two lines
on either side of a row that Vann-Nett only crosses (more than on_m from the
row's ends), but such a chain is only taken where the same plan redraws the
row it crosses. What ELVIS draws for the same stream is found through the
network's own topology, never by proximity alone:

- Detour: a chain's two ends anchor on one row and the row between them is
  off Vann-Nett. The row takes the chain -- a stream drawn twice between the
  same two points is one stream. Two anchors on one row that the
  component's lines join through its junctions whichever way they run make
  a detour too.
- Duplicate path: the anchors lie on different rows and the shortest path
  between them through the network is, apart from Vann-Nett line at either
  end, off Vann-Nett and close to the chain -- or runs the same way beside it
  nearly all along (a confluence Vann-Nett puts elsewhere). The rows on the
  path take the chain, each its own part of it. Where nothing joins the
  anchors, a chain anchored at the loose end of an earlier addition runs on
  along it to where it meets ELVIS (a node, or another row's interior), and
  the addition, drawn twice once the path takes the whole, goes.
- Duplicate head: one end anchors, the other is a free Vann-Nett end, and a
  path from the anchor runs off Vann-Nett beside the chain to a free ELVIS
  end, which goes where Vann-Nett's is -- or to the point of an ELVIS line
  nearest that end, or on from there along ELVIS where a connector from
  that point would run back along the line, which a straight connector
  joins to it (a walk ending at, or within JOG_M of, a node ELVIS runs on
  through must be close to the chain there: ELVIS drawing the stream
  elsewhere would reroute that stream through Vann-Nett's head). A walk to
  a free end longer than the line by more than max(JOG_M, half the line's
  length) means ELVIS runs on past Vann-Nett's end: its upper course is
  kept and a connector reaches it (or the line is added) -- unless ELVIS
  neither runs on past that end nor ends farther than that from it
  (meander), or the connector would turn back off the line by more than a
  right angle (a hook, or one running back along the line) and the walk
  runs the rows' way all along, longer than the line by at most twice that,
  ELVIS running on past Vann-Nett's end by at most HOOK_M: then the free
  end moves after all.

A node inside a duplicate moves onto the chain with every row at it: the rows
of the path take the chain's line, and the other rows there share one
straight connector (the row taking it is cut at the old node, where the rest
keep their ends) -- unless Vann-Nett draws a row's stream into a junction of
the chain, which the row then takes, junction and all, or a missing line at
the node carries the row's stream on (leaving it the way the row flows), in
which case the row keeps its end there and that line joins it. A row end
otherwise only leaves a node if that leaves no other row there without an
outflow or an inflow. No straight line longer than JOG_M is drawn except
those connectors, each at most CONNECT_M, crossing no line, and listed in the
QA. Whatever the rules do not cover is added as a new row, in Vann-Nett's
drawing, joined where it anchors, oriented to agree with the rows it meets,
and listed in the QA (``kept_beside``: every stretch of a row off Vann-Nett
beside a line another row, an addition or another piece of the same row
carries, with the reason the missing line there was not taken; listed once
over the passes, on the network as it ends up) -- unless the pass's own
lines already carry it within on_m (Vann-Nett drawing one stream twice) and
no other line meets it at a junction and no staying row relies on it; of
two additions joining the same two points within on_m of each other all
along, only the longer is drawn. A loop of more than 2 m leaving the
network at one point is added out and back, or once, flowing back in,
where that carries all of it within on_m; no addition is joined into a ring
(its second end lands on the line beside the first where it can, else it
is left where Vann-Nett has it, and listed). An addition crossing a line
within JOG_M of an end of either gets a node there; so does a row that
earlier additions ended on without one (a T), where the line now taking
them runs through it or other rows still end there. Additions in a directed
cycle are turned round where that breaks it -- failing all else, leaving a
source or sink at a free end of the rows or at the tip of a ring of two
additions, listed under the turned addition; ``flow_conflicts`` lists the
cycles left, additions that leave a node with only inflow or outflow either
way round, any other node the completion leaves so, and any row end left
meeting nothing where the input had others (kind loose); ``loose_ends``
lists ends that could not be joined and addition ends meeting nothing that
are no free end of Vann-Nett's.

A row nothing happens to comes out byte-identical, so a run made on the old
network can be redone for the rows that changed alone. Rows the completion
changes are tagged ``vn_completed`` (its line now follows Vann-Nett in part),
``vn_moved`` (an end carried onto a completed line by a connector) or
``vn_cut`` (cut where a line now joins); additions ``vn_added``.
"""

from __future__ import annotations

import heapq
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from loguru import logger
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import nearest_points, substring

from .vannnett_patch import (
    ADDED,
    PATCH,
    PATCH_VN,
    VN_ID,
    Params,
    Result,
    _ends,
    _matches,
    _node_ids,
    _restore_dtypes,
    _sample,
    _strahler,
    _tangent,
    vn_lines,
)

COMPLETED = (
    "vn_completed"  # the row now follows (part of) a missing Vann-Nett line
)
MOVED = "vn_moved"  # an end carried onto a completed line by a connector
CUT = "vn_cut"  # cut where a completed or added line now joins

JOG_M = 2.0  # the longest straight join drawn at an anchor
SNAP_M = 0.01  # a cut this close to a row's end or vertex goes there
CONNECT_M = 50.0  # the longest connector carrying a node or an end over
OVERLAP_M = 2.5  # two replacements this close on one row become one
BRIDGE_M = 10.0  # loose addition ends this close on one Vann-Nett line join
HOOK_M = 15.0  # ELVIS run-on a free head may lose rather than draw a hook

# reasons a later plan of the same pass can overtake (a line is not planned
# again within the pass): the QA says when they held
WHY_CROSS = (
    "when planned, it crossed an ELVIS line there: no place Vann-Nett leaves it"
)
WHY_NODE = "when planned, an end would have left a node without in- or outflow"
WHY_BACK = "the connector would run back along the line"
NOISE_M = 1e-3  # missing line reaching less than this beyond on_m is noise


# ── helpers ────────────────────────────────────────────────────────────────
def _join(parts):
    """Concatenate coordinate runs, dropping a point only where a seam
    repeats it exactly. (vannnett_patch._join compares with np.allclose's
    default relative tolerance, which at UTM northings is tens of metres:
    it drops a connector's far end.)"""
    parts = [
        np.asarray(p, dtype=float).reshape(-1, 2)
        for p in parts
        if p is not None and len(p)
    ]
    if not parts:
        return np.empty((0, 2))
    c = np.vstack(parts)
    keep = np.r_[True, np.any(np.abs(np.diff(c, axis=0)) > 1e-9, axis=1)]
    return c[keep]


def _sub(line, a, b):
    """Coordinates of line from measure a to b (reversed when a > b)."""
    lo, hi = min(a, b), max(a, b)
    c = shapely.get_coordinates(substring(line, lo, hi))
    if len(c) == 1:
        c = np.vstack([c, c])
    return c[::-1] if a > b else c


def _key(geoms):
    return shapely.to_wkb(shapely.normalize(shapely.set_precision(geoms, 0.01)))


def _pkey(xy):
    return tuple(np.round(np.asarray(xy, dtype=float) * 1000).astype(np.int64))


def _xy(p):
    return shapely.get_coordinates(p)[0]


def _share(a, b, d, step=2.0):
    """Share of line a lying within d of geometry b."""
    n = max(5, int(np.ceil(a.length / step)))
    pts = shapely.line_interpolate_point(a, np.linspace(0, a.length, n + 1))
    return float((shapely.distance(pts, b) <= d).mean())


def _dir_ok(pieces, downstream):
    """Every walked piece runs the way the walk does along Vann-Nett."""
    return all((b > a) == downstream for _, a, b in pieces)


# ── missing line ───────────────────────────────────────────────────────────
def missing(vg, ng, on_m: float) -> pd.DataFrame:
    """Pieces of Vann-Nett line farther than on_m from every network line.

    A part drawn like a row to the centimetre is skipped. The rest is screened
    with samples at most 1 m apart: the distance to the network changes by at
    most the distance moved, so a part none of whose samples is farther than
    on_m - 0.5 m lies within on_m everywhere. A part that fails the screen is
    cut by the on_m buffers of the rows near it, and each piece kept only if
    it really reaches beyond on_m (by more than NOISE_M, sampled every 2 cm:
    the buffer polygon lies up to 1.2 mm inside the true distance at curves).
    Returns part, m0, m1.
    """
    tree = shapely.STRtree(ng)
    have = set(_key(ng))
    rest = np.flatnonzero([k not in have for k in _key(vg)])
    cand = []
    for c0 in range(0, len(rest), 20_000):
        sub = rest[c0 : c0 + 20_000]
        idx, _, pts = _sample(vg[sub], 1.0)
        d = tree.query_nearest(pts, return_distance=True, all_matches=False)[1]
        cand.append(sub[np.unique(idx[d > on_m - 0.5])])
    cand = np.concatenate(cand) if cand else np.empty(0, dtype=np.int64)
    rows = []
    for p in cand:
        g = vg[p]
        near = tree.query(g, predicate="dwithin", distance=on_m + 0.01)
        if not len(near):
            rows.append((int(p), 0.0, float(g.length)))
            continue
        off = g.difference(
            shapely.union_all(shapely.buffer(ng[near], on_m, quad_segs=16))
        )
        for piece in shapely.get_parts(shapely.line_merge(off)):
            if piece.is_empty or piece.length < 1e-4:
                continue
            n = max(20, int(np.ceil(piece.length / 0.02)))
            pts = np.r_[
                shapely.line_interpolate_point(
                    piece, np.linspace(0, piece.length, n + 1)
                ),
                shapely.points(shapely.get_coordinates(piece)),
            ]
            dd = shapely.distance(pts[:, None], ng[near][None, :]).min(axis=1)
            if dd.max() <= on_m + NOISE_M:
                continue
            a = float(g.project(shapely.get_point(piece, 0)))
            b = float(g.project(shapely.get_point(piece, -1)))
            rows.append((int(p), min(a, b), max(a, b)))
    out = pd.DataFrame(rows, columns=["part", "m0", "m1"])
    out = out.sort_values(["part", "m0"]).reset_index(drop=True)
    merged: list[list] = []
    for r in out.itertuples(index=False):
        if merged and merged[-1][0] == r.part and r.m0 - merged[-1][2] < 0.01:
            merged[-1][2] = max(merged[-1][2], r.m1)
        else:
            merged.append([r.part, r.m0, r.m1])
    return pd.DataFrame(merged, columns=["part", "m0", "m1"])


# ── the missing lines ──────────────────────────────────────────────────────
@dataclass
class _Line:
    """One merged line of missing Vann-Nett line, in Vann-Nett's direction.

    ends[w]: ("anchor", row, measure on row, xy) | ("free", xy) | ("inner", xy)
    """

    geom: LineString
    vn_id: object
    ends: list
    comp: int


def _lines(miss, V, ng, on_m) -> list[_Line]:
    vg = V.geometry.values
    vn_ids = V[VN_ID].to_numpy() if VN_ID in V else np.full(len(V), None)
    a, b = _ends(vg)
    (vs, ve), vxy = _node_ids(a, b)
    vdeg = np.bincount(np.r_[vs, ve], minlength=len(vxy))
    L = shapely.length(vg)
    parent = list(range(len(miss)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    at_node = defaultdict(list)
    free_xy = set()
    for k, r in enumerate(miss.itertuples(index=False)):
        for at_end, node in (
            (r.m0 <= 1e-3, vs[r.part]),
            (r.m1 >= L[r.part] - 1e-3, ve[r.part]),
        ):
            if at_end:
                at_node[int(node)].append(k)
                if vdeg[node] == 1:
                    free_xy.add(_pkey(vxy[node]))
    for ks in at_node.values():
        for k in ks[1:]:
            parent[find(k)] = find(ks[0])
    comp = np.array([find(k) for k in range(len(miss))])
    pieces = np.array(
        [
            substring(vg[r.part], r.m0, r.m1)
            for r in miss.itertuples(index=False)
        ],
        dtype=object,
    )
    tree = shapely.STRtree(ng)
    out = []
    for c in np.unique(comp):
        ks = np.flatnonzero(comp == c)
        merged = shapely.line_merge(
            MultiLineString(list(pieces[ks])), directed=True
        )
        for g in shapely.get_parts(merged):
            if g.length < 1e-3:
                continue
            mid = g.interpolate(0.5, normalized=True)
            k = ks[np.argmin(shapely.distance(pieces[ks], mid))]
            ends = []
            for w in (0, -1):
                p = shapely.get_point(g, w)
                xy = _xy(p)
                (_, j), dist = tree.query_nearest(
                    [p], return_distance=True, all_matches=False
                )
                if len(j) and dist[0] <= on_m + 0.05:
                    row = int(j[0])
                    ends.append(("anchor", row, float(ng[row].project(p)), xy))
                elif _pkey(xy) in free_xy:
                    ends.append(("free", xy))
                else:
                    ends.append(("inner", xy))
            out.append(_Line(g, vn_ids[miss.part.iloc[k]], ends, int(c)))
    return out


@dataclass
class _Chain:
    """Lines following one another downstream through the component's own
    junctions, from an end that is no junction to another."""

    ks: list  # the lines, in Vann-Nett's direction
    geom: LineString
    ends: list  # the first line's start, the last line's end
    joints: list  # (measure along geom, junction key) between the lines
    vn_id: object  # the water body most of it belongs to
    n_ids: int  # how many water bodies it runs through


def _chains(lines, cap=64, cross=None, between=None) -> list[_Chain]:
    """Chains of lines downstream through the component's junctions -- and
    across a row a Vann-Nett line only crosses (``cross``: the line ending
    on that row to the one starting there), joined along Vann-Nett
    (``between``)."""
    cross = cross or {}
    starts = defaultdict(list)  # junction -> lines starting there
    for k, ln in enumerate(lines):
        if ln.ends[0][0] == "inner":
            starts[_pkey(ln.ends[0][1])].append(k)
    out = []
    for k0, ln in enumerate(lines):
        if ln.ends[0][0] == "inner" or (k0, 0) in cross:
            continue
        stack, n = [[k0]], 0
        while stack and n < cap:
            path = stack.pop()
            e = lines[path[-1]].ends[1]
            nxt = cross.get((path[-1], 1))
            if nxt is not None and nxt[1] == 0:
                if nxt[0] not in path:
                    stack.append([*path, nxt[0]])
                continue
            if e[0] != "inner":
                out.append(_chain(lines, path, between))
                n += 1
                continue
            for k in starts.get(_pkey(e[1]), ()):
                if k not in path:
                    stack.append([*path, k])
    return out


def _chain(lines, ks, between=None) -> _Chain:
    runs = [shapely.get_coordinates(lines[k].geom) for k in ks]
    parts = [runs[0]]
    for a_run, b_run in zip(runs[:-1], runs[1:], strict=True):
        if between is not None and np.hypot(*(a_run[-1] - b_run[0])) > 1e-9:
            parts.append(between(a_run[-1], b_run[0]))
        parts.append(b_run)
    geom = LineString(_join(parts))
    joints, m = [], 0.0
    for k, nxt in zip(ks[:-1], ks[1:], strict=True):
        m += lines[k].geom.length
        e = lines[k].ends[1]
        if e[0] == "inner":
            joints.append((m, _pkey(e[1])))
        else:  # across a row: what lies between counts too
            m += float(
                np.hypot(*(runs[ks.index(nxt)][0] - runs[ks.index(k)][-1]))
            )
    w8 = defaultdict(float)
    for k in ks:
        w8[lines[k].vn_id] += lines[k].geom.length
    return _Chain(
        list(ks),
        geom,
        [lines[ks[0]].ends[0], lines[ks[-1]].ends[1]],
        joints,
        max(w8, key=w8.get),
        len(w8),
    )


# ── the network as a graph ─────────────────────────────────────────────────
class _Graph:
    """Rows between end nodes (to the millimetre), digitising order = flow."""

    def __init__(self, ng):
        self.ng = ng
        self.L = shapely.length(ng)
        a, b = _ends(ng)
        (self.s, self.e), self.xy = _node_ids(a, b)
        self.at = defaultdict(set)  # node -> {(row, end)}, end 0 start 1 end
        for r in range(len(ng)):
            self.at[int(self.s[r])].add((r, 0))
            self.at[int(self.e[r])].add((r, 1))

    def node(self, r, w):
        return int(self.s[r] if w == 0 else self.e[r])

    def deg(self, n):
        return len(self.at[n])

    def can_leave(self, n, r, w):
        """May row r's end w leave node n? Not if the node would be left
        with an inflow and no outflow, or an outflow and no inflow."""
        rest = self.at[n] - {(r, w)}
        if not rest:
            return True
        ins = sum(1 for _, x in rest if x == 1)
        outs = len(rest) - ins
        return not (
            (w == 0 and outs == 0 and ins > 0)
            or (w == 1 and ins == 0 and outs > 0)
        )

    def leave(self, n, r, w):
        self.at[n].discard((r, w))

    def paths_from(self, rows, start, cutoff, exclude=()):
        """Dijkstra over the given rows from a position (row, measure), or
        from a node when start is an int: {node: (distance, previous node,
        piece)} where a piece is (row, measure from, measure to) in walking
        order. The start row and the rows in exclude are not walked
        through."""
        rows = {int(x) for x in rows}
        skip = {int(x) for x in exclude}
        if isinstance(start, (int, np.integer)):
            seeds = [(int(start), 0.0, None)]
        else:
            r0, t0 = start
            skip.add(int(r0))
            seeds = [
                (int(self.s[r0]), t0, (r0, t0, 0.0)),
                (int(self.e[r0]), self.L[r0] - t0, (r0, t0, float(self.L[r0]))),
            ]
        adj = defaultdict(list)
        for r in rows:
            if r in skip:
                continue
            u, v = int(self.s[r]), int(self.e[r])
            adj[u].append((v, self.L[r], (r, 0.0, float(self.L[r]))))
            adj[v].append((u, self.L[r], (r, float(self.L[r]), 0.0)))
        best = {}
        heap = []
        for n, d, piece in seeds:
            if d < best.get(n, (np.inf,))[0]:
                best[n] = (d, None, piece)
                heapq.heappush(heap, (d, n))
        while heap:
            d, u = heapq.heappop(heap)
            if d > best[u][0] or d > cutoff:
                continue
            for v, w, piece in adj[u]:
                nd = d + w
                if nd < best.get(v, (np.inf,))[0] and nd <= cutoff:
                    best[v] = (nd, u, piece)
                    heapq.heappush(heap, (nd, v))
        return best

    @staticmethod
    def trace(best, n):
        """Pieces from the start position to node n, in walking order."""
        out = []
        while n is not None:
            d, prev, piece = best[n]
            if piece is not None:
                out.append(piece)
            n = prev
        return out[::-1]


# ── one pass ───────────────────────────────────────────────────────────────
@dataclass
class _Rep:
    """Row measures [lo, hi] replaced by body, in the row's direction, joined
    to the row's own point at lo (p_lo) and at hi (p_hi) where those are kept
    (None: the row's end moves to the body's end)."""

    lo: float
    hi: float
    body: np.ndarray
    line: int
    how: str
    p_lo: np.ndarray | None = None
    p_hi: np.ndarray | None = None

    @property
    def coords(self):
        return _join([self.p_lo, self.body, self.p_hi])


@dataclass
class _Ctx:
    """What one pass knows while it plans."""

    lines: list
    ng: np.ndarray
    G: _Graph
    n_tree: shapely.STRtree
    P: Params
    off_frac: object
    reps: dict
    used: set = field(default_factory=set)  # lines taken
    moved: set = field(default_factory=set)  # nodes moved onto a line
    locked: set = field(default_factory=set)  # rows given a connector
    inner_at: dict = field(default_factory=dict)  # junction -> [(line, end)]
    reason: dict = field(default_factory=dict)  # line -> why not taken
    cross: dict = field(default_factory=dict)  # end -> the end across a row
    vg: np.ndarray | None = None
    v_tree: shapely.STRtree | None = None
    added: set = field(default_factory=set)  # rows an earlier patch added
    dropped: set = field(default_factory=set)  # rows whose line others carry
    line_tree: shapely.STRtree | None = None  # the missing lines
    anchor_tree: shapely.STRtree | None = None  # their anchored ends
    anchor_ends: list = field(default_factory=list)  # (line, end) per point
    carry: set = field(default_factory=set)  # lines a staying row relies on


def complete(
    net: gpd.GeoDataFrame,
    vn: gpd.GeoDataFrame,
    regine: gpd.GeoDataFrame | None = None,
    params: Params | None = None,
    max_passes: int = 4,
) -> Result:
    """The network with every Vann-Nett river line in it (see the module).

    Replacing a stretch of row can uncover a few metres of some other
    Vann-Nett line it also lay on; the next pass takes those in, so passes
    repeat until nothing is missing or a pass finds nothing it can take.
    """
    P = params or Params()
    V = vn_lines(vn.to_crs(net.crs))
    stats: dict[str, float] = {}
    qa: dict[str, list] = defaultdict(list)
    first = None
    net_in = net
    next_key = _next_vn_number(net)  # additions never reuse a key
    for k in range(1, max_passes + 1):
        net, q, st, next_key = _complete_pass(net, V, regine, P, next_key)
        first = first or st
        for name, g in q.items():
            if len(g):
                qa[name].append(g.assign(pass_no=k))
        for key, v in st.items():
            if key in ("rows_in", "rows_out", "missing_km_before"):
                continue
            stats[key] = (
                v if key == "missing_km_after" else stats.get(key, 0) + v
            )
        stats["passes"] = k
        if st["missing_km_after"] == 0.0 or not st.get("missing_lines", 0):
            break
    stats["rows_in"] = first["rows_in"]
    stats["missing_km_before"] = first["missing_km_before"]
    stats["rows_out"] = len(net)
    if stats["passes"] > 1 and qa.get("kept_beside"):
        qa["kept_beside"] = _kept_on(qa["kept_beside"], net, stats["passes"])
        stats["kept_beside"] = sum(len(g) for g in qa["kept_beside"])
        if not qa["kept_beside"]:
            del qa["kept_beside"]
    # every node the passes left with only inflow or only outflow that the
    # input does not have (where a pass moved an old one, say), and every
    # row end left meeting nothing where the input had others, is listed
    listed = {
        _pkey(xy)
        for g in qa.get("flow_conflicts", [])
        for xy in shapely.get_coordinates(g.geometry.values)
    }
    before = _flow_ends(net_in)
    now = _flow_ends(net)
    for key, v in _lone_ends(net_in, net).items():
        now.setdefault(key, v)
    vn_col = net[PATCH_VN] if PATCH_VN in net else None
    extra = [
        (
            Point(xy),
            _rowkey(net, r),
            None if vn_col is None else vn_col.iat[r],
            kind,
        )
        for key, (kind, r, xy) in now.items()
        # a lone row end is new wherever it is, an old sink's node or not
        if key not in listed and (key not in before or kind == "loose")
    ]
    if extra and first is not None:
        qa["flow_conflicts"].append(
            gpd.GeoDataFrame(
                extra,
                columns=["geometry", "stretch", VN_ID, "kind"],
                crs=net.crs,
            ).assign(pass_no=stats["passes"])
        )
    stats["flow_conflicts"] = stats.get("flow_conflicts", 0) + len(extra)
    out_qa = {
        name: gpd.GeoDataFrame(pd.concat(gs, ignore_index=True), crs=net.crs)
        for name, gs in qa.items()
    }
    if "still_missing" in out_qa:
        last = out_qa["still_missing"]
        out_qa["still_missing"] = last[last.pass_no == stats["passes"]]
    if stats["missing_km_after"] > 0:
        logger.warning(
            f"{stats['missing_km_after'] * 1000:.1f} m of Vann-Nett still more "
            f"than {P.on_m} m from the network after {stats['passes']} passes"
        )
    return Result(net, out_qa, stats)


def _kept_on(frames, net, last=None):
    """kept_beside of several passes, against the final network: a stretch
    an earlier pass listed keeps only what still lies on the network (within
    1 cm; a later pass may have redrawn its row), and a later pass lists
    only what no earlier one did for the same stretch (within 1 cm). Parts
    of 1 m or more are kept."""
    geoms = shapely.get_parts(net.geometry.values)
    tree = shapely.STRtree(geoms)
    last = last or max(int(f.pass_no.max()) for f in frames)
    seen = defaultdict(list)  # stretch -> geometries listed by earlier passes

    def parts(g):
        return [
            x
            for x in shapely.get_parts(shapely.line_merge(g))
            if x.geom_type == "LineString" and x.length >= 1.0
        ]

    out = []
    for f in sorted(frames, key=lambda f: int(f.pass_no.min())):
        keep, rows = [], []
        for idx, g in zip(f.index, f.geometry.values, strict=True):
            pieces = [g]
            earlier = [
                h for h in seen[f.at[idx, "stretch"]] if h.distance(g) <= 0.01
            ]
            if earlier:  # what an earlier pass lists goes
                cover = shapely.union_all([h.buffer(0.01) for h in earlier])
                if g.intersection(cover).length > 0.02 + 1e-6:  # an overlap
                    pieces = parts(g.difference(cover))
            if int(f.at[idx, "pass_no"]) < last:
                on = []
                for pc in pieces:
                    near = tree.query(pc, predicate="dwithin", distance=0.01)
                    if not len(near):
                        continue  # its row is drawn elsewhere now
                    buf = shapely.union_all(geoms[near]).buffer(0.01)
                    on += (
                        [pc] if pc.within(buf) else parts(pc.intersection(buf))
                    )
                pieces = on
            keep += [idx] * len(pieces)
            rows += pieces
        for idx, g in zip(keep, rows, strict=True):
            seen[f.at[idx, "stretch"]].append(g)
        if keep:
            sub = pd.DataFrame(f.drop(columns="geometry")).loc[keep]
            sub = sub.reset_index(drop=True)
            out.append(gpd.GeoDataFrame(sub, geometry=rows, crs=f.crs))
    return out


def _next_vn_number(net):
    """The first VN number no row of the network uses."""
    nums = [
        int(m.group(1))
        for m in (
            re.fullmatch(r"VN(\d+)", str(s))
            for s in net["strekningLnr"].dropna()
        )
        if m
    ]
    return max(nums, default=0) + 1


def _complete_pass(net, V, regine, P, next_key):
    """One pass of the completion: (network, qa layers, stats, the next
    addition number)."""
    crs = net.crs
    net = net.reset_index(drop=True)
    parts, pidx = shapely.get_parts(net.geometry.values, return_index=True)
    if (
        len(parts) != len(net)
        or (np.bincount(pidx, minlength=len(net)) != 1).any()
    ):
        raise ValueError("complete() needs one line part per row")
    ng = parts
    vg = V.geometry.values
    v_tree = shapely.STRtree(vg)
    stats: dict[str, float] = {"rows_in": len(net)}

    miss = missing(vg, ng, P.on_m)
    stats["missing_km_before"] = float((miss.m1 - miss.m0).sum() / 1000)
    logger.info(
        f"Vann-Nett: {len(miss):,} pieces, {stats['missing_km_before']:,.2f} km "
        f"more than {P.on_m} m from the network"
    )
    if miss.empty:
        stats["missing_km_after"] = 0.0
        return net.copy(), {}, stats, next_key
    lines = _lines(miss, V, ng, P.on_m)
    stats["missing_lines"] = len(lines)
    if not lines:
        stats["missing_km_after"] = stats["missing_km_before"]
        return net.copy(), {}, stats, next_key
    G = _Graph(ng)
    n_tree = shapely.STRtree(ng)

    def off_frac(coords):
        """Share of a stretch of network line lying on no Vann-Nett line."""
        g = LineString(coords) if not isinstance(coords, LineString) else coords
        if g.length < 1e-6:
            return 0.0
        n = max(4, int(np.ceil(g.length / P.step)))
        m = np.linspace(0, g.length, n + 1)[1:-1]
        pts = shapely.line_interpolate_point(g, m)
        tang = _tangent(np.array([g] * len(m), dtype=object), m)
        on = _matches(pts, tang, vg, v_tree, P.on_m, P.max_angle)[0] >= 0
        return float((~on).mean())

    reps: dict[int, list[_Rep]] = defaultdict(list)
    ctx = _Ctx(lines, ng, G, n_tree, P, off_frac, reps, vg=vg, v_tree=v_tree)
    if PATCH in net:
        tags_in = net[PATCH].fillna("").astype(str).str.split(",")
        ctx.added = {i for i, t in enumerate(tags_in) if ADDED in t}
    for k, ln in enumerate(lines):
        for w, e in enumerate(ln.ends):
            if e[0] == "inner":
                ctx.inner_at.setdefault(_pkey(e[1]), []).append((k, w))
            elif e[0] == "anchor":
                ctx.anchor_ends.append((k, w))
    ctx.line_tree = shapely.STRtree([ln.geom for ln in lines])
    ctx.anchor_tree = shapely.STRtree(
        shapely.points(
            np.array(
                [lines[k].ends[w][3] for k, w in ctx.anchor_ends], dtype=float
            ).reshape(-1, 2)  # a pass may have no anchored end at all
        ).reshape(-1)
    )
    crossing = ctx.cross = _crossings(lines, ng, P)

    def between(a, b):
        return _vn_path(vg, v_tree, a, b)

    chains = [
        c
        for c in _chains(lines, cross=crossing, between=between)
        if (c.ks[0], 0) not in crossing and (c.ks[-1], 1) not in crossing
    ]
    for k, _w in crossing:
        ctx.reason.setdefault(k, WHY_CROSS)
    attach = []  # (kind, owner, end, Vann-Nett point, row to prefer)
    conns = []  # (row, end, point on the line, the node it leaves)
    qa_rep, qa_kept, qa_conn = [], [], []

    # 1. detours: a chain's two ends anchor on one row
    def width(c):
        return abs(c.ends[1][2] - c.ends[0][2])

    detoured = []  # (row, chain) for the junction rule below

    for c in sorted(
        (
            c
            for c in chains
            if c.ends[0][0] == c.ends[1][0] == "anchor"
            and c.ends[0][1] == c.ends[1][1]
        ),
        key=lambda c: -width(c),
    ):
        if any(k in ctx.used for k in c.ks):
            continue
        r, t0, t1 = c.ends[0][1], c.ends[0][2], c.ends[1][2]
        lo, hi = min(t0, t1), max(t0, t1)
        if hi - lo < 1.0 and c.geom.length > 2.0:
            continue  # a loop off one point of the row: added beside it
        if hi - lo < 1.0 and any(
            x.lo - 1.0 < hi and lo < x.hi + 1.0 for x in reps[r]
        ):
            continue  # a second sub-metre detour at one point: added
        if hi - lo >= 1.0 and off_frac(_sub(ng[r], lo, hi)) < 0.5:
            continue  # the row between is Vann-Nett: a second channel, added
        if r in ctx.locked:
            continue
        if _crossed_rows(ctx, c.ks):
            continue  # the row it crosses would be left without a node
        body = shapely.get_coordinates(c.geom)
        rep = _Rep(
            lo,
            hi,
            body[::-1] if t0 > t1 else body,
            c.ks[0],
            "detour",
            _xy(ng[r].interpolate(lo)),
            _xy(ng[r].interpolate(hi)),
        )
        if not _fits(reps[r], rep):
            continue
        reps[r].append(rep)
        ctx.used.update(c.ks)
        qa_rep.append((LineString(_sub(ng[r], lo, hi)), r, c.ks[0], "detour"))
        detoured.append((r, lo, hi, c))
    _junction_rows(ctx, detoured, attach, qa_rep)

    # 2./3. duplicate paths and heads, found through the network; the
    # chains of one water body first, the longest first
    for c in sorted(chains, key=lambda c: (c.n_ids, -c.geom.length)):
        if any(k in ctx.used for k in c.ks):
            continue
        c.ends = [lines[c.ks[0]].ends[0], lines[c.ks[-1]].ends[1]]  # relinks
        kinds = [c.ends[0][0], c.ends[1][0]]
        if kinds == ["anchor", "anchor"]:
            if c.ends[0][1] == c.ends[1][1]:
                continue  # a second channel or loop beside its row: added
            status, payload = _duplicate_path(c, ctx)
        elif sorted(kinds) == ["anchor", "free"]:
            status, payload = _duplicate_head(c, ctx)
        else:
            continue
        if status == "ok" and _crossed_rows(ctx, c.ks) - {
            r for r, _ in payload["reps"]
        }:
            # a row the chain runs across that the plan does not redraw
            # would be crossed without a node
            status, payload = (
                "kept",
                (
                    payload["dgeom"],
                    WHY_CROSS,
                ),
            )
        if status != "ok":
            # the QA lists it only if an ELVIS line is still beside it in
            # the end (another line may yet replace this one)
            why = payload if status == "reject" else payload[1]
            for k in c.ks:
                ctx.reason.setdefault(k, why)
            continue
        for r, rep in payload["reps"]:
            reps[r].append(rep)
            qa_rep.append(
                (LineString(_sub(ng[r], rep.lo, rep.hi)), r, rep.line, rep.how)
            )
        for n, r, w in payload["leave"]:
            G.leave(n, r, w)
        ctx.moved |= payload["moved"]
        ctx.locked |= payload["locked"]
        ctx.used |= payload["used"]
        ctx.dropped |= payload.get("drop", set())
        ctx.carry |= payload.get("carry", set())
        attach += payload["attach"]
        conns += payload["conns"]
        qa_conn += payload["qa_conn"]
        for kk, ww, r_x, w_x in payload.get("relink", ()):
            # the line carrying a staying row's stream on joins that row's
            # end (it was anchored on the row the plan moves away)
            e = lines[kk].ends[ww]
            lines[kk].ends[ww] = (
                "anchor",
                r_x,
                0.0 if w_x == 0 else float(G.L[r_x]),
                e[-1],
            )

    # 3b. two anchors on one row that the component's lines join through
    # its junctions whichever way they run (Vann-Nett digitising two lines
    # out of one junction makes no chain): a detour, drawn along the row
    detoured = []
    for r, path, t0, t1 in _any_way_detours(lines, ctx):
        lo, hi = min(t0, t1), max(t0, t1)
        if hi - lo < 1.0 or r in ctx.locked:
            continue
        if any(k in ctx.used for k, _ in path):
            continue  # an earlier detour (or rule) took a line of it
        if off_frac(_sub(ng[r], lo, hi)) < 0.5:
            continue
        body = _join(
            [
                shapely.get_coordinates(lines[k].geom)[:: 1 if fwd else -1]
                for k, fwd in path
            ]
        )
        rep = _Rep(
            lo,
            hi,
            body[::-1] if t0 > t1 else body,
            path[0][0],
            "detour",
            _xy(ng[r].interpolate(lo)),
            _xy(ng[r].interpolate(hi)),
        )
        if not _fits(reps[r], rep):
            continue
        reps[r].append(rep)
        ctx.used.update(k for k, _ in path)
        qa_rep.append(
            (LineString(_sub(ng[r], lo, hi)), r, path[0][0], "detour")
        )
        joints, mm = [], 0.0
        for k, fwd in path[:-1]:
            mm += lines[k].geom.length
            joints.append((mm, _pkey(lines[k].ends[1 if fwd else 0][1])))
        chain = _Chain(
            [k for k, _ in path], LineString(body), [], joints, None, 1
        )
        detoured.append((r, lo, hi, chain))
    _junction_rows(ctx, detoured, attach, qa_rep)

    # 4. rows: replaced stretches put in, nodes moved
    new_geom: dict[int, np.ndarray] = {}
    for r, rl in reps.items():
        if not rl:
            continue
        rl = _merge_reps(sorted(rl, key=lambda x: x.lo), vg, v_tree)
        cs, cur = [], 0.0
        for rep in rl:
            if rep.lo > cur + 1e-9:
                cs.append(_sub(ng[r], cur, rep.lo))
            cs.append(rep.coords)
            cur = max(cur, rep.hi)
        if cur < G.L[r] - 1e-9:
            cs.append(_sub(ng[r], cur, G.L[r]))
        new_geom[r] = _join(cs)
    moved_rows = set()
    by_node = defaultdict(list)
    for r, w, x, n_xy in conns:
        by_node[_pkey(n_xy)].append((r, w, x, n_xy))
    for group in by_node.values():
        # one row takes the connector; the others keep their ends at the
        # old node, which that row is cut at
        r, w, x, n_xy = group[0]
        c = new_geom.get(r)
        c = shapely.get_coordinates(ng[r]) if c is None else c
        new_geom[r] = _join([np.vstack([x, c] if w == 0 else [c, x])])
        moved_rows.add(r)
        if len(group) > 1:  # at the row's own point (the node is rounded)
            attach.append(("cut", -1, 0, c[0] if w == 0 else c[-1], r))
    seen_conn, uniq = set(), []
    for item in qa_conn:
        key = tuple(np.round(np.asarray(item[0].coords).ravel(), 3))
        if key not in seen_conn:
            seen_conn.add(key)
            uniq.append(item)
    qa_conn = uniq

    # 5. what no row took is added -- but not a line that a line this pass
    # put in (a replaced stretch, or a longer addition: Vann-Nett drawing
    # one stream twice) already carries within on_m, unless other lines
    # meet it at a junction or a staying row relies on it
    covered = _covered(lines, ctx, new_geom)
    for i in covered:
        ctx.reason[i] = "carried within on_m by a line this pass put in"
    additions, tips = [], []
    for i, ln in enumerate(lines):
        if i in ctx.used or i in covered:
            continue
        if _is_loop(ln):
            # a loop off one point: its two halves, out and back -- or one,
            # flowing back in, where that carries all of it within on_m
            halves = _loop_halves(ln.geom, P.on_m)
            for half, w in halves:
                attach.append(
                    ("add", len(additions), w, ln.ends[w][3], ln.ends[w][1])
                )
                additions.append([half, i])
            if len(halves) == 1:
                half, w = halves[0]
                tips.append(("tip", np.asarray(half[-1 if w == 0 else 0])))
            continue
        k = len(additions)
        additions.append([shapely.get_coordinates(ln.geom), i])
        for w, e in enumerate(ln.ends):
            if e[0] == "anchor":
                attach.append(("add", k, w, e[3], e[1]))
            elif e[0] == "inner":
                attach.append(("inner", k, w, e[1], None))

    # 6. join loose ends to the final lines, cut those lines there
    final = [
        LineString()
        if r in ctx.dropped
        else LineString(new_geom[r])
        if r in new_geom
        else ng[r]
        for r in range(len(ng))
    ]
    additions, attach, freed = _bridge(additions, attach, final, vg, P)
    splits, loose = _attach_all(attach, final, new_geom, additions, P)
    additions = _node_crossings(additions, final, splits)
    additions, twice = _drawn_twice(additions, P)
    for i in twice:
        covered.add(i)
        ctx.reason[i] = "carried within on_m by a line this pass put in"
    loose = freed + loose + tips
    loose += _dangling(additions, final, loose, vg, v_tree)
    conflicts = _orient(additions, final, splits)
    conflicts += _cycles(additions, final, splits)
    pairs = _pairs(final, additions, new_geom, vg, v_tree, V, lines, ctx)

    out_coords, src, cut = [], [], []
    for r in range(len(ng)):
        if r in ctx.dropped:
            continue  # its line is carried by the rows that took it over
        c, pts = new_geom.get(r), splits.get(r, [])
        if c is None and not pts:
            out_coords.append(None)
            src.append(r)
            cut.append(False)
            continue
        c = shapely.get_coordinates(ng[r]) if c is None else c
        pieces = _cut_at(c, pts)
        if c is not new_geom.get(r) and len(pieces) == 1:
            out_coords.append(None)  # nothing changed after all
            src.append(r)
            cut.append(False)
            continue
        for piece in pieces:
            out_coords.append(piece)
            src.append(r)
            cut.append(len(pieces) > 1)
    for c, i in additions:
        out_coords.append(c)
        src.append(-1 - i)
        cut.append(False)

    tags = defaultdict(set)
    for r in new_geom:
        if reps.get(r):
            tags[r].add(COMPLETED)
    for r in moved_rows:
        tags[r].add(MOVED)
    patch_vn = {}
    for r, rl in reps.items():
        if rl:
            w8 = defaultdict(float)
            for x in rl:
                w8[lines[x.line].vn_id] += x.hi - x.lo
            patch_vn[r] = max(w8, key=w8.get)
    out = _assemble(
        net, out_coords, src, cut, tags, patch_vn, lines, regine, next_key
    )
    next_key += len(additions)
    add_keys = out.loc[np.asarray(src) < 0, ["elvID", "strekningLnr"]]
    add_keys = (add_keys.elvID + "_" + add_keys.strekningLnr).to_list()

    for g, r, vid, why in pairs:
        qa_kept.append((g, _rowkey(net, r), vid, why))

    stats["rows_completed"] = sum(1 for r in new_geom if reps.get(r))
    stats["rows_moved"] = len(moved_rows)
    stats["rows_dropped"] = len(ctx.dropped)
    stats["rows_cut"] = sum(1 for r in splits if splits[r])
    stats["added"] = len(additions)
    stats["added_km"] = float(
        sum(LineString(c).length for c, _ in additions) / 1000
    )
    stats["covered"] = len(covered)
    stats["replaced_row_km"] = float(sum(g.length for g, *_ in qa_rep) / 1000)
    stats["kept_beside"] = len(qa_kept)
    stats["connectors"] = len(qa_conn)
    stats["connector_km"] = float(sum(d for _, _, d, _ in qa_conn) / 1000)
    stats["flow_conflicts"] = len(conflicts)
    stats["loose_ends"] = len(loose)
    left = missing(vg, shapely.get_parts(out.geometry.values), P.on_m)
    stats["missing_km_after"] = float((left.m1 - left.m0).sum() / 1000)
    stats["rows_out"] = len(out)
    qa = {
        "missing_before": gpd.GeoDataFrame(
            {
                "comp": [ln.comp for ln in lines],
                VN_ID: [ln.vn_id for ln in lines],
                "end0": [ln.ends[0][0] for ln in lines],
                "end1": [ln.ends[1][0] for ln in lines],
                "taken": [i in ctx.used for i in range(len(lines))],
                "why": [
                    None if i in ctx.used else ctx.reason.get(i)
                    for i in range(len(lines))
                ],
            },
            geometry=[ln.geom for ln in lines],
            crs=crs,
        ),
        "replaced": gpd.GeoDataFrame(
            [
                (g, _rowkey(net, r), lines[i].vn_id, how)
                for g, r, i, how in qa_rep
            ],
            columns=["geometry", "stretch", VN_ID, "how"],
            crs=crs,
        ),
        "kept_beside": gpd.GeoDataFrame(
            qa_kept, columns=["geometry", "stretch", VN_ID, "why"], crs=crs
        ),
        "connectors": gpd.GeoDataFrame(
            [
                (g, _rowkey(net, r), round(d, 2), kind)
                for g, r, d, kind in qa_conn
            ],
            columns=["geometry", "stretch", "length_m", "kind"],
            crs=crs,
        ),
        "dropped": gpd.GeoDataFrame(
            [(ng[r], _rowkey(net, r)) for r in sorted(ctx.dropped)],
            columns=["geometry", "stretch"],
            crs=crs,
        ),
        "flow_conflicts": gpd.GeoDataFrame(
            [
                (Point(xy), add_keys[k], lines[additions[k][1]].vn_id, kind)
                for k, xy, kind in conflicts
            ],
            columns=["geometry", "stretch", VN_ID, "kind"],
            crs=crs,
        ),
        "loose_ends": gpd.GeoDataFrame(
            [(Point(xy), kind) for kind, xy in loose],
            columns=["geometry", "kind"],
            crs=crs,
        ),
        "still_missing": gpd.GeoDataFrame(
            geometry=[
                substring(vg[r.part], r.m0, r.m1) for r in left.itertuples()
            ],
            crs=crs,
        ),
    }
    logger.info(f"completion pass: {stats}")
    return out, qa, stats, next_key


def _crossed_rows(ctx, ks):
    """Rows a chain of the lines ks runs across (a line of it ends on the
    row and the next starts there)."""
    return {
        ctx.lines[k].ends[1][1]
        for k in ks[:-1]
        if ctx.lines[k].ends[1][0] == "anchor"
    }


def _snap_vertex(geom, m):
    """The point at measure m of geom, or the vertex there when one lies
    within 1e-6 (a junction of the lines joined into geom)."""
    x = _xy(geom.interpolate(m))
    vx = shapely.get_coordinates(geom)
    d = np.hypot(*(vx - x).T)
    return vx[int(np.argmin(d))].copy() if d.min() <= 1e-6 else x


def _junction_rows(ctx, detoured, attach, qa_rep):
    """A row ending at a detoured row's end node whose stream Vann-Nett draws
    into a junction of the detour takes that drawing, junction and all; the
    detoured row is cut there."""
    G, P, ng = ctx.G, ctx.P, ctx.ng
    for r, lo, hi, c in detoured:
        if not c.joints:
            continue
        for w_end, gap in ((0, lo), (1, G.L[r] - hi)):
            if gap > P.near_m:
                continue
            n = G.node(r, w_end)
            for r_x, w_x in sorted(G.at[n] - {(r, w_end)}):
                if r_x in ctx.locked or not G.can_leave(n, r_x, w_x):
                    continue
                hit = _junction_line(ctx, c.joints, r_x, w_x, {r})
                if hit is None:
                    continue
                m, jrep, lts = hit
                x = _snap_vertex(c.geom, m)
                body = jrep.body.copy()
                body[-1 if w_x == 1 else 0] = x
                jrep.body = body
                ctx.reps[r_x].append(jrep)
                ctx.used |= lts
                G.leave(n, r_x, w_x)
                attach.append(("cut", -1, 0, x, r))
                qa_rep.append(
                    (
                        LineString(_sub(ng[r_x], jrep.lo, jrep.hi)),
                        r_x,
                        jrep.line,
                        "junction",
                    )
                )


def _any_way_detours(lines, ctx):
    """(row, [(line, forward)], measure at start, measure at end) for two
    anchors on one row that lines not yet taken join through the
    component's junctions, either way; widest first."""
    by_comp = defaultdict(list)
    for i, ln in enumerate(lines):
        if i not in ctx.used:
            by_comp[ln.comp].append(i)
    out = []
    for ks in by_comp.values():
        at = defaultdict(list)  # junction -> (line, end)
        anchors = []
        for k in ks:
            for w, e in enumerate(lines[k].ends):
                if e[0] == "anchor" and (k, w) not in ctx.cross:
                    anchors.append((k, w, e[1], e[2]))
                elif e[0] == "inner":
                    at[_pkey(e[1])].append((k, w))
        for a in range(len(anchors)):
            for b in range(a + 1, len(anchors)):
                k1, w1, r1, t1 = anchors[a]
                k2, w2, r2, t2 = anchors[b]
                if r1 != r2 or (k1, w1) == (k2, w2):
                    continue
                path = _component_path(lines, at, k1, w1, k2, w2)
                if path is not None:
                    out.append((r1, path, t1, t2))
    out.sort(key=lambda x: -abs(x[3] - x[2]))
    return out


def _component_path(lines, at, k1, w1, k2, w2):
    """Lines walked from end w1 of line k1 to end w2 of line k2 through
    junctions, either way: [(line, walked forward)], or None."""
    start = (k1, w1 == 0)  # leaving line k1 from end w1
    stack = [(start, [start])]
    seen = {k1}
    while stack:
        (k, fwd), path = stack.pop()
        far = 1 if fwd else 0
        if k == k2 and far == w2:
            return path
        e = lines[k].ends[far]
        if e[0] != "inner":
            continue
        for kk, ww in at[_pkey(e[1])]:
            if kk in seen:
                continue
            seen.add(kk)
            nxt = (kk, ww == 0)
            stack.append((nxt, [*path, nxt]))
    return None


def _along(g, row, d, cos_min, step=2.0):
    """Share of line g running beside row: within d of it, the same way
    (either sense) within the angle whose cosine is cos_min."""
    n = max(5, int(np.ceil(g.length / step)))
    m = np.linspace(0, g.length, n + 1)
    pts = shapely.line_interpolate_point(g, m)
    tg = _tangent(np.array([g] * len(m), dtype=object), m)
    t = shapely.line_locate_point(row, pts)
    rt = _tangent(np.array([row] * len(m), dtype=object), t)
    ok = (shapely.distance(pts, row) <= d) & (
        np.abs(np.sum(tg * rt, axis=1)) >= cos_min
    )
    return float(ok.mean())


def _crossings(lines, ng, P):
    """Anchor ends where a Vann-Nett line only runs across a row: one line
    ends and the next starts there (within 2 on_m + 1 m along the row, and
    more than on_m from the row's ends -- a line meeting a row at its end
    meets it at a junction), both steeply to the row (more than max_angle),
    running on the same way, and neither of them a loop off that row or
    running beside it (as Vann-Nett zig-zagging about its own ELVIS line
    does). Vann-Nett does not leave the row there -- a tributary on its way
    to a junction drawn elsewhere -- so no duplicate starts or ends there."""
    cos_max = np.cos(np.radians(P.max_angle))
    by_row = defaultdict(list)
    for k, ln in enumerate(lines):
        g = ln.geom
        L = g.length
        for w, e in enumerate(ln.ends):
            if e[0] != "anchor":
                continue
            a, b = (0.0, min(5.0, L)) if w == 0 else (max(0.0, L - 5.0), L)
            d = _xy(g.interpolate(b)) - _xy(g.interpolate(a))
            n = float(np.hypot(*d))
            if n < 1e-6:
                continue
            other = ln.ends[1 - w]
            lobe = other[0] == "anchor" and other[1] == e[1]
            by_row[e[1]].append((e[2], k, w, d / n, lobe))
    out = {}  # (line, end) -> (line, end) across the row
    for r, lst in by_row.items():
        if len(lst) < 2:
            continue
        t = np.array([x[0] for x in lst])
        rt = _tangent(np.array([ng[r]] * len(lst), dtype=object), t)
        steep = [abs(float(x[3] @ rt[i])) < cos_max for i, x in enumerate(lst)]
        L_r = ng[r].length
        inside = [P.on_m + 0.05 < x[0] < L_r - P.on_m - 0.05 for x in lst]
        for i, (ti, ki, wi, di, li) in enumerate(lst):
            if wi != 1 or not steep[i] or not inside[i] or (ki, wi) in out:
                continue  # i: a line ending on the row, steeply, inside it
            best = None
            for j, (tj, kj, wj, dj, lj) in enumerate(lst):
                if (
                    wj == 0
                    and kj != ki
                    and (kj, wj) not in out
                    and steep[j]
                    and inside[j]
                    and abs(ti - tj) <= 2 * P.on_m + 1.0
                    and not (li or lj)
                    and float(di @ dj) > cos_max
                    and (best is None or abs(ti - tj) < best[0])
                ):
                    best = (abs(ti - tj), kj, wj)
            if best is not None and not any(
                _along(lines[k].geom, ng[r], P.near_m + P.slack_m, cos_max)
                >= 0.5
                for k in (ki, best[1])
            ):
                out[(ki, wi)] = (best[1], best[2])
                out[(best[1], best[2])] = (ki, wi)
    return out


def _fits(rl, rep):
    """Can rep go on a row that already has the replacements rl? Only beside
    them, or overlapping one by less than OVERLAP_M where both keep the
    row's point there (the two then become one)."""
    for x in rl:
        ov = min(rep.hi, x.hi) - max(rep.lo, x.lo)
        if ov <= 1e-6:
            continue
        if ov > OVERLAP_M:
            return False
        if x.lo <= rep.lo and x.hi < rep.hi:
            if x.p_hi is None or rep.p_lo is None:
                return False
        elif rep.lo <= x.lo and rep.hi < x.hi:
            if rep.p_hi is None or x.p_lo is None:
                return False
        else:
            return False  # one inside the other
    return True


def _merge_reps(rl, vg, v_tree):
    """Replacements on one row that meet or overlap by less than OVERLAP_M
    become one, the lines joined along Vann-Nett where it runs between them
    -- no notch back to the old line."""
    out = [rl[0]]
    for x in rl[1:]:
        cur = out[-1]
        if (
            x.lo - cur.hi <= OVERLAP_M
            and cur.p_hi is not None
            and x.p_lo is not None
        ):
            b1, a2 = cur.body[-1], x.body[0]
            gap = float(np.hypot(*(b1 - a2)))
            mid = _vn_path(vg, v_tree, b1, a2)
            if (
                mid is not None
                and len(mid) > 1
                and float(np.hypot(*np.diff(mid, axis=0).T).sum())
                > gap + 2 * JOG_M
            ):
                mid = None  # Vann-Nett goes round: not its way between them
            if mid is None and gap > JOG_M:
                out.append(x)  # no Vann-Nett line to follow: kept apart
                continue
            out[-1] = _Rep(
                cur.lo,
                max(cur.hi, x.hi),
                _join([cur.body, mid, x.body]),
                cur.line,
                cur.how,
                cur.p_lo,
                x.p_hi,
            )
        else:
            out.append(x)
    return out


def _vn_path(vg, v_tree, a, b):
    """Vann-Nett's own line from point a to point b: along one part that
    carries both, or two parts meeting at a node (the shortest such way);
    else None."""
    p, q = Point(a), Point(b)
    at_a = [int(j) for j in v_tree.query(p, predicate="dwithin", distance=1e-6)]
    best = None
    for j in at_a:
        if vg[j].distance(q) <= 1e-6:
            c = _sub(vg[j], vg[j].project(p), vg[j].project(q))
            n = float(np.hypot(*np.diff(c, axis=0).T).sum())
            if best is None or n < best[0]:
                best = (n, c)
    if best is not None:
        return best[1]
    at_b = [int(j) for j in v_tree.query(q, predicate="dwithin", distance=1e-6)]
    for i in at_a:
        for j in at_b:
            if i == j:
                continue
            # the node between: an end of either part lying on the other
            # (a part may end on another's interior), so a->b and b->a agree
            for ei in (
                shapely.get_point(vg[i], 0),
                shapely.get_point(vg[i], -1),
                shapely.get_point(vg[j], 0),
                shapely.get_point(vg[j], -1),
            ):
                if vg[j].distance(ei) <= 1e-6 and vg[i].distance(ei) <= 1e-6:
                    c = _join(
                        [
                            _sub(vg[i], vg[i].project(p), vg[i].project(ei)),
                            _sub(vg[j], vg[j].project(ei), vg[j].project(q)),
                        ]
                    )
                    n = float(np.hypot(*np.diff(c, axis=0).T).sum())
                    if best is None or n < best[0]:
                        best = (n, c)
    return None if best is None else best[1]


def _loop_halves(g, on_m):
    """A line leaving a row and coming back to the same point, as (coords,
    anchored end): out to its farthest point and back. A line that comes
    back the way it went is drawn once, flowing back in."""
    c = shapely.get_coordinates(g)
    mid = (c[0] + c[-1]) / 2
    tip = int(np.argmax(np.hypot(*(c - mid).T)))
    if tip == 0 or tip == len(c) - 1:
        return [(c, 0)]
    out_, back = c[: tip + 1], c[tip:]
    # samples 0.1 m apart within on_m - 0.1: all of it within on_m (a loop
    # this narrow, or this small, drawn twice would be a ring)
    if _share(LineString(out_), LineString(back), on_m - 0.1, step=0.1) == 1.0:
        return [(back, 1)]
    return [(out_, 0), (back, 1)]


def _is_loop(ln):
    """A line leaving the network and coming back to (within 1 m of) the
    point it left, longer than 2 m (shorter ones are sub-metre detours)."""
    e0, e1 = ln.ends
    return (
        e0[0] == e1[0] == "anchor"
        and np.hypot(*(np.asarray(e0[3]) - np.asarray(e1[3]))) < 1.0
        and ln.geom.length > 2.0
    )


def _covered(lines, ctx, new_geom):
    """Lines no row took that the pass's own lines already carry: every
    point within on_m (samples 0.1 m apart within on_m - 0.1) of the rows as
    they now run or of a longer line to be added (of a loop, the half that
    is drawn). Never a line other lines meet at a junction or one a staying
    row relies on: those would be left unjoined. {line}"""
    P, ng = ctx.P, ctx.ng
    todo = sorted(
        (i for i in range(len(lines)) if i not in ctx.used),
        key=lambda i: -lines[i].geom.length,
    )
    if not todo:
        return set()
    lim = P.on_m - 0.1
    rows = [
        r for r in new_geom if r not in ctx.dropped and len(new_geom[r]) > 1
    ]
    redrawn = [LineString(new_geom[r]) for r in rows]
    r_tree = shapely.STRtree(redrawn) if redrawn else None
    kept_g, kept_b = [], []  # lines to be added, and their bounds

    def near(g):
        out = [
            LineString(new_geom[j]) if j in new_geom else ng[j]
            for j in (
                int(x)
                for x in ctx.n_tree.query(
                    g, predicate="dwithin", distance=P.on_m
                )
            )
            if j not in ctx.dropped
        ]
        if r_tree is not None:
            out += [
                redrawn[int(x)]
                for x in r_tree.query(g, predicate="dwithin", distance=P.on_m)
            ]
        if kept_b:
            b = np.asarray(kept_b)
            x0, y0, x1, y1 = g.bounds
            hit = np.flatnonzero(
                (b[:, 0] <= x1 + P.on_m)
                & (b[:, 2] >= x0 - P.on_m)
                & (b[:, 1] <= y1 + P.on_m)
                & (b[:, 3] >= y0 - P.on_m)
            )
            out += [kept_g[h] for h in hit]
        return out

    def within(g, carriers, step):
        n = max(2, int(np.ceil(g.length / step)))
        pts = shapely.line_interpolate_point(g, np.linspace(0, g.length, n + 1))
        d = np.full(len(pts), np.inf)
        for h in carriers:
            d = np.minimum(d, shapely.distance(pts, h))
            if d.max() <= lim:
                return True
        return False

    out = set()
    for i in todo:
        ln = lines[i]
        g = ln.geom
        may = i not in ctx.carry and all(e[0] != "inner" for e in ln.ends)
        carriers = near(g) if may else []
        if carriers and within(g, carriers, 1.0) and within(g, carriers, 0.1):
            out.add(i)
            continue
        drawn = (
            [LineString(h) for h, _w in _loop_halves(g, P.on_m)]
            if _is_loop(ln)
            else [g]
        )
        for h in drawn:
            kept_g.append(h)
            kept_b.append(h.bounds)
    return out


def _dangling(additions, final, loose, vg, v_tree):
    """Addition ends that meet nothing -- no other addition, no line -- and
    are no free end of Vann-Nett's (the tip of a loop drawn once, say), not
    yet listed: [("tip", point)]."""
    tree = shapely.STRtree(np.array(final, dtype=object))
    idx = [k for k, (c, _) in enumerate(additions) if len(c) >= 2]
    if not idx:
        return []
    a_geoms = np.array([LineString(additions[k][0]) for k in idx], dtype=object)
    a_tree = shapely.STRtree(a_geoms)
    ends = defaultdict(int)
    for k in idx:
        c = additions[k][0]
        ends[_pkey(c[0])] += 1
        ends[_pkey(c[-1])] += 1
    listed = {_pkey(xy) for _, xy in loose}
    out = []
    for n_k, k in enumerate(idx):
        c = additions[k][0]
        for xy in (c[0], c[-1]):
            key = _pkey(xy)
            if ends[key] > 1 or key in listed:
                continue
            p = Point(xy)
            if len(tree.query(p, predicate="dwithin", distance=1e-3)):
                continue  # on a row (or within a millimetre of one)
            if any(
                int(j) != n_k
                for j in a_tree.query(p, predicate="dwithin", distance=1e-3)
            ):
                continue  # on another addition
            hits = v_tree.query(p, predicate="dwithin", distance=1e-6)
            n_end, through = 0, False
            for j in hits:
                m = vg[j].project(p)
                if m <= 1e-6 or m >= vg[j].length - 1e-6:
                    n_end += 1
                else:
                    through = True
            if n_end == 1 and not through:
                continue  # Vann-Nett's own free end
            listed.add(key)
            out.append(("tip", np.asarray(xy, dtype=float)))
    return out


# ── duplicates found through the network ───────────────────────────────────
def _classify(pieces, ng, off_frac):
    """Per walked piece: True where it runs off Vann-Nett."""
    return [off_frac(_sub(ng[r], a, b)) >= 0.5 for r, a, b in pieces]


def _map_pieces(i, g, pieces, k_at, G, ng, how, keep_first, keep_last):
    """Replacements for walked pieces: piece j of the walk takes the line
    from k_at[j] to k_at[j + 1] (line measures in walking order). A piece
    ending inside its row -- or at the walk's own ends where the row keeps
    its node -- joins the row's own point there."""
    out = []
    last = len(pieces) - 1
    for j, (r, a, b) in enumerate(pieces):
        body = _sub(g, k_at[j], k_at[j + 1])
        own_a = 1e-6 < a < G.L[r] - 1e-6 or (j == 0 and keep_first)
        own_b = 1e-6 < b < G.L[r] - 1e-6 or (j == last and keep_last)
        pa = _xy(ng[r].interpolate(a)) if own_a else None
        pb = _xy(ng[r].interpolate(b)) if own_b else None
        if a > b:  # walked against the row's direction
            out.append((r, _Rep(b, a, body[::-1], i, how, pb, pa)))
        else:
            out.append((r, _Rep(a, b, body, i, how, pa, pb)))
    return out


def _walk_split(pieces, flags):
    """(first off, last off) if the walk is on* off+ on*, else None."""
    offs = [k for k, f in enumerate(flags) if f]
    if not offs or not all(flags[offs[0] : offs[-1] + 1]):
        return None
    return offs[0], offs[-1]


def _ends_of(pieces, G):
    """Nodes between consecutive walked pieces."""
    nodes = []
    for (r, _a, b), _ in zip(pieces[:-1], pieces[1:], strict=True):
        nodes.append(G.node(r, 0 if b <= 1e-6 else 1))
    return nodes


def _alike(g, dgeom, dup, P, downstream, cap):
    """ELVIS drawing the same stream elsewhere: running the same way, beside
    the line nearly all along, no more than a confluence's worth longer."""
    return (
        _dir_ok(dup, downstream)
        and dgeom.length >= 0.4 * g.length
        and dgeom.length <= 2.5 * g.length + 2 * P.near_m
        and shapely.hausdorff_distance(dgeom, g) <= cap
        and _share(g, dgeom, P.near_m) >= 0.7
        and _share(dgeom, g, P.near_m) >= 0.5
    )


def _crosses(seg, ctx, skip, line, at, own=()):
    """Does the connector seg cross a network line (but those in skip), a
    missing line (but its own chain's), or the line it joins (anywhere but
    at)?"""
    e0 = np.asarray(seg.coords[0])
    e1 = np.asarray(seg.coords[-1])

    def away(x):
        return (
            len(x)
            and (
                np.minimum(np.hypot(*(x - e0).T), np.hypot(*(x - e1).T)) > 1e-6
            ).any()
        )

    for j in ctx.n_tree.query(seg, predicate="intersects"):
        if int(j) in skip:
            continue
        if away(shapely.get_coordinates(seg.intersection(ctx.ng[j]))):
            return True
    if ctx.line_tree is not None:
        for k in ctx.line_tree.query(seg, predicate="intersects"):
            if int(k) in own:
                continue
            x = shapely.get_coordinates(seg.intersection(ctx.lines[k].geom))
            if away(x):
                return True
    x = shapely.get_coordinates(seg.intersection(line))
    return bool(len(x) and (np.hypot(*(x - at).T) > 1e-6).any())


def _duplicate_path(c, ctx):
    """Both ends anchor on different rows: ("ok", plan), ("kept", (ELVIS
    geometry, why)) when a duplicate was found but cannot be replaced, or
    ("reject", why) when none was found.

    Where nothing joins the anchors that way and one of them lies at the
    loose end of a line an earlier patch added (Vann-Nett it put in beside
    an ELVIS line it could not redraw), the line goes on along that
    addition to where it meets ELVIS: the ELVIS line beside both then takes
    the whole, and the addition, now drawn twice, goes. An addition ending
    on another row's interior (a T the old network left without a node)
    meets ELVIS there, on that row; where that finds nothing, the addition
    is followed on as if the row were not there."""
    res = _path_core(c, ctx)
    if res[0] != "reject":
        return res
    kept = None
    for side in (1, 0):
        for stop_inside in (True, False):
            ext = _over_additions(c, ctx, side, stop_inside)
            if ext is None:
                continue
            c2, run, n_far, tees = ext
            tcuts = [("tcut", -1, 0, xy, row) for row, xy in tees]
            if (
                c2.ends[0][0] == c2.ends[1][0] == "anchor"
                and c2.ends[0][1] == c2.ends[1][1]
            ):
                plan = _detour_over(c2, run, ctx)
                if plan is not None:
                    plan["attach"] = plan["attach"] + tcuts
                    plan["locked"] = plan["locked"] | {r for r, _ in tees}
                    return "ok", plan
            res2 = _path_core(
                c2,
                ctx,
                start_node=n_far if side == 0 and n_far is not None else None,
                end_node=n_far if side == 1 and n_far is not None else None,
                exclude=run,
            )
            if res2[0] == "ok":
                res2[1]["drop"] = set(run)
                res2[1]["locked"] = (
                    res2[1]["locked"] | set(run) | {r for r, _ in tees}
                )
                res2[1]["leave"] += [
                    (ctx.G.node(r, w), r, w) for r in run for w in (0, 1)
                ]
                res2[1]["attach"] = res2[1]["attach"] + tcuts
                return res2
            if res2[0] == "kept" and kept is None:
                kept = res2
    return kept if kept is not None else res


def _over_additions(c, ctx, side, stop_inside=True):
    """The chain c carried on, at its end side (1: B, 0: A), along rows an
    earlier patch added that run on from a loose end where c anchors, end to
    end through plain nodes, to the first node where anything else meets
    them -- or, with stop_inside, to a node lying on another row's interior,
    which is then the far anchor: (the longer chain, those rows, that node
    or None, the Ts to cut) or None. A T is a row interior the run passes
    through (the line taking it crosses the row there) or ends at with
    other rows ending there too: (that row, the exact point)."""
    G, P, ng = ctx.G, ctx.P, ctx.ng
    _, q, t, _xy_ = c.ends[side]
    if q not in ctx.added or q in ctx.locked or ctx.reps[q]:
        return None
    if t <= P.on_m + 0.05:
        w_near = 0
    elif G.L[q] - t <= P.on_m + 0.05:
        w_near = 1
    else:
        return None
    if G.deg(G.node(q, w_near)) != 1:
        return None  # not a loose end: the path search had it
    run, parts = [q], [_sub(ng[q], t, G.L[q] if w_near == 0 else 0.0)]
    w_far = 1 - w_near
    inside = None
    tees = []  # (row, point): a T of the old network to be made a node
    while True:
        n = G.node(run[-1], w_far)
        on_row = _on_interior(ctx, n, run)
        if stop_inside and on_row is not None:
            inside = on_row  # it ends on another row's interior (no node)
            others = [
                (r_o, w_o)
                for r_o, w_o in sorted(G.at[n])
                if r_o not in run and r_o not in ctx.dropped
            ]
            if others:  # rows ending there too: that row is cut for them
                r_o, w_o = others[0]
                tees.append((on_row[0], _end_xy(ng[r_o], w_o)))
            break
        rest = [x for x in G.at[n] if x != (run[-1], w_far)]
        if len(rest) != 1:
            break
        r2, w2 = rest[0]
        if r2 not in ctx.added or r2 in run or r2 in ctx.locked or ctx.reps[r2]:
            break
        if on_row is not None:
            # the run goes on through a T: the line that takes it will
            # cross that row there, so the row is cut there (a node)
            tees.append((on_row[0], _end_xy(ng[run[-1]], w_far)))
        run.append(r2)
        parts.append(
            _sub(ng[r2], 0.0, G.L[r2])
            if w2 == 0
            else _sub(ng[r2], G.L[r2], 0.0)
        )
        w_far = 1 - w2
    n_far = G.node(run[-1], w_far)
    if inside is None:
        # the run stops at a node on a row's interior where other rows end
        # too (no row carries it on): a T, which the row is cut at
        on_row = _on_interior(ctx, n_far, run)
        others = [
            (r_o, w_o)
            for r_o, w_o in sorted(G.at[n_far])
            if r_o not in run and r_o not in ctx.dropped
        ]
        if on_row is not None and others:
            r_o, w_o = others[0]
            tees.append((on_row[0], _end_xy(ng[r_o], w_o)))
    t_far = 0.0 if w_far == 0 else float(G.L[run[-1]])
    run_c = _join(parts)
    far = ("anchor", run[-1], t_far, run_c[-1])
    ks = list(c.ks)
    if inside is not None:  # the far anchor: the row the run ends on
        far = ("anchor", inside[0], inside[1], run_c[-1])
        n_far = None
    elif G.deg(n_far) < 2:
        # loose at the far end too (a fragment): go on along the missing
        # line that starts there, to where that one anchors
        nxt = None
        for k, ln in enumerate(ctx.lines):
            if k in ctx.used or k in ks:
                continue
            for w, e in enumerate(ln.ends):
                if (
                    e[0] == "anchor"
                    and e[1] == run[-1]
                    and abs(e[2] - t_far) <= P.on_m + 0.05
                    and ln.ends[1 - w][0] == "anchor"
                    and ln.ends[1 - w][1] not in run
                    and w == (0 if side == 1 else 1)
                ):
                    nxt = (k, w)
                    break
            if nxt:
                break
        if nxt is None:
            return None  # nothing to join
        k, w = nxt
        lk = shapely.get_coordinates(ctx.lines[k].geom)
        run_c = _join([run_c, lk if w == 0 else lk[::-1]])
        far = ctx.lines[k].ends[1 - w]
        ks.append(k)
        n_far = None
    lc = shapely.get_coordinates(c.geom)
    if side == 1:
        geom = LineString(_join([lc, run_c]))
        ends = [c.ends[0], far]
        joints = list(c.joints)
    else:
        geom = LineString(_join([run_c[::-1], lc]))
        ends = [far, c.ends[1]]
        shift = geom.length - c.geom.length
        joints = [(m + shift, key) for m, key in c.joints]
    c2 = _Chain(ks, geom, ends, joints, c.vn_id, c.n_ids)
    return c2, run, n_far, tees


def _end_xy(g, w):
    """The exact first (w 0) or last (w 1) point of a line."""
    c = shapely.get_coordinates(g)
    return c[0 if w == 0 else -1].copy()


def _detour_over(c2, run, ctx):
    """The chain carried over earlier additions anchors twice on one row:
    that row takes the whole as a detour, the additions go. A plan, or
    None."""
    G, P, ng = ctx.G, ctx.P, ctx.ng
    r, t0, t1 = c2.ends[0][1], c2.ends[0][2], c2.ends[1][2]
    lo, hi = min(t0, t1), max(t0, t1)
    if hi - lo < 1.0 or r in ctx.locked:
        return None
    stretch = LineString(_sub(ng[r], lo, hi))
    if ctx.off_frac(stretch) < 0.5:
        return None
    dist = float(shapely.hausdorff_distance(stretch, c2.geom))
    if dist > max(P.near_m + P.slack_m, min(0.25 * c2.geom.length, 300.0)):
        return None  # the row is not beside what it would take (as walks)
    body = shapely.get_coordinates(c2.geom)
    rep_ = _Rep(
        lo,
        hi,
        body[::-1] if t0 > t1 else body,
        c2.ks[0],
        "detour",
        _xy(ng[r].interpolate(lo)),
        _xy(ng[r].interpolate(hi)),
    )
    if not _fits(ctx.reps[r], rep_):
        return None
    return {
        "dgeom": MultiLineString([stretch]),
        "reps": [(r, rep_)],
        "leave": [(G.node(x, w), x, w) for x in run for w in (0, 1)],
        "attach": [],
        "conns": [],
        "qa_conn": [],
        "moved": set(),
        "locked": set(run),
        "used": set(c2.ks),
        "drop": set(run),
        "dist": dist,
        "conn_len": 0.0,
        "extra_m": 0.0,
    }


def _on_interior(ctx, n, skip):
    """The row (not in skip) whose interior node n lies on, within 1 mm, and
    the measure there: (row, measure) or None."""
    p = Point(ctx.G.xy[n])
    best = None
    for j in ctx.n_tree.query(p, predicate="dwithin", distance=1e-3):
        j = int(j)
        if j in skip or j in ctx.dropped:
            continue
        m = float(ctx.ng[j].project(p))
        if SNAP_M < m < ctx.G.L[j] - SNAP_M:
            d = float(ctx.ng[j].distance(p))
            if best is None or d < best[0]:
                best = (d, j, m)
    return None if best is None else (best[1], best[2])


def _path_core(c, ctx, start_node=None, end_node=None, exclude=()):
    """The duplicate path search between the chain's two ends (or from or
    to a given node instead of an anchor)."""
    G, P = ctx.G, ctx.P
    (_, q0, t0, _A), (_, q1, t1, _B) = c.ends
    g = c.geom
    cutoff = 2.5 * g.length + 2 * P.near_m
    near = ctx.n_tree.query(
        g, predicate="dwithin", distance=P.near_m + P.node_r
    )
    start = start_node if start_node is not None else (q0, t0)
    skip = set(exclude) | ({q1} if end_node is None else set())
    best = G.paths_from(near, start, cutoff, exclude=skip)
    cands = []
    if end_node is not None:
        if end_node in best:
            cands.append((best[end_node][0], end_node))
    else:
        for n, off in ((G.node(q1, 0), t1), (G.node(q1, 1), G.L[q1] - t1)):
            if n in best:
                cands.append((best[n][0] + off, n))
    if not cands:
        return "reject", "no network path between the anchors"
    _d, n = min(cands)
    pieces = _Graph.trace(best, n)
    if end_node is None:
        pieces.append((q1, 0.0 if n == G.node(q1, 0) else float(G.L[q1]), t1))
    pieces = [p for p in pieces if abs(p[2] - p[1]) > 1e-9]
    if not pieces:
        return "reject", "no network path between the anchors"
    if len({p[0] for p in pieces}) < len(pieces):
        # both anchors on one row: out to a node and back along that row
        return "reject", "the network path runs back along a row"
    split = _walk_split(pieces, _classify(pieces, ctx.ng, ctx.off_frac))
    if split is None:
        return "reject", "the network path between is on Vann-Nett"
    f0, f1 = split
    if (start_node is not None and f0 != 0) or (
        end_node is not None and f1 != len(pieces) - 1
    ):
        return "reject", "no network path between the anchors"
    return _replace_walk(
        c,
        g,
        c.joints,
        pieces,
        f0,
        f1,
        ctx,
        "path",
        True,
        (
            None if start_node is not None else q0,
            None if end_node is not None else q1,
        ),
    )


def _duplicate_head(c, ctx):
    """One end anchors, the other is a free Vann-Nett end: the walk from the
    anchor ends at a free ELVIS end near it, or at the point of an ELVIS
    line running on past it that lies nearest it."""
    G, P, ng = ctx.G, ctx.P, ctx.ng
    w_a = 0 if c.ends[0][0] == "anchor" else 1
    _, q, t, _A = c.ends[w_a]
    F = np.asarray(c.ends[1 - w_a][1], dtype=float)
    if w_a == 0:
        g, joints = c.geom, c.joints
    else:  # walk from the anchor: against Vann-Nett's direction
        g = LineString(shapely.get_coordinates(c.geom)[::-1])
        joints = [(c.geom.length - m, key) for m, key in reversed(c.joints)]
    cutoff = 2.5 * g.length + 2 * P.near_m
    near = {
        int(x)
        for x in ctx.n_tree.query(
            g, predicate="dwithin", distance=P.near_m + P.node_r
        )
    }
    best = G.paths_from(near, (q, t), cutoff)
    targets = []
    for n in best:
        if G.deg(n) == 1 and np.hypot(*(G.xy[n] - F)) <= P.near_m + P.slack_m:
            pieces = [
                p for p in _Graph.trace(best, n) if abs(p[2] - p[1]) > 1e-9
            ]
            if pieces:
                targets.append(("free", pieces))
    Fp = Point(F)
    Lg = g.length
    into = _xy(g.interpolate(max(0.0, Lg - min(10.0, Lg)))) - F
    rows = {q} | {r for n in best for r, _ in G.at[n] if r in near}
    for r in sorted(rows):
        if ng[r].distance(Fp) > CONNECT_M:
            continue
        tP = float(ng[r].project(Fp))
        w_end = 0 if tP <= SNAP_M else 1 if tP >= G.L[r] - SNAP_M else None
        if w_end is not None and G.deg(G.node(r, w_end)) >= 3:
            continue  # the row's end at a junction: ELVIS does not run on
        if r == q:
            a0 = t
            pieces = [(q, t, tP)]
        else:
            opts = []
            for w in (0, 1):
                n = G.node(r, w)
                if n not in best:
                    continue
                walk = _Graph.trace(best, n)
                if any(p[0] == r for p in walk):
                    continue  # reached through the row itself
                d = best[n][0] + (tP if w == 0 else G.L[r] - tP)
                opts.append((d, w, walk))
            if not opts:
                continue
            _d, w, walk = min(opts, key=lambda o: o[0])
            a0 = 0.0 if w == 0 else float(G.L[r])
            pieces = [*walk, (r, a0, tP)]
        if sum(abs(b - a) for _, a, b in pieces) >= 0.4 * Lg:
            # ELVIS runs the line's way (not a spur off the row): the
            # connector may join it a little on
            pieces[-1] = (r, a0, _reach_on(ng[r], a0, tP, F, into, P))
        pieces = [p for p in pieces if abs(p[2] - p[1]) > 1e-9]
        if pieces:
            # the walk ends at a node ELVIS runs on through (one row in, one
            # out): taking it reroutes that stream, so it must be the line's
            r_l, _a_l, b_l = pieces[-1]
            d0_, d1_ = b_l, float(G.L[r_l]) - b_l
            w_l = (0 if d0_ <= d1_ else 1) if min(d0_, d1_) <= JOG_M else None
            through = w_l is not None and sorted(
                x for _, x in G.at[G.node(r_l, w_l)]
            ) == [0, 1]
            targets.append(("reach", pieces, through))
    results, modes, seen = [], [], set()
    for mode, pieces, *strict in targets:
        key = (mode, tuple((r, round(a, 3), round(b, 3)) for r, a, b in pieces))
        if key in seen:
            continue
        seen.add(key)
        split = _walk_split(pieces, _classify(pieces, ng, ctx.off_frac))
        if split is None or split[1] != len(pieces) - 1:
            continue  # the walk must be off Vann-Nett up to its end
        res = _replace_walk(
            c,
            g,
            joints,
            pieces,
            split[0],
            split[1],
            ctx,
            mode,
            w_a == 0,
            (q, None),
            strict=bool(strict and strict[0]),
        )
        results.append(res)
        modes.append(mode)
    if not results:
        return "reject", "no ELVIS line from the anchor towards the free end"
    frees = [
        r
        for r, m in zip(results, modes, strict=True)
        if r[0] == "ok" and m == "free"
    ]
    plain = [r for r in frees if not r[1]["capped"]]
    if plain:  # a free ELVIS end within reach moves to Vann-Nett's
        return min(plain, key=lambda r: r[1]["dist"])
    oks = [
        r
        for r, m in zip(results, modes, strict=True)
        if r[0] == "ok" and m == "reach"
    ]

    def cost(r):
        return r[1]["dist"] + r[1]["conn_len"] + 0.5 * r[1]["extra_m"]

    lim2 = 2.0 * max(JOG_M, 0.5 * c.geom.length)
    capped = [
        r
        for r in frees
        if r[1]["dir_ok"]
        and r[1]["extra_m"] <= lim2
        and r[1]["run_on"] <= HOOK_M
    ]
    if oks:
        best = min(oks, key=cost)
        if best[1]["hook"] and capped:
            return min(capped, key=lambda r: r[1]["dist"])
        return best
    runback = any(
        m == "reach" and r == ("reject", WHY_BACK)
        for r, m in zip(results, modes, strict=True)
    )
    if runback and capped:  # K2: the most hooked connector of all
        return min(capped, key=lambda r: r[1]["dist"])
    # nothing reaches ELVIS's upper course: a walk past it is refused
    results = [
        ("reject", "the free ELVIS end lies far beyond the line's end")
        if r[0] == "ok" and m == "free"
        else r
        for r, m in zip(results, modes, strict=True)
    ]
    kept = [r for r in results if r[0] == "kept"]
    return kept[0] if kept else results[0]


def _reach_on(row, a, tP, B, into, P):
    """Where a head connector from Vann-Nett's free end B joins row (walked
    from measure a): tP, its point nearest B -- or, where a connector there
    would run back along the line (into: from B back into the line's last
    10 m), the first point on along the row, every metre and no farther
    than CONNECT_M from B, where it does not. tP when there is none."""
    n_into = float(np.hypot(*into))
    if n_into < 1e-9:
        return tP
    cos_ma = np.cos(np.radians(P.max_angle))

    def back(s):
        p = _xy(row.interpolate(s))
        d = float(np.hypot(*(p - B)))
        return d > JOG_M and float((p - B) @ into) >= cos_ma * d * n_into, d

    if not back(tP)[0]:
        return tP
    sgn = 1.0 if tP >= a else -1.0
    s = tP
    while True:
        s += sgn
        if not 0.0 < s < row.length:
            return tP
        bk, d = back(s)
        if d > CONNECT_M:
            return tP
        if not bk:
            return s


def _junction_line(ctx, joints, r_x, w_x, dup_rows):
    """A Vann-Nett line, not yet taken, from row r_x into a junction of the
    chain: (measure of the junction along the chain, the row's replacement,
    the lines) or None. Row r_x's end w_x is at an ELVIS node the chain's
    duplicate runs through; the row takes the line from where it anchors on
    it, the node going to the junction. On its way the line may cross rows
    of that duplicate (which it splits into a line on either side)."""
    ng, G = ctx.ng, ctx.G
    for m, key in joints:
        for lt, wt in ctx.inner_at.get(key, ()):
            if lt in ctx.used:
                continue
            seq = [(lt, wt)]  # lines from the junction back towards r_x
            e = ctx.lines[lt].ends[1 - wt]
            while (
                e[0] == "anchor"
                and e[1] != r_x
                and e[1] in dup_rows
                and (seq[-1][0], 1 - seq[-1][1]) in ctx.cross
            ):
                k2, w2 = ctx.cross[(seq[-1][0], 1 - seq[-1][1])]
                if k2 in ctx.used or any(k2 == k for k, _ in seq):
                    break
                seq.append((k2, w2))
                e = ctx.lines[k2].ends[1 - w2]
            if e[0] != "anchor" or e[1] != r_x:
                continue
            t_a = e[2]
            lo, hi = (t_a, float(G.L[r_x])) if w_x == 1 else (0.0, t_a)
            if hi - lo < 1e-6 or r_x in ctx.locked:
                continue
            if any(
                min(hi, x.hi) - max(lo, x.lo) > -1e-6 for x in ctx.reps[r_x]
            ):
                continue
            tail = _sub(ng[r_x], lo, hi)
            if LineString(tail).length >= 2.0 and ctx.off_frac(tail) < 0.5:
                continue
            # from the anchor on r_x to the junction, along Vann-Nett
            runs = []
            for k, w in reversed(seq):  # w: the end towards the junction
                lc = shapely.get_coordinates(ctx.lines[k].geom)
                runs.append(lc if w == 1 else lc[::-1])
            parts = [runs[0]]
            for a_run, b_run in zip(runs[:-1], runs[1:], strict=True):
                mid = _vn_path(ctx.vg, ctx.v_tree, a_run[-1], b_run[0])
                parts += [mid, b_run]
            body = _join(parts)
            p_a = _xy(ng[r_x].interpolate(t_a))
            if w_x == 1:  # the row flows into the node: anchor -> junction
                rep = _Rep(lo, hi, body, lt, "junction", p_a, None)
            else:  # the row flows out of it: junction -> anchor
                rep = _Rep(lo, hi, body[::-1], lt, "junction", None, p_a)
            return m, rep, {k for k, _ in seq}
    return None


def _passed_over(ctx, g, r_x, w_x):
    """The line g passes over row r_x (within on_m), whose end w_x is at a
    node the duplicate moves: the row ends there instead, its tail beyond
    -- off Vann-Nett -- dropped. (measure along g, the row's replacement,
    no lines) or None."""
    P, ng, G = ctx.P, ctx.ng, ctx.G
    row = ng[r_x]
    if row.distance(g) > P.on_m + 0.05:
        return None
    pa, pb = nearest_points(row, g)
    t_x = float(row.project(pa))
    k = float(g.project(pb))
    lo, hi = (t_x, float(G.L[r_x])) if w_x == 1 else (0.0, t_x)
    if hi - lo < 1.0 or hi - lo > CONNECT_M or r_x in ctx.locked:
        return None  # (a tail of at most CONNECT_M goes)
    if any(min(hi, x.hi) - max(lo, x.lo) > -1e-6 for x in ctx.reps[r_x]):
        return None
    if ctx.off_frac(_sub(row, lo, hi)) < 0.5:
        return None  # the tail carries Vann-Nett: not this
    p_a = _xy(row.interpolate(t_x))
    end = _xy(g.interpolate(k)).reshape(1, 2)
    if w_x == 1:
        rep_ = _Rep(lo, hi, end, -1, "junction", p_a, None)
    else:
        rep_ = _Rep(lo, hi, end, -1, "junction", None, p_a)
    return k, rep_, set()


def _replace_walk(
    c, g, joints, pieces, f0, f1, ctx, mode, downstream, anchors, strict=False
):
    """Replace the off-Vann-Nett stretch pieces[f0:f1 + 1] of a walk from the
    line's first end A to its last end B. mode: "path" (B anchors), "free"
    (the walk ends at a free ELVIS end, which goes to B) or "reach" (the walk
    ends at the ELVIS point nearest B, which a connector joins to B). strict:
    the walk must be close to the line (ELVIS drawing the same stream
    elsewhere, _alike, is not enough)."""
    G, P, ng = ctx.G, ctx.P, ctx.ng
    dup = pieces[f0 : f1 + 1]
    Lk = g.length
    A = _xy(shapely.get_point(g, 0))
    B = _xy(shapely.get_point(g, -1))
    dgeom = MultiLineString([LineString(_sub(ng[r], a, b)) for r, a, b in dup])
    dist = float(shapely.hausdorff_distance(dgeom, g))
    dlen = float(dgeom.length)
    tol = max(P.near_m + P.slack_m, min(0.25 * Lk, 300.0))
    ratio_ok = 0.4 <= dlen / max(Lk, 1e-9) <= 2.5 or (
        mode != "reach"
        and Lk < P.near_m
        and abs(dlen - Lk) <= P.near_m
        and _dir_ok(dup, downstream)
    )
    close = dist <= tol and ratio_ok
    if mode == "reach":
        close = close and _dir_ok(dup, downstream)
    cap = 300.0
    if not close:
        if not _alike(g, dgeom, dup, P, downstream, cap):
            return (
                "reject",
                f"not the same stream: Hausdorff {dist:.0f} m (tolerance "
                f"{tol:.0f} m), {dlen:.0f} m of ELVIS against {Lk:.0f} m",
            )
        if strict:
            return (
                "reject",
                f"ELVIS draws the stream elsewhere ({dlen:.0f} m of ELVIS "
                f"against {Lk:.0f} m): taking it would reroute the stream "
                "ELVIS runs on through a node",
            )
    why = None
    leave, attach, conns, qa_conn = [], [], [], []
    extra_reps, extra_used = [], set()
    conn_len = 0.0
    capped = hook = False
    run_on = 0.0
    nodes = _ends_of(pieces, G)
    last = len(pieces) - 1
    dup_rows = {r for r, _, _ in dup}
    # the duplicate's first end: inside the anchor row, or at a node where it
    # leaves Vann-Nett line kept as it is
    keep_first = keep_last = False
    if f0 > 0:
        n = nodes[f0 - 1]
        r, a, _ = dup[0]
        w = 0 if a <= 1e-6 else 1
        if np.hypot(*(A - G.xy[n])) <= JOG_M:
            keep_first = True  # the row keeps its node, a jog reaches A
        elif n in ctx.moved or not G.can_leave(n, r, w):
            why = WHY_NODE
        else:
            leave.append((n, r, w))
            attach.append(("row", r, w, A, anchors[0]))
    else:
        keep_first = True
    if f1 < last:
        n = nodes[f1]
        r, _, b = dup[-1]
        w = 0 if b <= 1e-6 else 1
        if np.hypot(*(B - G.xy[n])) <= JOG_M:
            keep_last = True
        elif n in ctx.moved or not G.can_leave(n, r, w):
            why = why or WHY_NODE
        else:
            leave.append((n, r, w))
            attach.append(("row", r, w, B, anchors[1]))
    elif mode == "free":
        # a walk longer than the line by more than lim: ELVIS runs on past
        # Vann-Nett's end, its upper course is kept -- a connector reaches
        # it, or the line is added (but see _duplicate_head on hooks) --
        # unless it neither runs on past B (beyond its point nearest B) nor
        # ends farther than lim from B: then the excess is meander
        lim = max(JOG_M, 0.5 * Lk)
        walk = LineString(_join([_sub(ng[r_], a_, b_) for r_, a_, b_ in dup]))
        run_on = float(walk.length - walk.project(Point(B)))
        r, _, b = dup[-1]
        gap = float(np.hypot(*(B - _xy(ng[r].interpolate(b)))))
        capped = dlen - Lk > lim and not (run_on <= lim and gap <= lim)
        w = 0 if b <= 1e-6 else 1
        n = G.node(r, w)
        if n in ctx.moved:
            why = why or "a junction already moved"
        leave.append((n, r, w))
    elif mode == "reach":
        keep_last = True
        r, _, b = dup[-1]
        p_end = _xy(ng[r].interpolate(b))
        d = float(np.hypot(*(B - p_end)))
        if d > JOG_M:
            # a connector running back along the line draws it out and back
            into = _xy(g.interpolate(max(0.0, Lk - min(10.0, Lk)))) - B
            n_into = float(np.hypot(*into))
            if (
                n_into > 1e-9
                and float((p_end - B) @ into)
                >= np.cos(np.radians(P.max_angle)) * d * n_into
            ):
                return "reject", WHY_BACK
            # it turns more than a right angle off the line: a hook
            hook = n_into > 1e-9 and float((p_end - B) @ into) > 0.0
            seg = LineString([B, p_end])
            if d > CONNECT_M:
                why = why or f"a connector of {d:.0f} m"
            elif _crosses(seg, ctx, dup_rows, g, B, own=set(c.ks)):
                why = why or "the connector would cross a line"
            conn_len += d
            qa_conn.append((seg, r, d, "head"))
    else:
        keep_last = True
    # nodes inside the duplicate move onto the line, in walking order, with
    # every row at them: each to its nearest point on the line (or to the
    # junction Vann-Nett draws its other row into), else -- where those run
    # out of order, as where ELVIS runs on past the line's end to a
    # confluence elsewhere -- in proportion along the duplicate
    inner_nodes = nodes[f0:f1]
    extras, hits = [], []
    for idx, n in enumerate(inner_nodes):
        (r_in, a_in, b_in), (r_out, a_out, _b) = dup[idx], dup[idx + 1]
        own = {
            (r_in, 1 if b_in > a_in else 0),
            (r_out, 0 if a_out <= 1e-6 else 1),
        }
        extra = sorted(x for x in G.at[n] if x not in own)
        if extra and n in ctx.moved:
            why = why or "a junction inside already moved"
        kj, taken = None, {}
        if extra and joints:
            for r_x, w_x in extra:
                hit = _junction_line(ctx, joints, r_x, w_x, dup_rows)
                if hit is not None and (kj is None or abs(hit[0] - kj) < 1e-6):
                    kj = hit[0]
                    taken[(r_x, w_x)] = hit
        for r_x, w_x in extra:
            if (r_x, w_x) in taken or kj is not None:
                continue
            hit = _passed_over(ctx, g, r_x, w_x)
            if hit is not None:
                hit[1].line = c.ks[0]
                kj = hit[0]
                taken[(r_x, w_x)] = hit
        extras.append(extra)
        hits.append((kj, taken))
    run = np.cumsum([abs(b - a) for _r, a, b in dup])
    ks = None
    for way in ("nearest", "along"):
        cand = []
        for idx, n in enumerate(inner_nodes):
            kj = hits[idx][0]
            if kj is not None:
                cand.append(kj)
            elif way == "nearest":
                cand.append(float(g.project(Point(G.xy[n]))))
            else:
                cand.append(float(Lk * run[idx] / run[-1]))
        seq = [0.0, *cand, Lk]
        if all(b > a + 1e-6 for a, b in zip(seq[:-1], seq[1:], strict=True)):
            ks = cand
            break
    if ks is None:
        why = why or "nodes out of order along the line"
        ks = []
    k_at = [0.0]
    stays = set()  # (node, row, end): the row's stream goes on on Vann-Nett
    carry = set()  # the lines those rows rely on: always added
    relink = []  # (line, end, row, row end): that line now joins the row
    for idx, k in enumerate(ks):
        if why:
            break
        n = inner_nodes[idx]
        extra, (_kj, taken) = extras[idx], hits[idx]
        k_at.append(k)
        x = _xy(g.interpolate(k))
        at_n = {rr for rr, _ in G.at[n]}
        # Vann-Nett's own line for a row's stream, leaving this node the way
        # the row flows on: the row keeps its end, the line joins it there
        goes_on, near_lines = {}, []
        for j in ctx.anchor_tree.query(
            Point(G.xy[n]), predicate="dwithin", distance=P.on_m + 0.05
        ):
            kk, ww = ctx.anchor_ends[int(j)]
            if kk not in ctx.used and kk not in c.ks:
                goes_on[(ww, ctx.lines[kk].ends[ww][1])] = kk
                near_lines.append((ww, kk))
        for (r_x, w_x), (_m, rep, lts) in taken.items():
            body = rep.body.copy()
            body[-1 if w_x == 1 else 0] = x
            rep.body = body
            extra_reps.append((r_x, rep))
            extra_used |= lts
        for r_x, w_x in extra:
            if (r_x, w_x) in taken:
                continue
            if r_x in dup_rows:
                why = "a row at a junction inside is on the duplicate too"
                break
            # flowing in (w 1), a line starting here (end 0) carries it on;
            # flowing out (w 0), a line ending here (end 1) feeds it -- one
            # anchored on the row, or one at the node running the row's way
            # (anchored on whichever row lay nearest), which then joins it
            kk = goes_on.get((1 - w_x, r_x))
            if kk is None:
                kk = _continues(ctx, r_x, w_x, near_lines)
                if kk is not None:
                    relink.append((kk, 1 - w_x, r_x, w_x))
            if kk is not None:
                stays.add((n, r_x, w_x))
                carry.add(kk)
                continue
            d = float(np.hypot(*(G.xy[n] - x)))
            if d > CONNECT_M:
                why = f"a junction inside lies {d:.0f} m from the line"
                break
            seg = LineString([G.xy[n], x])
            if d > SNAP_M and _crosses(
                seg, ctx, at_n | dup_rows, g, x, own=set(c.ks)
            ):
                why = "a connector would cross a line"
                break
            conns.append((r_x, w_x, x, G.xy[n].copy()))
            if d > JOG_M:
                qa_conn.append((seg, r_x, d, "junction"))
                conn_len += d
    k_at.append(Lk)

    def kept(why):
        if capped:  # only a plan in hand may stand in for a hooked reach
            return "reject", "the free ELVIS end lies far beyond the line's end"
        return "kept", (dgeom, why)

    if why:
        return kept(why)
    if len(k_at) != len(dup) + 1:
        return kept("nodes out of order along the line")
    how = "path" if mode == "path" else "head"
    plan = _map_pieces(c.ks[0], g, dup, k_at, G, ng, how, keep_first, keep_last)
    for r, rep in plan:
        if r in ctx.locked or not _fits(ctx.reps[r], rep):
            return kept("row already changed")
    inner = set(nodes[f0:f1])
    for n in inner:
        for r, w in list(G.at[n]):
            if (n, r, w) not in stays:
                leave.append((n, r, w))
    return "ok", {
        "dgeom": dgeom,
        "reps": plan + extra_reps,
        "leave": leave,
        "attach": attach,
        "conns": conns,
        "qa_conn": qa_conn,
        "moved": inner,
        "locked": {r for r, *_ in conns},
        "used": set(c.ks) | extra_used,
        "dist": dist,
        "conn_len": conn_len,
        "extra_m": max(0.0, dlen - Lk),
        "relink": relink,
        "carry": carry,
        "capped": capped,
        "run_on": run_on,
        "hook": hook,
        "dir_ok": _dir_ok(dup, downstream),
    }


def _continues(ctx, r_x, w_x, near_lines):
    """A missing line at row r_x's end w_x that carries the row's stream on
    (r_x flowing in: a line starting there; flowing out: a line ending
    there), running the row's way within max_angle over their first (last)
    5 m: its index, or None."""
    row, L = ctx.ng[r_x], float(ctx.G.L[r_x])
    if w_x == 1:
        a, b = row.interpolate(max(0.0, L - 5.0)), row.interpolate(L)
    else:
        a, b = row.interpolate(0.0), row.interpolate(min(5.0, L))
    rd = _xy(b) - _xy(a)
    nr = float(np.hypot(*rd))
    if nr < 1e-9:
        return None
    cos_min = np.cos(np.radians(ctx.P.max_angle))
    for ww, kk in near_lines:
        if ww != 1 - w_x:
            continue
        g = ctx.lines[kk].geom
        Lg = g.length
        if ww == 0:
            p, q = g.interpolate(0.0), g.interpolate(min(5.0, Lg))
        else:
            p, q = g.interpolate(max(0.0, Lg - 5.0)), g.interpolate(Lg)
        d = _xy(q) - _xy(p)
        nd = float(np.hypot(*d))
        if nd > 1e-9 and float(rd @ d) >= cos_min * nr * nd:
            return kk
    return None


# ── joining, cutting, orienting ────────────────────────────────────────────
def _bridge(additions, attach, final, vg, P):
    """Additions whose anchors lay on a row that has since taken another line
    (a Vann-Nett tributary crossing an ELVIS line now redrawn) face each
    other across the stretch that row used to cover: they become one, joined
    by Vann-Nett's own line between them."""
    tree = shapely.STRtree(np.array(final, dtype=object))
    v_tree = shapely.STRtree(vg)
    orphan = []  # (addition, end, point) with no line left under the anchor
    for kind, k, w, xy, _ in attach:
        if kind != "add":
            continue
        if not len(tree.query(Point(xy), predicate="dwithin", distance=JOG_M)):
            orphan.append((k, w, np.asarray(xy, dtype=float)))
    if not orphan:
        return additions, attach, []
    ends_tree = shapely.STRtree(
        shapely.points(np.array([o[2] for o in orphan]))
    )
    parent = {}  # addition -> the one it went into

    def root(k):
        while k in parent:
            k = parent[k]
        return k

    done = set()
    for a_i, (k1, w1, p) in enumerate(orphan):
        if a_i in done or w1 != 1:
            continue
        for b_i in ends_tree.query(
            Point(p), predicate="dwithin", distance=BRIDGE_M
        ):
            k2, w2, q = orphan[b_i]
            if b_i in done or w2 != 0 or root(k1) == root(k2) or k2 in parent:
                continue
            if additions[root(k1)][1] == additions[k2][1]:
                continue  # the two halves of one loop
            on = [
                j
                for j in v_tree.query(
                    Point(p), predicate="dwithin", distance=1e-6
                )
                if vg[j].distance(Point(q)) <= 1e-6
            ]
            if not on:
                continue
            j = on[0]
            mid = _sub(vg[j], vg[j].project(Point(p)), vg[j].project(Point(q)))
            r1 = root(k1)
            additions[r1][0] = _join([additions[r1][0], mid, additions[k2][0]])
            parent[k2] = r1
            done.update((a_i, b_i))
            break
    # a loose end left over goes on along its own Vann-Nett line, back to
    # the network (within on_m) or to that line's end, at most BRIDGE_M
    moved = {}  # (addition, end) -> its new end, or None: a free end now
    freed = []
    for o_i, (k, w, p) in enumerate(orphan):
        if o_i in done:
            continue
        kk = root(k)
        c = additions[kk][0]
        if len(c) < 2:
            continue
        inside = Point(c[1] if w == 0 else c[-2])
        on = v_tree.query(Point(p), predicate="dwithin", distance=1e-6)
        if not len(on):
            continue
        line = vg[on[0]]
        m_p = float(line.project(Point(p)))
        sign = 1.0 if float(line.project(inside)) <= m_p else -1.0
        stop = None
        for s in np.arange(0.5, BRIDGE_M + 1e-9, 0.5):
            m = m_p + sign * s
            if m >= line.length or m <= 0.0:
                stop = line.length if sign > 0 else 0.0
                break
            q = line.interpolate(m)
            if len(tree.query(q, predicate="dwithin", distance=P.on_m)):
                stop = m
                break
        if stop is None:
            continue
        ext = _sub(line, m_p, stop)
        additions[kk][0] = _join([ext[::-1], c]) if w == 0 else _join([c, ext])
        end = ext[-1]
        near = tree.query(Point(end), predicate="dwithin", distance=JOG_M)
        moved[(k, w)] = end if len(near) else None
        if not len(near):
            freed.append(("freed", np.asarray(end, dtype=float)))
    if not parent and not moved:
        return additions, attach, freed
    keep = [k for k in range(len(additions)) if k not in parent]
    new_index = {k: n for n, k in enumerate(keep)}
    gone = {(orphan[n][0], orphan[n][1]) for n in done}
    out_attach = []
    for kind, k, w, xy, pref in attach:
        if kind in ("add", "inner"):
            if (k, w) in gone:
                continue
            if kind == "add" and (k, w) in moved:
                if moved[(k, w)] is None:
                    continue
                xy, pref = moved[(k, w)], None
            out_attach.append((kind, new_index[root(k)], w, xy, pref))
        else:
            out_attach.append((kind, k, w, xy, pref))
    return [additions[k] for k in keep], out_attach, freed


def _attach_all(attach, final, new_geom, additions, P):
    """Join every loose end to the line it lands on (within JOG_M of its
    Vann-Nett point, on the row it anchored on where that is still there),
    with a straight jog of at most JOG_M; that line is cut there. Returns
    ({row: [cut points]}, [(kind, point) left loose])."""
    tree = shapely.STRtree(np.array(final, dtype=object))
    splits: dict[int, list] = defaultdict(list)
    loose = []
    add_ends = defaultdict(list)  # rounded point -> additions ending there
    for k, (c, _) in enumerate(additions):
        for w in (0, 1):
            add_ends[tuple(np.round(c[0 if w == 0 else -1], 6))].append((k, w))

    def place(xy, skip=None, prefer=None, snap_ends=True):
        """The point on the final network nearest xy, snapped to a row end
        (a junction just beside one goes onto it rather than cutting a
        sliver off the row), a cut already made or a vertex: (row, point)
        or None."""
        p = Point(xy)
        cand = [
            int(j)
            for j in tree.query(p, predicate="dwithin", distance=JOG_M)
            if j != skip
        ]
        if not cand:
            return None
        if prefer is not None and prefer in cand:
            j = prefer
        else:
            j = min(cand, key=lambda j: final[j].distance(p))
        gj = final[j]
        x = _xy(gj.interpolate(gj.project(p)))
        vx = shapely.get_coordinates(gj)
        if snap_ends:
            for e in (vx[0], vx[-1]):
                if (
                    np.hypot(*(e - x)) <= P.on_m + 0.5
                    and np.hypot(*(e - xy)) <= JOG_M
                ):
                    return j, e
        for q in splits.get(j, ()):
            if np.hypot(*(q - x)) <= SNAP_M:
                return j, q
        dv = np.hypot(*(vx - x).T)
        if dv.min() <= SNAP_M:
            x = vx[int(np.argmin(dv))]
        return j, x

    def cut_all(j, x):
        """Cut row j at x, and every other row passing exactly there."""
        _cut_request(splits, final, j, x)
        for jj in tree.query(Point(x), predicate="dwithin", distance=1e-6):
            if int(jj) != j and final[int(jj)].distance(Point(x)) <= 1e-6:
                _cut_request(splits, final, int(jj), x)

    for kind, owner, w, xy, prefer in attach:
        if kind == "tcut":
            # a T the old network left without a node, which a plan's line
            # now ends at or runs through: row `prefer` is cut at that very
            # point, and so is every row passing exactly through it
            p = Point(xy)
            if final[prefer].distance(p) <= 1e-3:
                _cut_request(splits, final, prefer, xy)
            for jj in tree.query(p, predicate="dwithin", distance=1e-6):
                if int(jj) != prefer and final[int(jj)].distance(p) <= 1e-6:
                    _cut_request(splits, final, int(jj), xy)
            continue
        if kind in ("inner", "cut"):
            # a junction on a line cuts it there (other additions ending
            # there join each other); 'cut': a row a moved end now meets
            key = tuple(np.round(xy, 6))
            hit = place(xy, prefer=prefer, snap_ends=False)
            if hit is not None and np.hypot(*(hit[1] - xy)) <= 1e-6:
                cut_all(hit[0], hit[1])
            elif kind == "inner" and len(add_ends.get(key, [])) < 2:
                loose.append(("inner", np.asarray(xy, dtype=float)))
            continue
        hit = place(xy, skip=owner if kind == "row" else None, prefer=prefer)
        if hit is None or np.hypot(*(hit[1] - xy)) > JOG_M:
            loose.append((kind, np.asarray(xy, dtype=float)))
            continue
        j, x = hit
        if kind == "row":
            c = new_geom[owner]
            c = np.vstack([x, c] if w == 0 else [c, x])
            new_geom[owner] = _join([c])
            final[owner] = LineString(new_geom[owner])
        else:
            c = additions[owner][0]
            other = c[-1] if w == 0 else c[0]
            if np.hypot(*(other - x)) <= 1e-6:
                # its other end already landed there: joined, it would be a
                # ring. It lands on the line itself instead (no snap to the
                # row end) where that is somewhere else; failing that, this
                # end is left where Vann-Nett has it, and listed
                hit2 = place(xy, prefer=prefer, snap_ends=False)
                if (
                    hit2 is None
                    or np.hypot(*(hit2[1] - xy)) > JOG_M
                    or np.hypot(*(hit2[1] - other)) <= 0.5
                ):
                    loose.append(("tip", np.asarray(xy, dtype=float)))
                    continue
                j, x = hit2
            additions[owner][0] = _join(
                [np.vstack([x, c] if w == 0 else [c, x])]
            )
        cut_all(j, x)
    return splits, loose


def _node_crossings(additions, final, splits):
    """An addition crossing a line -- a row as it now runs, or another
    addition -- within JOG_M of an end of either (Vann-Nett and the network
    meeting a little short of or past a junction) gets a node there. All
    crossings are found first. An addition whose stretch beyond a crossing
    (at most JOG_M) runs on to the other line anyway drops that stretch -- at
    the crossing farthest from that end, where there are several; of two
    additions crossing, only the one dropping more; crossings on a dropped
    stretch cut nothing, and the cut it was joined at goes where nothing
    else meets the line there -- otherwise both lines are cut there (onto a
    cut the row already has within SNAP_M). A drop whose crossing, or whose
    far end (which must stay met), lies on a stretch the other addition
    drops gives way to the next smaller drop of that end, or to none.
    Crossings away from ends (a culvert, say) are left. Returns the
    additions; rows are cut through splits."""
    if not additions:
        return additions
    geoms = [LineString(c) for c, _ in additions]
    r_tree = shapely.STRtree(np.array(final, dtype=object))
    a_tree = shapely.STRtree(geoms)

    def snap(h, p):
        """p, or h's end where p lies within SNAP_M of it."""
        m = float(h.project(Point(p)))
        if m <= SNAP_M:
            return _end_xy(h, 0), m
        if m >= h.length - SNAP_M:
            return _end_xy(h, 1), m
        return p, m

    def tail(g, m, h):
        """g's end beyond m (0 or 1) if that stretch is at most JOG_M and
        ends on h, else None."""
        for w, (n, e) in enumerate(
            ((m, _end_xy(g, 0)), (g.length - m, _end_xy(g, 1)))
        ):
            if n <= JOG_M and h.distance(Point(e)) <= 1e-6:
                return w
        return None

    xs = []  # (addition, "row"/"add", the other line, point, m on g, m on h)
    for k, g in enumerate(geoms):
        hits = [("row", int(j)) for j in r_tree.query(g, predicate="crosses")]
        hits += [
            ("add", int(j))
            for j in a_tree.query(g, predicate="crosses")
            if int(j) > k
        ]
        for kind, j in hits:
            h = final[j] if kind == "row" else geoms[j]
            for pt in shapely.get_parts(shapely.intersection(g, h)):
                if pt.geom_type != "Point":
                    continue
                p, m_h = snap(h, np.asarray(pt.coords[0], dtype=float))
                if kind == "row":  # onto a cut the row already has
                    for q_ in splits.get(j, ()):
                        if np.hypot(*(np.asarray(q_) - p)) <= SNAP_M:
                            p = np.asarray(q_, dtype=float)
                            m_h = float(h.project(Point(p)))
                            break
                m_g = float(g.project(Point(p)))
                if not SNAP_M < m_g < g.length - SNAP_M:
                    continue  # at the addition's end: attached already
                if min(m_g, g.length - m_g, m_h, h.length - m_h) > JOG_M:
                    continue  # streams crossing away from any end
                xs.append((k, kind, j, p, m_g, m_h))
    if not xs:
        return additions
    # per addition end: (measure, point, length dropped, the other line)
    trunc = defaultdict(dict)
    alts = defaultdict(list)  # (addition, end) -> every candidate
    for k, kind, j, p, m_g, m_h in xs:
        h = final[j] if kind == "row" else geoms[j]
        sides = [(k, geoms[k], m_g, h, (kind, j))]
        if kind == "add":
            sides.append((j, geoms[j], m_h, geoms[k], ("add", k)))
        cand = []
        for kk, g, m, hh, oth in sides:
            w = tail(g, m, hh)
            if w is None:
                continue
            cand.append((m if w == 0 else g.length - m, kk, w, m, oth))
        if len(cand) == 2:  # K1: one of the two drops a stretch here
            cand = [max(cand, key=lambda x: x[0])]
        for drop, kk, w, m, oth in cand:
            alts[(kk, w)].append((m, p, drop, oth))
            cur = trunc[kk].get(w)
            if cur is None or drop > cur[2]:
                trunc[kk][w] = (m, p, drop, oth)
    for kk in list(trunc):  # a line cannot lose all of itself
        t = trunc[kk]
        if 0 in t and 1 in t and t[1][0] - t[0][0] <= SNAP_M:
            del t[0 if t[0][2] <= t[1][2] else 1]

    def dropped(kk, m):
        t = trunc.get(kk, {})
        return (0 in t and m < t[0][0] - 1e-9) or (
            1 in t and m > t[1][0] + 1e-9
        )

    def on_gone(x, kk, w):
        # the truncation point, or the end of the stretch it drops (which
        # must stay met), lies on a stretch the other addition drops
        _m, p_, _d, (kind_, j_) = x
        if kind_ != "add":
            return False
        e_ = _end_xy(geoms[kk], w)
        return any(
            dropped(j_, float(geoms[j_].project(Point(q_)))) for q_ in (p_, e_)
        )

    changed = True
    while changed:  # K1b, falling back to a smaller truncation of the end
        changed = False
        for kk in list(trunc):
            for w in list(trunc[kk]):
                if on_gone(trunc[kk][w], kk, w):
                    d0 = trunc[kk][w][2]
                    nxt = sorted(
                        (
                            x
                            for x in alts[(kk, w)]
                            if x[2] < d0 and not on_gone(x, kk, w)
                        ),
                        key=lambda x: -x[2],
                    )
                    if nxt:
                        trunc[kk][w] = nxt[0]
                    else:
                        del trunc[kk][w]
                    changed = True

    def at_end(kk, m):
        return any(abs(x[0] - m) <= 1e-9 for x in trunc.get(kk, {}).values())

    cuts = defaultdict(list)
    for k, kind, j, p, m_g, m_h in xs:
        if dropped(k, m_g) or (kind == "add" and dropped(j, m_h)):
            continue  # on a stretch that goes
        if not at_end(k, m_g):
            cuts[k].append(p)
        if kind == "row":
            _cut_request(splits, final, j, p)
        elif not at_end(j, m_h):
            cuts[j].append(p)
    out, gone = [], []  # gone: (row, point) where a dropped stretch ended
    for k, (c, i) in enumerate(additions):
        t = trunc.get(k)
        if t:
            g = geoms[k]
            m0 = t[0][0] if 0 in t else 0.0
            m1 = t[1][0] if 1 in t else g.length
            c = _sub(g, m0, m1)
            for w, (_m, p, _d, (kind, j)) in t.items():
                c[0 if w == 0 else -1] = p
                if kind == "row":
                    gone.append((j, _end_xy(g, w)))
            c = _join([c])
        out += [[pc, i] for pc in _cut_at(c, cuts.get(k, []))]
    if gone:
        # the row a dropped stretch was joined at keeps no cut there for it,
        # unless another addition or a row ends there too
        a_end = {_pkey(pc[x]) for pc, _ in out for x in (0, -1)}
        for j, e in gone:
            if _pkey(e) in a_end:
                continue
            ep = Point(e)
            if any(
                int(jj) != j
                and min(
                    np.hypot(*(_end_xy(final[int(jj)], 0) - e)),
                    np.hypot(*(_end_xy(final[int(jj)], 1) - e)),
                )
                <= 1e-6
                for jj in r_tree.query(ep, predicate="dwithin", distance=1e-6)
                if not final[int(jj)].is_empty
            ):
                continue
            splits[j] = [
                q_ for q_ in splits.get(j, []) if np.hypot(*(q_ - e)) > 1e-6
            ]
    return out


def _drawn_twice(additions, P):
    """Two additions joining the same two points and lying within on_m of
    each other all along are one stream drawn twice (Vann-Nett digitising a
    narrow lens): the shorter goes. Returns the additions and the lines no
    addition draws any more."""
    by_ends = defaultdict(list)
    for k, (c, _i) in enumerate(additions):
        if len(c) < 2:
            continue
        a, b = _pkey(c[0]), _pkey(c[-1])
        if a != b:
            by_ends[tuple(sorted((a, b)))].append(k)
    drop = set()
    for ks in by_ends.values():
        if len(ks) < 2:
            continue
        geoms = {k: LineString(additions[k][0]) for k in ks}
        ks = sorted(ks, key=lambda k: (-geoms[k].length, k))
        for n, x in enumerate(ks[1:], start=1):
            if any(  # a longer one it lies along is kept
                y not in drop
                and shapely.hausdorff_distance(geoms[x], geoms[y])
                <= P.on_m - 0.1
                for y in ks[:n]
            ):
                drop.add(x)
    if not drop:
        return additions, set()
    out = [a for k, a in enumerate(additions) if k not in drop]
    left = {i for _c, i in out}
    return out, {additions[k][1] for k in drop} - left


def _cut_request(splits, final, j, x):
    gj = final[j]
    m = gj.project(Point(x))
    if SNAP_M < m < gj.length - SNAP_M:
        splits[j].append(np.asarray(x, dtype=float))


def _cut_at(coords, points):
    """The line cut at points lying on it; each piece ends exactly there."""
    if not len(points):
        return [coords]
    line = LineString(coords)
    at = sorted(
        {
            round(float(line.project(Point(p))), 6): np.asarray(p)
            for p in points
        }.items()
    )
    at = [(m, p) for m, p in at if SNAP_M < m < line.length - SNAP_M]
    merged = []
    for m, p in at:  # cuts closer than SNAP_M are one
        if merged and m - merged[-1][0] <= SNAP_M:
            continue
        merged.append((m, p))
    if not merged:
        return [coords]
    edges = [(0.0, None), *merged, (line.length, None)]
    out = []
    for (a, pa), (b, pb) in zip(edges[:-1], edges[1:], strict=True):
        c = shapely.get_coordinates(substring(line, a, b))
        if pa is not None:
            c[0] = pa
        if pb is not None:
            c[-1] = pb
        c = _join([c])  # substring repeats a cut that falls on a vertex
        if len(c) >= 2:
            # on the cut point itself, also where _join took it for a repeat
            # of a vertex within 1e-9 m: both pieces then share it exactly
            if pa is not None:
                c[0] = pa
            if pb is not None:
                c[-1] = pb
            out.append(c)
    return out


def _orient(additions, final, splits):
    """An addition keeps Vann-Nett's direction unless turning it round
    takes away a node where nothing flows in or nothing flows out. Returns
    [(addition, node)] where either way leaves one."""
    if not additions:
        return []
    ends = defaultdict(lambda: [0, 0])  # node -> [inflow, outflow] of rows
    for g in final:
        c = shapely.get_coordinates(g)
        if not len(c):
            continue  # a dropped row
        ends[tuple(np.round(c[-1], 3))][0] += 1
        ends[tuple(np.round(c[0], 3))][1] += 1
    for pts in splits.values():
        for x in pts:  # a cut point is a through node: one in, one out
            k = tuple(np.round(x, 3))
            ends[k][0] += 1
            ends[k][1] += 1
    a_ends = defaultdict(list)
    for k, (c, _) in enumerate(additions):
        a_ends[tuple(np.round(c[0], 3))].append((k, 0))
        a_ends[tuple(np.round(c[-1], 3))].append((k, 1))
    out = []
    for k, (c, _) in enumerate(additions):
        s, e = tuple(np.round(c[0], 3)), tuple(np.round(c[-1], 3))

        def conflicts(start, end, k=k):
            bad = []
            for node, out_ in ((start, True), (end, False)):
                i_n, o_n = ends.get(node, [0, 0])
                for kk, w in a_ends.get(node, []):
                    if kk == k:
                        continue
                    if w == 0:
                        o_n += 1
                    else:
                        i_n += 1
                if i_n + o_n == 0:
                    continue  # a free end
                if out_:
                    o_n += 1
                else:
                    i_n += 1
                if i_n == 0 or o_n == 0:
                    bad.append(node)
            return bad

        keep, turn = conflicts(s, e), conflicts(e, s)
        if len(turn) < len(keep):
            additions[k][0] = c[::-1]
            a_ends[s] = [(kk, 1 - w if kk == k else w) for kk, w in a_ends[s]]
            a_ends[e] = [(kk, 1 - w if kk == k else w) for kk, w in a_ends[e]]
            keep = turn
        out += [(k, np.asarray(node, dtype=float), "divide") for node in keep]
    return out


def _cycles(additions, final, splits):
    """Water going round: directed cycles through additions. An addition in
    one is turned round where that breaks the cycle and leaves no node with
    only inflow or only outflow -- failing that, where the only such node
    is the tip of a ring of two additions, or a free end of the rows (a
    head or mouth: a cycle traps the water, a source there does not).
    Returns [(addition, node, kind)]: each node such a turn left a sink or
    source (kind "sink"/"source", under the turned addition), and the
    cycles left ("cycle")."""
    if not additions:
        return []
    starts, ends_, kinds = [], [], []  # edges: row pieces, then additions
    for j, g in enumerate(final):
        c = shapely.get_coordinates(g)
        if not len(c):
            continue
        pts = [c[0]]
        if splits.get(j):
            line = LineString(c)
            ms = sorted(line.project(Point(p)) for p in splits[j])
            pts += [_xy(line.interpolate(m)) for m in ms]
        pts.append(c[-1])
        for a, b in zip(pts[:-1], pts[1:], strict=True):
            starts.append(a)
            ends_.append(b)
            kinds.append(-1)
    first_add = len(starts)
    for k, (c, _) in enumerate(additions):
        starts.append(c[0])
        ends_.append(c[-1])
        kinds.append(k)
    q = np.round(np.vstack([starts, ends_]) * 1000).astype(np.int64)
    uq, inv = np.unique(q, axis=0, return_inverse=True)
    inv = inv.ravel()
    ne = len(starts)
    u, v = inv[:ne].copy(), inv[ne:].copy()
    nn = int(inv.max()) + 1

    def scc():
        m = sparse.csr_matrix(
            (np.ones(ne, dtype=np.int8), (u, v)), shape=(nn, nn)
        )
        return connected_components(m, directed=True, connection="strong")[1]

    def flow(node):
        i_n = int((v == node).sum())
        o_n = int((u == node).sum())
        return i_n, o_n

    rdeg = np.bincount(np.r_[u[:first_add], v[:first_add]], minlength=nn)

    def end_ok(n, level):
        """May node n be left with only inflow or only outflow? At level 1
        the tip of two additions hanging off one node (a ring of two); at
        level 2 also a free end of the rows (a stream head or mouth that an
        addition now shares). Either is then listed, as a sink or source."""
        if level >= 1 and rdeg[n] == 0:
            es = np.flatnonzero((u == n) | (v == n))
            far = {int(v[x]) if u[x] == n else int(u[x]) for x in es}
            if len(es) == 2 and len(far) == 1 and n not in far:
                return True
        return level >= 2 and rdeg[n] == 1

    lab = scc()
    out = []
    relaxed = []  # (addition, node): ends of additions turned at level 1/2
    for level in (0, 1, 2):
        for _round in range(3):
            same = lab[u] == lab[v]
            cyc_adds = [
                e for e in range(first_add, ne) if same[e] and u[e] != v[e]
            ]
            if not cyc_adds:
                break
            changed = False
            for e in cyc_adds:
                a, b = u[e], v[e]
                ia, oa = flow(a)
                ib, ob = flow(b)
                # turned round: a loses an outflow and gains an inflow, b
                # the other way; not if that leaves either without one
                if (oa - 1 == 0 and ia + oa >= 2 and not end_ok(a, level)) or (
                    ib - 1 == 0 and ib + ob >= 2 and not end_ok(b, level)
                ):
                    continue
                u[e], v[e] = b, a
                lab2 = scc()
                if lab2[u[e]] != lab2[v[e]]:
                    k = kinds[e]
                    additions[k][0] = additions[k][0][::-1]
                    if level:
                        relaxed += [(k, int(a)), (k, int(b))]
                    lab = lab2
                    changed = True
                else:
                    u[e], v[e] = a, b
            if not changed:
                break
    # a node such a turn left with only inflow or only outflow is listed
    # under the turned addition (by this pass, not the after-pass check)
    seen = set()
    for k, n in reversed(relaxed):
        i_n, o_n = flow(n)
        if n not in seen and i_n + o_n >= 2 and (i_n == 0 or o_n == 0):
            seen.add(n)
            out.append((k, uq[n] / 1000.0, "sink" if o_n == 0 else "source"))
    same = lab[u] == lab[v]
    for e in range(first_add, ne):
        if same[e]:  # in a cycle, or a ring (start = end)
            k = kinds[e]
            out.append(
                (k, np.asarray(additions[k][0][0], dtype=float), "cycle")
            )
    return out


def _pairs(final, additions, new_geom, vg, v_tree, V, lines, ctx):
    """ELVIS drawing a stream twice, near what this pass changed: stretches of
    a row off Vann-Nett (more than on_m) that run beside a Vann-Nett line
    another row, an addition or another piece of the same row carries --
    within near_m, the same way -- each from half a sample spacing before
    its first such sample to half one after its last (gaps of up to 5
    samples bridged). [(stretch, row, Vann-Nett id, why)]; why is the reason
    the missing line beside it was not taken."""
    P = ctx.P
    fa = np.array(final, dtype=object)
    adds = [LineString(c) for c, _ in additions if len(c) >= 2]
    changed = [fa[r] for r in new_geom] + adds
    if not changed:
        return []
    tree = shapely.STRtree(fa)
    _, hit = tree.query(
        np.array(changed, dtype=object), predicate="dwithin", distance=P.near_m
    )
    rows = np.array(
        sorted(
            int(j)
            for j in set(hit.tolist())
            if not fa[j].is_empty and fa[j].length >= 2.0
        ),
        dtype=np.int64,
    )
    if not len(rows):
        return []
    L = shapely.length(fa[rows])
    nper = np.maximum(4, np.ceil(L / 2.0).astype(int))
    rid = np.repeat(np.arange(len(rows)), nper - 1)
    m = np.concatenate(
        [np.linspace(0, ll, k + 1)[1:-1] for ll, k in zip(L, nper, strict=True)]
    )
    geoms = fa[rows][rid]
    pts = shapely.line_interpolate_point(geoms, m)
    tg = _tangent(geoms, m)
    on = _matches(pts, tg, vg, v_tree, P.on_m, P.max_angle)[0] >= 0
    vi, _vd, vm = _matches(pts, tg, vg, v_tree, P.near_m, P.max_angle)
    cand = np.flatnonzero(~on & (vi >= 0))
    beside = np.zeros(len(m), dtype=bool)
    near_k = np.full(len(m), -1)  # the missing line under the matched point
    if len(cand):
        carriers = np.array(list(fa) + adds, dtype=object)
        c_tree = shapely.STRtree(carriers)
        vps = shapely.line_interpolate_point(vg[vi[cand]], vm[cand])
        a_i, c_i = c_tree.query(vps, predicate="dwithin", distance=P.on_m)
        own = rows[rid[cand[a_i]]]
        other = c_i != own
        same = np.flatnonzero(~other)
        if len(same):
            # the row itself carries the point -- farther along it than the
            # gap between them: another piece of it, beside this one, where
            # the sample lies on its old ELVIS course (not on a connector)
            s_i = cand[a_i[same]]
            mp = shapely.line_locate_point(fa[own[same]], vps[a_i[same]])
            gap = shapely.distance(pts[s_i], vps[a_i[same]])
            old = shapely.distance(pts[s_i], ctx.ng[own[same]]) <= 1e-6
            # beside the Vann-Nett line, not beyond its end (an ELVIS head
            # running on past Vann-Nett's, a centreline in a lake)
            vlen = shapely.length(vg[vi[s_i]])
            mid = (vm[s_i] > 1e-6) & (vm[s_i] < vlen - 1e-6)
            far = np.abs(mp - m[s_i]) > gap + 2 * P.on_m
            other[same[far & old & mid]] = True
        beside[cand[a_i[other]]] = True
        (qi, qk), qd = ctx.line_tree.query_nearest(
            vps, return_distance=True, all_matches=False
        )
        on_l = qd <= 1e-6
        near_k[cand[qi[on_l]]] = qk[on_l]
    vn_ids = V[VN_ID].to_numpy() if VN_ID in V else np.full(len(V), None)
    out = []
    bounds = np.r_[0, np.cumsum(nper - 1)]
    for i, j in enumerate(rows):
        lo_, hi_ = bounds[i], bounds[i + 1]
        on_r, bes_r, m_r = on[lo_:hi_], beside[lo_:hi_], m[lo_:hi_]
        if on_r.all() or not bes_r.any():
            continue
        half = 0.5 * L[i] / nper[i]  # half a sample spacing
        for s_, e_ in _beside_runs(on_r, bes_r):
            if m_r[e_] - m_r[s_] < 1.0:
                continue
            run = slice(s_, e_ + 1)
            sel = bes_r[run]
            ids = [vn_ids[x] for x in vi[lo_:hi_][run][sel]]
            # the most frequent; a tie goes to the first along the row (a
            # set's order would change with the string hash seed per run)
            vid = Counter(ids).most_common(1)[0][0] if ids else None
            stretch = LineString(
                _sub(fa[j], max(0.0, m_r[s_] - half), min(L[i], m_r[e_] + half))
            )
            kr = near_k[lo_:hi_][run][sel]
            kr = kr[kr >= 0]
            why = "beside a Vann-Nett line another row carries"
            if len(kr):  # the reason of the missing line it runs beside
                k = int(Counter(kr.tolist()).most_common(1)[0][0])
                if k not in ctx.used:
                    why = ctx.reason.get(
                        k, "an ELVIS line beside it: no rule applies"
                    )
            out.append((stretch, int(j), vid, why))
    return out


def _beside_runs(on_r, bes_r, gap=5):
    """(first, last) sample of each stretch of a row's samples that lies
    off Vann-Nett and beside a carried line, gaps of up to gap samples
    (about 2 m each) bridged within one off run."""
    out = []
    s_, n = 0, len(on_r)
    while s_ < n:
        if on_r[s_]:
            s_ += 1
            continue
        e_ = s_
        while e_ + 1 < n and not on_r[e_ + 1]:
            e_ += 1
        hit = np.flatnonzero(bes_r[s_ : e_ + 1]) + s_
        if len(hit):
            a = b = int(hit[0])
            for h in hit[1:]:
                if h - b - 1 <= gap:
                    b = int(h)
                else:
                    out.append((a, b))
                    a = b = int(h)
            out.append((a, b))
        s_ = e_ + 1
    return out


def _rowkey(net, r):
    eid, lnr = net.at[r, "elvID"], net.at[r, "strekningLnr"]
    return f"{eid}" if pd.isna(lnr) else f"{eid}_{lnr}"


def _node_table(net):
    """Row ends to the millimetre: node keys (mm), inflows, outflows and a
    row there, per node."""
    g = shapely.get_geometry(net.geometry.values, 0)
    g = np.where(shapely.get_type_id(g) == 1, g, net.geometry.values)
    a, b = _ends(g)
    n = len(g)
    q = np.round(np.vstack([a, b]) * 1000).astype(np.int64)
    uniq, inv = np.unique(q, axis=0, return_inverse=True)
    inv = inv.ravel()
    outs = np.bincount(inv[:n], minlength=len(uniq))
    ins = np.bincount(inv[n:], minlength=len(uniq))
    first = np.full(len(uniq), n)
    np.minimum.at(first, inv, np.r_[np.arange(n), np.arange(n)])
    return uniq, ins, outs, first


def _flow_ends(net):
    """Nodes (to the millimetre) where two or more rows meet and water only
    flows in ("sink") or only out ("source"): {key: (kind, a row there,
    point)}."""
    if not len(net):
        return {}
    uniq, ins, outs, first = _node_table(net)
    bad = np.flatnonzero((ins + outs >= 2) & ((ins == 0) | (outs == 0)))
    return {
        tuple(uniq[i]): (
            "sink" if outs[i] == 0 else "source",
            int(first[i]),
            uniq[i] / 1000.0,
        )
        for i in bad
    }


def _lone_ends(net_in, net):
    """Row ends that met other row ends in the input and meet none in the
    output -- a row left joined to nothing: {key: ("loose", the row,
    point)}."""
    if not len(net) or not len(net_in):
        return {}
    u0, i0, o0, _ = _node_table(net_in)
    u1, i1, o1, f1 = _node_table(net)
    was = pd.DataFrame({"x": u0[:, 0], "y": u0[:, 1], "d": i0 + o0})
    one = np.flatnonzero(i1 + o1 == 1)
    now = pd.DataFrame({"x": u1[one, 0], "y": u1[one, 1], "i": one})
    hit = now.merge(was[was.d >= 2], on=["x", "y"])["i"].to_numpy()
    return {tuple(u1[x]): ("loose", int(f1[x]), u1[x] / 1000.0) for x in hit}


def _assemble(net, out_coords, src, cut, tags, patch_vn, lines, regine, start):
    """The output rows: untouched rows as they were, changed and cut rows
    with their tags, additions as new rows after them."""
    crs = net.crs
    cols = [c for c in net.columns if c != "geometry"]
    src = np.asarray(src)
    old = src >= 0
    out = pd.DataFrame(index=range(len(src)), columns=cols, dtype=object)
    out.loc[old, cols] = net.iloc[src[old]][cols].to_numpy()
    geoms = np.array(
        [
            net.geometry.values[s]
            if c is None
            else MultiLineString([LineString(c)])
            for c, s in zip(out_coords, src, strict=True)
        ],
        dtype=object,
    )
    if PATCH not in out:
        out[PATCH] = None
    if PATCH_VN not in out:
        out[PATCH_VN] = None
    changed = np.array([c is not None for c in out_coords])
    for k, s in enumerate(src):
        if s < 0 or out_coords[k] is None:
            continue
        new = set(tags.get(int(s), ()))
        if cut[k]:
            new.add(CUT)
        old_t = out.at[k, PATCH]
        t = set(old_t.split(",")) if isinstance(old_t, str) and old_t else set()
        out.at[k, PATCH] = ",".join(sorted(t | new)) or None
        if int(s) in patch_vn:
            out.at[k, PATCH_VN] = patch_vn[int(s)]

    add = np.flatnonzero(~old)
    if len(add):
        vnid = [lines[i].vn_id for i in -1 - src[add]]
        for c in cols:
            out[c] = out[c].astype("object")
        out.loc[add, cols] = None
        out.loc[add, "objektType"] = "ElvBekk"
        out.loc[add, "elvID"] = [f"VN-{v}" for v in vnid]
        out.loc[add, "strekningLnr"] = [
            f"VN{k:07d}" for k in range(start, start + len(add))
        ]
        out.loc[add, PATCH] = ADDED
        out.loc[add, PATCH_VN] = vnid
        if regine is not None and len(regine):
            mid = shapely.line_interpolate_point(
                np.array(
                    [shapely.get_parts(g)[0] for g in geoms[add]], dtype=object
                ),
                0.5,
                normalized=True,
            )
            hit = gpd.sjoin(
                gpd.GeoDataFrame(geometry=mid, crs=crs),
                regine.to_crs(crs)[
                    ["vassdragsnummer", "elvNavnHierarki", "geometry"]
                ],
                how="left",
                predicate="within",
            )
            hit = hit[~hit.index.duplicated()]
            out.loc[add, "vassdragsnummer"] = hit["vassdragsnummer"].to_numpy()
            out.loc[add, "elvNavnHierarki"] = hit["elvNavnHierarki"].to_numpy()
        fixed = np.array(
            pd.to_numeric(out["elveordenStrahler"], errors="coerce"),
            dtype=float,
        )
        fixed[add] = np.nan
        order = _strahler(
            np.array([shapely.get_parts(g)[0] for g in geoms], dtype=object),
            fixed,
        )
        out.loc[add, "elveordenStrahler"] = order[add]
    L = np.array([g.length for g in geoms])
    if "SHAPE_Length" in out:
        out.loc[changed, "SHAPE_Length"] = L[changed]
    if "elvelengde_m" in out:
        out.loc[changed, "elvelengde_m"] = np.round(L[changed], 2)
    out = _restore_dtypes(out, net)
    return gpd.GeoDataFrame(out, geometry=list(geoms), crs=crs)
