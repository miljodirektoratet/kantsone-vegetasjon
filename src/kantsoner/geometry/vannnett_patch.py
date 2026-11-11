"""ELVIS with Vann-Nett's river lines patched in: the ones ELVIS lacks, and
Vann-Nett's drawing wherever both have a stream but draw it differently.

The pipeline walks ELVIS (stages._elvis), so a Vann-Nett water body with no
ELVIS line under it gets no segments, and a segment is coupled to its water
body only where it lies within 1 m of the Vann-Nett line. Measured
2026-10-26: Vann-Nett's river layer is ELVIS stretches copied and grouped
into water bodies -- 1.49 M of its 1.51 M parts lie on ELVIS to the metre,
split at the same junctions -- but about 400 km has no ELVIS at all
(Ranfjorden, Fosen, Jæren ...) and about 6,000 km is the same stream drawn
again, tens of metres off, with junctions and sources in other places
(Nea-Nidelva throughout, closer to FKB-Vann than ELVIS there; Voss, Vefsn,
Sogn in parts).

ELVIS stays the backbone: every stretch keeps its row, identity and
attributes, and the network stays noded the way ELVIS is (stretches meet at
end points; ELVIS itself has only 51 ends lying on another line's interior).

- Redrawn: an ELVIS stretch that has a Vann-Nett counterpart off its line
  takes the Vann-Nett line. Its end nodes are mapped onto Vann-Nett first,
  once per node so every stretch meeting there agrees: to the Vann-Nett node
  the counterparts run into (followed the way the ELVIS lines run, and no
  Vann-Nett node taken twice), else to the nearest point of the main
  counterpart. The stretch then follows the Vann-Nett path between its
  mapped ends, or keeps its ELVIS line if that path strays from it. A head
  or tail Vann-Nett doesn't draw (ELVIS reaching further up a valley) keeps
  the ELVIS line, joined to the Vann-Nett path where the two part.
- ELVIS stretches meeting a moved node with no counterpart of their own
  (lake midlines, streams Vann-Nett lacks) have their end pulled to the new
  position over the last few tens of metres; the rest of the line stays put.
- Vann-Nett lines left over -- not on ELVIS, not taken by a redrawn stretch,
  not alongside an ELVIS line that kept its own drawing -- are added as new
  stretches. Where one ends on another line within 1 m it is snapped there,
  and a line it ends on mid-way is split so the junction is a node.
- ELVIS's exact duplicate features (same key, same line) are dropped, as
  elvis._stretches does anyway; multipart features come out one row per part.
  So is a short ELVIS-only line whose two nodes land on the same place (a
  pond Vann-Nett runs straight through): 17 nationally, listed in the QA.

Additions are ``objektType`` ElvBekk with ``elvID`` "VN-<vannforekomstid>",
``strekningLnr`` "VN<number>", Strahler order worked out from the patched
network, and vassdragsnummer / elvNavnHierarki from the REGINE unit they lie
in. The ``patch`` column says what happened to each row; ``patch_vn`` is the
water body whose line a row now follows.
"""

from __future__ import annotations

import heapq
from collections import defaultdict
from dataclasses import dataclass, field

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from loguru import logger
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import substring

VN_ID = "vannforekomstid"
PATCH = "patch"
PATCH_VN = "patch_vn"

# patch tags, comma-joined in the PATCH column
REDRAWN = "vn_geometry"  # (part of) the line now follows Vann-Nett
MOVED = "end_moved"  # an end pulled to a node that moved
SPLIT = "split"  # cut where another line now joins mid-way
ADDED = "vn_added"  # a Vann-Nett line ELVIS doesn't have


@dataclass(frozen=True)
class Params:
    step: float = 2.0  # sample spacing along lines
    on_m: float = 1.0  # this close, the two datasets draw the same line
    near_m: float = 50.0  # a redrawn stream stays this close to ELVIS's
    max_angle: float = 45.0  # ... running the same way (either sense)
    min_run: float = 20.0  # shorter covered/uncovered runs are absorbed
    end_m: float = 100.0  # a stretch end's counterpart is read this far in
    node_r: float = 150.0  # how far a junction may move to Vann-Nett's
    sliver_m: float = 10.0  # left-over Vann-Nett shorter than this is dropped
    parallel_frac: float = 0.5  # share alongside a kept ELVIS line = its stream
    slack_m: float = 30.0  # path check tolerance on top of near_m


@dataclass
class Result:
    network: gpd.GeoDataFrame
    qa: dict[str, gpd.GeoDataFrame] = field(default_factory=dict)
    stats: dict[str, float] = field(default_factory=dict)


# ── small geometry helpers ─────────────────────────────────────────────────
def _ends(geoms):
    return (
        shapely.get_coordinates(shapely.get_point(geoms, 0)),
        shapely.get_coordinates(shapely.get_point(geoms, -1)),
    )


def _node_ids(*xys):
    """One id per distinct coordinate (to the millimetre) across the arrays."""
    q = np.round(np.vstack(xys) * 1000).astype(np.int64)
    uniq, inv = np.unique(q, axis=0, return_inverse=True)
    inv = inv.ravel()
    out, i = [], 0
    for xy in xys:
        out.append(inv[i : i + len(xy)])
        i += len(xy)
    return out, uniq / 1000.0


def _nearest_dist(tree, pts):
    """Distance from each point to the nearest geometry in tree (inf if empty)."""
    if not len(tree.geometries):
        return np.full(len(pts), np.inf)
    _, d = tree.query_nearest(pts, return_distance=True, all_matches=False)
    return d


def _groups(idx):
    """Slices of the runs of equal values in a sorted index array."""
    if not len(idx):
        return []
    b = np.flatnonzero(np.r_[True, idx[1:] != idx[:-1], True])
    return [slice(b[k], b[k + 1]) for k in range(len(b) - 1)]


def _sample(geoms, step):
    """(line index, measure, point) about every ``step`` m, both ends included."""
    L = shapely.length(geoms)
    n = np.maximum(1, np.ceil(L / step)).astype(np.int64)
    cnt = n + 1
    idx = np.repeat(np.arange(len(geoms)), cnt)
    k = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    m = k * np.repeat(L / n, cnt)
    return idx, m, shapely.line_interpolate_point(geoms[idx], m)


def _tangent(geoms, m, h=2.0):
    """Unit direction of each line at measure m."""
    L = shapely.length(geoms)
    a = shapely.get_coordinates(
        shapely.line_interpolate_point(geoms, np.clip(m - h, 0, L))
    )
    b = shapely.get_coordinates(
        shapely.line_interpolate_point(geoms, np.clip(m + h, 0, L))
    )
    d = b - a
    n = np.hypot(d[:, 0], d[:, 1])
    return d / np.where(n > 0, n, 1)[:, None]


def _piece(line, m0, m1):
    """Coordinates of line from measure m0 to m1 (reversed when m0 > m1)."""
    lo, hi = min(m0, m1), max(m0, m1)
    if hi - lo < 1e-6:
        return np.empty((0, 2))
    c = shapely.get_coordinates(substring(line, lo, hi))
    return c[::-1] if m0 > m1 else c


def _join(parts):
    """Concatenate coordinate runs, dropping repeated points at the seams."""
    parts = [p for p in parts if len(p)]
    if not parts:
        return np.empty((0, 2))
    out = [parts[0]]
    for p in parts[1:]:
        if np.allclose(out[-1][-1], p[0], atol=1e-6):
            p = p[1:]
        out.append(p)
    c = np.vstack(out)
    keep = np.r_[True, np.any(np.abs(np.diff(c, axis=0)) > 1e-9, axis=1)]
    return c[keep]


def _move_end(coords, at_start, xy, min_d=30.0):
    """Pull one end of a line to xy, the displacement fading out along the
    line over max(5x the move, min_d) metres; the other end stays put."""
    c = np.array(coords[::-1] if not at_start else coords, dtype=float)
    d = np.asarray(xy, dtype=float) - c[0]
    dist = float(np.hypot(*d))
    if dist == 0:
        return coords
    s = np.r_[0.0, np.cumsum(np.hypot(*np.diff(c, axis=0).T))]
    D = min(max(5 * dist, min_d), s[-1]) if s[-1] > 0 else 1.0
    w = np.clip(1 - s / D, 0, 1)
    w[0] = 1.0
    c = c + w[:, None] * d
    return c if at_start else c[::-1]


def _runs(m, flag, min_len):
    """Runs [(m0, m1, flag)] along one line; runs shorter than min_len are
    absorbed by their neighbours, shortest first."""
    cut = np.flatnonzero(flag[1:] != flag[:-1])
    edges = np.r_[m[0], (m[cut] + m[cut + 1]) / 2, m[-1]]
    firsts = np.r_[0, cut + 1]
    runs = [
        [edges[i], edges[i + 1], bool(flag[f])] for i, f in enumerate(firsts)
    ]
    while len(runs) > 1:
        lens = [r[1] - r[0] for r in runs]
        i = int(np.argmin(lens))
        if lens[i] >= min_len:
            break
        runs[i][2] = not runs[i][2]
        merged = []
        for r in runs:
            if merged and merged[-1][2] == r[2]:
                merged[-1][1] = r[1]
            else:
                merged.append(r)
        runs = merged
    return [tuple(r) for r in runs]


def _matches(pts, tang, lines, line_tree, dist, max_angle):
    """For each point, the nearest line within dist running the same way
    (either sense): (line index or -1, distance, measure on that line)."""
    best = np.full(len(pts), -1)
    bd = np.full(len(pts), np.inf)
    bm = np.zeros(len(pts))
    cos_min = np.cos(np.radians(max_angle))
    CH = 400_000
    for a in range(0, len(pts), CH):
        p, j = line_tree.query(
            pts[a : a + CH], predicate="dwithin", distance=dist
        )
        if not len(p):
            continue
        d = shapely.distance(pts[a + p], lines[j])
        mm = shapely.line_locate_point(lines[j], pts[a + p])
        t = _tangent(lines[j], mm)
        ok = np.abs(np.sum(t * tang[a + p], axis=1)) >= cos_min
        p, j, d, mm = p[ok], j[ok], d[ok], mm[ok]
        o = np.lexsort((d, p))
        p, j, d, mm = p[o], j[o], d[o], mm[o]
        first = np.r_[True, p[1:] != p[:-1]][: len(p)]
        p, j, d, mm = p[first], j[first], d[first], mm[first]
        best[a + p], bd[a + p], bm[a + p] = j, d, mm
    return best, bd, bm


# ── inputs ─────────────────────────────────────────────────────────────────
def elvis_lines(elvis: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """One row per distinct line part, exact copies dropped."""
    ex = elvis.explode(index_parts=False)
    ex = ex[ex.geometry.notna() & (ex.geometry.length > 0)].reset_index(
        drop=True
    )
    key = ex["elvID"].astype(str) + "_" + ex["strekningLnr"].astype(str)
    wkb = shapely.to_wkb(shapely.normalize(ex.geometry.values))
    dup = pd.DataFrame({"k": key.to_numpy(), "w": wkb}).duplicated().to_numpy()
    if dup.any():
        logger.info(f"ELVIS: {int(dup.sum()):,} exact duplicate parts dropped")
    return ex[~dup].reset_index(drop=True)


def vn_lines(vn: gpd.GeoDataFrame, tol: float = 0.5) -> gpd.GeoDataFrame:
    """Vann-Nett river parts, one row per line, noded at every junction.

    The layer comes split at junctions, apart from parts that end on
    another's interior rather than at its end -- a branch leaving mid-line;
    103 nationally to the centimetre, more within ``tol``. Such a part is cut
    there and the end snapped onto the cut, so the junction is a node.
    """
    cols = [c for c in (VN_ID, "vannforekomstnavn") if c in vn.columns]
    ex = vn[cols + ["geometry"]].explode(index_parts=False)
    ex = ex[ex.geometry.notna() & (ex.geometry.length > 0)].reset_index(
        drop=True
    )
    wkb = shapely.to_wkb(shapely.normalize(ex.geometry.values))
    ex = ex[~pd.Series(wkb).duplicated().to_numpy()].reset_index(drop=True)

    g = ex.geometry.values
    n = len(g)
    a, b = _ends(g)
    pts = shapely.points(np.vstack([a, b]))
    owner = np.r_[np.arange(n), np.arange(n)]
    (ia, ib), _ = _node_ids(a, b)
    ids = np.r_[ia, ib]
    dangling = np.bincount(ids)[ids] == 1  # meets no other part's end
    p, j = shapely.STRtree(g).query(pts, predicate="dwithin", distance=tol)
    keep = (owner[p] != j) & dangling[p]
    p, j = p[keep], j[keep]
    inner = (shapely.distance(pts[p], shapely.get_point(g[j], 0)) > tol) & (
        shapely.distance(pts[p], shapely.get_point(g[j], -1)) > tol
    )
    p, j = p[inner], j[inner]
    o = np.lexsort((shapely.distance(pts[p], g[j]), p))
    p, j = p[o], j[o]
    first = np.r_[True, p[1:] != p[:-1]][: len(p)]
    cuts = defaultdict(list)  # part -> [(measure, point)]
    snap = {}  # end (index into pts) -> point
    for pi, ji in zip(p[first], j[first], strict=True):
        m = float(shapely.line_locate_point(g[ji], pts[pi]))
        xy = shapely.get_coordinates(shapely.line_interpolate_point(g[ji], m))[
            0
        ]
        cuts[ji].append((m, xy))
        snap[pi] = xy
    if not cuts:
        return ex
    rows, geoms = [], []
    for i in range(n):
        ms = sorted(cuts.get(i, []), key=lambda x: x[0])
        edges = [(0.0, None), *ms, (g[i].length, None)]
        pieces = []
        for (m0, p0), (m1, p1) in zip(edges[:-1], edges[1:], strict=True):
            if m1 - m0 <= 1e-6:
                continue
            c = shapely.get_coordinates(substring(g[i], m0, m1))
            if p0 is not None:
                c[0] = p0
            if p1 is not None:
                c[-1] = p1
            pieces.append(c)
        if i in snap:
            pieces[0][0] = snap[i]
        if i + n in snap:
            pieces[-1][-1] = snap[i + n]
        rows += [i] * len(pieces)
        geoms += [LineString(c) for c in pieces]
    out = ex.iloc[rows].reset_index(drop=True)
    out = out.set_geometry(np.array(geoms, dtype=object), crs=ex.crs)
    logger.info(
        f"Vann-Nett: {len(snap)} part ends lying on another part's interior "
        f"(within {tol} m) joined there"
    )
    return out


# ── the Vann-Nett network ──────────────────────────────────────────────────
class _Net:
    """Directed network over a subset of Vann-Nett parts (digitising order).

    A position on it is ("n", node) or ("p", part, measure); a position at
    either end of a part is always given as the node.
    """

    def __init__(self, geoms, idx):
        self.g = geoms[idx]
        self.idx = np.asarray(idx)
        self.local = {int(gi): i for i, gi in enumerate(self.idx)}
        a, b = _ends(self.g)
        (self.s, self.e), self.xy = _node_ids(a, b)
        self.L = shapely.length(self.g)
        self.deg = np.bincount(np.r_[self.s, self.e], minlength=len(self.xy))
        self.out = defaultdict(list)
        self.inc = defaultdict(list)
        for i, (u, v) in enumerate(zip(self.s, self.e, strict=True)):
            self.out[u].append(i)
            self.inc[v].append(i)

    def pos(self, i, m):
        if m <= 1e-6:
            return ("n", int(self.s[i]))
        if m >= self.L[i] - 1e-6:
            return ("n", int(self.e[i]))
        return ("p", int(i), float(m))

    def xy_of(self, pos):
        if pos[0] == "n":
            return self.xy[pos[1]]
        return shapely.get_coordinates(
            shapely.line_interpolate_point(self.g[pos[1]], pos[2])
        )[0]

    def walk(self, i, forward, limit):
        """From part i along the line, through nodes where it just continues,
        to the first junction or end: (parts passed, that node)."""
        parts, dist = [i], 0.0
        node = self.e[i] if forward else self.s[i]
        while self.deg[node] == 2 and dist < limit:
            nxt = self.out[node] if forward else self.inc[node]
            if len(nxt) != 1 or nxt[0] in parts:
                break
            i = nxt[0]
            parts.append(i)
            dist += self.L[i]
            node = self.e[i] if forward else self.s[i]
        return parts, int(node)

    def path(self, a, b, directed, cutoff):
        """Pieces [(part, m_from, m_to)] from position a to b, or None."""
        pa = a[1] if a[0] == "p" else None
        pb = b[1] if b[0] == "p" else None
        if a == b:
            return []
        if pa is not None and pa == pb and (a[2] <= b[2] or not directed):
            return [(pa, a[2], b[2])]
        if a[0] == "n":
            starts = [(0.0, a[1], None)]
        else:
            starts = [(self.L[pa] - a[2], self.e[pa], (pa, a[2], self.L[pa]))]
            if not directed:
                starts.append((a[2], self.s[pa], (pa, a[2], 0.0)))
        if b[0] == "n":
            goals = {b[1]: (0.0, None)}
        else:
            goals = {self.s[pb]: (b[2], (pb, 0.0, b[2]))}
            if not directed:
                g2 = (self.L[pb] - b[2], (pb, self.L[pb], b[2]))
                if g2[0] < goals.get(self.e[pb], (np.inf,))[0]:
                    goals[self.e[pb]] = g2
        dist, prev, heap = {}, {}, []
        for c, node, piece in starts:
            if c < dist.get(node, np.inf):
                dist[node], prev[node] = c, (None, piece)
                heapq.heappush(heap, (c, node))
        best, best_node = np.inf, None
        while heap:
            c, u = heapq.heappop(heap)
            if c > dist.get(u, np.inf):
                continue
            if c > cutoff or c >= best:
                break
            if u in goals and c + goals[u][0] < best:
                best, best_node = c + goals[u][0], u
            steps = [(i, self.e[i], (i, 0.0, self.L[i])) for i in self.out[u]]
            if not directed:
                steps += [
                    (i, self.s[i], (i, self.L[i], 0.0)) for i in self.inc[u]
                ]
            for i, v, piece in steps:
                if i in (pa, pb):
                    continue
                nc = c + self.L[i]
                if nc < dist.get(v, np.inf):
                    dist[v], prev[v] = nc, (u, piece)
                    heapq.heappush(heap, (nc, v))
        if best_node is None:
            return None
        pieces, u = [goals[best_node][1]], best_node
        while True:
            pu, piece = prev[u]
            pieces.append(piece)
            if pu is None:
                break
            u = pu
        return [p for p in pieces[::-1] if p is not None]

    def reach(self, a, forward, cutoff):
        """Parts passed following the lines from position a downstream
        (forward) or upstream, within cutoff metres; a's own part included."""
        parts = {a[1]} if a[0] == "p" else set()
        if a[0] == "n":
            heap = [(0.0, a[1])]
        else:
            i, m = a[1], a[2]
            heap = [(self.L[i] - m, self.e[i]) if forward else (m, self.s[i])]
        dist = {heap[0][1]: heap[0][0]}
        while heap:
            c, u = heapq.heappop(heap)
            if c > dist.get(u, np.inf) or c > cutoff:
                continue
            for i in self.out[u] if forward else self.inc[u]:
                parts.add(i)
                v = self.e[i] if forward else self.s[i]
                if c + self.L[i] < dist.get(v, np.inf):
                    dist[v] = c + self.L[i]
                    heapq.heappush(heap, (c + self.L[i], v))
        return parts

    def reach_nearest(self, a, xy, forward, cutoff):
        """Following the lines from position a downstream (forward) or
        upstream, the node nearest xy: (node, pieces in flow order) or None."""
        if a[0] == "n":
            starts = [(0.0, a[1], None)]
        else:
            i, m = a[1], a[2]
            starts = [
                (self.L[i] - m, self.e[i], (i, m, self.L[i]))
                if forward
                else (m, self.s[i], (i, 0.0, m))
            ]
        dist, prev, heap = {}, {}, []
        for c, node, piece in starts:
            dist[node], prev[node] = c, (None, piece)
            heapq.heappush(heap, (c, node))
        while heap:
            c, u = heapq.heappop(heap)
            if c > dist.get(u, np.inf) or c > cutoff:
                continue
            for i in self.out[u] if forward else self.inc[u]:
                v = self.e[i] if forward else self.s[i]
                nc = c + self.L[i]
                if nc < dist.get(v, np.inf):
                    dist[v], prev[v] = nc, (u, (i, 0.0, self.L[i]))
                    heapq.heappush(heap, (nc, v))
        if not dist:
            return None
        nodes = np.fromiter(dist.keys(), dtype=np.int64)
        d = np.hypot(*(self.xy[nodes] - xy).T)
        best = int(nodes[np.lexsort(([dist[int(x)] for x in nodes], d))[0]])
        pieces, u = [], best
        while True:
            pu, piece = prev[u]
            pieces.append(piece)
            if pu is None:
                break
            u = pu
        pieces = [p for p in pieces if p is not None]
        return best, pieces[::-1] if forward else pieces

    def coords(self, pieces):
        return _join([_piece(self.g[i], m0, m1) for i, m0, m1 in pieces])

    def nearest_on(self, parts, xy):
        """Nearest position on any of the parts to xy."""
        pt = Point(xy)
        parts = list(parts)
        i = parts[int(np.argmin(shapely.distance(self.g[parts], pt)))]
        return self.pos(i, float(shapely.line_locate_point(self.g[i], pt)))


class _Consumed:
    """Which stretches of which Vann-Nett parts a redrawn ELVIS line took."""

    def __init__(self):
        self.iv = defaultdict(list)

    def overlap(self, pieces):
        """(metres of these pieces already taken, by whom)."""
        tot, who = 0.0, set()
        for i, m0, m1 in pieces:
            lo, hi = min(m0, m1), max(m0, m1)
            for a, b, owner in self.iv[i]:
                o = max(0.0, min(hi, b) - max(lo, a))
                if o > 0:
                    tot += o
                    who.add(owner)
        return tot, who

    def add(self, pieces, owner=None):
        for i, m0, m1 in pieces:
            self.iv[i].append((min(m0, m1), max(m0, m1), owner))

    def free(self, i, lo, hi):
        """Sub-intervals of [lo, hi] on part i nobody took."""
        out, cur = [], lo
        for a, b, _ in sorted(self.iv[i], key=lambda x: x[:2]):
            if b <= cur:
                continue
            if a >= hi:
                break
            if a > cur:
                out.append((cur, a))
            cur = max(cur, b)
        if cur < hi:
            out.append((cur, hi))
        return out


# ── the patch ──────────────────────────────────────────────────────────────
@dataclass
class _End:
    """One end of an ELVIS stretch, as the node mapping sees it."""

    state: str  # "V" redrawn up to here, "X" on Vann-Nett here, "E" neither
    parts: list = field(default_factory=list)  # Vann-Nett line it runs into
    own: int = -1  # the first Vann-Nett junction or end along that line


def patch(
    elvis: gpd.GeoDataFrame,
    vn: gpd.GeoDataFrame,
    regine: gpd.GeoDataFrame | None = None,
    params: Params | None = None,
) -> Result:
    P = params or Params()
    crs = elvis.crs
    E = elvis_lines(elvis)
    V = vn_lines(vn.to_crs(crs))
    eg, vg = E.geometry.values, V.geometry.values
    e_tree, v_tree = shapely.STRtree(eg), shapely.STRtree(vg)
    vn_ids = V[VN_ID].to_numpy() if VN_ID in V else np.full(len(V), None)
    stats: dict[str, float] = {
        "elvis_rows": len(elvis),
        "elvis_lines": len(E),
        "vn_parts": len(V),
    }

    # 1. Vann-Nett off ELVIS. A coarse pass finds the parts with anything
    # more than on_m from ELVIS; only those are sampled finely.
    idx, m, pts = _sample(vg, 10.0)
    cand_v = np.unique(idx[_nearest_dist(e_tree, pts) > P.on_m])
    vf_rows = []  # (part, m0, m1): runs of Vann-Nett off ELVIS
    if len(cand_v):
        idx, m, pts = _sample(vg[cand_v], P.step)
        d = _nearest_dist(e_tree, pts)
        for sl in _groups(idx):
            for m0, m1, off in _runs(m[sl], d[sl] > P.on_m, P.min_run):
                if off and m1 > m0:
                    vf_rows.append((int(cand_v[idx[sl][0]]), m0, m1))
    del idx, m, pts
    vf = pd.DataFrame(vf_rows, columns=["part", "m0", "m1"])
    vf_g = np.array(
        [substring(vg[p], a, b) for p, a, b in vf.itertuples(index=False)],
        dtype=object,
    )
    stats["vn_off_elvis_km"] = float((vf.m1 - vf.m0).sum() / 1000)
    logger.info(
        f"Vann-Nett: {len(cand_v):,} parts leave ELVIS somewhere, "
        f"{stats['vn_off_elvis_km']:,.0f} km of line more than {P.on_m} m off"
    )
    if not len(vf):
        return Result(
            _assemble(E, {}, {}, [], [], {}, crs, regine, stats), {}, stats
        )

    # 2. ELVIS stretches near that line, sampled: on Vann-Nett (exact), on a
    # same-way off-ELVIS Vann-Nett line within near_m (matched), or neither.
    # Each covered sample gets its counterpart position on Vann-Nett.
    cand_e = np.unique(
        e_tree.query(vf_g, predicate="dwithin", distance=P.near_m)[1]
    )
    sidx, sm, spts = _sample(eg[cand_e], P.step)
    stang = _tangent(eg[cand_e][sidx], sm)
    # exact: on a Vann-Nett line running the same way -- at a junction the
    # other lines meeting there, and lines crossing, don't count
    xj, _, xm = _matches(spts, stang, vg, v_tree, P.on_m, P.max_angle)
    exact = xj >= 0
    mj, md, mm = _matches(
        spts, stang, vf_g, shapely.STRtree(vf_g), P.near_m, P.max_angle
    )
    matched = (mj >= 0) & ~exact
    cp_part = np.where(matched, vf.part.to_numpy()[np.maximum(mj, 0)], -1)
    cp_m = np.where(matched, vf.m0.to_numpy()[np.maximum(mj, 0)] + mm, 0.0)
    cp_part[exact] = xj[exact]
    cp_m[exact] = xm[exact]
    covered = cp_part >= 0
    # does the counterpart run the same way as the ELVIS line?
    same_way = np.zeros(len(sm), dtype=bool)
    same_way[covered] = (
        np.sum(
            _tangent(vg[cp_part[covered]], cp_m[covered]) * stang[covered],
            axis=1,
        )
        >= 0
    )

    # 3. per stretch: the covered stretch [a, b] to take from Vann-Nett, the
    # uncovered head/tail kept as ELVIS line
    span: dict[int, tuple[float, float]] = {}
    samples: dict[int, slice] = {}
    exact_end: dict[tuple[int, int], bool] = {}
    matched_m: dict[int, float] = {}  # length with an off-line counterpart
    for sl in _groups(sidx):
        s = int(cand_e[sidx[sl][0]])
        samples[s] = sl
        mm_s, L = sm[sl], sm[sl][-1]
        exact_end[(s, 0)] = bool(exact[sl][mm_s <= min(10.0, L)].all())
        exact_end[(s, 1)] = bool(exact[sl][mm_s >= max(0.0, L - 10.0)].all())
        runs = _runs(mm_s, covered[sl], P.min_run)
        cov = [(m0, m1) for m0, m1, f in runs if f]
        if not cov:
            continue
        a, b = cov[0][0], cov[-1][1]
        inside = (mm_s >= a) & (mm_s <= b)
        step = mm_s[1] - mm_s[0] if len(mm_s) > 1 else 0.0
        matched_m[s] = matched[sl][inside].sum() * step
        # (a short connector counts when half of it has a counterpart)
        if matched_m[s] >= min(P.sliver_m, L / 2):
            span[s] = (a, b)
    redrawn = set(span)
    logger.info(
        f"ELVIS: {len(cand_e):,} stretches near off-ELVIS Vann-Nett line, "
        f"{len(redrawn):,} with a redrawn counterpart"
    )

    # 4. the local Vann-Nett network those stretches can be routed on
    zone = np.unique(
        v_tree.query(
            eg[cand_e], predicate="dwithin", distance=P.near_m + P.node_r
        )[1]
    )
    net = _Net(vg, zone)
    in_net = np.zeros(len(vg), dtype=bool)
    in_net[zone] = True

    ea, eb = _ends(eg)
    (e_s, e_e), e_xy = _node_ids(ea, eb)
    e_deg = np.bincount(np.r_[e_s, e_e], minlength=len(e_xy))
    strahler = (
        pd.to_numeric(E["elveordenStrahler"], errors="coerce")
        .fillna(0)
        .to_numpy()
        if "elveordenStrahler" in E
        else np.zeros(len(E))
    )

    def end_info(s, w):
        if s not in samples:  # no off-ELVIS Vann-Nett anywhere near it
            line = eg[s]
            m_in = (
                min(10.0, line.length)
                if w == 0
                else max(0.0, line.length - 10.0)
            )
            probe = [
                shapely.get_point(line, -w),
                shapely.line_interpolate_point(line, m_in),
            ]
            _, dd = v_tree.query_nearest(probe, return_distance=True)
            return _End("X" if np.all(dd <= P.on_m) else "E")
        sl = samples[s]
        mm_s, L = sm[sl], sm[sl][-1]
        if exact_end[(s, w)]:
            state = "X"
        elif s in span and (span[s][0] <= 0 if w == 0 else span[s][1] >= L):
            state = "V"
        else:
            return _End("E")
        # The Vann-Nett line this end runs into: the one the stretch follows
        # most over its last end_m -- not just the nearest to the end, which
        # at a junction is as often the other tributary's -- traced on along
        # Vann-Nett the way the ELVIS line runs, to the next junction or end.
        ok = covered[sl] & in_net[np.maximum(cp_part[sl], 0)]
        to_end = mm_s if w == 0 else L - mm_s
        win = np.flatnonzero(ok & (to_end <= P.end_m))
        if not len(win):
            win = np.flatnonzero(ok)
        if not len(win):
            return _End("E")
        parts = np.array([net.local[int(p)] for p in cp_part[sl][win]])
        vals, cnt = np.unique(parts, return_counts=True)
        top = vals[cnt == cnt.max()]
        p = top[0] if len(top) == 1 else parts[np.argmin(to_end[win])]
        aligned = same_way[sl][win][parts == p].mean() >= 0.5
        walked, own = net.walk(
            int(p), forward=aligned == (w == 1), limit=P.node_r
        )
        return _End(state, walked, own)

    # 5. place every node a redrawn end reaches on Vann-Nett, once per node
    ends_at = defaultdict(list)
    for s in redrawn:
        for node in (e_s[s], e_e[s]):
            ends_at[int(node)] = []
    nodes = np.fromiter(ends_at.keys(), dtype=np.int64)
    for arr, w in ((e_s, 0), (e_e, 1)):
        for s in np.flatnonzero(np.isin(arr, nodes)):
            ends_at[int(arr[s])].append((int(s), w))
    info = {se: end_info(*se) for ends in ends_at.values() for se in ends}

    node_pos: dict[int, tuple] = {}  # moved node -> its Vann-Nett position
    end_pos: dict[tuple[int, int], tuple] = {}  # (s, w) -> (position, join)
    moved: dict[int, np.ndarray] = {}
    node_rows, picks = [], []
    for node, ends in ends_at.items():
        nxy = e_xy[node]
        vends = [se for se in ends if info[se].state == "V" and info[se].parts]
        if any(info[se].state == "X" for se in ends):
            # An ELVIS stretch lies on Vann-Nett here: the node stays. A
            # stretch redrawn into it ends where its counterpart meets the
            # network, if that is a Vann-Nett junction or end within reach
            # on an ELVIS line (attached there later), else at its
            # counterpart's nearest point with a connector to the node.
            for se in ends:
                i = info[se]
                if not i.parts or i.state == "E":
                    continue
                pos = net.nearest_on(i.parts, nxy)
                if (
                    i.state == "V"
                    and np.hypot(*(net.xy_of(pos) - nxy)) > P.on_m
                ):
                    # its counterpart doesn't pass the node: end where that
                    # line meets the network instead, if it does nearby
                    cxy = net.xy[i.own]
                    on_net = len(
                        e_tree.query(
                            Point(cxy), predicate="dwithin", distance=P.on_m
                        )
                    )
                    if on_net and np.hypot(*(cxy - nxy)) <= P.node_r:
                        end_pos[se] = (("n", i.own), None)
                        continue
                end_pos[se] = (pos, nxy)
            continue
        if not vends:
            continue
        # Otherwise the node moves: to the Vann-Nett node the redrawn ends
        # run into if they agree on it -- a junction for a junction, an end
        # for a source or outlet, never a junction for a node between two
        # stretches (a tributary's bank node, a lake inlet: that junction is
        # the next ELVIS node's) -- else onto the main counterpart.
        k, deg = len(vends), int(e_deg[node])
        votes = defaultdict(int)
        for se in vends:
            votes[info[se].own] += 1
        cands = []
        for c, n_own in votes.items():
            dd = float(np.hypot(*(net.xy[c] - nxy)))
            if dd > P.node_r:
                continue
            if deg == 1:
                ok = net.deg[c] == 1
            elif deg == 2:
                ok = net.deg[c] <= 2 and dd <= P.near_m and n_own == k
            elif n_own == k:
                ok = k >= 2 or dd <= P.near_m
            else:
                ok = n_own >= 2 and net.deg[c] >= 3
            if ok:
                cands.append((c, (-n_own, dd)))
        best = min(cands, key=lambda x: x[1]) if cands else (-1, (0, np.inf))
        picks.append((node, best[0], best[1]))

    used = set()
    for node, c, _ in sorted(picks, key=lambda x: (x[1] < 0, x[2])):
        nxy = e_xy[node]
        vends = [
            se
            for se in ends_at[node]
            if info[se].state == "V" and info[se].parts
        ]
        if c >= 0 and c not in used:
            used.add(c)
            pos, how = ("n", c), "vn_node"
        else:
            # the outgoing stretch with the highest order, else the incoming
            main = max(vends, key=lambda se: (se[1] == 0, strahler[se[0]]))
            pos = net.nearest_on(info[main].parts, nxy)
            pos = _pull_back(net, pos, nxy, used, junctions=e_deg[node] < 3)
            how = "projection"
        if pos[0] == "n":
            used.add(pos[1])
        node_pos[node] = pos
        moved[node] = net.xy_of(pos)
        node_rows.append(
            (
                Point(moved[node]),
                float(nxy[0]),
                float(nxy[1]),
                float(np.hypot(*(moved[node] - nxy))),
                how,
                len(ends_at[node]),
                len(vends),
            )
        )

    # 6. redrawn stretches follow Vann-Nett between their ends, best matched
    # first, so a line two stretches claim goes to the closer one
    consumed = _Consumed()
    new_geom: dict[int, np.ndarray] = {}
    tags: dict[int, set] = defaultdict(set)
    patch_vn: dict[int, str] = {}
    kept: list = []  # ELVIS line kept where Vann-Nett had a counterpart
    failed_rows, before_rows = [], []
    loose: list[tuple[int, int]] = []

    def quality(s):
        sl = samples[s]
        dd = md[sl][matched[sl]]
        return float(dd.mean()) if len(dd) else np.inf

    def counterpart(s, at):
        sl = samples[s]
        ok = np.flatnonzero(covered[sl] & in_net[np.maximum(cp_part[sl], 0)])
        if not len(ok):
            return None
        o = ok[np.argmin(np.abs(sm[sl][ok] - at))]
        return net.pos(net.local[int(cp_part[sl][o])], float(cp_m[sl][o]))

    def counterpart_within(s, part_set, from_end):
        """(measure, Vann-Nett position) of the covered sample furthest from
        the kept end whose run of counterparts on ``part_set`` reaches that
        end unbroken -- where the stretch's own Vann-Nett line begins."""
        sl = samples[s]
        ok = covered[sl] & in_net[np.maximum(cp_part[sl], 0)]
        loc = np.array([net.local.get(int(p), -1) for p in cp_part[sl]])
        inside = ok & np.isin(loc, list(part_set))
        order = np.flatnonzero(ok)[::-1] if from_end else np.flatnonzero(ok)
        last = None
        for o in order:
            if not inside[o]:
                break
            last = o
        if last is None:
            return None
        return float(sm[sl][last]), net.pos(
            int(loc[last]), float(cp_m[sl][last])
        )

    def end_of_route(s, w, on_route):
        """(Vann-Nett position the route starts/ends at, join): join is the
        point the line must reach, None to attach it later."""
        node = int(e_s[s] if w == 0 else e_e[s])
        if on_route and (s, w) in end_pos:
            return end_pos[(s, w)]
        if on_route and node in node_pos:
            return node_pos[node], moved[node]
        return None, moved.get(node, e_xy[node])

    for s in sorted(redrawn, key=quality):
        line = eg[s]
        a, b = span[s]
        L = sm[samples[s]][-1]
        (A, jA), (B, jB) = (
            end_of_route(s, 0, a <= 0),
            end_of_route(s, 1, b >= L),
        )
        A = A or counterpart(s, a)
        B = B or counterpart(s, b)
        seg = _piece(line, a, b)
        route, why = None, "no anchor"
        if A is not None and B is not None and len(seg) >= 2:
            # Candidate routes, first good one wins. Straight between the
            # mapped ends; failing that -- Vann-Nett meeting at a junction
            # ELVIS doesn't have, two tributaries joining just above the
            # main one, say -- ending (or starting) at the junction the
            # counterpart itself runs into, or at the node nearest the mapped
            # one that the line does reach; last, against digitising order.
            # An end moved off its mapped place is attached where it lands.
            cutoff = 2.5 * (b - a) + 2 * P.node_r
            own_a = info.get((s, 0), _End("E")).own if a <= 0 else -1
            own_b = info.get((s, 1), _End("E")).own if b >= L else -1
            tries = [(A, B, True, a, b)]
            if own_b >= 0 and ("n", own_b) != B:
                tries.append((A, ("n", own_b), True, a, b))
            if own_a >= 0 and ("n", own_a) != A:
                tries.append((("n", own_a), B, True, a, b))
                if own_b >= 0:
                    tries.append((("n", own_a), ("n", own_b), True, a, b))
            fw = net.reach_nearest(A, net.xy_of(B), True, cutoff)
            if fw and np.hypot(*(net.xy[fw[0]] - net.xy_of(B))) <= P.node_r:
                tries.append((A, ("n", fw[0]), True, a, b))
            bw = net.reach_nearest(B, net.xy_of(A), False, cutoff)
            if bw and np.hypot(*(net.xy[bw[0]] - net.xy_of(A))) <= P.node_r:
                tries.append((("n", bw[0]), B, True, a, b))
            # The counterpart can be two Vann-Nett streams: the head of the
            # ELVIS line drawn as part of a neighbour that turns away. Keep
            # to the Vann-Nett line that drains to the end (or is fed from
            # the start), leaving the rest as an ELVIS head (or tail).
            for part_set, keep_end in (
                (net.reach(B, False, cutoff), 1),
                (net.reach(A, True, cutoff), 0),
            ):
                cut = counterpart_within(s, part_set, from_end=keep_end == 1)
                if cut is None:
                    continue
                at, pos = cut
                if keep_end == 1 and a < at < b:
                    tries.append((pos, B, True, at, b))
                elif keep_end == 0 and a < at < b:
                    tries.append((A, pos, True, a, at))
            tries.append((A, B, False, a, b))
            why = "no route"
            for A2, B2, directed, a2, b2 in tries:
                r = net.path(A2, B2, directed, cutoff)
                rc = net.coords(r) if r else np.empty((0, 2))
                seg2 = _piece(line, a2, b2) if (a2, b2) != (a, b) else seg
                if len(rc) < 2 or len(seg2) < 2:
                    continue
                if not _close(seg2, rc, P):
                    why = "strays"
                    continue
                took, who = consumed.overlap(r)
                if took > P.sliver_m:
                    why = "taken by " + ",".join(
                        _key(E, o) for o in sorted(who)
                    )
                    continue
                route, why = r, None
                if A2 != A and a2 <= 0:
                    jA = None
                if B2 != B and b2 >= L:
                    jB = None
                if a2 > 0 and a <= 0:  # an ELVIS head now: it goes to the node
                    jA = moved.get(int(e_s[s]), e_xy[e_s[s]])
                if b2 < L and b >= L:
                    jB = moved.get(int(e_e[s]), e_xy[e_e[s]])
                A, B, a, b = A2, B2, a2, b2
                break
        if why:
            failed_rows.append((line, _key(E, s), b - a, why))
            kept.append(shapely.get_coordinates(line))
            continue
        rc = net.coords(route)
        consumed.add(route, s)
        head, tail = _piece(line, 0, a), _piece(line, b, L)
        for piece in (head, tail):
            if len(piece):
                kept.append(piece)
        coords = _join([head, rc, tail])
        # Each end reaches its node: the Vann-Nett position it was given if
        # the node moved, else where the node stands -- a straight connector
        # from the route, or an ELVIS head/tail pulled there. An end given
        # its counterpart's own junction at a node that stays is attached to
        # whatever line it lands on, later.
        for w, j, on_route in ((0, jA, not len(head)), (1, jB, not len(tail))):
            if on_route and j is None:
                loose.append((s, w))
                continue
            end = coords[0] if w == 0 else coords[-1]
            if np.hypot(*(end - j)) <= 1e-6:
                continue
            if on_route:
                coords = np.vstack([j, coords] if w == 0 else [coords, j])
            else:
                coords = _move_end(coords, w == 0, j)
        new_geom[s] = coords
        tags[s].add(REDRAWN)
        ids = [vn_ids[net.idx[i]] for i, _, _ in route]
        patch_vn[s] = pd.Series(ids).mode().iloc[0] if ids else None
        before_rows.append((line, _key(E, s)))

    # ELVIS lines meeting a moved node that were not redrawn: pull the end
    for node, xy in moved.items():
        for s, w in ends_at[node]:
            if REDRAWN in tags.get(s, ()):
                continue
            c = new_geom.get(s, shapely.get_coordinates(eg[s]))
            new_geom[s] = _move_end(c, w == 0, xy)
            tags[s].add(MOVED)
    collapsed = [s for s, c in new_geom.items() if LineString(c).length < 0.5]
    stats["redrawn_stretches"] = sum(REDRAWN in t for t in tags.values())
    stats["ends_moved_stretches"] = sum(t == {MOVED} for t in tags.values())
    stats["nodes_moved"] = len(moved)
    stats["redraw_failed"] = len(failed_rows)
    # A short ELVIS-only line both of whose nodes landed on the same place
    # (a pond Vann-Nett runs straight through) has no extent left; it is
    # dropped rather than handed on as a line of no length, and listed.
    stats["collapsed"] = len(collapsed)
    collapsed_rows = [(eg[s], _key(E, s), eg[s].length) for s in collapsed]
    if collapsed:
        logger.warning(
            f"{len(collapsed)} stretches shrank to nothing and are dropped: "
            f"{[_key(E, s) for s in collapsed[:10]]}"
        )

    # 7. what is left of Vann-Nett off ELVIS becomes new stretches, unless an
    # ELVIS line that kept its own drawing runs alongside
    # (only ELVIS lines Vann-Nett draws elsewhere: one Vann-Nett copies
    # exactly has its counterpart already, and a Vann-Nett-only tributary
    # running into it is no copy of it)
    # ELVIS lines pulled to a moved node count too: the Vann-Nett line
    # between the same two nodes is the same stretch.
    keep_g = [LineString(c) for c in kept if len(c) >= 2]
    keep_g += [
        LineString(new_geom[s]) if s in new_geom else eg[s]
        for s in cand_e.tolist()
        if s not in redrawn and matched_m.get(s, 0.0) >= P.sliver_m
    ]
    keep_g += [
        LineString(c)
        for s, c in new_geom.items()
        if tags.get(s) == {MOVED}
        and (s in redrawn or matched_m.get(s, 0.0) < P.sliver_m)
    ]
    keep_g = np.array(keep_g, dtype=object)
    add_rows = []  # (global part, m0, m1)
    for p, m0, m1 in vf.itertuples(index=False):
        loc = net.local.get(int(p))
        for lo, hi in (
            consumed.free(loc, m0, m1) if loc is not None else [(m0, m1)]
        ):
            if hi - lo >= P.sliver_m:
                add_rows.append((int(p), lo, hi))
    add_g = np.array(
        [substring(vg[p], a, b) for p, a, b in add_rows], dtype=object
    )
    not_added = np.zeros(len(add_rows), dtype=bool)
    if len(add_rows) and len(keep_g):
        aidx, am, apts = _sample(add_g, P.step)
        atang = _tangent(add_g[aidx], am)
        kj, _, _ = _matches(
            apts, atang, keep_g, shapely.STRtree(keep_g), P.near_m, P.max_angle
        )
        n_s = np.bincount(aidx, minlength=len(add_g))
        along = np.bincount(aidx, weights=kj >= 0, minlength=len(add_g))
        not_added = along / n_s >= P.parallel_frac
    additions = [(add_g[i], add_rows[i][0]) for i in np.flatnonzero(~not_added)]
    stats["vn_added_km"] = float(sum(g.length for g, _ in additions) / 1000)
    stats["vn_not_added_km"] = float(
        shapely.length(add_g[not_added]).sum() / 1000
    )
    logger.info(
        f"Vann-Nett: {len(additions):,} lines added ({stats['vn_added_km']:,.0f} km); "
        f"{stats['vn_not_added_km']:,.0f} km left out alongside ELVIS lines that "
        "kept their own drawing"
    )

    net_out = _assemble(
        E,
        new_geom,
        tags,
        loose,
        additions,
        patch_vn,
        crs,
        regine,
        stats,
        vn_meta=V,
        drop=set(collapsed),
        collapsed=collapsed_rows,
    )
    stats["collapsed"] = len(collapsed_rows)
    qa = {
        "moved_nodes": gpd.GeoDataFrame(
            node_rows,
            columns=[
                "geometry",
                "from_x",
                "from_y",
                "moved_m",
                "how",
                "n_ends",
                "n_redrawn",
            ],
            crs=crs,
        ),
        "redrawn_before": gpd.GeoDataFrame(
            before_rows, columns=["geometry", "stretch"], crs=crs
        ),
        "redraw_failed": gpd.GeoDataFrame(
            failed_rows,
            columns=["geometry", "stretch", "length_m", "why"],
            crs=crs,
        ),
        "collapsed": gpd.GeoDataFrame(
            collapsed_rows, columns=["geometry", "stretch", "length_m"], crs=crs
        ),
        "vn_not_added": gpd.GeoDataFrame(
            {
                VN_ID: [
                    vn_ids[add_rows[i][0]] for i in np.flatnonzero(not_added)
                ]
            },
            geometry=list(add_g[not_added]),
            crs=crs,
        ),
    }
    return Result(net_out, qa, stats)


def _pull_back(net, pos, nxy, taken, junctions):
    """A projection that landed on a Vann-Nett node another ELVIS node
    already has -- or on a junction, for a node that isn't one -- pulled back
    along a part ending there by the node's own distance from it (at most
    half the part), so two ELVIS nodes never share a place and the stretch
    between them keeps a length."""
    if pos[0] == "p":  # within a metre of a node counts as on it
        i, m = pos[1], pos[2]
        if min(m, net.L[i] - m) > 1.0:
            return pos
        pos = ("n", int(net.s[i] if m < net.L[i] - m else net.e[i]))
    c = pos[1]
    if c not in taken and not (junctions and net.deg[c] >= 3):
        return pos
    back = float(np.hypot(*(net.xy[c] - nxy)))
    best = None
    for i in net.inc[c] + net.out[c]:
        step = min(max(back, 1.0), net.L[i] / 2)
        cand = net.pos(i, net.L[i] - step if net.e[i] == c else step)
        dd = float(np.hypot(*(net.xy_of(cand) - nxy)))
        if best is None or dd < best[0]:
            best = (dd, cand)
    return best[1] if best else pos


def _key(E, s):
    return f"{E.at[s, 'elvID']}_{E.at[s, 'strekningLnr']}"


def _close(seg, route, P):
    """Does the Vann-Nett route follow the ELVIS line it replaces?

    Every point of the ELVIS line within near_m + slack of the route, and
    every point of the route that close to the ELVIS line or to where the
    route's ends moved to; the route no shorter than half the line.
    """
    ls, lr = LineString(seg), LineString(route)
    tol = P.near_m + P.slack_m
    moved = max(
        float(np.hypot(*(seg[0] - route[0]))),
        float(np.hypot(*(seg[-1] - route[-1]))),
    )
    if not (0.5 * ls.length <= lr.length <= 2.0 * ls.length + 2 * moved):
        return False
    _, _, pe = _sample(np.array([ls], dtype=object), 5.0)
    _, _, pr = _sample(np.array([lr], dtype=object), 5.0)
    return bool(
        shapely.distance(pe, lr).max() <= tol
        and shapely.distance(pr, ls).max() <= max(tol, moved + P.slack_m)
    )


# ── assembly: attach loose ends, split, attributes ─────────────────────────
def _mm(xy):
    return tuple(np.round(np.asarray(xy) * 1000).astype(np.int64))


def _attach(lines, ends, on_m, snap_m):
    """Close the network at the ends the patch made or moved.

    An end that meets no other line's end goes, in this order: onto another
    line's end within on_m; onto another line within on_m, which is split
    there; onto another line's end within snap_m (a gap the node mapping
    left). Ends meeting nothing that close stay free -- sources, outlets.
    Returns (lines, {line: [points to split it at]}).
    """
    lines = [np.array(c, dtype=float) for c in lines]
    geoms = np.array([LineString(c) for c in lines], dtype=object)
    tree = shapely.STRtree(geoms)
    count = defaultdict(int)  # line ends at each point
    for c in lines:
        count[_mm(c[0])] += 1
        count[_mm(c[-1])] += 1
    splits = defaultdict(list)
    for i, w in ends:
        c = lines[i]
        xy = c[0] if w == 0 else c[-1]
        if count[_mm(xy)] > 1:
            continue
        p = Point(xy)
        cand = [
            j
            for j in tree.query(p, predicate="dwithin", distance=snap_m)
            if j != i
        ]
        best_end = None
        for j in cand:
            for e in (lines[j][0], lines[j][-1]):
                dd = float(np.hypot(*(e - xy)))
                if dd <= snap_m and (best_end is None or dd < best_end[0]):
                    best_end = (dd, e.copy())
        on_line = [j for j in cand if geoms[j].distance(p) <= on_m]
        if best_end is not None and best_end[0] <= on_m:
            new = best_end[1]
        elif on_line:
            j = min(on_line, key=lambda j: geoms[j].distance(p))
            new = np.array(geoms[j].interpolate(geoms[j].project(p)).coords[0])
            splits[j].append(new)
        elif best_end is not None:
            new = best_end[1]
        else:
            continue
        count[_mm(xy)] -= 1
        count[_mm(new)] += 1
        if w == 0:
            c[0] = new
        else:
            c[-1] = new
    return lines, splits


def _cut(coords, points):
    """The line cut at the given points, each piece meeting the next exactly
    at that point."""
    line = LineString(coords)
    at = sorted(
        (float(line.project(Point(xy))), tuple(xy))
        for xy in {tuple(x) for x in points}
    )
    at = [(m, xy) for m, xy in at if 1e-6 < m < line.length - 1e-6]
    if not at:
        return [coords]
    edges = [(0.0, None), *at, (line.length, None)]
    out = []
    for (a, pa), (b, pb) in zip(edges[:-1], edges[1:], strict=True):
        if b - a <= 1e-6:
            continue
        piece = shapely.get_coordinates(substring(line, a, b))
        if pa is not None:
            piece[0] = pa
        if pb is not None:
            piece[-1] = pb
        if len(piece) >= 2:
            out.append(piece)
    return out


def _strahler(geoms, fixed):
    """Strahler order for rows where ``fixed`` is NaN, from the network.

    Upstream of a line = lines ending where it starts. Rows already ordered
    (ELVIS) are taken as given; the rest are ordered upstream first. A line
    caught in a loop, or fed by one, takes the largest order known upstream,
    or 1.
    """
    order = np.array(fixed, dtype=float)
    todo = np.flatnonzero(np.isnan(order))
    if not len(todo):
        return order
    a, b = _ends(geoms)
    (s, e), _ = _node_ids(a, b)
    ending = defaultdict(list)
    for i, v in enumerate(e):
        ending[v].append(i)
    starting = defaultdict(list)
    for i in todo:
        starting[s[i]].append(i)
    up = {int(i): [j for j in ending[s[i]] if j != i] for i in todo}
    waiting = {i: sum(np.isnan(order[j]) for j in ups) for i, ups in up.items()}
    ready = [i for i, n in waiting.items() if n == 0]
    while ready:
        i = ready.pop()
        vals = [order[j] for j in up[i]]
        mx = max(vals, default=0.0)
        order[i] = 1.0 if not vals else mx + 1 if vals.count(mx) >= 2 else mx
        for k in starting[e[i]]:
            if k != i and np.isnan(order[k]):
                waiting[k] -= 1
                if waiting[k] == 0:
                    ready.append(int(k))
    for i in np.flatnonzero(np.isnan(order)):  # loops
        known = [order[j] for j in up.get(int(i), []) if not np.isnan(order[j])]
        order[i] = max(known, default=1.0)
    return order


def _assemble(
    E,
    new_geom,
    tags,
    loose,
    additions,
    patch_vn,
    crs,
    regine,
    stats,
    vn_meta=None,
    on_m=1.0,
    snap_m=5.0,
    drop=frozenset(),
    collapsed=None,
):
    nE = len(E)
    collapsed = [] if collapsed is None else collapsed
    lines = [
        new_geom[s] if s in new_geom else shapely.get_coordinates(g)
        for s, g in enumerate(E.geometry.values)
    ]
    lines += [shapely.get_coordinates(g) for g, _ in additions]
    # every end the patch made or moved: loose route ends, both ends of each
    # changed ELVIS line and of each added one
    ends = list(loose) + [(i, w) for i in new_geom for w in (0, 1)]
    ends += [(nE + k, w) for k in range(len(additions)) for w in (0, 1)]
    lines, splits = _attach(lines, ends, on_m, snap_m)
    stats["lines_split"] = len(splits)

    geoms, kinds, src = [], [], []
    for i, c in enumerate(lines):
        if i in drop:
            continue
        if i < nE and i in new_geom and LineString(c).length < 0.5:
            # shrank to nothing in the attaching (see patch: collapsed)
            collapsed.append(
                (
                    E.geometry.values[i],
                    _key(E, i),
                    float(E.geometry.values[i].length),
                )
            )
            continue
        pieces = _cut(c, splits.get(i, []))
        for piece in pieces:
            geoms.append(LineString(piece))
            src.append(i)
            kinds.append(len(pieces) > 1)

    geoms = np.array(geoms, dtype=object)
    src = np.array(src)
    is_add = src >= nE

    # ELVIS rows, attributes untouched
    cols = [c for c in E.columns if c != "geometry"]
    out = pd.DataFrame(index=range(len(geoms)), columns=cols, dtype=object)
    if (~is_add).any():
        out.loc[~is_add, cols] = E.iloc[src[~is_add]][cols].to_numpy()
    tag = [set(tags.get(int(i), ())) if i < nE else {ADDED} for i in src]
    for t, split in zip(tag, kinds, strict=True):
        if split:
            t.add(SPLIT)

    # additions: identity, order, REGINE unit
    n_add = int(is_add.sum())
    if n_add:
        part = np.array([additions[i - nE][1] for i in src[is_add]])
        vnid = vn_meta[VN_ID].to_numpy()[part] if vn_meta is not None else None
        for c in cols:
            out[c] = out[c].astype("object")
        add_rows = np.flatnonzero(is_add)
        out.loc[add_rows, cols] = None
        out.loc[add_rows, "objektType"] = "ElvBekk"
        out.loc[add_rows, "elvID"] = [f"VN-{v}" for v in vnid]
        out.loc[add_rows, "strekningLnr"] = [
            f"VN{k:07d}" for k in range(1, n_add + 1)
        ]
        if regine is not None and len(regine):
            mid = shapely.line_interpolate_point(
                geoms[is_add], 0.5, normalized=True
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
            out.loc[add_rows, "vassdragsnummer"] = hit[
                "vassdragsnummer"
            ].to_numpy()
            out.loc[add_rows, "elvNavnHierarki"] = hit[
                "elvNavnHierarki"
            ].to_numpy()
        fixed = np.array(
            pd.to_numeric(out.get("elveordenStrahler"), errors="coerce"),
            dtype=float,
        )
        fixed[is_add] = np.nan
        out["elveordenStrahler"] = _strahler(geoms, fixed)
        out.loc[add_rows, PATCH_VN] = vnid
    if PATCH_VN not in out:
        out[PATCH_VN] = None
    for s, v in patch_vn.items():
        out.loc[(src == s) & ~is_add, PATCH_VN] = v
    out[PATCH] = [",".join(sorted(t)) if t else None for t in tag]
    # ELVIS's two length fields are the line's length; keep them so on every
    # row whose line changed, and leave the rest as delivered
    changed = np.array([bool(t) for t in tag])
    L = shapely.length(geoms[changed])
    if "SHAPE_Length" in out:
        out.loc[changed, "SHAPE_Length"] = L
    if "elvelengde_m" in out:
        out.loc[changed, "elvelengde_m"] = np.round(L, 2)
    out = _restore_dtypes(out, E)
    return gpd.GeoDataFrame(
        out,
        geometry=[MultiLineString([g]) for g in geoms],
        crs=crs,
    )


def _restore_dtypes(out, E):
    """ELVIS's own column types back, after rows were built as objects."""
    if not len(E):  # nothing to take them from
        return out
    for c in out.columns:
        if c not in E.columns:
            continue
        dt = E[c].dtype
        if pd.api.types.is_integer_dtype(dt):
            out[c] = (
                pd.to_numeric(out[c], errors="coerce").round().astype("Int64")
            )
        elif pd.api.types.is_float_dtype(dt):
            out[c] = pd.to_numeric(out[c], errors="coerce")
        elif pd.api.types.is_datetime64_any_dtype(dt):
            out[c] = pd.to_datetime(out[c], errors="coerce")
    return out
