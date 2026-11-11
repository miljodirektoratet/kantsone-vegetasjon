"""Kantsoner via a catchment-wide partition of the water between ELVIS stretches.

ELVIS (Vann-Nett) is the ownership and topology system; the water polygons are
used only to find where the physical channel and its banks are. Per stretch:

    anchor    both ends onto the local medial axis inside the water, not the
              nearest bank -- a bank anchor drags the route along that side
    route     least-cost path between the anchors, cheapest along the medial
              axis, expensive on land; routes keep their full length, so
              two stretches meeting at a confluence both run a little way
              into the shared water below it
    cell      catchment-wide Voronoi over the routes, each cell grown by
              CELL_OVERLAP_M; one cell per stretch, and neighbouring
              stretches may overlap up to that margin at their boundary
    channel   the water polygons the route actually enters, inside the cell
    band      that water buffered out to BAND_M, minus the water itself
    claim     every point of the route assigned to the nearest 10 m ELVIS
              segment -- the topological link between route and delivery key
    slice     Voronoi over the claims n band, split left/right by the route

The 10 m segments are cut on the ELVIS line and carry the delivery keys;
everything else is derived from the route. `build` is pure computation --
reading sources and writing results belongs to the stage.
"""

from __future__ import annotations

import multiprocessing as mp
import os
from dataclasses import dataclass, field

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from loguru import logger
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.windows import from_bounds
from scipy.ndimage import distance_transform_edt, zoom
from shapely.geometry import LineString, MultiPoint, Point
from shapely.ops import (
    nearest_points,
    substring,
    unary_union,
    voronoi_diagram,
)
from skimage.graph import route_through_array

from .. import membudget
from ..spec import constants as C
from . import flowacc

RES_M = 1.0  # cost-surface pixel; sets the route's precision only
CORRIDOR_M = 250.0  # how far off the ELVIS line the route may wander
# ── elevation (see _elevation_walls) ──────────────────────────────────────
# Water does not climb, and a route that climbs has left its channel: over a
# ridge into the next creek, or up a tributary from a confluence. Between two
# fixed anchors every route must lose the same height, so charging for total
# elevation change -- up and down alike -- charges exactly the climbing and
# nothing else, and needs neither the flow direction nor the path history.
# Realised as contour walls: the DTM is read coarse, quantised into
# DTM_LEVEL_M steps, and the pixels where the level changes cost extra. A
# path pays per level crossed; one running along a hillside crosses none.
# (A slope-magnitude cost would instead charge the valley side under a creek
# -- 19% median at 10 m in Lærdal against the creek's own 5-10% fall.)
DTM_COARSE_M = 10.0  # read resolution; served from the COG's overviews
DTM_LEVEL_M = 2.0  # contour interval; above the DTM's water-surface noise
# Cost of one wall pixel. Crossing a level runs through ~2 wall pixels, so
# one DTM_LEVEL_M step costs ~2x this. Off: measured on Lærdal at 500 (a
# metre of climb = a metre of land) the walls moved 247 routes by > 20 m,
# 92 for the worse -- in-water sidesteps around contours the smoothed field
# paints across river channels -- against 88 for the better, and the detour
# statistics did not move. Kept for the pure-climb case; the topographic
# position term below carries the terrain now.
DTM_WALL_COST = 0.0
# ── DTM channels (see _flow_channels) ────────────────────────────────────
# Where the DTM says water runs but nothing is mapped, the route may treat
# the ground as water: a D8 flow accumulation over the routing window (plus
# a margin, so the window edge does not cut a channel's catchment short) at
# the overview's own resolution, and every dry pixel above the threshold
# priced like the shallowest mapped water. Window-local accumulation is a
# floor on the true value, so the threshold is small and speaks of what
# gathers *nearby*: a creek, a gully, the reach inside the window.
FLOWACC_RES_M = (
    8.0  # the DTM's 8x overview, read natively (no resampling stripes)
)
FLOWACC_PAD_M = 50.0
FLOWACC_CHANNEL_HA = (
    4.0  # the 8 m D8 web on a planar hillside sits well below this
)
FLOWACC_CHANNEL_COST = 20.0  # what a mapped-water pixel at the bank costs
# Only channels this close to the stretch's own ELVIS line count. The term
# is there to find the stretch's *own* unmapped channel; priced everywhere
# in the corridor it became bait, and a headwater stretch whose gully
# gathers too little to clear the threshold looped 500-850 m through a
# neighbouring gully or lakeshore rather than walk 300 m of land along
# its line (45 routes on Lærdal, 3x their length at worst).
FLOWACC_NEAR_LINE_M = 50.0
ANCHOR_SEARCH_M = 60.0  # radius to look in for a medial-axis anchor
ANCHOR_DEPTH_WEIGHT = 2.0  # metres of detour a metre of extra depth is worth
SMOOTH_WINDOW_M = 7.0  # moving-average window taking the lattice off the route
# Per pixel of dry ground, against ~1 for the deepest water and
# FLOWACC_CHANNEL_COST where the DTM drains a channel. Halving it (with the
# contour walls at the same price) moved 79 routes on Lærdal, 7 of them by
# more than 20 m; there is no reason for it now the walls are off.
LAND_COST = 1000.0
SEGMENT_M = 10.0
# A zone is one segment's share of one bank, so it should be about as long as
# the segment's claim on the route and BAND_M wide. Anything past this
# multiple of that is not a bank slice -- it is a cell that reached into a
# pool or across a confluence. Measured against the claim rather than a flat
# SEGMENT_M: nearest-segment assignment hands out claims from under a metre
# to 74 m, and a fixed cap either throttles the long ones or lets the short
# ones run away.
MAX_ZONE_FACTOR = 5.0
CLAIM_SAMPLE_M = 1.0  # spacing at which the route is assigned + cells sampled
# A Voronoi assigns everything, so without a cap a lake far from any ELVIS
# line is adopted by whichever route happens to be nearest and gets a band.
MAX_OWN_DIST_M = 1000.0
# A zone that comes out in several disconnected parts is one bank slice
# plus, sometimes, debris the cell picked up elsewhere -- across a
# confluence, behind an islet. A part carrying at least this fraction of
# the claim's full band area (claim length x BAND_M) is a candidate bank
# slice; anything under it is debris outright. A fraction rather than a
# flat area for the same reason as MAX_ZONE_FACTOR: claims run from under
# a metre to 74 m, and a fixed floor either swallows the short claims'
# zones whole or lets debris through beside the long ones.
ZONE_PART_MIN_FACTOR = 0.2
# Parts that clear the floor are then judged on what separates them -- see
# _merge_parts. Keeping only the nearest, as this did before, threw away
# the far half of every bank a mapped side-stream happened to cross:
# measured over 150 Lærdal stretches, 3.3% of segment-sides split this way
# and 60% of their area was being discarded, usually the larger half.
#
# The connector between two parts must be water over at least this much of
# its length for them to count as one bank interrupted, rather than two
# places the cell reached separately.
WATER_BRIDGE_MIN_FRAC = 0.5
# Width of the probe along that connector. Thin on purpose: it should pick
# up the water actually blocking the gap, not the main channel running
# alongside it.
BRIDGE_PROBE_M = 0.5
# Each stretch's cell grown past the strict Voronoi boundary. A hard
# boundary cut zones off mid-band wherever it fell within BAND_M of the
# water -- most visibly along a wide lake, where the cell edge lands in
# open water and the far shore's band was clipped to nothing. The width of
# the widest band: enough that a band reaching the boundary can complete,
# and neighbouring stretches' zones overlap at most that far.
CELL_OVERLAP_M = 20.0
# ── the no-mapped-water fallback (see _zone_fallback) ─────────────────────
# Half-width of the channel assumed where a real watercourse has no mapped
# water to measure a bank from. Measured over Lærdal: of the real-watercourse
# segments that end up with no usable zone, 91% are Strahler order 1-2 and
# 61% have no mapped water within 20 m of the line at all -- headwater
# streams grunnkart never drew. 1 m each side (a 2 m channel) is the scale
# of those. It would be wrong for a large river, but large rivers are not
# what lands here.
FALLBACK_HALF_WIDTH_M = 1.0
# Only real watercourses get the fallback. A lake midline does not: it has
# no riparian zone within the widest band because it is open water, and the
# kantsone there belongs to the shoreline, which the lake track measures.
# The taxonomy, and why ElvelinjeFiktiv is now included and ElvBekkRegulert
# is not, is in spec/constants.py.
FALLBACK_OBJEKTTYPER = C.RIVER_OBJEKTTYPER
# ── route capture (see _captured) ─────────────────────────────────────────
# Two routes this close are on the same water. Some of that is by design --
# stretches meeting at a confluence both run a little way into the shared
# channel below it -- so coincidence alone means nothing; over Lærdal the
# median route shares 10 m with a neighbour.
ROUTE_COINCIDENT_M = 5.0
# Past this much shared length it is not a confluence any more. p90 of the
# shared length is 52 m, p99 393 m, and the worst observed was a Strahler-2
# stretch whose 1,859 m line produced a 2,442 m route running 537 m down
# the neighbouring river, delivering its kantsone a median 60 m -- and up
# to 228 m -- from its own line.
CAPTURE_MIN_M = 50.0
# Half-baseline the route's direction is measured over when deciding which
# bank is which. The route is a pixel path, so a short baseline reads lattice
# jitter rather than the channel; too long and it stops following real bends.
# Swept against the flip rate on Hurumlandet: 5 m gave 4.23%, 35 m 3.09%,
# 120 m 5.17%.
TANGENT_M = 35.0
# Bearings tried for a route's extension ray before giving up on finding
# one that clears the route (see _ray): 16 is 22.5 degree steps.
RAY_STEPS = 16
# Spacing of the side votes cast along the route (see _sides). Zones are
# ~10 m long, so 5 m guarantees every piece that borders the route votes.
BANK_VOTE_M = 5.0


# ── route ─────────────────────────────────────────────────────────────────
def _anchor(point, half_width, minx, maxy, shape, lake=None):
    """Pixel near `point` that trades nearness off against depth.

    The nearest water pixel is usually on the bank, and a route starting
    there hugs it, so it is worth stepping inwards towards the local medial
    axis. Ranking on depth alone is too strong a pull: it takes the deepest
    pixel anywhere in the window, so an endpoint already sitting in its own
    creek gets dragged up to 85 m onto whatever neighbouring channel happens
    to be half a metre wider. Scoring `distance - WEIGHT * depth` keeps the
    step inwards and rejects the jump: crossing to another channel has to
    buy 1/WEIGHT metres of depth per metre travelled. Where no water is in
    reach the endpoint itself is returned and the route pays the land cost.

    `lake` masks pixels inside lake polygons: those earn no depth reward,
    so an endpoint at or near a lake snaps to the *closest* lake pixel
    instead of being dragged toward deeper open water. Stepping to a river
    channel's medial axis is the point of the reward; a lake has no channel
    to centre on.
    """
    h, w = shape
    row = int((maxy - point.y) / RES_M)
    col = int((point.x - minx) / RES_M)
    r = int(ANCHOR_SEARCH_M / RES_M)
    r0, r1 = max(row - r, 0), min(row + r + 1, h)
    c0, c1 = max(col - r, 0), min(col + r + 1, w)
    win = half_width[r0:r1, c0:c1]
    if win.size and win.any():
        rr, cc = np.nonzero(win)
        reach = np.hypot(rr + r0 - row, cc + c0 - col) * RES_M
        depth = win[rr, cc]
        if lake is not None:
            depth = np.where(lake[r0:r1, c0:c1][rr, cc], 0.0, depth)
        best = np.argmin(reach - ANCHOR_DEPTH_WEIGHT * depth)
        return int(rr[best] + r0), int(cc[best] + c0)
    return (row, col) if 0 <= row < h and 0 <= col < w else None


def _smooth(line):
    """Moving average along the path, both ends pinned.

    route_through_array steps pixel to pixel, so the path arrives as a
    staircase of 1 m moves in eight directions. The Voronoi generators sit
    on that line and their cell boundaries stand square to it, so every
    lattice step turns into a kink in a zone boundary. Averaging over a
    window a few metres long takes the staircase out without moving the line
    off the channel: the window is short against the corridor and the ends
    are held exactly where the anchors put them.
    """
    n = int(SMOOTH_WINDOW_M / RES_M) | 1
    xy = np.asarray(line.coords)
    if len(xy) <= n:
        return line
    pad, k = n // 2, np.ones(n) / n
    out = np.empty_like(xy)
    for d in (0, 1):
        held = np.r_[np.full(pad, xy[0, d]), xy[:, d], np.full(pad, xy[-1, d])]
        out[:, d] = np.convolve(held, k, mode="valid")
    out[0], out[-1] = xy[0], xy[-1]
    return LineString(out)


def _coarse_elevation(dtm, minx, miny, maxx, maxy, h, w):
    """The window's elevation at DTM_COARSE_M, or None where it is all
    nodata.

    Read through the raster's overviews, so no full-resolution pixels move.
    Nodata -- sea, foreign ground -- is flattened to the window mean, so it
    neither walls off nor attracts.
    """
    hc, wc = (
        max(2, int(np.ceil(h * RES_M / DTM_COARSE_M))),
        max(2, int(np.ceil(w * RES_M / DTM_COARSE_M))),
    )
    z = dtm.read(
        1,
        window=from_bounds(minx, miny, maxx, maxy, dtm.transform),
        out_shape=(hc, wc),
        resampling=Resampling.average,
        boundless=True,
        fill_value=dtm.nodata,
    ).astype("float64")
    bad = ~np.isfinite(z)
    if dtm.nodata is not None:
        bad |= z <= dtm.nodata + 1.0
    if bad.all():
        return None
    z[bad] = z[~bad].mean()
    return z


def _upsample(a, h, w):
    """A coarse window array interpolated onto the (h, w) cost grid."""
    hc, wc = a.shape
    out = zoom(a, (h / hc, w / wc), order=1)[:h, :w]
    if out.shape != (h, w):
        out = np.pad(
            out, ((0, h - out.shape[0]), (0, w - out.shape[1])), mode="edge"
        )
    return out


def _flow_channels(dtm, minx, miny, maxx, maxy, h, w):
    """Boolean (h, w) mask of dry pixels the DTM drains enough water
    through to count as a channel -- see the FLOWACC_* constants."""
    pad = FLOWACC_PAD_M
    px0, py0, px1, py1 = minx - pad, miny - pad, maxx + pad, maxy + pad
    hc = max(3, int(round((py1 - py0) / FLOWACC_RES_M)))
    wc = max(3, int(round((px1 - px0) / FLOWACC_RES_M)))
    z = dtm.read(
        1,
        window=from_bounds(px0, py0, px1, py1, dtm.transform),
        out_shape=(hc, wc),
        resampling=Resampling.average,
        boundless=True,
        fill_value=dtm.nodata,
    ).astype("float64")
    bad = ~np.isfinite(z)
    if dtm.nodata is not None:
        bad |= z <= dtm.nodata + 1.0
    if bad.all():
        return None
    # Nodata is sea or foreign ground: an outlet, so it sits below every
    # valid cell. A high wall would be wrong twice over -- the priority
    # flood fills depressions relative to the window edge, so walls around
    # the edge turn the whole window into one flat lake that drains in a
    # spiral to the border, and no channel survives.
    z[bad] = z[~bad].min() - 1000.0
    acc = flowacc.accumulation(z, FLOWACC_RES_M)
    channel = acc >= FLOWACC_CHANNEL_HA * 1e4
    # Back to the unpadded window on the cost grid: pad cells off, then
    # nearest-neighbour up to RES_M.
    r0 = int(round(pad / FLOWACC_RES_M))
    c0 = int(round(pad / FLOWACC_RES_M))
    core = channel[
        r0 : r0 + int(round(h * RES_M / FLOWACC_RES_M)) + 1,
        c0 : c0 + int(round(w * RES_M / FLOWACC_RES_M)) + 1,
    ]
    up = (
        zoom(
            core.astype("float32"),
            (h / core.shape[0], w / core.shape[1]),
            order=0,
        )[:h, :w]
        > 0.5
    )
    if up.shape != (h, w):
        up = np.pad(
            up, ((0, h - up.shape[0]), (0, w - up.shape[1])), mode="edge"
        )
    return up


def _elevation_walls(zc, h, w):
    """Boolean (h, w) mask of the pixels where the quantised elevation
    changes -- the contour walls the route pays to cross. See the
    DTM_* constants for why walls rather than slope; off by default.
    """
    z = _upsample(zc, h, w)
    level = np.floor(z / DTM_LEVEL_M).astype(np.int32)
    walls = np.zeros((h, w), dtype=bool)
    dy = level[1:, :] != level[:-1, :]
    dx = level[:, 1:] != level[:, :-1]
    walls[1:, :] |= dy
    walls[:-1, :] |= dy
    walls[:, 1:] |= dx
    walls[:, :-1] |= dx
    return walls


def _route(line, water_parts, lake_parts=None, dtm=None):
    """Least-cost path down the channel between the line's two ends.

    `lake_parts` are the lake polygons near this stretch: pixels inside
    them get no depth reward when anchoring the endpoints (see _anchor).
    The cost surface itself is untouched -- a route crossing a lake still
    follows the medial axis like it always did. `dtm` is an open elevation
    raster; with one, dry ground the DTM drains a channel through is priced
    as water (_flow_channels) and, if DTM_WALL_COST is on, the contour walls
    of _elevation_walls are added.
    """
    minx, miny, maxx, maxy = line.buffer(CORRIDOR_M).bounds
    # Snap the window to a global RES_M grid. Every stretch then samples the
    # same pixels, so two stretches meeting at an ELVIS node see an identical
    # depth field there and anchor on the identical pixel -- with per-window
    # origins the same node landed 10-35 m apart from one stretch to the
    # next.
    minx = np.floor(minx / RES_M) * RES_M
    miny = np.floor(miny / RES_M) * RES_M
    maxx = np.ceil(maxx / RES_M) * RES_M
    maxy = np.ceil(maxy / RES_M) * RES_M
    w = int(round((maxx - minx) / RES_M))
    h = int(round((maxy - miny) / RES_M))
    if w < 3 or h < 3 or w * h > 60_000_000:
        return None
    burn = {
        "out_shape": (h, w),
        "transform": from_origin(minx, maxy, RES_M, RES_M),
        "fill": 0,
        "all_touched": True,
    }
    # Clip each polygon to the window before burning. The wet mask only ever
    # holds in-window pixels, so geometry beyond the window is pure cost --
    # a main-stem polygon touching thousands of corridors was re-rasterized
    # in full for each one. The clip cannot change which pixels burn.
    clipped = [
        g
        for geom in shapely.clip_by_rect(water_parts, minx, miny, maxx, maxy)
        for g in shapely.get_parts(shapely.make_valid(geom))
        if g.geom_type == "Polygon" and not g.is_empty
    ]
    if not clipped:
        return None
    wet = rasterize(clipped, default_value=1, **burn).astype(bool)
    if not wet.any():
        return None
    corridor = rasterize(
        [line.buffer(CORRIDOR_M)], default_value=1, **burn
    ).astype(bool)

    # Cheapest along the medial axis. With every wet pixel costing the same
    # the path was free to meander inside the channel and came out up to 43%
    # longer than the ELVIS line it should track.
    half_width = distance_transform_edt(wet) * RES_M
    depth = half_width / max(half_width.max(), 1e-6)
    cost = np.where(wet, 1.0 / (0.05 + depth), LAND_COST)
    if dtm is not None:
        zc = _coarse_elevation(dtm, minx, miny, maxx, maxy, h, w)
        if zc is not None:
            if DTM_WALL_COST > 0.0:
                cost[_elevation_walls(zc, h, w)] += DTM_WALL_COST
        if FLOWACC_CHANNEL_HA > 0.0:
            channel = _flow_channels(dtm, minx, miny, maxx, maxy, h, w)
            if channel is not None:
                near_line = rasterize(
                    [line.buffer(FLOWACC_NEAR_LINE_M)], default_value=1, **burn
                ).astype(bool)
                # Only dry ground is repriced; mapped water keeps its depth
                # cost, and anchoring (half_width) never sees this.
                cost[channel & near_line & ~wet] = FLOWACC_CHANNEL_COST
    cost[~corridor] = LAND_COST * 50

    lake = None
    if lake_parts is not None and len(lake_parts):
        lclip = [
            g
            for geom in shapely.clip_by_rect(lake_parts, minx, miny, maxx, maxy)
            for g in shapely.get_parts(shapely.make_valid(geom))
            if g.geom_type == "Polygon" and not g.is_empty
        ]
        if lclip:
            lake = rasterize(lclip, default_value=1, **burn).astype(bool)

    ends = [Point(line.coords[0]), Point(line.coords[-1])]
    rc = [_anchor(p, half_width, minx, maxy, (h, w), lake=lake) for p in ends]
    if any(a is None for a in rc) or rc[0] == rc[1]:
        return None
    path, _ = route_through_array(
        cost, rc[0], rc[1], fully_connected=True, geometric=True
    )
    if len(path) < 2:
        return None
    pts = [
        (minx + (c + 0.5) * RES_M, maxy - (r + 0.5) * RES_M) for r, c in path
    ]
    return _smooth(LineString(pts)).simplify(0.2, preserve_topology=False)


def _areal(geom):
    """The polygonal part of `geom`, and nothing else.

    difference and intersection leave a stray line or point wherever two
    boundaries graze each other. One such row is enough to make a written
    layer mixed-type, and gpd.overlay -- which every extractor runs against
    the rings -- refuses that outright with "df2 contains mixed geometry
    types". Osen and Laerdal both got through geometry and died in
    attributes on a single bad row in 3,988.

    Also the validity gate: difference and intersection can emit a polygon
    with a broken ring ('Ring edge missing'), and one written invalid zone
    makes every GEOS op in attributes throw. Repaired here, at the source,
    so everything downstream can assume valid input.
    """
    if geom.is_empty:
        return geom
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    if geom.geom_type in ("Polygon", "MultiPolygon"):
        return geom
    parts = [
        g
        for g in shapely.get_parts(geom)
        if g.geom_type in ("Polygon", "MultiPolygon")
    ]
    if not parts:
        return shapely.Polygon()
    return parts[0] if len(parts) == 1 else unary_union(parts)


def _cells(points, envelope):
    """Voronoi cells over `points`, repaired.

    GEOS emits the occasional self-intersecting cell -- 6 of 67,039 on
    Tysnes -- and one of those is enough to break every union downstream.
    Seven of thirteen catchments died on it, reporting "non-noded
    intersection", "side location conflict" and "unable to assign free hole
    to a shell": four GEOS messages, one cause, and the reported coordinate
    is the self-intersection itself. Hurumlandet passed only by being small
    enough that the diagram came out clean.

    Not duplicate sites: the inputs were checked and are all distinct.
    """
    out = []
    for cell in voronoi_diagram(MultiPoint(points), envelope=envelope).geoms:
        if not cell.is_valid:
            parts = [
                g
                for g in shapely.get_parts(shapely.make_valid(cell))
                if g.geom_type in ("Polygon", "MultiPolygon")
            ]
            if not parts:
                continue
            cell = parts[0] if len(parts) == 1 else unary_union(parts)
        out.append(cell)
    return out


def _dissolve_cells(cells):
    """Merge the Voronoi cells belonging to one owner.

    They tile without overlapping and share their edges exactly -- a
    polygonal coverage -- so they dissolve by dropping shared edges, with no
    intersection computed anywhere. That is both faster than a general union
    and far harder to break. Repairs above can leave a cell overlapping its
    neighbour by a sliver, so a general union stays as the fallback.
    """
    if len(cells) == 1:
        return cells[0]
    try:
        return shapely.coverage_union_all(np.array(cells, dtype=object))
    except shapely.errors.GEOSException:
        return unary_union(cells)


def _partition(routes, extent):
    """{stretch_id: cell}, the plane divided by nearest route.

    Voronoi over every route vertex, then the cells belonging to one route
    dissolved back together and grown by CELL_OVERLAP_M. The growth means
    the cells are no longer disjoint: two stretches' zones may overlap up
    to that far into each other's territory, which is accepted so a band
    reaching the boundary completes instead of being clipped mid-band.
    Euclidean nearest: where two channels run parallel with land between,
    a point can fall to the wrong one, which is worth measuring but rare
    enough to start here.
    """
    pts, owner = [], []
    for sid, route in routes.items():
        for x, y in shapely.get_coordinates(route):
            pts.append(shapely.Point(x, y))
            owner.append(sid)
    if len(pts) < 3:
        return {}
    env = shapely.box(*extent).buffer(MAX_OWN_DIST_M)
    cells = _cells(pts, env)
    tree = shapely.STRtree(np.array(pts, dtype=object))
    groups = {}
    for cell in cells:
        # By nearest generator, not contains: a degenerate diagram from
        # near-collinear generators puts a point on its own cell's boundary,
        # and a contains test silently drops that cell.
        groups.setdefault(
            owner[int(tree.nearest(cell.representative_point()))], []
        ).append(cell)
    return {
        k: _dissolve_cells(v).buffer(CELL_OVERLAP_M) for k, v in groups.items()
    }


def _valid_retry(fn, *geoms):
    """fn(*geoms); on a GEOS topology error, repair the inputs and retry.

    Source polygons are repaired at read, but GEOS overlay can still throw
    'side location conflict' on rare near-degenerate combinations (seen in
    Sira og Sokna). make_valid is a no-op on valid geometry, so retrying
    with repaired inputs changes nothing on the paths that already worked.
    """
    try:
        return fn(*geoms)
    except shapely.errors.GEOSException:
        logger.warning("GEOS topology error; retrying with repaired inputs")
        return fn(*(shapely.make_valid(g) for g in geoms))


# ── the water this stretch owns ───────────────────────────────────────────
def _channel(route, cell, parts, tree):
    """The water this stretch owns.

    Three conditions, and all are needed:

      the route enters the polygon -- a side creek joined at a confluence is
          part of the same FKB polygon but the route never runs down it, and
          without this the cell clips a wedge out of it and the band follows
          that wedge, leaving spikes off the side of the ribbon
      it lies inside this stretch's cell -- which is what stops two stretches
          claiming the same water where their routes meet
      it lies within MAX_OWN_DIST_M of the route -- a Voronoi assigns
          everything, so without a cap a distant lake is adopted by whichever
          route happens to be nearest

    And one thing deliberately absent: no invented water. An earlier
    fallback buffered the route itself into a channel when less than 20% of
    it was wet, meant for short unmapped gaps -- but it fired per stretch,
    so a route 83% on dry land grew a full-length fabricated channel and
    1.5 ha of kantsoner through forest. Bands come from mapped water --
    river polygons, creek polygons, lakes -- or not at all; a route section
    over land yields no kantsone.
    """
    hit = tree.query(route.buffer(RES_M * 2), predicate="intersects")
    if len(hit) == 0:
        return None
    local = _valid_retry(
        shapely.intersection, _valid_retry(unary_union, parts[hit]), cell
    )
    if local.is_empty:
        return None
    local = _valid_retry(
        shapely.intersection, local, route.buffer(MAX_OWN_DIST_M)
    )
    if local.is_empty:
        return None
    return _areal(local)


# ── segments, ownership, sides ────────────────────────────────────────────
def _segments(line, stretch_id):
    n = max(1, round(line.length / SEGMENT_M))
    out = []
    for i in range(n):
        piece = substring(line, i / n, (i + 1) / n, normalized=True)
        if not piece.is_empty and piece.length > 0:
            out.append((f"{stretch_id}_{i:05d}", piece))
    return out


def _claim(route, line, segs):
    """{segment_id: the piece(s) of route nearest to it}.

    The segments carry the delivery keys and so must stay ELVIS's, but the
    band is built around the route, and the two lines diverge. Cutting the
    band with cells generated from the ELVIS segments put the slicing points
    beside the thing being sliced, and 4% of segments came out with only one
    bank.

    So every point of the route goes to the segment it is nearest to. The
    nearest segment is the one holding the nearest point on the whole line,
    which `project` gives directly. Claims come out unequal in length -- a
    10 m segment takes a median 10 m of route but up to 74 m where the route
    swings away from the ELVIS line and one segment is nearest for a long
    run -- and the boundaries between them sit midway between samples so
    they still tile the route without seams.
    """
    n = len(segs)
    if n < 2:
        return {segs[0][0]: route} if segs else {}
    k = max(2, int(route.length / CLAIM_SAMPLE_M) + 1)
    at = np.linspace(0.0, route.length, k)
    seg_len = line.length / n
    idx = np.clip(
        [int(line.project(route.interpolate(d)) / seg_len) for d in at],
        0,
        n - 1,
    )
    cut = np.flatnonzero(idx[:-1] != idx[1:])
    edges = np.r_[0.0, (at[cut] + at[cut + 1]) / 2, route.length]
    runs = {}
    for e, (a, b) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        if b - a > 1e-9:
            runs.setdefault(
                idx[0] if e == 0 else idx[cut[e - 1] + 1], []
            ).append(substring(route, a, b))
    out = {}
    for i, (seg_id, _) in enumerate(segs):
        if i not in runs:
            continue
        piece = unary_union(runs[i])
        if not piece.is_empty:
            out[seg_id] = piece
    return out


def _voronoi(claims, envelope):
    """{segment_id: cell}, the ground nearest each segment's piece of route.

    Sampling every claim at a fixed spacing and unioning the cells that
    share a key, rather than generating from one point per claim. A point
    each makes the cell boundaries perpendicular bisectors between
    generators, so claims of unequal length -- which is what nearest-segment
    assignment produces -- splay them into visible fans. Sampling instead
    approximates the line-segment Voronoi, whose boundaries stand square to
    the channel however long the claim is.
    """
    pts, keys = [], []
    for seg_id, piece in claims.items():
        parts = (
            [piece]
            if piece.geom_type in ("Point", "LineString")
            else list(piece.geoms)
        )
        for part in parts:
            if part.geom_type == "Point":
                sampled = [part]
            else:
                m = max(2, int(part.length / CLAIM_SAMPLE_M) + 1)
                sampled = [
                    part.interpolate(t, normalized=True)
                    for t in np.linspace(0.0, 1.0, m)
                ]
            pts.extend(sampled)
            keys.extend([seg_id] * len(sampled))
    if len(pts) < 2:
        return {keys[0]: envelope} if keys else {}
    tree = shapely.STRtree(np.array(pts, dtype=object))
    out = {}
    for cell in _cells(pts, envelope):
        # By nearest generator, not by which one the cell contains: samples
        # 1 m apart down a straight reach are all but collinear, the diagram
        # comes back degenerate, and a strict contains test threw away the
        # cells whose point landed on their own boundary -- as white gaps
        # between finished zones. A cell's nearest generator is its own.
        out.setdefault(
            keys[int(tree.nearest(cell.representative_point()))], []
        ).append(cell)
    return {k: _dissolve_cells(v) for k, v in out.items()}


def _ray(end, toward, avoid, reach):
    """Straight extension from `end`, pointing away from `toward`, out to
    `reach`, on the first bearing that does not cross `avoid`.

    A route that arrives at its end from the wrong direction -- an oxbow,
    or a confluence it could only reach by looping through a side channel
    -- has an end bearing pointing back over its own length, and a ray
    along it slices the whole stretch in two. So the bearing is only the
    first candidate: RAY_STEPS bearings are tried in order of angular
    distance from it, and the first that exits without crossing the route
    (or the other ray) wins. Where every bearing crosses, the one that
    gets furthest before doing so is the least bad.
    """
    vx, vy = end[0] - toward[0], end[1] - toward[1]
    base = np.arctan2(vy, vx) if (vx or vy) else 0.0
    origin = Point(end)
    best = None
    for k in range(RAY_STEPS):
        step = (k + 1) // 2 * (1 if k % 2 else -1)
        ang = base + step * 2 * np.pi / RAY_STEPS
        far = (end[0] + np.cos(ang) * reach, end[1] + np.sin(ang) * reach)
        seg = LineString([end, far])
        crossings = [
            g
            for g in shapely.get_parts(seg.intersection(avoid))
            if g.distance(origin) > 1.0
        ]
        if not crossings:
            return far
        first = min(g.distance(origin) for g in crossings)
        if best is None or first > best[0]:
            best = (first, far)
    return best[1]


def _sides(pieces, route):
    """+1 (left) or -1 (right) for each piece, voted along the route.

    Every BANK_VOTE_M along the route, a point 1 m to the left and one 1 m
    to the right; whichever piece each lands in gets that vote. A piece's
    side is then the side of the route it actually lies on, however many
    pieces the split produced and whatever shape they have. A single
    representative point decided it before, and on a piece that wrapped
    round a bend the point could sit on the far side of the route from
    the ground the piece mostly covers. A piece that never borders the
    route -- the wedge between the two rays -- has no votes and falls back
    to that representative-point test.
    """
    n = max(2, int(route.length / BANK_VOTE_M) + 1)
    d = np.linspace(0.0, route.length, n)
    p = shapely.line_interpolate_point(route, d)
    a = shapely.line_interpolate_point(
        route, np.clip(d - 2.0, 0.0, route.length)
    )
    b = shapely.line_interpolate_point(
        route, np.clip(d + 2.0, 0.0, route.length)
    )
    tx = shapely.get_x(b) - shapely.get_x(a)
    ty = shapely.get_y(b) - shapely.get_y(a)
    norm = np.hypot(tx, ty)
    ok = norm > 0
    nx, ny = -ty[ok] / norm[ok], tx[ok] / norm[ok]
    px, py = shapely.get_x(p)[ok], shapely.get_y(p)[ok]
    probes = shapely.points(np.r_[px + nx, px - nx], np.r_[py + ny, py - ny])
    sign = np.r_[np.ones(int(ok.sum())), -np.ones(int(ok.sum()))]
    tree = shapely.STRtree(np.array(pieces, dtype=object))
    qi, ci = tree.query(probes, predicate="within")
    votes = np.zeros(len(pieces))
    np.add.at(votes, ci, sign[qi])
    out = np.sign(votes)
    for i in np.flatnonzero(out == 0):
        rp = pieces[i].representative_point()
        at = route.project(rp)
        near = route.interpolate(at)
        e0 = route.interpolate(max(at - TANGENT_M, 0.0))
        e1 = route.interpolate(min(at + TANGENT_M, route.length))
        cross = (e1.x - e0.x) * (rp.y - near.y) - (e1.y - e0.y) * (
            rp.x - near.x
        )
        out[i] = 1 if cross > 0 else -1
    return out


def _banks(route, extent):
    """(left, right) halves of `extent`, split by the extended route.

    The route is the channel, so extending it out of the extent divides
    the extent into one bank and the other, and a point's side is a
    property of which region it falls in -- decided once per stretch, not
    per cell. Cutting each cell against the route's local tangent instead
    meant every cell decided its own side, and 3% of them disagreed with
    their neighbour.

    The extension rays are chosen so they cannot cross the route (_ray),
    and the pieces are sided by where they lie along it (_sides); with a
    ray cutting back over the stretch, the split came out in five pieces
    and both banks of a 1.5 km stretch were written as "right".
    """
    coords = list(route.coords)
    if len(coords) < 2:
        return None, None
    # The box must contain the route, not just the band. A wide lake can
    # leave a stretch's only band ground on the far shore, and a box drawn
    # around that band alone sits entirely beside the route -- the extended
    # route then never crosses it, split returns one piece, and the whole
    # stretch was silently discarded.
    rb = route.bounds
    extent = (
        min(extent[0], rb[0]),
        min(extent[1], rb[1]),
        max(extent[2], rb[2]),
        max(extent[3], rb[3]),
    )
    box = shapely.box(*extent)
    reach = 2 * np.hypot(extent[2] - extent[0], extent[3] - extent[1])

    # Take the bearing over TANGENT_M of line, not one pixel step, so the
    # extension follows the channel rather than lattice jitter.
    head = route.interpolate(min(TANGENT_M, route.length))
    tail = route.interpolate(max(route.length - TANGENT_M, 0.0))
    p0 = _ray(coords[0], (head.x, head.y), route, reach)
    p1 = _ray(
        coords[-1],
        (tail.x, tail.y),
        unary_union([route, LineString([coords[0], p0])]),
        reach,
    )
    extended = LineString([p0] + coords + [p1])
    try:
        pieces = list(shapely.ops.split(box, extended).geoms)
    except Exception:
        return None, None
    if len(pieces) < 2:
        return None, None

    left, right = [], []
    for piece, side in zip(pieces, _sides(pieces, route), strict=True):
        (left if side > 0 else right).append(piece)
    return (
        unary_union(left) if left else None,
        unary_union(right) if right else None,
    )


# ── pass-2 parallelism ────────────────────────────────────────────────────
# The per-stretch zone construction is independent between stretches and
# entirely GEOS-bound, so it forks. The heavy read-only state — the water
# polygons and their STRtree — is placed in module globals
# before the fork and inherited copy-on-write: nothing large is ever
# pickled. The per-task payload is one stretch's own geometry; the result
# is its rows. Worker budget: the machine's cores divided by the number of
# catchments running beside us (KANTSONER_OUTER_WORKERS, exported by the
# CLI), so eight regions with four inner workers each fill 32 cores once.
_SHARED: dict = {}


# Extra memory while the zoning pool runs, GB per 1,000 stretches: forked
# workers copy the parent's pages as they touch them and the parent gathers
# every result. Glomma on 8 workers went from 43 to 83.6 GB; a parallel
# run's regions zone on 32 // workers, far fewer.
ZONING_GB_PER_1000 = 0.3


def _stretch_workers() -> int:
    outer = int(os.environ.get("KANTSONER_OUTER_WORKERS", "1") or "1")
    return max(1, (os.cpu_count() or 8) // max(outer, 1))


def _captured(sid, route, claims) -> set:
    """Segments whose claim sits on a watercourse that is not theirs.

    A route that finds nothing of its own stretch's channel inside
    CORRIDOR_M will take whatever water is in reach, and the neighbouring
    river usually is. Nothing downstream notices: the zones come out well
    formed, on the wrong watercourse, and are delivered with Beregnet
    "yes".

    Three things have to hold. The route runs on another route's water for
    more than a confluence's worth; this segment's claim is on that shared
    part; and some other stretch's elvenett line is nearer to the claim
    than its own. The last is what decides which of the two is the
    interloper, without having to arbitrate on size -- and it is why a
    genuinely wide river, whose line may be drawn well off its own
    channel, is never caught: no competing line is nearer.

    Captured segments are left unmeasured here, which hands them to the
    elvenett-line fallback in pass 3 (Usikkerhet "Ja").
    """
    route_tree = _SHARED.get("route_tree")
    lines = _SHARED.get("lines_by_sid")
    if route_tree is None or not lines:
        return set()
    routes_arr = _SHARED["routes_arr"]
    route_sids = _SHARED["route_sids"]
    near = [
        j
        for j in route_tree.query(
            route.buffer(ROUTE_COINCIDENT_M), predicate="intersects"
        )
        if route_sids[j] != sid and route_sids[j] in lines
    ]
    if not near:
        return set()
    shared = _valid_retry(
        unary_union, [routes_arr[j].buffer(ROUTE_COINCIDENT_M) for j in near]
    )
    overlap = _valid_retry(shapely.intersection, route, shared)
    if overlap.is_empty or overlap.length <= CAPTURE_MIN_M:
        return set()

    own = lines.get(sid)
    if own is None:
        return set()
    rival = [lines[route_sids[j]] for j in near]
    out = set()
    for seg_id, claim in claims.items():
        mid = claim.interpolate(0.5, normalized=True)
        if mid.distance(shared) > 0:
            continue
        d_own = mid.distance(own)
        if any(mid.distance(line) < d_own for line in rival):
            out.add(seg_id)
    return out


def _cut_by_unrouted_water(near, far, water, route_tree) -> bool:
    """Is what lies between these two parts water that nothing routes down?

    The gap has to be water over most of the connector -- otherwise the two
    parts are places the cell reached separately, not one bank cut in half
    -- and no route may run down that water. A watercourse with a route has
    its own cell and its own zones, so annexing its bank would hand the
    same ground to two segments; one with no route has no zones to take.
    """
    if water is None:
        return False
    a, b = nearest_points(near, far)
    conn = shapely.LineString([a, b])
    if conn.length <= 0:
        return True  # parts touch at a point; nothing separates them
    probe = conn.buffer(BRIDGE_PROBE_M)
    gap = _valid_retry(shapely.intersection, water, probe)
    if gap.is_empty:
        return False
    wet = _valid_retry(shapely.intersection, water, conn)
    if wet.is_empty or wet.length < conn.length * WATER_BRIDGE_MIN_FRAC:
        return False
    if route_tree is None:
        return True
    return len(route_tree.query(gap, predicate="intersects")) == 0


def _merge_parts(piece, claim, band_m, water, route_tree):
    """One bank slice out of a multi-part intersection.

    The nearest substantial part is the bank; any other substantial part
    joins it when only unrouted water separates the two (see
    _cut_by_unrouted_water). Everything else is dropped, which is what
    keeps a cell that reached across a confluence from delivering the
    neighbouring watercourse's bank as this segment's.
    """
    if piece.geom_type != "MultiPolygon":
        return piece
    parts = list(piece.geoms)
    floor = claim.length * band_m * ZONE_PART_MIN_FACTOR
    big = [g for g in parts if g.area > floor]
    if not big:
        return min(parts, key=claim.distance)
    near = min(big, key=claim.distance)
    keep = [near] + [
        g
        for g in big
        if g is not near and _cut_by_unrouted_water(near, g, water, route_tree)
    ]
    if len(keep) == 1:
        return near
    return _areal(_valid_retry(unary_union, keep))


def _zone_stretch(task):
    """Everything pass 2 does for one stretch. Runs in a forked worker."""
    sid, route, cell, line = task
    parts = _SHARED["parts"]
    tree = _SHARED["tree"]
    bands = _SHARED["bands"]
    band_m = float(max(bands))

    out = {
        "sid": sid,
        "skipped": False,
        "zones": [],
        "claims": [],
        "rings": {},
        "route_row": None,
        "captured": set(),
    }
    channel = _channel(route, cell, parts, tree) if cell is not None else None
    if channel is None:
        out["skipped"] = True
        return out
    # Subtract *all* water, not just this stretch's, or the band runs into
    # the neighbour's channel. Clip to the (CELL_OVERLAP_M-grown) cell, so
    # a stretch claims at most that far past its Voronoi boundary. "All"
    # water is bounded
    # to the polygons touching the widest band: anything further away cannot
    # intersect a ring, and differencing against the catchment-wide union
    # made this scale with region size instead of stretch size.
    near = tree.query(channel.buffer(band_m), predicate="intersects")
    water = _valid_retry(unary_union, parts[near]) if len(near) else None
    rings = {}
    for b in bands:
        ring = channel.buffer(b)
        if water is not None:
            ring = _valid_retry(shapely.difference, ring, water)
        ring = _areal(_valid_retry(shapely.intersection, ring, cell))
        if not ring.is_empty:
            rings[b] = ring
    out["rings"] = rings
    band = rings.get(int(band_m))
    if band is None:
        out["skipped"] = True
        return out

    segs = _segments(line, sid)
    claims = _claim(route, line, segs)
    # A captured claim is on the neighbouring river, so it is not this
    # segment's channel and nothing may be measured from it -- not the
    # zones below, and not bekk_bredde or avstand_beregning either, which
    # would otherwise report the other river's width and a matching
    # distance. Dropping the claim leaves both NaN, which is honest, and
    # pass 3 rebuilds the segment on its own elvenett line.
    captured = _captured(sid, route, claims)
    # Pass 3 needs to know: a captured stretch's route is on the wrong
    # watercourse, so it is the one case where the ELVIS line is the
    # better guess and the route must not be used.
    out["captured"] = set(captured)
    out["claims"] = [
        {"segment_id": seg_id, "stretch_id": sid, "geometry": piece}
        for seg_id, piece in claims.items()
        if seg_id not in captured
    ]
    seg_cells = _voronoi(
        claims, shapely.box(*channel.buffer(band_m * 2).bounds)
    )
    left, right = _banks(route, band.bounds)
    if left is None or right is None:
        out["skipped"] = True
        return out
    for seg_id, _seg in segs:
        cell_i = seg_cells.get(seg_id)
        if cell_i is None or seg_id in captured:
            continue
        cap = claims[seg_id].length * band_m * MAX_ZONE_FACTOR
        # A degenerate Voronoi cell can be invalid; the retry repairs it.
        slice_ = _valid_retry(shapely.intersection, cell_i, band)
        if slice_.is_empty or slice_.area < 1:
            continue
        for side, bank in (("L", left), ("R", right)):
            piece = _areal(_valid_retry(shapely.intersection, slice_, bank))
            if piece.is_empty:
                continue
            # Down to the ground that belongs to this bank (see
            # _merge_parts). Before the size tests, not after -- the cap
            # asks whether *this bank slice* is implausibly large, and
            # summing it with debris across the river threw away zones
            # that are fine once the debris is gone.
            piece = _merge_parts(
                piece, claims[seg_id], band_m, water, _SHARED.get("route_tree")
            )
            if piece.area < 1 or piece.area > cap:
                continue
            out["zones"].append(
                {
                    "segment_id": seg_id,
                    "stretch_id": sid,
                    "side": side,
                    "area_m2": round(piece.area, 2),
                    "geometry": piece,
                }
            )
    out["route_row"] = {
        "stretch_id": sid,
        "length_m": round(route.length, 1),
        "elvis_length_m": round(line.length, 1),
        "geometry": route,
    }
    return out


def _zone_fallback(
    sid, line, pieces, parts, tree, bands, real_rings, route_tree=None
):
    """Zones for segments the normal pass could not measure, built on the
    elvenett line itself.

    Pass 2 leaves a segment empty when there is no mapped water to take a
    bank from -- either grunnkart never drew this headwater stream, or the
    router, finding nothing of this stream's own inside CORRIDOR_M, took the
    larger watercourse running beside it and produced a band that the size
    checks then threw out. Either way the segment is a real watercourse that
    the delivery has to carry, and the ground beside it is measurable even
    though the water's edge is not.

    So the water is assumed: FALLBACK_HALF_WIDTH_M either side of the
    elvenett line stands in for the channel, and the bands are built off
    that exactly as they are off a real one. Real water is still subtracted,
    so where the line does turn out to run through mapped water the rings
    collapse and the segment stays unmeasured -- the fallback cannot invent
    a kantsone inside a lake.

    `pieces` is [(segment_id, geometry)] for the segments to rescue, and
    `real_rings` the rings pass 2 already produced for this stretch. The
    returned rings have those subtracted, so a segment pass 2 did measure
    keeps exactly the area it had: its zone lies inside the real ring, and
    what this adds is disjoint from it.
    """
    out = {"sid": sid, "zones": [], "rings": {}}
    if not pieces:
        return out
    band_m = float(max(bands))
    run = unary_union([p for _, p in pieces])
    channel = _areal(run.buffer(FALLBACK_HALF_WIDTH_M, cap_style="flat"))
    if channel.is_empty:
        return out

    near = tree.query(run.buffer(band_m), predicate="intersects")
    water = _valid_retry(unary_union, parts[near]) if len(near) else None
    rings = {}
    for b in bands:
        ring = _valid_retry(shapely.difference, run.buffer(b), channel)
        if water is not None:
            ring = _valid_retry(shapely.difference, ring, water)
        prior = real_rings.get(b)
        if prior is not None and not prior.is_empty:
            ring = _valid_retry(shapely.difference, ring, prior)
        ring = _areal(ring)
        if not ring.is_empty:
            rings[b] = ring
    out["rings"] = rings
    band = rings.get(int(band_m))
    if band is None:
        return out

    # Generated from EVERY segment of the stretch, not just the rescued
    # ones. A rescued segment sitting alone between two that pass 2 did
    # measure has no rescued neighbour to bound its cell, so generating
    # from the rescued subset gave it everything up to the next rescued
    # segment -- 53 segments away in one observed case -- and its zone
    # sprawled the length of the band. Zones are still only emitted for
    # the rescued ones; the rest are here to hold the boundaries.
    claims = dict(_segments(line, sid))
    seg_cells = _voronoi(claims, shapely.box(*line.buffer(band_m * 2).bounds))
    # Sides come off the whole elvenett line, not the rescued run, so a
    # segment's bank is the same bank its neighbours in the stretch got.
    left, right = _banks(line, band.bounds)
    if left is None or right is None:
        return out
    for seg_id, piece_line in pieces:
        cell_i = seg_cells.get(seg_id)
        if cell_i is None:
            continue
        cap = piece_line.length * band_m * MAX_ZONE_FACTOR
        slice_ = _valid_retry(shapely.intersection, cell_i, band)
        if slice_.is_empty or slice_.area < 1:
            continue
        for side, bank in (("L", left), ("R", right)):
            part = _areal(_valid_retry(shapely.intersection, slice_, bank))
            if part.is_empty:
                continue
            part = _merge_parts(part, piece_line, band_m, water, route_tree)
            if part.area < 1 or part.area > cap:
                continue
            out["zones"].append(
                {
                    "segment_id": seg_id,
                    "stretch_id": sid,
                    "side": side,
                    "area_m2": round(part.area, 2),
                    "geometry": part,
                }
            )
    return out


# ── orchestration ─────────────────────────────────────────────────────────
@dataclass
class Result:
    """Everything the geometry stage persists, one GeoDataFrame each."""

    segments: gpd.GeoDataFrame  # 10 m ELVIS pieces -- the delivery keys
    claims: gpd.GeoDataFrame  # each segment's piece of route
    zones: gpd.GeoDataFrame  # the kantsoner, one row per (segment, side)
    routes: gpd.GeoDataFrame  # reconstructed channel line per stretch
    buffers: dict[int, gpd.GeoDataFrame] = field(default_factory=dict)


def _stretches(elvis):
    """[(stretch_id, line, feature)], one per distinct piece of line.

    Keyed on elvID + strekningLnr: NVE's own identifiers, stable across
    downloads, where OBJECTID is only the row number of one export (it
    equals the GeoPackage fid) and the dataframe index the order the driver
    returned rows in.

    ELVIS repeats that key, though: nationally 69,783 ids sit on more than
    one feature. All but 260 are exact copies -- up to 11 of one line, never
    differing in any attribute -- and collapse here to one. The rest are
    different pieces under one id, usually a stretch stored as consecutive
    pieces, a few with gaps between them. Keyed on the id alone, each piece
    overwrote the last in build() and only the one read last was routed and
    delivered: 48 km in an earlier run, 6.3 km of Nidelva's main stem
    through Trondheim among it.

    So every piece of a repeated key (the parts of a multipart feature
    included), and every piece with no strekningLnr, gets a suffix _p0,
    _p1, ... numbered by position, not OBJECTID, so the ids survive a
    re-download. Pieces are not joined: the graph links segments by shared
    end points, not ids, and each piece keeps its own objektType and
    drawing direction -- two of the 260 ids hold a piece drawn the other
    way.
    """
    pieces: dict[str, list] = {}
    numbered = set()  # keys suffixed even when they hold a single piece
    seen = set()
    copies = 0
    for _, feat in elvis.iterrows():
        lnr = feat["strekningLnr"]
        if pd.isna(lnr) or not str(lnr).strip():
            # elvID alone would lump every such stretch of a river into one
            # key, so these are always numbered.
            key = str(feat["elvID"])
            numbered.add(key)
        else:
            key = f"{feat['elvID']}_{lnr}"
        for line in shapely.get_parts(feat.geometry):
            # normalize: a copy drawn the other way is still a copy
            copy = (key, shapely.to_wkb(shapely.normalize(line)))
            if copy in seen:
                copies += 1
                continue
            seen.add(copy)
            pieces.setdefault(key, []).append((line, feat))

    out, numbered_pieces = [], 0
    for key, group in pieces.items():
        if len(group) == 1 and key not in numbered:
            out.append((key, *group[0]))
            continue
        group.sort(key=lambda piece: _position(piece[0]))
        out += [
            (f"{key}_p{n}", line, feat) for n, (line, feat) in enumerate(group)
        ]
        numbered_pieces += len(group)
    if copies or numbered_pieces:
        logger.info(
            f"ELVIS: {copies} exact copies dropped, {numbered_pieces} pieces "
            "under a repeated or missing strekningLnr numbered _pN"
        )
    return out


def _position(line):
    """Sort key placing a piece by its end points, to the millimetre."""
    c = shapely.get_coordinates(line)
    return (*np.round(c[0], 3), *np.round(c[-1], 3), round(line.length, 3))


def build(
    elvis: gpd.GeoDataFrame,
    zone_water: gpd.GeoDataFrame,
    route_water_extra,
    *,
    bands: tuple[int, ...] = (3, 10, 20),
    lakes=None,
    dtm=None,
) -> Result:
    """The whole derivation for one region.

    `zone_water` is what banks are measured from (river polygons + trimmed
    creek polygons). `route_water_extra` is geometry that exists only for
    routing connectivity -- lakes, canals, untrimmed creeks -- and never
    grows a band. `lakes` is the lake polygons on their own: pixels inside
    them earn no depth reward when anchoring route endpoints, so endpoints
    snap to the closest lake pixel instead of deeper open water (_anchor).
    `dtm` is an open elevation raster (or None): dry ground the DTM drains a
    channel through is then routed like water, see _flow_channels. The one
    source read inside build, because the windows are only known per stretch.
    """
    crs = zone_water.crs
    parts = zone_water.explode(index_parts=False).geometry.values
    tree = shapely.STRtree(parts)

    # The routing surface, exploded and indexed: intersecting one stretch's
    # corridor against a catchment-sized union is quadratic in region size.
    route_parts = np.array(
        list(parts)
        + [
            p
            for g in route_water_extra
            for p in (g.geoms if hasattr(g, "geoms") else [g])
        ],
        dtype=object,
    )
    route_tree = shapely.STRtree(route_parts)

    stretches = _stretches(elvis)
    logger.info(
        f"{len(stretches)} ELVIS stretches, {len(parts)} water polygons"
    )

    # Lake polygons on their own index: _route burns the ones near a
    # stretch into a mask so anchoring can withhold their depth reward.
    lake_parts = np.array(
        [g for g in lakes if g is not None and not g.is_empty]
        if lakes is not None
        else [],
        dtype=object,
    )
    lake_tree = shapely.STRtree(lake_parts) if len(lake_parts) else None

    # ── pass 1: a route per stretch ───────────────────────────────────────
    routes, lines_by_sid, feat_by_sid = {}, {}, {}
    seg_rows = []
    for n, (sid, line, feat) in enumerate(stretches, 1):
        lines_by_sid[sid] = line
        feat_by_sid[sid] = feat
        if line.length < SEGMENT_M:
            continue
        corridor = line.buffer(CORRIDOR_M)
        near_idx = route_tree.query(corridor, predicate="intersects")
        route = None
        if len(near_idx):
            # All water touching the window: the depth field must not depend
            # on where this stretch's corridor happens to end, or
            # neighbouring stretches disagree about the same channel. _route
            # clips to its raster window; the corridor still bounds the
            # *route* through the cost surface.
            lk = (
                lake_parts[lake_tree.query(corridor)]
                if lake_tree is not None
                else None
            )
            route = _route(line, route_parts[near_idx], lake_parts=lk, dtm=dtm)
        if route is not None and route.length > 0:
            routes[sid] = route
        if n % 25 == 0:
            logger.info(f"routed {n}/{len(stretches)} stretches")

    logger.info(f"{len(routes)} routes, partitioning the catchment")
    cells = _partition(routes, zone_water.total_bounds)

    # ── pass 2: water, band and zones inside each cell ────────────────────
    # One task per stretch, forked workers, heavy state shared copy-on-write
    # (see _zone_stretch). Results are assembled in sorted-stretch order so
    # the written files are identical run to run regardless of completion
    # order; the serial path exists for one-core budgets and non-fork
    # platforms and produces byte-identical output.
    tasks = [
        (sid, route, cells.get(sid), lines_by_sid[sid])
        for sid, route in routes.items()
    ]
    # Shared with the forked workers: _merge_parts asks whether the water
    # separating two parts of a bank slice carries a route of its own.
    route_sids = list(routes)
    routes_arr = np.array([routes[s] for s in route_sids], dtype=object)
    route_tree = shapely.STRtree(routes_arr) if len(routes_arr) else None
    _SHARED.update(
        parts=parts,
        tree=tree,
        bands=tuple(bands),
        route_tree=route_tree,
        routes_arr=routes_arr,
        route_sids=route_sids,
        lines_by_sid=lines_by_sid,
    )
    n_workers = _stretch_workers()
    results = []
    if n_workers > 1 and "fork" in mp.get_all_start_methods():
        logger.info(f"zoning {len(tasks)} stretches on {n_workers} workers")
        with (
            membudget.surcharge(ZONING_GB_PER_1000, "zoning"),
            mp.get_context("fork").Pool(processes=n_workers) as pool,
        ):
            for k, res in enumerate(
                pool.imap_unordered(_zone_stretch, tasks, chunksize=8), 1
            ):
                results.append(res)
                if k % 200 == 0:
                    logger.info(f"zoned {k}/{len(tasks)} stretches")
    else:
        for k, task in enumerate(tasks, 1):
            results.append(_zone_stretch(task))
            if k % 200 == 0:
                logger.info(f"zoned {k}/{len(tasks)} stretches")
    _SHARED.clear()

    zones, out_routes, claim_rows, skipped = [], [], [], 0
    captured_all: set = set()
    rings_by_sid: dict[str, dict] = {}
    for res in sorted(results, key=lambda r: r["sid"]):
        rings_by_sid[res["sid"]] = res["rings"]
        if res["skipped"]:
            skipped += 1
            continue
        zones += res["zones"]
        claim_rows += res["claims"]
        captured_all |= res.get("captured") or set()
        if res["route_row"] is not None:
            out_routes.append(res["route_row"])

    # ── pass 3: the no-mapped-water fallback ──────────────────────────────
    # Every real-watercourse segment pass 2 left without a zone on both
    # banks gets one built on the elvenett line instead. See _zone_fallback
    # for why, and FALLBACK_OBJEKTTYPER for what is deliberately left out.
    both_sides: dict[str, set] = {}
    for z in zones:
        both_sides.setdefault(z["segment_id"], set()).add(z["side"])
    rescued: set = set()
    fallback_zones = []
    n_stretches = 0
    for sid, line in lines_by_sid.items():
        feat = feat_by_sid[sid]
        if str(feat.get("objektType")) not in FALLBACK_OBJEKTTYPER:
            continue
        # The channel is assumed along the ELVIS line, not along the
        # route. Following the route was tried and is wrong here: the
        # route runs *on* the mapped water, so the assumed channel lands
        # on water, and the band -- which subtracts all real water --
        # comes out empty. 525 of Orreaana's 5,922 segments lost both
        # bands that way, having had 90 m2 and 190 m2 before.
        #
        # The case that motivated it is narrower: a route that left the
        # water and followed a DTM depression, where the line is beside
        # the real channel and the route is on it. Telling the two apart
        # needs the route's own relation to mapped water, not just its
        # existence, so it is left undone rather than done wrongly.
        missing = [
            (seg_id, piece)
            for seg_id, piece in _segments(line, sid)
            if both_sides.get(seg_id, set()) != {"L", "R"}
        ]
        if not missing:
            continue
        res = _zone_fallback(
            sid,
            line,
            missing,
            parts,
            tree,
            tuple(bands),
            rings_by_sid.get(sid, {}),
            route_tree,
        )
        if not res["zones"]:
            continue
        n_stretches += 1
        fallback_zones += res["zones"]
        rescued.update(z["segment_id"] for z in res["zones"])
        for b, ring in res["rings"].items():
            prior = rings_by_sid.setdefault(sid, {}).get(b)
            rings_by_sid[sid][b] = (
                _valid_retry(unary_union, [prior, ring])
                if prior is not None and not prior.is_empty
                else ring
            )
    # Only segments that came out with both banks are really rescued; a
    # half-measured one is no better than before and is dropped, so
    # Beregnet cannot go "yes" on the strength of one side.
    have = {}
    for z in fallback_zones:
        have.setdefault(z["segment_id"], set()).add(z["side"])
    rescued = {s for s, sides in have.items() if sides == {"L", "R"}}
    # A segment reaches the fallback when it is short of a bank, which
    # includes having one. That one came off the route, and the two are
    # not the same measurement -- keeping both would deliver a segment
    # twice over, once from each. The fallback pair replaces it.
    zones = [z for z in zones if z["segment_id"] not in rescued]
    zones += [z for z in fallback_zones if z["segment_id"] in rescued]
    logger.info(
        f"fallback: {len(rescued):,} segments rescued on the elvenett line "
        f"across {n_stretches:,} stretches"
    )

    buffer_rows: dict[int, list] = {b: [] for b in bands}
    for sid in sorted(rings_by_sid):
        for b, ring in rings_by_sid[sid].items():
            buffer_rows[b].append(
                {"id": sid, "stretch_id": sid, "geometry": ring}
            )

    # ── the delivery keys: every stretch, routed or not ───────────────────
    # A stretch without a route still exists officially and is delivered
    # with empty attributes rather than silently missing.
    carry = [
        c
        for c in (
            "elvID",
            "elvenavn",
            "objektType",
            "elveordenStrahler",
            "vassdragsnummer",
        )
        if c in elvis.columns
    ]
    for sid, line in lines_by_sid.items():
        feat = feat_by_sid[sid]
        for seg_id, piece in _segments(line, sid):
            row = {
                "segment_id": seg_id,
                "stretch_id": sid,
                "routed": sid in routes,
                # Measured off an assumed channel on the elvenett line
                # rather than mapped water -- see _zone_fallback. Carried
                # through to Beregnet and Usikkerhet in stages.py.
                "fallback": seg_id in rescued,
                "geometry": piece,
            }
            for c in carry:
                row[c] = feat[c]
            seg_rows.append(row)

    def _gdf(rows, columns):
        if rows:
            return gpd.GeoDataFrame(rows, geometry="geometry", crs=crs)
        return gpd.GeoDataFrame(
            {c: [] for c in columns}, geometry="geometry", crs=crs
        )

    result = Result(
        segments=_gdf(seg_rows, ["segment_id", "stretch_id", "geometry"]),
        claims=_gdf(claim_rows, ["segment_id", "stretch_id", "geometry"]),
        zones=_gdf(
            zones, ["segment_id", "stretch_id", "side", "area_m2", "geometry"]
        ),
        routes=_gdf(
            out_routes,
            ["stretch_id", "length_m", "elvis_length_m", "geometry"],
        ),
        buffers={
            b: _gdf(rows, ["id", "stretch_id", "geometry"])
            for b, rows in buffer_rows.items()
        },
    )
    logger.info(
        f"{len(result.zones):,} zones over "
        f"{result.zones['segment_id'].nunique():,} of "
        f"{len(result.segments):,} segments, {skipped} stretches skipped"
    )
    return result
