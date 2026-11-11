"""Kantsoner around lakes, measured from the shoreline.

Elvenett draws a line down the middle of every lake, and the river track
measures outward from it. On a narrow lake that happens to work; on a wide
one the band never reaches land and the segment comes out unclassified --
7,725 km of it in an earlier national run. The kantsone of a lake belongs to its
shore, not to a line in the middle of the water, so this measures it there.

Three things decide the geometry:

* **Which water is a lake.** A Vann-Nett *innsjøvannforekomst* says so, and
  it is also the unit the result aggregates to. But its outline is
  generalised -- perimeters differ from the mapped shoreline by up to 3.4x
  -- so it is used to *select and name* a water surface, never to buffer.
  The shoreline buffered is grunnkart/FKB's, matched by largest overlap.

* **Islands count.** An island in a lake has a shore and a riparian zone
  like any other, so interior rings are segmented alongside the exterior.

* **The band is land.** It is the ground within the band width of the
  shore and outside *all* water, not just this lake's -- otherwise a
  neighbouring pond or an inflowing river would be measured as kantsone.

The band is split between shore segments by nearest-segment allocation, the
same way the river track splits a bank between the segments claiming it.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import shapely
from loguru import logger

from ..spec import constants as C

LAKE_LAYER = "VannforekomstInnsjø"
# Vann-Nett id, and the two fields worth carrying into the delivery.
ID_FIELD = "identifikasjon_lokalId"
EU_FIELD = "eUVannforekomstId"
NAME_FIELD = "navn"
# A mapped water surface must share at least this much of its area with a
# Vann-Nett lake to be called that lake's shoreline. Generous, because the
# two outlines are drawn to different generalisations; the point is only to
# reject a river that happens to clip a lake's smoothed outline.
MIN_OVERLAP = 0.30
# A water surface smaller than this is not part of the lake even when it
# sits wholly inside the lake's smoothed outline -- it is a pond beside it,
# and the overlap test alone cannot tell the two apart, because a pond
# fully inside the outline scores 1.0. In Orreaana five such slivers,
# 0.0-0.4 ha, added 771 m of shoreline that belongs to no lake.
MIN_WATER_M2 = 5_000.0
# Shore pieces this short are a rounding artefact of segmentize on a closed
# ring, not a stretch of shore.
MIN_PIECE_M = 1.0
SEGMENT_M = C.SEGMENT_M if hasattr(C, "SEGMENT_M") else 10.0
# A ring narrower than this on average (2 x area / perimeter) is a sliver
# the union left, not shore. Where two of a lake's polygons overlap along
# almost the same edge, the gap between them survives as a hole a few
# nanometres wide, or a scrap of one survives as a part of its own. Its
# boundary runs out and back along one line, so it is cut into segments
# lying on top of each other: 87 km of shore that is not there on Glomma's
# lakes (174 holes, 11 parts), and two such segments sharing a midpoint
# stopped the Voronoi split at Skasen. Real islands are far wider; between
# 5 and 50 cm Glomma has 3 rings in all.
SLIVER_M = 0.1

# Where Grunnkart cuts a lake at a kommune border, the two polygons stop
# centimetres short of each other (36 cm on Njakkajavri), and their union
# keeps a hairline slit along the border. Its two edges came out as shore:
# open water drawn twice, 97 km in 290 lakes nationally (Storevatnet
# 064-26215-L: 1.3 km of 7.6 km). Closing the union by this much (grow, then
# shrink) joins polygons up to twice as far apart. It also fills an island or
# inlet narrower than twice this, which no lake has but a needle of land
# 15 cm wide on 002-182-L.
CLOSE_GAP_M = 0.5


def _sliver(ring) -> bool:
    return 2.0 * shapely.Polygon(ring).area <= SLIVER_M * ring.length


def _surface(geoms):
    """A lake's mapped polygons merged into one surface, slivers dropped."""
    merged = shapely.make_valid(shapely.union_all(shapely.make_valid(geoms)))
    if len(geoms) > 1:
        # Only a lake mapped as several polygons can have a seam to close.
        merged = shapely.make_valid(
            merged.buffer(CLOSE_GAP_M).buffer(-CLOSE_GAP_M)
        )
    parts = []
    for part in shapely.get_parts(merged):
        if part.geom_type != "Polygon" or _sliver(part.exterior):
            continue
        holes = [r for r in part.interiors if not _sliver(r)]
        parts.append(shapely.Polygon(part.exterior, holes))
    return parts[0] if len(parts) == 1 else shapely.MultiPolygon(parts)


def _boundaries(poly) -> list:
    """Exterior and interior rings of a polygon or multipolygon."""
    out = []
    for part in shapely.get_parts(poly):
        if part.geom_type != "Polygon":
            continue
        out.append(shapely.LineString(shapely.get_exterior_ring(part).coords))
        for i in range(shapely.get_num_interior_rings(part)):
            ring = shapely.get_interior_ring(part, i)
            out.append(shapely.LineString(ring.coords))
    return out


def _cut(line, step: float) -> list:
    """A line chopped into pieces of `step`, the last one taking the rest.

    The pieces shapely.ops.substring cuts, in one pass over the line.
    substring copies all of the line's coordinates for every piece, so a
    shore took time in the square of its length: three quarters of an hour
    for Glomma's lakes, a third of it Mjøsa's 31,000 pieces.
    """
    n = max(1, int(round(line.length / step)))
    if n == 1:
        return [line]
    edges = np.linspace(0.0, line.length, n + 1)
    xy = shapely.get_coordinates(line)
    # Distance along the line at each vertex, summed the way substring sums
    # it, so a vertex beside a cut falls on the same side of it. The last
    # vertex is never inside a piece; it ends the last one.
    pts = xy.tolist()
    run, dist = [0.0], 0.0
    for (x1, y1), (x2, y2) in zip(pts[:-2], pts[1:-1], strict=True):
        dist += ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
        run.append(dist)
    first = np.searchsorted(run, edges[:-1], side="right")
    last = np.searchsorted(run, edges[1:], side="left")
    ends = shapely.get_coordinates(shapely.line_interpolate_point(line, edges))
    out = []
    for k in range(n):
        piece = shapely.LineString(
            np.vstack([ends[k], xy[first[k] : last[k]], ends[k + 1]])
        )
        if piece.length >= MIN_PIECE_M:
            out.append(piece)
    return out


def select(
    water: gpd.GeoDataFrame, lakes: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Mapped water surfaces that are a Vann-Nett lake, carrying its id.

    Matched by largest overlap, not by containment: neither outline
    contains the other reliably, because they are generalised differently.
    """
    if water.empty or lakes.empty:
        return water.iloc[0:0].assign(lake_id=[], lake_navn=[], lake_eu=[])
    wg = shapely.make_valid(water.geometry.values)
    lg = shapely.make_valid(lakes.geometry.values)
    tree = shapely.STRtree(lg)
    hit_w, hit_l = tree.query(wg, predicate="intersects")
    if not len(hit_w):
        return water.iloc[0:0].assign(lake_id=[], lake_navn=[], lake_eu=[])
    share = shapely.area(
        shapely.intersection(wg[hit_w], lg[hit_l])
    ) / np.maximum(shapely.area(wg[hit_w]), 1.0)
    keep = (share >= MIN_OVERLAP) & (shapely.area(wg[hit_w]) >= MIN_WATER_M2)
    hit_w, hit_l, share = hit_w[keep], hit_l[keep], share[keep]
    # One water surface can clip several lakes; the biggest overlap wins.
    order = np.lexsort((-share, hit_w))
    hit_w, hit_l = hit_w[order], hit_l[order]
    first = np.ones(len(hit_w), bool)
    first[1:] = hit_w[1:] != hit_w[:-1]
    hit_w, hit_l = hit_w[first], hit_l[first]

    out = water.iloc[hit_w].copy()
    for col, field in (
        ("lake_id", ID_FIELD),
        ("lake_navn", NAME_FIELD),
        ("lake_eu", EU_FIELD),
    ):
        out[col] = (
            lakes[field].to_numpy()[hit_l] if field in lakes.columns else None
        )
    return out.reset_index(drop=True)


def polygons(selected: gpd.GeoDataFrame, lakes_src: gpd.GeoDataFrame):
    """The reporting unit: one row per matched Vann-Nett lake.

    Vann-Nett geometry, not the mapped surfaces the shoreline was cut
    from. The two differ -- the outline is generalised, perimeters by up
    to 3.4x -- and this is the layer the result is aggregated to and
    reported on, so it has to be the water body as the forvaltning draws
    it. The mapped area is carried alongside so the difference is visible
    rather than hidden.
    """
    if selected.empty:
        return lakes_src.iloc[0:0]
    area = (
        selected.assign(_a=selected.geometry.area)
        .groupby("lake_id")["_a"]
        .sum()
    )
    out = lakes_src[lakes_src[ID_FIELD].isin(area.index)].copy()
    out = out.rename(
        columns={
            ID_FIELD: "lake_id",
            NAME_FIELD: "lake_navn",
            EU_FIELD: "lake_eu",
        }
    )
    keep = [c for c in ("lake_id", "lake_navn", "lake_eu") if c in out.columns]
    out = out[keep + ["geometry"]]
    out["areal_vannflate_m2"] = out["lake_id"].map(area).round(1)
    out["areal_vannforekomst_m2"] = out.geometry.area.round(1)
    return out.reset_index(drop=True)


def shoreline(selected: gpd.GeoDataFrame, step: float = SEGMENT_M):
    """10 m shore segments for every selected lake, exterior and islands."""
    rows = []
    # Grouped by lake, not by water surface: one Vann-Nett lake is often
    # mapped as several polygons -- Froeylandsvatnet is seven -- and a
    # counter restarting per polygon gave every one of them a segment
    # called _000000. 3,599 of Orreaana's 9,567 ids collided that way.
    for lake_id, grp in selected.groupby("lake_id", sort=False):
        # The union first, then its boundary. Taking each polygon's
        # boundary separately puts a "shoreline" wherever two of a lake's
        # polygons meet -- which is open water, not shore. On
        # Froeylandsvatnet that was 702 of 4,113 segments, every one of
        # them inside the lake with no land within 20 m.
        merged = _surface(grp.geometry.values)
        pieces = []
        for ring in _boundaries(merged):
            pieces += _cut(ring, step)
        head = grp.iloc[0]
        for i, piece in enumerate(pieces):
            rows.append(
                {
                    "segment_id": f"{lake_id}_{i:06d}",
                    "stretch_id": str(lake_id),
                    "lake_id": lake_id,
                    "lake_navn": head.lake_navn,
                    "lake_eu": head.lake_eu,
                    "objektType": "Innsjøstrand",
                    "geometry": piece,
                }
            )
    gdf = gpd.GeoDataFrame(rows, crs=selected.crs)
    # One bank only: the other side of the shoreline is the lake. Marked so
    # that every both-bank formula reduces to the one that exists, via the
    # same path a border river takes (see kantsoner.sides).
    gdf["K_L_Ingen_bredd"] = "Nei"
    gdf["K_R_Ingen_bredd"] = "Ja"
    return gdf


def bands(
    shore: gpd.GeoDataFrame,
    surfaces: dict,
    wet,
    cells: dict,
    band: tuple[int, int],
):
    """The land ring of one band, split between the shore segments.

    Outside *all* mapped water (`wet`, the parts of its union), not just
    this lake's: an inflowing river or a neighbouring pond is water, and
    measuring it as kantsone would read as bare ground that can never be
    vegetated. `surfaces` holds each lake's surface and `cells` each
    segment's Voronoi cell; all three are the same for every band, so
    build() makes them once.
    """
    # Only the water near a lake is taken off its ring. Off the whole union
    # it cost a second a lake on Tana, and the better part of an hour a band
    # on Glomma; the parts within the ring's box give the same ring to the
    # bit, in milliseconds.
    tree = shapely.STRtree(wet)
    rows = []
    for lake_id, grp in shore.groupby("stretch_id", sort=False):
        poly = surfaces[lake_id]
        ring = shapely.difference(poly.buffer(float(band[1])), poly)
        if ring.is_empty:
            continue
        near = wet[tree.query(shapely.box(*ring.bounds))]
        ring = shapely.difference(ring, shapely.multipolygons(near))
        if ring.is_empty:
            continue
        for idx in grp.index:
            part = shapely.intersection(ring, cells[idx])
            if part.is_empty or shapely.area(part) <= 0:
                continue
            rows.append(
                {
                    "segment_id": shore.at[idx, "segment_id"],
                    "stretch_id": lake_id,
                    "area_m2": float(shapely.area(part)),
                    "geometry": part,
                }
            )
    return gpd.GeoDataFrame(rows, crs=shore.crs)


def _nearest_cells(points, selected) -> dict:
    """Voronoi cells of the shore midpoints, so each piece of ring goes to
    the segment it is nearest -- the shoreline equivalent of the river
    track's claim cells.

    Built on the distinct midpoints: GEOS refuses a repeated one
    ("Multiple input coordinates in cell"), which a shoreline doubling back
    on itself gives (see SLIVER_M). Segments sharing a midpoint share its
    cell.
    """
    env = shapely.buffer(
        shapely.union_all(shapely.make_valid(selected.geometry.values)), 200.0
    )
    xy = shapely.get_coordinates(points)
    sites, which = np.unique(xy, axis=0, return_inverse=True)
    if len(sites) < len(xy):
        logger.warning(
            f"Lakes: {len(xy) - len(sites)} shore segments share a midpoint "
            "with another, and its zone"
        )
    cells = shapely.get_parts(
        shapely.voronoi_polygons(
            shapely.multipoints(sites), extend_to=env, ordered=True
        )
    )
    # GEOS can return a cell that crosses itself, and no overlay takes one
    # as input ("found non-noded intersection" -- a segment on Sandvatnet,
    # 002-32712-L). Rebuilt from its rings it is the same area, valid.
    bad = ~shapely.is_valid(cells)
    if bad.any():
        logger.debug(f"Lakes: {int(bad.sum())} invalid Voronoi cells repaired")
        cells[bad] = shapely.make_valid(
            cells[bad], method="structure", keep_collapsed=False
        )
    return dict(enumerate(cells[which.reshape(-1)]))


def tag_inside(segments, polys):
    """Mark river segments that lie inside a lake with that lake's id.

    The link the inheritance needs. It cannot go by water-body id: 98.7%
    of lake midline carries no Vann-Nett id at all, and the 1.3% that does
    carries a *river* body's, leaked across the shoreline by the 1 m
    tagging tolerance. Position is the only honest join.
    """
    out = segments.copy()
    out["lake_id"] = out.get("lake_id")
    if polys is None or polys.empty or segments.empty:
        return out
    mids = shapely.line_interpolate_point(
        segments.geometry.values, 0.5, normalized=True
    )
    tree = shapely.STRtree(shapely.make_valid(polys.geometry.values))
    hit_m, hit_p = tree.query(mids, predicate="within")
    if len(hit_m):
        ids = polys["lake_id"].to_numpy()
        col = out["lake_id"].to_numpy(dtype=object).copy()
        col[hit_m] = ids[hit_p]
        out["lake_id"] = col
    return out


def build(water, lakes_src, crs, bands_=None):
    """The whole lake track's geometry, in the river track's own shapes.

    Returns ``(segments, zones, rings, polys)`` where the first three are
    byte-compatible with what :func:`kantsoner.stages.geometry` writes for
    rivers -- same columns, same keys -- so the attributes, derived and
    aggregation stages read lake shore segments without knowing they are
    not riverbank. The fourth is the reporting unit.

    The one difference is declared in the data rather than in code: the
    right bank is marked absent, and every both-bank formula reduces to
    the left through :mod:`kantsoner.sides`, exactly as on a border river.
    """
    bands_ = tuple(bands_ or C.BANDS)
    empty = gpd.GeoDataFrame(geometry=[], crs=crs)
    sel = select(water, lakes_src)
    if sel.empty:
        return empty, empty, dict.fromkeys(bands_, empty), empty

    shore = shoreline(sel)
    segments = gpd.GeoDataFrame(
        {
            "segment_id": shore["segment_id"],
            "stretch_id": shore["stretch_id"],
            "routed": True,
            "fallback": False,
            "elvID": None,
            "elvenavn": shore["lake_navn"],
            "objektType": "Innsjøstrand",
            "elveordenStrahler": None,
            "vassdragsnummer": None,
            "K_L_Utenfor_Norge": "Nei",
            "K_R_Utenfor_Norge": "Nei",
            "K_S_Grenseelv": "Nei",
            "K_L_Ingen_bredd": "Nei",
            "K_R_Ingen_bredd": "Ja",
            C.KILDE_COL: C.KILDE_LAKE,
            # The lake is the water body; its id is the segment's too, so
            # the shore aggregates to the same unit the forvaltning reports.
            # A "string" column, as the river segments' is: concatenated
            # with an object column, their missing ids became pd.NA objects
            # that GeoPackage stores as the text "<NA>" -- 5,294 Alna
            # segments delivered under a water body called "<NA>".
            "vannforekomstid": shore["lake_eu"]
            .fillna(shore["lake_id"])
            .astype("string"),
            "lake_id": shore["lake_id"],
        },
        geometry=shore.geometry.values,
        crs=crs,
    )

    # Per-segment zone at the widest band, plus one ring per band keyed by
    # lake. That is the shape the attributes stage expects: it cuts each
    # segment's zone against the ring of the band it is measuring, so the
    # zone has to reach as far as the widest ring does.
    ids = sorted(shore["stretch_id"].unique())
    # The same for every band, so made once: the union of all the water
    # alone takes ~20 minutes on Glomma.
    surfaces = {
        str(lake_id): _surface(grp.geometry.values)
        for lake_id, grp in sel.groupby("lake_id", sort=False)
    }
    wet = shapely.get_parts(
        shapely.union_all(shapely.make_valid(water.geometry.values))
    )
    mids = shapely.line_interpolate_point(
        shore.geometry.values, 0.5, normalized=True
    )
    cells = _nearest_cells(mids, sel)
    zones, rings = None, {}
    for band in bands_:
        z = bands(shore, surfaces, wet, cells, band)
        by_lake = dict(tuple(z.groupby("stretch_id")))
        rings[band] = gpd.GeoDataFrame(
            {"id": ids},
            geometry=[
                shapely.union_all(by_lake[s].geometry.values)
                if s in by_lake
                else shapely.Polygon()
                for s in ids
            ],
            crs=crs,
        )
        if band == C.AREAL_BAND:
            zones = z
    logger.info(
        f"Lakes: {sel['lake_id'].nunique()} water bodies, "
        f"{len(shore):,} shore segments, {len(zones):,} zones"
    )
    return segments, zones, rings, polygons(sel, lakes_src)
