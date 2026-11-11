"""Stage orchestration — what each CLI command actually runs.

Four stages, boundaries on disk so a multi-hour run survives a crash and every
intermediate opens in QGIS:

    geometry     ELVIS-anchored kantsoner: routes, the catchment partition,
                 the 10 m ELVIS segments and their zones, and the 3/10/20 m
                 buffers the attribute stage intersects against
    attributes   every measured attribute, computed per batch of stretches,
                 keyed directly on the 10 m ELVIS segments
    derived      system A + B formulas, plus the topological attributes --
                 the graph is built straight on the ELVIS segments
    aggregation  roll-ups along Vann-Nett, unchanged

Topology is settled in `geometry`: each kantsone knows its 10 m ELVIS segment
and each segment its stretch, so nothing downstream ever re-derives ownership
— the transfer step that used to move attributes from a synthetic centerline
onto Vann-Nett no longer exists.
"""

from __future__ import annotations

import os

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from loguru import logger
from shapely.geometry import box
from shapely.ops import unary_union

from . import derived as derived_layer
from . import sides as sides_mod
from .config import Config, RegionPaths
from .geometry import elvis, landegrense
from .geometry import graph as graph_mod
from .geometry import lakes as lakes_mod
from .geometry.creek_polygons import CREEK_TRIM_M, _dissolve, _trim
from .io import (
    open_raster,
    read_all_layers,
    read_gpkg,
    read_vector,
    save_gpkg,
    save_gpkg_incremental,
)
from .network import barriers, connectivity, priority
from .primary import (
    canopy,
    channel,
    grofter,
    infrastructure,
    landcover,
    laser_metadata,
    rorlagt,
    singletree,
    terrain,
    vegetation,
)
from .spec import constants as C
from .spec.attributes import band_token

# Grunnkart (Arealregnskap) filters, pushed into the driver by pyogrio.
# The årsversjon2025 national FGDB carries the full arealdekkeNiva1/Niva2
# hierarchy the specification workbook describes -- niva1 for the classes
# that split further at niva2 (Skog into 5 species, Jordbruk into 5 land
# uses) where the old flat schema had one value, niva2 everywhere the old
# schema's granularity already matched it one for one.
RIVER_WHERE = "arealdekkeniva2 IN ('ferskvannElvBekk', 'ferskvannKanal')"
SAMFERDSEL_WHERE = "arealdekkeniva2 = 'samferdsel'"
# Lakes yield kantsoner like any other mapped water. (They were routing-only
# for a while; that meant a stretch running along a lakeshore either invented
# a channel or delivered nothing, and the invention was worse.)
LAKE_WHERE = (
    "arealdekkeniva2 IN ('ferskvannInnsjoTjern', "
    "'ferskvannInnsjoTjernRegulert')"
)
# The ditch prediction layer carries other line types in the same table.
GROFT_WHERE = "type = 'ditch'"
# "Where vegetation could grow at all" — the reference state for system A.
POTVEG_WHERE = (
    "arealdekkeniva1 = 'skog' OR arealdekkeniva2 LIKE 'jordbruk%' "
    "OR (arealdekkeniva2 = 'snaumarkUspesifisert' "
    "AND grunnforhold = 'jorddekt')"
)
# The same filter as a pandas query, applied to the batch's grunnkart frame
# instead of a second driver read.
POTVEG_QUERY = (
    "arealdekkeniva1 == 'skog' or arealdekkeniva2.str.startswith('jordbruk') "
    "or (arealdekkeniva2 == 'snaumarkUspesifisert' "
    "and grunnforhold == 'jorddekt')"
)
CANAL_TYPE = "KanalGrøft"
CREEK_TYPE = "ElvBekk"
CREEK_HALF_WIDTH_M = 0.5  # FKB carries vannbredde=1.0 for every creek here
# Islands are filled out of a water polygon depending on the Strahler order
# of the stretch running past them -- see _fill_holes. elvenett's own
# stream-order field.
STRAHLER_FIELD = "elveordenStrahler"
STRAHLER_FILL_MIN = 5  # at/above this order, every island is filled
# Below STRAHLER_FILL_MIN, an island is filled only up to a ceiling that
# grows with the local Strahler order: ISLAND_FILL_BASE_M2 * order ** 2.
# order=1 -> 150 m^2, order=2 -> 600 m^2, order=3 -> 1350 m^2, order=4 ->
# 2400 m^2. An island with no nearby Strahler order is treated as order 1,
# the most conservative (smallest-ceiling) assumption.
ISLAND_FILL_BASE_M2 = 150.0
ISLAND_FILL_EXPONENT = 2
# Each 10 m segment is tagged with the Vann-Nett water body whose line its
# buffer touches -- see _with_vannforekomst.
VANNFOREKOMST_ID = "vannforekomstid"
VANNFOREKOMST_MATCH_M = 1.0


def _aoi(paths: RegionPaths) -> gpd.GeoDataFrame:
    return read_gpkg(paths.aoi)


def _source(root, aoi, *, crs, **kw) -> gpd.GeoDataFrame:
    """AOI subset of a national dataset, in memory.

    Nothing is written: the L1/L2 algorithms take their inputs through
    ``io.as_gdf``, which passes a GeoDataFrame straight through. Only features
    intersecting the AOI are ever read.
    """
    return read_vector(root, aoi, crs=crs, **kw)


# ── Source subsets, read on demand ─────────────────────────────────────────
# These are filtered reads of national data, not computed products. Caching
# them under geometry/ made them look like results and cost ~35 MB per region
# to save ~60 s of driver-side filtering. Each stage that needs one reads it.


def _fkb_aoi(cfg: Config, aoi):
    """The AOI grown by the tree-detection margin.

    FKB only matters within reach of the riparian zones, and a river whose
    zones cross the catchment boundary still needs the buildings on the far
    side.
    """
    out = aoi.copy()
    out["geometry"] = aoi.buffer(C.DETECTION_MARGIN_M)
    return out


def _buildings(
    cfg: Config, aoi, *, exact_mask: bool = False
) -> gpd.GeoDataFrame:
    """FKB Bygning + BygnAnlegg footprints. Polygon (FLATE) layer only."""
    fkb_aoi = _fkb_aoi(cfg, aoi)
    opts = {
        "crs": cfg.crs,
        "layer": "polygons",
        "clip": True,
        "exact_mask": exact_mask,
    }
    return gpd.GeoDataFrame(
        pd.concat(
            [
                _source(cfg.paths.fkb.bygning_dir, fkb_aoi, **opts),
                _source(cfg.paths.fkb.bygnanlegg_dir, fkb_aoi, **opts),
            ],
            ignore_index=True,
        ),
        geometry="geometry",
        crs=cfg.crs,
    )


def _grunnkart(cfg: Config, aoi, where: str) -> gpd.GeoDataFrame:
    return _source(
        cfg.paths.grunnkart_dir, aoi, crs=cfg.crs, where=where, clip=True
    )


def _fkb_vann_lines(cfg: Config, aoi) -> gpd.GeoDataFrame:
    lines = _source(
        cfg.paths.fkb.vann_dir,
        _fkb_aoi(cfg, aoi),
        crs=cfg.crs,
        layer="lines",
        clip=True,
    )
    # No FKB-Vann lines in the AOI (or empty tiles, as in data/testomraade):
    # the read has no columns, and the type filters below need this one.
    if "objekttypenavn" not in lines.columns:
        lines["objekttypenavn"] = pd.Series(dtype="object", index=lines.index)
    return lines


def _elvis(cfg: Config, aoi) -> gpd.GeoDataFrame:
    """NVE ELVIS river-network stretches intersecting the AOI, whole.

    elvenett rather than elvevannforekomster: the water-body layer stops at
    every lakeshore (the lake is its own innsjøvannforekomst), so a
    watercourse through a lake had no owner and its shoreline no kantsone.
    The network layer runs InnsjøMidtlinje through every lake and carries
    real topology (strekningLnr, Strahler order). (elvID, strekningLnr) is
    the stable identity -- though ELVIS repeats it on some features; see
    elvis._stretches.

    Not clipped: a stretch crossing the catchment boundary keeps its full
    geometry, because the 10 m segments cut from it are the delivery keys
    and must be the same segments a neighbouring region would cut. The
    intersects test is against the AOI polygon itself — the driver-side
    filter is only a bbox, which admits stretches from the bbox corners
    that never touch the region.
    """
    gdf = _source(cfg.paths.elvenett, aoi, crs=cfg.crs, clip=False)
    gdf = gdf[gdf.geometry.notna()]
    return gdf[gdf.intersects(aoi.union_all())]


# ── L1 geometry ────────────────────────────────────────────────────────────
def _island_ceiling_m2(order: np.ndarray) -> np.ndarray:
    """Max area filled at each (already order>=1, NaN-free) Strahler order."""
    return ISLAND_FILL_BASE_M2 * order**ISLAND_FILL_EXPONENT


def _fill_holes(
    gdf: gpd.GeoDataFrame, elvis_net: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Water polygons with some or all of their islands filled in.

    An interior ring is an island; treating it as land grew kantsoner on
    island shores and stopped bands mid-river at every islet. Filled, the
    pipeline sees the water body as one sheet there: no band is measured
    from an island shore, and the water subtraction in zoning removes the
    full sheet. A pond on a large island is its own grunnkart feature, so it
    survives the fill on its own row.

    Each island is judged on its own, by the Strahler order of the nearest
    elvenett stretch to *that island* (not the whole enclosing water
    polygon, which can be one connected sheet spanning several orders and
    would otherwise borrow a neighbouring reach's order for a small
    tributary's island). At/above STRAHLER_FILL_MIN every island is filled
    regardless of size -- on a big river an interior ring is a mid-channel
    bar, an artefact of where the two arms were drawn. Below that, only
    islands up to _island_ceiling_m2(order) are filled; a real landform
    beside a channel narrow enough that its shore is a genuine kantsone is
    left alone. A stretch with no Strahler order is treated as order 1.
    """
    out = gdf.explode(index_parts=False, ignore_index=True)
    geoms = out.geometry.values
    holed = shapely.get_num_interior_rings(geoms) > 0
    if not holed.any():
        return out

    if STRAHLER_FIELD in elvis_net.columns and len(elvis_net):
        order = pd.to_numeric(
            elvis_net[STRAHLER_FIELD], errors="coerce"
        ).to_numpy()
    else:
        order = np.empty(0)
    order = (
        np.where(np.isnan(order) | (order < 1), 1, order)
        if len(order)
        else order
    )
    tree = shapely.STRtree(elvis_net.geometry.values) if len(order) else None

    geoms = geoms.copy()
    n_islands = n_filled = n_broken = 0
    for i in np.flatnonzero(holed):
        poly = geoms[i]
        try:
            rings = list(poly.interiors)
            keep = []
            for ring in rings:
                n_islands += 1
                island = shapely.polygons(ring)
                local_order = 1.0
                if tree is not None:
                    hit = tree.query(island, predicate="intersects")
                    if len(hit):
                        local_order = order[hit].max()
                    else:
                        local_order = order[tree.nearest(island)]
                if local_order >= STRAHLER_FILL_MIN or (
                    island.area <= _island_ceiling_m2(local_order)
                ):
                    n_filled += 1
                    continue
                keep.append(ring)
            if len(keep) != len(rings):
                # shapely.polygons(shell, holes=[]) raises (a numpy shape
                # mismatch inside its linearrings() coercion, not a real "no
                # holes" case) where omitting holes entirely does not -- and
                # keep=[] is the single most common outcome here, every
                # island in the polygon qualifying to be filled.
                geoms[i] = (
                    shapely.polygons(shapely.get_exterior_ring(poly))
                    if not keep
                    else shapely.polygons(
                        shapely.get_exterior_ring(poly), holes=keep
                    )
                )
        except Exception:
            # Belt and braces: some other malformed-geometry case a national
            # grunnkart merge could still produce. Leave this one polygon's
            # islands as-is rather than fail the whole region's geometry
            # stage over it.
            n_broken += 1
    out["geometry"] = geoms
    logger.info(
        f"islands filled: {n_filled} of {n_islands}"
        + (
            f" ({n_broken} water polygons skipped, unfillable)"
            if n_broken
            else ""
        )
    )
    return out


def _creek_polygons(creek_lines: gpd.GeoDataFrame) -> list:
    """FKB creek lines as bank polygons: dissolve, trim, flat-cap buffer.

    Dissolving first means only the true ends of a creek are pulled back;
    the 5 m trim plus the flat cap keep invented water out of the band where
    a creek line runs to the centre of the water body it joins rather than
    to its bank.
    """
    runs = _dissolve(list(creek_lines.geometry))
    return [
        _trim(g, CREEK_TRIM_M).buffer(CREEK_HALF_WIDTH_M, cap_style="flat")
        for g in runs
        if g is not None and not g.is_empty
    ]


def _vannforekomst_names(cfg: Config | None) -> pd.Series:
    """Vann-Nett's own name per id, for its rivers and its lakes.

    Empty without a config: merge_delivery.py runs ``aggregation`` with
    ``cfg=None`` to get the lake-shore table, where names play no part.
    """
    paths = getattr(cfg, "paths", None)
    frames = []
    for path, key, name in (
        (
            getattr(paths, "elvevannforekomster", None),
            VANNFOREKOMST_ID,
            "vannforekomstnavn",
        ),
        (
            getattr(paths, "vannforekomster", None),
            "identifikasjon_lokalId",
            "navn",
        ),
    ):
        if path:
            table = read_gpkg(path, columns=[key, name], read_geometry=False)
            frames.append(table.rename(columns={key: "id", name: "name"}))
    if not frames:
        return pd.Series(dtype="string")
    names = pd.concat(frames).drop_duplicates("id").set_index("id")["name"]
    return names.astype("string")


def _with_vannforekomst(cfg: Config, segments: gpd.GeoDataFrame):
    """`segments` with a VANNFOREKOMST_ID column.

    The elvevannforekomst whose line the segment, grown by
    VANNFOREKOMST_MATCH_M, intersects. Where several do (a segment sitting
    on the node between two water bodies) the first in file order is
    taken; where none does the id is empty. The column is the link from
    the ELVIS delivery keys to Vann-Nett's reporting units.
    """
    if not cfg.paths.elvevannforekomster or segments.empty:
        return segments
    aoi = gpd.GeoDataFrame(
        geometry=[box(*segments.total_bounds).buffer(VANNFOREKOMST_MATCH_M)],
        crs=segments.crs,
    )
    vf = _source(cfg.paths.elvevannforekomster, aoi, crs=cfg.crs, clip=False)
    ids = pd.Series(pd.NA, index=segments.index, dtype="string")
    if len(vf):
        tree = shapely.STRtree(vf.geometry.values)
        si, vi = tree.query(
            segments.geometry.buffer(VANNFOREKOMST_MATCH_M).values,
            predicate="intersects",
        )
        first = pd.Series(vi, index=si).groupby(level=0).min()
        ids.iloc[first.index.to_numpy()] = (
            vf[VANNFOREKOMST_ID].astype("string").to_numpy()[first.to_numpy()]
        )
    out = segments.copy()
    out[VANNFOREKOMST_ID] = ids
    logger.info(
        f"{int(ids.notna().sum()):,} of {len(segments):,} segments matched "
        "an elvevannforekomst"
    )
    return out


def geometry(cfg: Config, paths: RegionPaths) -> None:
    aoi = _aoi(paths)

    logger.info("Reading grunnkart water and FKB-Vann for the AOI")
    # Read before the water: which islands get filled is decided per-island
    # by the Strahler order of the nearest elvenett stretch.
    elvis_net = _elvis(cfg, aoi)
    river_polys = _fill_holes(_grunnkart(cfg, aoi, RIVER_WHERE), elvis_net)
    fkb_vann = _fkb_vann_lines(cfg, aoi)
    creek_lines = fkb_vann[fkb_vann["objekttypenavn"] == CREEK_TYPE]

    # The water banks are measured from — and the only water kantsoner can
    # come from: river polygons, creek polygons, lakes. A route section over
    # anything else yields no kantsone.
    lakes = _fill_holes(_grunnkart(cfg, aoi, LAKE_WHERE), elvis_net)
    geoms = list(river_polys.geometry) + list(lakes.geometry)
    if cfg.run.include_creeks and len(creek_lines):
        geoms += _creek_polygons(creek_lines)
    zone_water = gpd.GeoDataFrame(
        {"id": range(len(geoms))}, geometry=geoms, crs=cfg.crs
    )
    save_gpkg(zone_water, paths.rivers)

    # Routing-only water: canals/ditches carry the network but never grow a
    # band, and the raw creek lines are here untrimmed and round-capped on
    # purpose — the 5 m trim exists to keep invented water out of the *band*,
    # and applying it here would open gaps in the very connectivity these are
    # added for: a 10 m gap at a junction is a rounding error in length and
    # an absolute barrier to a least-cost path.
    extras = []
    for lines in (
        fkb_vann[fkb_vann["objekttypenavn"] == CANAL_TYPE],
        creek_lines,
    ):
        if len(lines):
            extras.append(
                unary_union(lines.geometry.values).buffer(CREEK_HALF_WIDTH_M)
            )

    logger.info("Routing every ELVIS stretch and building the partition")
    result = elvis.build(
        elvis_net,
        zone_water,
        extras,
        bands=tuple(b[1] for b in C.BANDS),
        lakes=lakes.explode(index_parts=False).geometry.values,
        # Routing pays for climbing (elvis._elevation_walls); windows are
        # read through the DTM's overviews, so this is cheap even nationally.
        dtm=open_raster(str(cfg.paths.dtm)) if cfg.paths.dtm else None,
    )

    # A river that *is* the national border has grunnkart on one bank only;
    # the buffer on the other lands on ground no Norwegian source covers.
    # Clip both zones and bands to Norwegian territory so every area, and
    # therefore every percentage derived from one, is measured over data
    # that exists. See geometry/landegrense.py.
    zones, buffers, border = landegrense.clip(
        result.zones,
        _source(cfg.paths.landegrense, aoi, crs=cfg.crs)
        if cfg.paths.landegrense
        else None,
        result.buffers,
    )

    segments = _with_vannforekomst(
        cfg, landegrense.mark(result.segments, border)
    )
    segments[C.KILDE_COL] = C.KILDE_RIVER
    left = zones[zones["side"] == "L"].drop(columns="side")
    right = zones[zones["side"] == "R"].drop(columns="side")
    rings = {(0, outer): gdf for outer, gdf in buffers.items()}

    # ── the lake track ────────────────────────────────────────────────────
    # Elvenett's line down the middle of a lake has no kantsone within the
    # widest band, so 7,725 km of it across an earlier national run came out
    # unclassified. The kantsone of a lake is on its shore, and that is
    # what this measures -- as 10 m shore segments carrying the same
    # columns a riverbank segment does, so everything downstream reads
    # them without knowing the difference. See geometry/lakes.py.
    if cfg.paths.vannforekomster:
        vn = _source(
            cfg.paths.vannforekomster,
            aoi,
            crs=cfg.crs,
            layer=lakes_mod.LAKE_LAYER,
        )
        l_seg, l_zone, l_rings, l_poly = lakes_mod.build(
            zone_water, vn, cfg.crs, bands_=C.BANDS
        )
        if len(l_seg):
            # Which river segments sit inside which lake -- the link the
            # derived stage inherits Barrier/Oppstuving/Fragmentering over.
            segments = lakes_mod.tag_inside(segments, l_poly)
            # A region with lakes but no rivers (Røst) has an empty,
            # all-object river frame; concatenated, it turns the lakes'
            # bool flags into objects that GeoPackage stores as "True".
            segments = (
                pd.concat([segments, l_seg], ignore_index=True)
                if len(segments)
                else l_seg
            )
            left = pd.concat([left, l_zone], ignore_index=True)
            for band, ring in l_rings.items():
                rings[band] = pd.concat([rings[band], ring], ignore_index=True)
            save_gpkg(l_poly, paths.geometry_dir / "lake_polygons.gpkg")

    # A lake midline keeps its segment -- it carries the network topology
    # the graph and the lake's own Barrier/Oppstuving need. Inside a lake
    # that is a Vann-Nett innsjøvannforekomst it gets no zone: the lake
    # track measures that shore, and a band built outward from a line in
    # the middle of the water is the artefact it replaces (2,460 on
    # Orreaana). A lake that is NOT its own vannforekomst is part of the
    # river water body it lies in, so its midline keeps the river routine
    # -- zones from the mapped water edge -- and is reported with that
    # river. Dropping those too left 60,000 segments on Bjerkreimvassdraget
    # with no kantsone and no lake to carry it.
    midline = segments["objektType"].isin(C.LAKE_OBJEKTTYPER)
    if "lake_id" in segments:
        midline &= segments["lake_id"].notna()
    else:
        midline &= False
    if midline.any():
        drop_seg = set(segments.loc[midline, "segment_id"])
        # rings are per stretch: drop only stretches lying wholly in a lake
        kept_str = set(segments.loc[~midline, "stretch_id"].dropna())
        drop_str = set(segments.loc[midline, "stretch_id"]) - kept_str
        before = len(left) + len(right)
        left = left[~left["segment_id"].isin(drop_seg)]
        right = right[~right["segment_id"].isin(drop_seg)]
        rings = {b: g[~g["id"].isin(drop_str)] for b, g in rings.items()}
        logger.info(
            f"Lake midlines: dropped {before - len(left) - len(right):,} "
            f"zones measured from the middle of a Vann-Nett lake"
        )

    save_gpkg(segments, paths.segments)
    save_gpkg(result.claims, paths.claims)
    save_gpkg(result.routes, paths.routes)
    for band, gdf in rings.items():
        save_gpkg(gdf, paths.zones(band))
    save_gpkg(left, paths.zone_side("left"))
    save_gpkg(right, paths.zone_side("right"))


# ── L2 attributes ──────────────────────────────────────────────────────────
def _batch_name(ids) -> str:
    """A batch's file-safe name: its first stretch id."""
    return ids[0].replace("/", "_").replace(" ", "_")


def _batch_plan(paths, budget):
    """Groups of whole stretches that neighbour each other on the ground.

    Ordering stretches by id scattered each batch across the catchment, so
    its tree detection mask spanned the whole region and singletree walked
    every CHM tile again. A Hilbert ordering keeps a batch's mask compact,
    so each tile is touched about once.

    Whole stretches stay together because a segment's zone is derived from
    its stretch's channel, so splitting one would change nothing but cost
    the batch two reads of the same water.
    """
    counts = read_gpkg(
        paths.segments, columns=["stretch_id"], read_geometry=False
    )["stretch_id"].value_counts()
    ring = read_gpkg(paths.zones(C.AREAL_BAND), columns=["id"])
    mid = ring.geometry.centroid
    pts = gpd.GeoSeries(
        gpd.points_from_xy(mid.x, mid.y), index=ring["id"], crs=ring.crs
    )

    out, batch, n = [], [], 0
    for sid in pts.hilbert_distance().sort_values().index:
        c = int(counts.get(sid, 0))
        if not c:
            continue
        if batch and n + c > budget:
            out.append((batch, n))
            batch, n = [], 0
        batch.append(str(sid))
        n += c
    # Stretches with no band at all (never routed, or all water) still carry
    # delivery keys; they ride in the last batch so their rows exist with
    # empty attributes rather than going missing from the output.
    banded = set(pts.index)  # once: inside the loop it made this quadratic
    rest = [str(s) for s in counts.index if s not in banded]
    if rest:
        batch += rest
        n += int(counts[rest].sum())
    if batch:
        out.append((batch, n))
    return out


def _in(ids):
    joined = ",".join("'" + i.replace("'", "''") + "'" for i in ids)
    return f"IN ({joined})"


def _attributes_batch(cfg, paths, ids):
    """Every measured attribute for one batch of stretches."""
    where = f"stretch_id {_in(ids)}"
    seg_b = read_gpkg(paths.segments, where=where)
    left = read_gpkg(paths.zone_side("left"), where=where)
    right = read_gpkg(paths.zone_side("right"), where=where)
    claims_b = read_gpkg(paths.claims, where=where)
    rings = {
        band: read_gpkg(paths.zones(band), where=f"id {_in(ids)}")
        for band in C.BANDS
    }

    # Grunnkart, the ditch layer and the water are read for this batch's
    # own ground, not its bounding box: a Hilbert batch of sprawling
    # northern stretches spans a huge bbox, and a rectangle read pulled
    # millions of land-cover polygons (25 GB) only to clip them away. The
    # exact-mask read materialises only features touching the rings.
    bounds = box(*rings[C.AREAL_BAND].total_bounds).buffer(C.CROWN_MARGIN_M)
    batch_aoi = gpd.GeoDataFrame(geometry=[bounds], crs=seg_b.crs)
    footprint = gpd.GeoDataFrame(
        geometry=[rings[C.AREAL_BAND].union_all().buffer(C.CROWN_MARGIN_M)],
        crs=seg_b.crs,
    )
    grunnkart = _source(
        cfg.paths.grunnkart_dir,
        footprint,
        crs=cfg.crs,
        clip=True,
        exact_mask=True,
    )
    # Sources that were read region-wide before the batching: on a giant
    # catchment each was gigabytes held for the whole stage. Per batch and
    # exact-masked they are a few hundred rows each.
    buildings = _buildings(cfg, footprint, exact_mask=True)
    samferdsel = grunnkart[grunnkart["arealdekkeniva2"] == "samferdsel"]
    potveg = grunnkart.query(POTVEG_QUERY)
    water_b = read_gpkg(paths.rivers, bbox=bounds.bounds)
    grofter_b = (
        _source(
            cfg.paths.grofter,
            batch_aoi,
            crs=cfg.crs,
            where=GROFT_WHERE,
            clip=True,
        )
        if cfg.paths.grofter
        else None
    )
    prosjektdekning_b = (
        _source(
            cfg.paths.nhm_prosjektdekning, batch_aoi, crs=cfg.crs, clip=True
        )
        if cfg.paths.nhm_prosjektdekning
        else None
    )

    # Trees are detected per batch, over that batch's zones only. A crown
    # straddling two batches is found twice, which is harmless: each batch
    # only measures the part overlapping its own slices.
    try:
        mask = left.geometry.union_all().union(right.geometry.union_all())
    except shapely.errors.GEOSException:
        sides_geom = np.concatenate(
            [left.geometry.values, right.geometry.values]
        )
        # Near-coincident zone edges can defeat a float union (Reisa:
        # 'non-noded intersection'). Snapping to a 1 cm grid always nodes,
        # and the crown margin below dwarfs the shift.
        logger.warning(
            "GEOS topology error in the detection mask; snapping to 1 cm"
        )
        mask = shapely.union_all(sides_geom, grid_size=0.01)
    mask = mask.buffer(C.CROWN_MARGIN_M)
    if len(buildings):
        mask = mask.difference(buildings.union_all())
    safe = _batch_name(ids)
    mask_path = paths.attributes_dir / f"_detection_mask_{safe}.gpkg"
    save_gpkg(gpd.GeoDataFrame(geometry=[mask], crs=seg_b.crs), mask_path)
    try:
        _, trees = singletree.run_singletree(
            cfg.paths.chm,
            cfg.paths.dtm,
            None,
            None,
            None,
            clip_mask_path=mask_path,
            centerlines_path=paths.routes,
            rivers_path=paths.rivers,
            write_heights_raster=False,
            min_height=C.MIN_TREE_HEIGHT_M,
            gaussian_sigma=C.TREE_GAUSSIAN_SIGMA,
            median_size=C.TREE_MEDIAN_SIZE,
            compactness=C.TREE_COMPACTNESS,
        )
    finally:
        mask_path.unlink(missing_ok=True)
    # The per-batch tree-crown polygons (delivered as single_trees), otherwise
    # only held in memory for vegetation.extract below and then discarded. Opt-in --
    # a national run has too many batches, and crown count scales with
    # canopy cover, to pay this write by default.
    if os.environ.get("KANTSONER_EXPORT_SINGLETREE") and len(trees):
        save_gpkg(trees, paths.attributes_dir / f"_singletree_{safe}.gpkg")

    cols = {}
    for band in C.BANDS:
        cols.update(
            vegetation.extract(
                seg_b,
                (left, right),
                band,
                rings[band],
                trees,
                buildings,
                samferdsel,
            )
        )

    cols.update(
        canopy.extract(
            seg_b,
            (left, right),
            C.PRIMARY_BAND,
            rings[C.PRIMARY_BAND],
            cfg.paths.chm,
        )
    )

    ch = channel.extract(seg_b, claims_b, water_b)
    cols.update(ch)
    cols.update(laser_metadata.extract(seg_b, prosjektdekning_b))
    cols.update(
        grofter.extract(
            seg_b,
            (left, right),
            C.PRIMARY_BAND,
            rings[C.PRIMARY_BAND],
            grofter_b,
            ch["K_S_bekk_azimut"],
        )
    )

    near = (left, right, rings[C.PRIMARY_BAND])
    wide = (left, right, rings[C.AREAL_BAND])
    nve = cfg.paths.nve
    cols.update(terrain.extract(seg_b, near[:2], C.PRIMARY_BAND, cfg.paths.dtm))
    cols.update(landcover.extract(seg_b, near, C.PRIMARY_BAND, potveg))
    cols.update(
        landcover.classes(
            seg_b, near[:2], C.PRIMARY_BAND, rings[C.PRIMARY_BAND], grunnkart
        )
    )
    cols.update(
        landcover.causes(
            seg_b, near[:2], C.PRIMARY_BAND, rings[C.PRIMARY_BAND], grunnkart
        )
    )
    cols.update(rorlagt.extract(seg_b, water_b, grunnkart))
    cols.update(
        landcover.veg_threshold(
            seg_b, near[:2], C.PRIMARY_BAND, rings[C.PRIMARY_BAND], grunnkart
        )
    )
    # Expected cover for the graded Kantvegetasjon (derived/kantveg_grad.py).
    cols.update(
        landcover.forventet_dekning(
            seg_b, near[:2], C.PRIMARY_BAND, rings[C.PRIMARY_BAND], grunnkart
        )
    )
    cols.update(
        landcover.areas(
            seg_b, near[:2], C.AREAL_BAND, rings[C.AREAL_BAND], grunnkart
        )
    )
    # system_a also classifies at 0-20 m, for stretches whose channel is wide
    # enough that beregningsbredde picks that band as the official one --
    # only the inputs it actually reads need a 0-20 pass; classes/terrain/
    # canopy/grofter/infrastructure stay 0-10-only, unaffected by the choice.
    cols.update(landcover.extract(seg_b, wide, C.AREAL_BAND, potveg))
    cols.update(
        landcover.causes(
            seg_b, wide[:2], C.AREAL_BAND, rings[C.AREAL_BAND], grunnkart
        )
    )
    cols.update(
        landcover.veg_threshold(
            seg_b, wide[:2], C.AREAL_BAND, rings[C.AREAL_BAND], grunnkart
        )
    )
    # PlastringErosjonssikri splits its features across a polygon layer and
    # a line layer of the same GeoPackage; gpd.read_file on a bare path
    # would silently read only the default one (~47% of its features live
    # in the other).
    vegvesen = (
        read_all_layers(cfg.paths.vegvesen) if cfg.paths.vegvesen else None
    )
    cols.update(
        infrastructure.extract(
            seg_b,
            near,
            wide,
            flomsikring=(nve.flomsikring_line, nve.flomsikring_punkt),
            erosjonssikring=(
                nve.erosjonssikring_line,
                nve.erosjonssikring_punkt,
                vegvesen,
            ),
        )
    )
    cols.update(
        infrastructure.barrier(seg_b, claims_b, nve.dam_line, nve.dam_punkt)
    )
    return seg_b.assign(**cols)


def attributes(cfg: Config, paths: RegionPaths) -> None:
    batches = list(_batch_plan(paths, cfg.run.batch_segments))
    total = sum(n for _, n in batches)
    logger.info(
        f"{total:,} segments in {len(batches)} batch(es) "
        f"of <= {cfg.run.batch_segments:,}"
    )

    def frames():
        for i, (ids, n) in enumerate(batches, 1):
            logger.info(
                f"batch {i}/{len(batches)}: "
                f"{len(ids):,} stretches, {n:,} segments"
            )
            yield _attributes_batch(cfg, paths, ids)

    save_gpkg_incremental(frames(), paths.attributes_file)


# ── L3 derived ─────────────────────────────────────────────────────────────
# Treslag measured again from arealdekkeniva2 for a region whose attributes
# stage read it from a Grunnkart field årsversjon 2025 does not have. When
# this file sits next to attributes.gpkg, derived takes Treslag from it; a
# fresh attributes stage measures it right and writes no such file.
TRESLAG_PATCH = "treslag_patch.csv.gz"
TRESLAG_COLUMNS = ("K_L_B0-10_Treslag", "K_R_B0-10_Treslag")


def _apply_treslag_patch(df: pd.DataFrame, paths: RegionPaths) -> pd.DataFrame:
    """Replace the Treslag columns by the patch file's, when it is there."""
    file = paths.attributes_dir / TRESLAG_PATCH
    if not file.exists():
        return df
    patch = pd.read_csv(
        file, dtype=str, keep_default_na=False, na_values=[""]
    ).set_index("segment_id")
    if patch.index.has_duplicates:
        raise ValueError(f"{file} has the same segment_id twice")
    df = df.copy()
    for col in TRESLAG_COLUMNS:
        df[col] = df["segment_id"].map(patch[col])
    logger.info(
        f"Treslag from {file.name}: "
        + ", ".join(
            f"{c} {df[c].notna().mean():.0%} filled" for c in TRESLAG_COLUMNS
        )
    )
    return df


def _fremmedart_observations(cfg: Config, paths: RegionPaths):
    """Artsdatabanken occurrences matched to the segment whose kantsone they
    lie in: DataFrame(segment_id, art), or None without a source.

    Matched against the kantsoner, not the ELVIS line — the line can sit
    200 m from the channel, the zones are where the plant was recorded.
    Points are dropped when flagged erroneous or unsure, or when their stated
    coordinate uncertainty could put them on a different reach.
    """
    if not cfg.paths.fremmedart:
        return None
    zones = pd.concat(
        [
            read_gpkg(paths.zone_side("left"), columns=["segment_id"]),
            read_gpkg(paths.zone_side("right"), columns=["segment_id"]),
        ],
        ignore_index=True,
    )
    if zones.empty:
        return None
    aoi = gpd.GeoDataFrame(
        geometry=[box(*zones.total_bounds).buffer(C.FREMMEDART_MATCH_M)],
        crs=zones.crs,
    )
    pts = _source(cfg.paths.fremmedart, aoi, crs=cfg.crs, clip=False)
    if pts.empty:
        return None
    flags = pts.reindex(
        columns=["absent", "hasErrors", "unsureIdentification"]
    ).astype(str)
    bad = flags.eq("True").any(axis=1)
    unc = pd.to_numeric(
        pts.get("coordinateUncertaintyInMeters"), errors="coerce"
    )
    pts = pts[~bad & ~(unc > C.FREMMEDART_MAX_UNCERTAINTY_M)]
    if pts.empty:
        return None
    import shapely

    tree = shapely.STRtree(zones.geometry.values)
    hit = tree.query_nearest(
        pts.geometry.values,
        max_distance=C.FREMMEDART_MATCH_M,
        all_matches=False,
    )
    if hit.shape[1] == 0:
        return None
    name_col = (
        "preferredPopularName"
        if "preferredPopularName" in pts
        else "validScientificName"
    )
    out = pd.DataFrame(
        {
            "segment_id": zones["segment_id"].values[hit[1]],
            "art": pts[name_col].values[hit[0]],
        }
    )
    logger.info(
        f"fremmedart: {len(pts)} usable occurrences, {len(out)} within "
        f"{C.FREMMEDART_MATCH_M:.0f} m of a kantsone on "
        f"{out['segment_id'].nunique()} segments"
    )
    return out


def _lake_network_inheritance(out: pd.DataFrame) -> pd.DataFrame:
    """Give a lake's shore the network attributes of its own midline.

    Barrier, Oppstuving and Fragmentering are properties of the river
    graph -- where a segment sits relative to dams and to the reaches above
    and below. A shoreline sits nowhere in that graph: it is a closed ring
    touching the lake, not the network. But the lake's midline *is* in the
    graph, so the answer already exists one join away.

    Worst-of per lake, not the mode: a dam on the outlet is a dam for the
    whole lake, and a shore that is half unaffected is not half-passable.
    """
    from .columns import is_lake

    if "lake_id" not in out.columns:
        return out
    lake = is_lake(out)
    if not lake.any():
        return out
    inherit = ("E_S_Barrier", "E_S_Oppstuving", "E_S_Fragmentering")
    src = out[~lake & out["lake_id"].notna()]
    if src.empty:
        return out
    order = {
        lab: i
        for i, lab in enumerate(
            ("Svært dårlig", "Dårlig", "Moderat", "God", "Svært god")
        )
    }
    for col in inherit:
        if col not in out.columns:
            continue
        if col == "E_S_Fragmentering":
            rank = src[col].map(order)
            pick = (
                src.assign(_r=rank)
                .sort_values("_r")
                .groupby("lake_id")[col]
                .first()
            )
        else:
            # Ja beats Nei: any barrier on the lake is the lake's barrier.
            pick = src.groupby("lake_id")[col].apply(
                lambda v: (
                    "Ja"
                    if (v.astype("string") == "Ja").any()
                    else ("Nei" if v.notna().any() else None)
                )
            )
        out.loc[lake, col] = out.loc[lake, "lake_id"].map(pick)
    return out


def derived(
    cfg: Config, paths: RegionPaths, *, allow_missing: bool = False
) -> None:
    """System A + B formulas, then the topological attributes.

    The graph is built directly on the 10 m ELVIS segments — they carry the
    official topology, so there is nothing to transfer first.
    """
    gdf = read_gpkg(paths.attributes_file)
    geoms, crs = gdf.geometry.values, gdf.crs
    df = pd.DataFrame(gdf.drop(columns="geometry"))
    del gdf
    df = _apply_treslag_patch(df, paths)
    logger.info(f"Deriving system A + B attributes for {len(df)} segments")

    # Every formula is row-wise, so the frame is derived in chunks: the
    # producers copy their input, and on a 600k-segment giant those copies
    # stacked to ~25 GB — the one memory spike batching the attributes stage
    # did not remove. The graph walk below is the only cross-row computation
    # and runs on the reassembled frame.
    chunk = 150_000
    if len(df) > chunk:
        out = pd.concat(
            [
                derived_layer.compute_all(
                    df.iloc[i : i + chunk], allow_missing=allow_missing
                )
                for i in range(0, len(df), chunk)
            ]
        )
    else:
        out = derived_layer.compute_all(df, allow_missing=allow_missing)

    logger.info("Building the river network graph")
    # The graph and the network walks need the geometry and a handful of
    # columns, not the ~200-column frame; a slim frame keeps the giants flat.
    official = gpd.GeoDataFrame(
        out[["segment_id", "stretch_id", "elvID"]], geometry=geoms, crs=crs
    )
    # Chains group by Vann-Nett's own reporting unit (delivery/aggregate.py
    # aggregates on it, not elvID -- see stages.py:aggregation), falling
    # back to elvID per segment where elvevannforekomster is not configured
    # or nothing matched within VANNFOREKOMST_MATCH_M.
    vannforekomst = out.get(VANNFOREKOMST_ID)
    official["river_key"] = (
        vannforekomst.fillna(out["elvID"])
        if vannforekomst is not None
        else out["elvID"]
    )
    # Flow direction is the elvenett digitising order, not a DTM guess.
    g = graph_mod.build(official)
    # Along-flow order for the block aggregation: blocks may only ever
    # contain consecutive water, and the segment ids alone cannot say which
    # stretch follows which.
    out = pd.concat(
        [out, graph_mod.chains(official, g, river_col="river_key")], axis=1
    )
    out = barriers.compute(out, g)
    out = _lake_network_inheritance(out)
    out = connectivity.compute(out, g, official)
    out = priority.compute(out, g, _fremmedart_observations(cfg, paths))

    save_gpkg(
        gpd.GeoDataFrame(out, geometry=geoms, crs=crs),
        paths.derived,
    )
    logger.info(f"Added {len(out.columns) - len(df.columns)} derived columns")


# ── L4 aggregation ─────────────────────────────────────────────────────────
# What the centerline is: the rivers Vann-Nett reports on, and nothing
# else. Two kinds go:
#
# * A lake line. Elvenett draws one down the middle of every lake, but
#   Vann-Nett's river water bodies stop at the shore -- checked: their
#   own lines put 0 m inside a lake polygon. A lake is reported as a
#   lake, through lake_polygon.gpkg, and its shore through
#   lake_shore.gpkg. Matched on objektType rather than on position,
#   because the Vann-Nett outlines are generalised and would swallow
#   genuine river segments near the shore.
#
# * A segment with no vannforekomstid. It belongs to no water body, so
#   there is nothing to report it against.
#
# Both are excluded before aggregation, not just from the delivery:
# left in, they sit in the denominator of every share at every level.
# A midline counts as lake line only inside a Vann-Nett lake, whose shore
# the lake track reports. A lake that is not its own innsjøvannforekomst
# is part of the river water body it lies in, and its midline is
# delivered as river (see the geometry stage). No lake_id column means
# the lake track found no Vann-Nett lake -- none in the region, or no
# lake layer configured -- so every midline is river, as the geometry
# stage already treated it. Krossåa: 964 of its 1,067 midline segments
# carry their river's vannforekomstid.


def _delivered(gdf: pd.DataFrame) -> pd.Series:
    """Which centerline segments go to the delivery and the aggregation."""
    from .columns import is_lake

    in_vn_lake = (
        gdf["lake_id"].notna()
        if "lake_id" in gdf
        else pd.Series(False, index=gdf.index)
    )
    keep = ~(
        gdf["objektType"].isin(C.LAKE_OBJEKTTYPER) & in_vn_lake
    ) & ~is_lake(gdf)
    if VANNFOREKOMST_ID in gdf:
        keep &= gdf[VANNFOREKOMST_ID].notna()
    return keep


def aggregation(cfg: Config, paths: RegionPaths) -> None:
    from .delivery import aggregate

    gdf = read_gpkg(paths.derived)
    # The segments are natively Vann-Nett, so the delivery names are aliases
    # of what geometry already settled, not the product of a transfer step.
    gdf["segment_id_vann-nett"] = gdf["segment_id"]
    # Aggregation (100/500/1000 m + river blocks) groups by Vann-Nett's own
    # reporting unit, not ELVIS's elvID -- a vannforekomst is the actual
    # water body Vann-Nett reports against, and can span or split from an
    # elvID boundary. Falls back to elvID per segment where a segment has no
    # vannforekomstid (elvevannforekomster not configured, or nothing within
    # VANNFOREKOMST_MATCH_M matched) so a gap in that optional source does
    # not drop the segment out of every aggregated level.
    vannforekomst = gdf.get(VANNFOREKOMST_ID)
    gdf["river_id_vann-nett"] = (
        vannforekomst.fillna(gdf["elvID"])
        if vannforekomst is not None
        else gdf["elvID"]
    )
    # The water body's name in Vann-Nett. ELVIS's elvenavn names the river a
    # stretch belongs to, which is often not the water body: eight bodies in
    # Figgjo read "Figgjo" that way, and most have none at all.
    names = _vannforekomst_names(cfg)
    gdf["vannforekomstnavn"] = gdf.get("elvenavn")
    if vannforekomst is not None and len(names):
        named = vannforekomst.map(names)
        gdf["vannforekomstnavn"] = named.where(
            named.notna(), gdf["vannforekomstnavn"]
        )
    gdf["segment_length_vann-nett"] = gdf.geometry.length.round(2)
    if "routed" in gdf:
        # Matched is pass-1 only: a channel could be reconstructed for this
        # stretch at all. It says nothing about whether pass 2 (channel,
        # band, banks) succeeded for this particular segment and side.
        gdf["Matched"] = gdf["routed"].map({True: "yes", False: "no"})

        # Beregnet is pass 2: did this segment actually end up with a real
        # zone on both sides. Pass 2 can still fail per segment, and per
        # side within a segment -- a stretch's first/last segment often
        # loses one side to the extended-tangent bank split, and a whole
        # stretch is dropped if no channel was found near it at all.
        # Downstream ratios (aggregate.py) and the Kantvegetasjon_A
        # override below key on this, not on Matched, since a routed
        # stretch can still deliver nothing for a given segment.
        # A fallback segment has no route -- its bands were built on the
        # elvenett line itself (geometry/elvis.py:_zone_fallback) -- but it
        # does have real zones and real measurements, so it counts as
        # computed. Usikkerhet below is what says the channel was assumed.
        beregnet = gdf["routed"]
        if "fallback" in gdf:
            beregnet = beregnet | gdf["fallback"].fillna(False).astype(bool)
        measured = []
        for side in ("L", "R"):
            col = f"K_{side}_{band_token(C.AREAL_BAND)}_Atot_m2"
            if col in gdf:
                # A bank abroad has no zone by design, not by failure --
                # its half of the river was never Norway's to measure. A
                # border segment with a good Norwegian bank is computed;
                # Grenseelv is what says only one bank went into it.
                ok = gdf[col] > 0
                abroad = sides_mod.absent(gdf, side)
                # Evidence only from a bank that has area *and* is Norway's
                # to measure. Area alone is not enough: a bank is called
                # abroad when it keeps less than GRENSE_KEEP_MIN_FRAC of
                # itself, not when it keeps nothing, so it can be excused
                # from the formulas and still carry a sliver -- 0.12 m2 on
                # one of these. Counting that sliver as evidence let 13
                # segments whose banks are both abroad through as computed,
                # and every both-bank formula skips both banks, so they were
                # delivered with an empty Kantvegetasjon_A: inside
                # length_matched, in none of the class buckets.
                measured.append(ok & ~abroad)
                beregnet = beregnet & (ok | abroad)
        if measured:
            # ...because excusing both banks would call a segment lying
            # wholly outside Norway computed on the strength of no
            # measurement at all. ELVIS draws the network across the
            # border, so a Finnish tributary's last stretches sit in the
            # catchment with both banks abroad -- 344 of them in Tana.
            # Without a usable bank the segment is Uklassifisert, which the
            # override below applies.
            beregnet = beregnet & np.logical_or.reduce(measured)
        gdf["Beregnet"] = beregnet.map({True: "yes", False: "no"})

        # system_a.py classifies before Beregnet exists, so it never sees
        # the "we didn't compute this" case -- a segment with no zone reads
        # "Naturlig fravær av kantskog" by construction (NaN compares below
        # threshold). Overwrite those here to Uklassifisert instead of
        # leaving a data gap mislabelled as a natural absence.
        unberegnet = ~beregnet
        if unberegnet.any():
            # Band-free is what aggregate.py and delivery actually read; the
            # raw B0-10/B0-20 variants are overwritten too so a QA export
            # never shows a "we didn't compute this" segment as a natural
            # absence in either one.
            for b in ("", "B0-10", "B0-20"):
                tok = f"{b}_" if b else ""
                for side in ("L", "R", "S"):
                    status_col = f"E_{side}_{tok}Kantvegetasjon_A"
                    arsak_col = f"E_{side}_{tok}Kantvegetasjon_A_arsak"
                    flag_col = f"E_{side}_{tok}Menneskelig_inngrep"
                    if status_col in gdf:
                        gdf.loc[unberegnet, status_col] = (
                            C.UKLASSIFISERT_KANTVEG
                        )
                    if arsak_col in gdf:
                        gdf.loc[unberegnet, arsak_col] = ""
                    if flag_col in gdf:
                        gdf.loc[unberegnet, flag_col] = ""
            for col in (
                "E_L_Kantvegetasjon_fjernet",
                "E_R_Kantvegetasjon_fjernet",
                "E_S_Kantvegetasjon_fjernet",
            ):
                if col in gdf:
                    gdf.loc[unberegnet, col] = np.nan
            if "E_S_Kantvegetasjon_grad" in gdf:
                gdf.loc[unberegnet, "E_S_Kantvegetasjon_grad"] = (
                    C.KANTVEG_GRAD_IKKE_BEREGNET
                )
        # A bekkelukking is not assessed at all: no class, so it falls out
        # of every share (CLASS_COUNTS has no bucket for it) rather than
        # counting as removed vegetation along a stream that is underground.
        if "K_S_Rorlagt_mistanke" in gdf:
            lukket = (
                gdf["K_S_Rorlagt_mistanke"].astype(str).str.strip().eq("Ja")
            )
            if lukket.any():
                for b in ("", "B0-10", "B0-20"):
                    tok = f"{b}_" if b else ""
                    for side in ("L", "R", "S"):
                        for col, val in (
                            (
                                f"E_{side}_{tok}Kantvegetasjon_A",
                                C.IKKE_VURDERT_BEKKELUKKING,
                            ),
                            (f"E_{side}_{tok}Kantvegetasjon_A_arsak", ""),
                            (f"E_{side}_{tok}Menneskelig_inngrep", ""),
                        ):
                            if col in gdf:
                                gdf.loc[lukket, col] = val
                for col in (
                    "E_L_Kantvegetasjon_fjernet",
                    "E_R_Kantvegetasjon_fjernet",
                    "E_S_Kantvegetasjon_fjernet",
                ):
                    if col in gdf:
                        gdf.loc[lukket, col] = np.nan
                if "E_S_Kantvegetasjon_grad" in gdf:
                    gdf.loc[lukket, "E_S_Kantvegetasjon_grad"] = (
                        C.IKKE_VURDERT_BEKKELUKKING
                    )
                logger.info(
                    f"Bekkelukking: {int(lukket.sum()):,} segments left unassessed"
                )
    # ── the lake deliverables ─────────────────────────────────────────────
    # Before the centerline is cut down, because this is the only place the
    # shore segments still exist. Same roll-up, grouped on the lake instead
    # of a block of consecutive segments, so a lake's share means what a
    # river's does.
    from .delivery import lakes as lakes_del

    poly_path = paths.geometry_dir / "lake_polygons.gpkg"
    if poly_path.exists():
        polys = read_gpkg(poly_path)
        agg = lakes_del.aggregate_lakes(gdf, polys)
        if agg is not None:
            save_gpkg(agg, paths.aggregation_dir / "lake_polygon.gpkg")
        zones = read_gpkg(paths.zone_side("left"))
        shore = lakes_del.shore_delivery(gdf, zones)
        if shore is not None:
            save_gpkg(shore, paths.aggregation_dir / "lake_shore.gpkg")

    keep = _delivered(gdf)
    dropped = int((~keep).sum())
    if dropped:
        logger.info(
            f"Centerline: {dropped:,} of {len(gdf):,} segments are lake line "
            f"or carry no vannforekomstid -- left out of the delivery"
        )
    gdf = gdf[keep].reset_index(drop=True)

    save_gpkg(gdf, paths.official_10m)
    aggregate.aggregate_all(gdf, paths)

    # The trimmed client tables, written here rather than by a later pass:
    # the two corrections that pass also carried are fixed in the pipeline,
    # so nothing is left for it to do. See delivery/trim.py.
    from .delivery import trim as trim_mod

    trim_mod.run(paths)
    if cfg is not None and os.environ.get("KANTSONER_EXPORT_SINGLETREE"):
        _join_single_trees(cfg, paths)


def _apex_keys(trees) -> np.ndarray:
    """A crown on the seam of two batches is found by both, at the same apex
    pixel, so the apex identifies it across the region."""
    x = np.rint(trees["top_x"].to_numpy() * 100).astype(np.int64)
    y = np.rint(trees["top_y"].to_numpy() * 100).astype(np.int64)
    return x * (1 << 31) + y


def _join_single_trees(cfg: Config, paths: RegionPaths) -> None:
    """delivery/single_trees.gpkg from the attributes stage's per-batch crowns.

    Batches are joined in plan order; a crown already kept from an earlier
    batch is dropped, crowns under SINGLE_TREE_MIN_HEIGHT_M are left out, and
    tree_id numbers the rest from 1. scripts/merge_delivery.py reads this
    file for the national single_trees layer.
    """
    seen = np.array([], dtype=np.int64)
    parts, next_id = [], 1
    for ids, _ in _batch_plan(paths, cfg.run.batch_segments):
        f = paths.attributes_dir / f"_singletree_{_batch_name(ids)}.gpkg"
        if not f.exists():
            continue
        trees = read_gpkg(f)
        trees = trees[trees["tree_height"] >= C.SINGLE_TREE_MIN_HEIGHT_M]
        key = _apex_keys(trees)
        _, first = np.unique(key, return_index=True)
        fresh = np.zeros(len(key), dtype=bool)
        fresh[first] = True
        fresh &= ~np.isin(key, seen)
        trees = trees[fresh].copy()
        seen = np.concatenate([seen, key[fresh]])
        trees["tree_id"] = np.arange(next_id, next_id + len(trees))
        next_id += len(trees)
        parts.append(trees.reindex(columns=list(C.SINGLE_TREE_FIELDS)))
    if not parts or next_id == 1:
        logger.info("single_trees: no crowns to deliver")
        return
    out = gpd.GeoDataFrame(
        pd.concat(parts, ignore_index=True), geometry="geometry", crs=cfg.crs
    )
    save_gpkg(
        out, paths.root / "delivery" / "single_trees.gpkg", layer="single_trees"
    )
    logger.info(
        f"single_trees: {len(out)} crowns -> delivery/single_trees.gpkg"
    )
