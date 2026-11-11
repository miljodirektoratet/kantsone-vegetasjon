"""Keep the kantsoner inside Norwegian territory.

Where a river *is* the national border the border runs down the middle of
the channel, and grunnkart -- like every other Norwegian source -- stops
there. The far bank is then foreign ground that nothing in this pipeline
has data for: no water polygon (so the unmapped half of the channel is
counted as land), no land cover, no laser. Left alone it is measured as a
bank with no vegetation and no potential for any, which is not a finding
about the river but an artefact of where the data ends.

So every zone is clipped to the Norway polygon as it is built. Two things
follow from doing it here rather than by a rule further down:

* A bank that is *partly* foreign keeps its Norwegian part and its
  percentages are computed over that, because Atot is clipped with
  everything else. There is no threshold to choose and no arm of the logic
  where an 80 %-Norwegian bank gets thrown away or kept whole.
* A bank that is *entirely* foreign clips to nothing, which the rest of
  the pipeline already understands as "this side has no zone".

The second case still needs saying out loud, because "no zone" otherwise
means "we failed to build one" and the combination rules treat the two
differently -- a bank with no data must be skipped, not scored as bare.
:data:`OUTSIDE_COLS` carries that: a bank keeping less than
``C.GRENSE_KEEP_MIN_FRAC`` of its area is marked, and the both-bank
formulas fall back to the surviving side (see :mod:`kantsoner.sides`).

Coastal water is not affected. The administrative outline includes
Norway's maritime areas -- 470,107 km2 against ~324,000 km2 of land -- and
measured over Hurumlandet and Lærdal not one zone of the 393,340 built
came within 200 m of it. Only the land borders with Sweden, Finland and
Russia can clip anything.
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd
import shapely
from loguru import logger

from ..sides import BORDER_COL, NO, OUTSIDE_COLS, YES
from ..spec import constants as C


def _foreign(norge: gpd.GeoDataFrame, bounds, crs):
    """Everything within ``bounds`` that is *not* Norway, or None.

    Built as the complement rather than used directly, because the clip is
    then a difference against a few slivers along the border instead of an
    intersection against a national multipolygon with ~10^6 vertices. Inside
    the bounding box of one catchment the two are the same operation.
    """
    land = shapely.union_all(norge.to_crs(crs).geometry.values)
    outside = shapely.box(*bounds).difference(land)
    return None if outside.is_empty else outside


_POLYGONAL = ("Polygon", "MultiPolygon")


def _polygonal(geom):
    """Only the 2-D part of a difference result.

    A bank lying along the border can come back from ``difference`` as a
    bare line, or as a collection with one in it, where the sliver that was
    left collapsed to its own edge. Those are not empty, so filtering on
    ``is_empty`` keeps them, and one non-polygonal row is enough for
    ``gpd.overlay`` to refuse the entire frame -- it raises "df1 contains
    mixed geometry types" rather than dropping the row. That is how the
    first Tana run died four hours into the attributes stage, on 4 such
    rows out of 1,160,991.
    """
    if geom.geom_type in _POLYGONAL:
        return geom
    parts = [g for g in shapely.get_parts(geom) if g.geom_type in _POLYGONAL]
    return shapely.union_all(parts) if parts else shapely.Polygon()


def clip(zones: gpd.GeoDataFrame, norge, buffers=None):
    """Clip zones (and bands) to Norway; report the banks that vanished.

    Returns ``(zones, buffers, flags)`` where ``flags`` is a DataFrame of
    ``segment_id, side`` for every bank left with less than
    ``GRENSE_KEEP_MIN_FRAC`` of the area it was built with. Zones that clip
    away entirely are dropped, so a wholly foreign bank simply has no row.
    """
    empty = pd.DataFrame({"segment_id": [], "side": []})
    if norge is None or not len(norge) or zones.empty:
        return zones, buffers, empty

    # Padded so a zone touching the edge is clipped against real border
    # geometry rather than against the box.
    foreign = _foreign(norge, zones.total_bounds, zones.crs)
    if foreign is None:
        logger.info("landegrense: the catchment is wholly inside Norway")
        return zones, buffers, empty

    parts = shapely.get_parts(shapely.normalize(foreign))
    tree = shapely.STRtree(parts)

    def _cut(gdf):
        """Difference only the rows that actually reach foreign ground."""
        geoms = gdf.geometry.values
        hit = tree.query(geoms, predicate="intersects")
        if not len(hit[0]):
            return gdf, {}
        out = geoms.copy()
        was = {}
        for i in sorted(set(hit[0].tolist())):
            near = shapely.union_all(parts[hit[1][hit[0] == i]])
            was[i] = shapely.area(geoms[i])
            out[i] = _polygonal(shapely.difference(geoms[i], near))
        return gdf.set_geometry(out, crs=gdf.crs), was

    zones, was = _cut(zones)
    if buffers is not None:
        buffers = {b: _cut(g)[0] for b, g in buffers.items()}

    area = zones.geometry.area
    gone = pd.Series(False, index=zones.index)
    for i, before in was.items():
        gone.iat[i] = (
            before > 0 and area.iat[i] < before * C.GRENSE_KEEP_MIN_FRAC
        )
    if "area_m2" in zones.columns:
        zones["area_m2"] = area.round(2)

    flags = zones.loc[gone, ["segment_id", "side"]].reset_index(drop=True)
    zones = zones[~(zones.geometry.isna() | zones.geometry.is_empty)]
    logger.info(
        f"landegrense: {len(was):,} zone(s) clipped at the national border, "
        f"{len(flags):,} bank(s) left with no Norwegian ground"
    )
    return zones.reset_index(drop=True), buffers, flags


def mark(segments: gpd.GeoDataFrame, flags: pd.DataFrame) -> gpd.GeoDataFrame:
    """Add the per-side and combined border columns to the segments layer."""
    out = segments.copy()
    border = pd.Series(False, index=out.index)
    for side, col in OUTSIDE_COLS.items():
        ids = (
            set(flags.loc[flags["side"] == side, "segment_id"])
            if len(flags)
            else set()
        )
        is_out = out["segment_id"].isin(ids)
        out[col] = is_out.map({True: YES, False: NO})
        border |= is_out
    out[BORDER_COL] = border.map({True: YES, False: NO})
    return out
