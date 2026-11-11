"""Rolling the shoreline up to the lake, and the two delivered files.

The shore is measured in 10 m segments like a riverbank, so the roll-up is
the same one the river levels use -- additive contributions per segment,
summed per group, divided in `_ratios` -- with the lake as the group
instead of a block of consecutive segments. Reusing those two functions is
deliberate: a lake's Kantvegetasjon share has to mean exactly what a
river's does, or the two cannot appear in the same report.
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd
from loguru import logger

from ..columns import is_lake
from . import aggregate

LEVEL = "innsjo"
LAKE_KEY = "lake_id"


def aggregate_lakes(gdf: gpd.GeoDataFrame, polys: gpd.GeoDataFrame):
    """One row per lake, carrying the shoreline's aggregated values."""
    if polys is None or polys.empty or LAKE_KEY not in gdf.columns:
        return None
    shore = gdf[is_lake(gdf) & gdf[LAKE_KEY].notna()]
    if shore.empty:
        return None

    base = aggregate._base_counts(shore)
    base[LAKE_KEY] = shore[LAKE_KEY].to_numpy()
    block = base.groupby(LAKE_KEY, dropna=True).sum(numeric_only=True)
    ratios = aggregate._ratios(block, LEVEL)
    stats = pd.concat([block, ratios], axis=1).reset_index()

    out = polys.merge(stats, on=LAKE_KEY, how="left")
    out[f"{LEVEL}_strandlinje_m"] = (
        out[LAKE_KEY]
        .map(shore.groupby(LAKE_KEY)[aggregate.LENGTH].sum())
        .round(1)
    )
    out[f"{LEVEL}_segmenter"] = out[LAKE_KEY].map(
        shore.groupby(LAKE_KEY).size()
    )
    logger.info(
        f"Lakes: aggregated {len(shore):,} shore segments -> {len(out)} lakes"
    )
    return out


def shore_delivery(gdf: gpd.GeoDataFrame, zones: gpd.GeoDataFrame):
    """The per-segment file: every measured value on the zone it was
    measured over.

    The geometry is the 0-20 m zone rather than the shoreline, because
    that is the area the numbers describe -- a line would invite reading
    them as properties of the water's edge.
    """
    shore = gdf[is_lake(gdf)]
    if shore.empty or zones is None or zones.empty:
        return None
    attrs = pd.DataFrame(shore.drop(columns="geometry"))
    out = zones[["segment_id", "area_m2", "geometry"]].merge(
        attrs, on="segment_id", how="inner"
    )
    return gpd.GeoDataFrame(out, geometry="geometry", crs=zones.crs)
