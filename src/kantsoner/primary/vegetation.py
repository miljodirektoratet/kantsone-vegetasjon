"""Coverage of the riparian zones by canopy, buildings and roads.

One of four extractors; each measures the same per-segment, per-side slices
from `geometry.segment_zones` against a different source.

Every percentage is `area of the feature inside the slice / Atot_m2`, where
Atot_m2 is the slice's own area — computed once per band and side, and reused
by every column. Deriving the denominator separately per column let two
`gpd.overlay` calls disagree by a hundredth of a square metre, which on a
2 m² slice reads as 111 %.
"""

import geopandas as gpd
import numpy as np
import pandas as pd

from ..geometry import segment_zones
from ..io import as_gdf
from ..spec.attributes import band_token
from ..spec.constants import (
    AREAL_BAND,
    HEIGHT_CLASSES,
    PRIMARY_BAND,
    VTOT_CLASSES,
)


def _area_by_segment(slices, features, by=None):
    """Area of `features` inside each slice, summed per segment.

    With `by`, returns one series per distinct value of that column.
    """
    cols = ["segment_id", "geometry"]
    right = features[[by, "geometry"]] if by else features[["geometry"]]
    ov = gpd.overlay(
        slices[cols], right, how="intersection", keep_geom_type=True
    )
    if ov.empty:
        return {} if by else pd.Series(dtype=float)
    ov["_a"] = ov.geometry.area
    if by:
        return {
            v: part.groupby("segment_id")["_a"].sum()
            for v, part in ov.groupby(by)
        }
    return ov.groupby("segment_id")["_a"].sum()


def extract(
    segments, zones, band, ring, trees, buildings=None, samferdsel=None
):
    """Coverage columns for one band, both sides, keyed on segment_id.

    The 0-3 band exists only to supply Vtot for the width ladder, so the
    height-class breakdown and the built-up columns are skipped there. 0-10
    and 0-20 both get the full set: system_a's Kantvegetasjon_A needs Bygg/
    Samferdsel at 0-20 too, for stretches wide enough to use that band as
    their official one (see primary/channel.py:beregningsbredde).
    """
    segments_gdf = as_gdf(segments)
    ids = segments_gdf["segment_id"]
    crs = zones[0].crs
    trees_gdf = as_gdf(trees).to_crs(crs)
    full = tuple(band) in (PRIMARY_BAND, AREAL_BAND)
    buildings_gdf = (
        as_gdf(buildings, crs=crs) if full and buildings is not None else None
    )
    samferdsel_gdf = (
        as_gdf(samferdsel, crs=crs) if full and samferdsel is not None else None
    )

    out = {}
    for side, cells in (("L", zones[0]), ("R", zones[1])):
        prefix = f"K_{side}_{band_token(band)}"
        sl = segment_zones.cut(cells, ring)

        # The one denominator every percentage on this band and side uses.
        atot = ids.map(
            sl.assign(_a=sl.geometry.area).groupby("segment_id")["_a"].sum()
        )
        atot = atot.fillna(0.0)
        out[f"{prefix}_Atot_m2"] = atot.round(2)

        def pct(area, atot=atot):
            share = ids.map(area).fillna(0.0) / atot.replace(0, np.nan) * 100
            return share.clip(upper=100).fillna(0.0).round(2)

        by_class = _area_by_segment(sl, trees_gdf, by="height_class")
        if full:
            for hc in HEIGHT_CLASSES:
                out[f"{prefix}_V{hc}"] = pct(
                    by_class.get(hc, pd.Series(dtype=float))
                )

        # Index-aligned addition would give NaN for any segment missing one
        # class; concat+groupby treats absent classes as zero.
        parts = [by_class[hc] for hc in VTOT_CLASSES if hc in by_class]
        vtot = (
            pd.concat(parts).groupby(level=0).sum()
            if parts
            else pd.Series(dtype=float)
        )
        out[f"{prefix}_Vtot"] = pct(vtot)
        if full:
            out[f"{prefix}_Vtot_m2"] = ids.map(vtot).fillna(0.0).round(2)

        if buildings_gdf is not None:
            out[f"{prefix}_Bygg"] = pct(_area_by_segment(sl, buildings_gdf))
        if samferdsel_gdf is not None:
            out[f"{prefix}_Samferdsel"] = pct(
                _area_by_segment(sl, samferdsel_gdf)
            )

    return out
