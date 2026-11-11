"""Suspected piped/culverted stream ("rørlagt bekk") detection.

A stream can be piped underground through farmland while the official
network (ELVIS) still draws a continuous line across the gap -- it is a
schematic network, not a surface-mapped one, so it does not stop just
because grunnkart/FKB-Vann has no polygon there. `elvis.py`'s router
doesn't stop either: its cost surface charges LAND_COST for a non-water
pixel rather than refusing it, so a route silently crosses the gap on
terrain alone and nothing else records that it happened.

This flags candidates, it does not detect anything: a real gap can also be
rock outcrop, dense canopy occlusion, or just stale mapping. The signal
combines three things that are each weak alone -- no mapped water nearby,
farmland on both banks, and a run of several such segments in a row rather
than one -- because together they are the specific pattern an agricultural
culvert produces and the others mostly don't.

Land cover is read from grunnkart directly, against strips either side of
the segment's own elvenett line, and NOT from the K_*_Jordbruk columns the
other extractors produce. Those are measured inside the segment's buffer
zone, and a segment in the middle of the gap has no zone -- so its Jordbruk
came back 0% and the flag could never fire on the very case it is for.
Measured over Hurumlandet + Lærdal: of 31,656 segments with no zone, the
number with farmland on both banks by the zone-derived columns was zero.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from ..io import as_gdf
from ..spec import constants as C

# How far a segment's own ELVIS line may sit from the nearest real mapped
# water polygon before counting as "no visible water here". ELVIS is a
# schematic digitisation, not a survey of the bank, so it sits some way off
# the mapped water even where the mapping is perfect -- too tight a radius
# reads that routine offset as a gap.
#
# Measured over the 216,358 segments of the two test catchments: distance
# from the segment line to the nearest mapped water is 0 m at the median
# and 10.4 m at p90 over the whole network, and tighter still in farmland
# (mean 2.9 m, p90 7.1 m, p95 13.0 m) because flat, accessible ground is
# better surveyed. 20 m clears that p95 with margin.
GAP_SEARCH_M = 20.0
# Width of the strip taken either side of the line to ask what the ground
# is. The primary band's own width: wide enough to be about the field
# rather than the verge, narrow enough not to reach across it.
JORDBRUK_SIDE_M = 10.0
# Share of each side strip that must be grunnkart jordbruk. Full coverage,
# not a majority: a stream is piped to make a field workable, so the case
# this is for runs through the middle of cultivated ground with farmland
# hard against both banks. Anything less than complete -- a verge, a
# treeline, a yard -- is a stream running past a field, and those keep
# their open channel. Strict, deliberately: this flag is for review, and a
# short list that is right beats a long one that is not.
JORDBRUK_MIN_PCT = 100.0
# Floating-point slack only, not a coverage allowance: an exact 100.0 %
# never survives polygon intersection arithmetic.
_PCT_EPS = 1e-6
# A single flagged segment is often a digitising wobble at a bend, not a
# real buried reach. Require a run of consecutive flagged segments (same
# stretch) whose combined length reaches this before calling it suspected
# piping -- length rather than a segment count, since segments are not a
# fixed 10 m (geometry/elvis.py:_segments divides each stretch into
# round(length / 10) *equal* pieces).
MIN_RUN_LENGTH_M = 30.0
# A fiktiv line is held to a lower bar, and to no run length at all.
#
# The two tests above are proxies for "elvenett says there is a watercourse
# here and nothing has mapped its surface". On an ElvelinjeFiktiv that is
# not a proxy -- it is what the objektType means. ELVIS drew a line without
# a surface to draw it from, so the evidence the 20 m search and the 30 m
# run were standing in for is already on the record.
#
# They were also failing on exactly these: of 12,166 fiktiv segments in
# Glomma sør, 41 were flagged -- 0.3%, against 3,329 ElvBekk -- because a
# verge, a treeline or a field road breaks a 100% strip. A majority of each
# strip is enough to say the line crosses cultivated ground rather than
# running past it.
FIKTIV_JORDBRUK_MIN_PCT = 50.0
JORDBRUK_CLASS = "jordbruk"  # grunnkart arealdekkeniva1
LEVEL1 = "arealdekkeniva1"
# The urban twin of the farmland test. A stream piped through a town leaves
# the same trace as one piped through a field -- elvenett draws the line,
# no water is mapped under it -- but with buildings and roads on both
# banks instead of fields. Elvenett does not mark these as fiktiv: Alna
# has 0 ElvelinjeFiktiv, and 3,205 of its segments with no mapped water
# within 20 m sit between banks that are >= 50 % built; none were flagged
# and 2,368 of them were classed Nedbygd kantsone. A majority of each
# strip, not all of it: a city block always has a yard or a verge.
BEBYGD_CLASS = "bebygdSamferdsel"  # grunnkart arealdekkeniva1
BEBYGD_MIN_PCT = 50.0
BEBYGD_ENABLED = True  # the built-ground test above
FIKTIV_BEBYGD_ENABLED = True  # and the fiktiv shortcut through built ground


def _fully_farmland(
    seg: gpd.GeoDataFrame,
    grunnkart,
    min_pct: float = JORDBRUK_MIN_PCT,
    land_class: str = JORDBRUK_CLASS,
) -> pd.Series:
    """True where both side strips of the segment's line are jordbruk.

    Both strips rather than either: a stream with a field on one side and
    forest or bedrock on the other runs along the edge of farmland, and
    those keep their channel. The one that gets buried crosses the middle.

    ``min_pct`` is the share each strip must reach. The default demands all
    of it; a fiktiv line is held to FIKTIV_JORDBRUK_MIN_PCT instead, for
    the reason given there.
    """
    gk = as_gdf(grunnkart)
    if gk is None or gk.empty or LEVEL1 not in gk.columns:
        return pd.Series(False, index=seg.index)
    jb = gk[gk[LEVEL1].astype("string").str.strip() == land_class]
    if jb.empty:
        return pd.Series(False, index=seg.index)
    jb = jb.to_crs(seg.crs)[["geometry"]]

    covered = pd.Series(True, index=seg.index)
    for side in (1.0, -1.0):
        strip = shapely.buffer(
            seg.geometry.values, side * JORDBRUK_SIDE_M, single_sided=True
        )
        strips = gpd.GeoDataFrame(
            {"_i": np.arange(len(seg))}, geometry=strip, crs=seg.crs
        )
        strips = strips[shapely.area(strip) > 0]
        if strips.empty:
            return pd.Series(False, index=seg.index)
        ov = gpd.overlay(strips, jb, how="intersection", keep_geom_type=True)
        got = pd.Series(0.0, index=np.arange(len(seg)))
        if not ov.empty:
            got.update(ov.assign(_a=ov.geometry.area).groupby("_i")["_a"].sum())
        whole = pd.Series(shapely.area(strip), index=np.arange(len(seg)))
        pct = (got / whole.replace(0, np.nan) * 100).fillna(0.0)
        covered &= (pct >= min_pct - _PCT_EPS).to_numpy()
    return covered


def extract(segments, water, grunnkart) -> dict:
    """K_S_Rorlagt_mistanke: 'Ja' for a run of segments with no nearby
    mapped water, with jordbruk -- or built ground -- on both sides of the
    line, or for any fiktiv line crossing farmland, which needs neither
    test.
    """
    seg = as_gdf(segments)
    if seg.empty:
        return {"K_S_Rorlagt_mistanke": pd.Series(dtype="object")}
    water_gdf = as_gdf(water)

    no_water = pd.Series(True, index=seg.index)
    if len(water_gdf):
        tree = shapely.STRtree(water_gdf.geometry.values)
        near = tree.query(
            seg.geometry.values, predicate="dwithin", distance=GAP_SEARCH_M
        )
        no_water.iloc[np.unique(near[0])] = False

    candidate = no_water & _fully_farmland(seg, grunnkart)
    if BEBYGD_ENABLED:
        candidate |= no_water & _fully_farmland(
            seg, grunnkart, BEBYGD_MIN_PCT, BEBYGD_CLASS
        )

    # A fiktiv line through farmland -- or built ground -- is flagged on its
    # own, without the water search or the run length; see
    # FIKTIV_JORDBRUK_MIN_PCT and BEBYGD_MIN_PCT.
    fiktiv = (
        seg["objektType"].isin(C.FIKTIV_OBJEKTTYPER)
        if "objektType" in seg.columns
        else pd.Series(False, index=seg.index)
    )
    direct = pd.Series(False, index=seg.index)
    if fiktiv.any():
        direct = fiktiv & _fully_farmland(
            seg, grunnkart, min_pct=FIKTIV_JORDBRUK_MIN_PCT
        )
        # The same shortcut for a fiktiv line through built ground: the
        # objektType already says the surface is unmapped, and a majority
        # of buildings and roads on both strips says why.
        if FIKTIV_BEBYGD_ENABLED:
            direct |= fiktiv & _fully_farmland(
                seg, grunnkart, min_pct=BEBYGD_MIN_PCT, land_class=BEBYGD_CLASS
            )

    order = np.lexsort(
        (seg["segment_id"].to_numpy(), seg["stretch_id"].to_numpy())
    )
    stretch_o = seg["stretch_id"].to_numpy()[order]
    cand_o = candidate.to_numpy()[order]
    length_o = seg.geometry.length.to_numpy()[order]

    # One id per maximal run of same-stretch, same-candidate-value segments.
    new_run = np.empty(len(order), dtype=bool)
    new_run[0] = True
    new_run[1:] = (stretch_o[1:] != stretch_o[:-1]) | (
        cand_o[1:] != cand_o[:-1]
    )
    run_id = np.cumsum(new_run)

    run_length = pd.Series(length_o).groupby(run_id).transform("sum").to_numpy()
    flagged_o = cand_o & (run_length >= MIN_RUN_LENGTH_M)

    flag = np.empty(len(order), dtype=bool)
    flag[order] = flagged_o
    flag = flag | direct.to_numpy()

    return {
        "K_S_Rorlagt_mistanke": pd.Series(
            np.where(flag, "Ja", "Nei"), index=seg.index
        )
    }
