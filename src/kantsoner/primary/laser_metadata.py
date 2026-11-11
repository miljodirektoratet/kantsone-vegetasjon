"""Laser survey year per segment, from NHM_ProsjektDekning.

NHM_ProsjektDekning (Nasjonal database for høydedata) is a polygon-per-
laser-survey-project coverage layer, each carrying the acquisition year in
AARSTALL. Projects abut and overlap -- a segment straddling two of them
takes the year of whichever polygon covers the greater length of it, not an
average or the first match, the same "dominance" rule landcover.py uses for
Arealbruk/Treslag.
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd

from ..io import as_gdf


def extract(segments, prosjektdekning) -> dict:
    """{'K_S_Laser_aar': Series} keyed to segments' own index, one value per
    segment.

    ``segments`` are the 10 m ELVIS pieces (LineStrings); ``prosjektdekning``
    is NHM_ProsjektDekning already read for the segments' own footprint, or
    None if the source is not configured -- an optional attribute like
    Grofte_tetthet, absent rather than wrong when the source is missing.
    """
    segments_gdf = as_gdf(segments)
    ids = segments_gdf["segment_id"]
    col = "K_S_Laser_aar"

    if prosjektdekning is None:
        return {col: pd.Series(pd.NA, index=ids.index, dtype="Int64")}
    pd_gdf = as_gdf(prosjektdekning)

    if pd_gdf.empty:
        return {col: pd.Series(pd.NA, index=ids.index, dtype="Int64")}
    pd_gdf = pd_gdf.to_crs(segments_gdf.crs)

    ov = gpd.overlay(
        segments_gdf[["segment_id", "geometry"]],
        pd_gdf[["AARSTALL", "geometry"]],
        how="intersection",
        keep_geom_type=True,
    )
    if ov.empty:
        return {col: pd.Series(pd.NA, index=ids.index, dtype="Int64")}

    ov["_len"] = ov.geometry.length
    top = (
        ov.groupby(["segment_id", "AARSTALL"])["_len"]
        .sum()
        .reset_index()
        .sort_values("_len")
        .groupby("segment_id")
        .last()
    )
    laser_aar = ids.map(top["AARSTALL"]).astype("Int64")
    return {col: laser_aar}
