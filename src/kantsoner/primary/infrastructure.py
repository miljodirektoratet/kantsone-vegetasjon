"""Flood and erosion protection on the banks, from NVE and Statens vegvesen.

Both are presence tests rather than coverage: a bank either carries an
engineering structure or it does not.

  Flomvern    "Ja"             flood protection intersects the 0-10 m zone
              "Tilbaketrukket" only the 0-20 m zone — a set-back embankment
              "Nei"            neither
  Forbygning  "Ja" / "Nei"     erosion protection intersects the 0-20 m zone

The combined `E_S_*` column scores each bank 0.5 and sums them, so a segment
protected on both sides scores 1.0.

  Barrier     "Ja" / "Nei"     an NVE dam in operation crosses the segment's
                               channel — tested against the claim (the
                               segment's piece of the routed channel line),
                               because the ELVIS line can sit well off the
                               water while the dam is on the water
"""

import geopandas as gpd
import pandas as pd
import shapely

from .. import sides
from ..geometry import segment_zones
from ..io import as_gdf

JA, NEI, TILBAKETRUKKET = "Ja", "Nei", "Tilbaketrukket"
# NVE dam status: D = i drift. V, U (ute av drift) and P1/P2 (planned) are not
# barriers today.
DAM_STATUS_ACTIVE = ("D",)
DAM_LINE_REACH_M = 10.0  # a dam line this close to the channel line crosses it
DAM_POINT_REACH_M = 15.0  # a dam point this close to the channel line is on it


def _load(*paths, crs):
    """Union of the supplied NVE/SVV layers, reprojected. None entries skipped."""
    frames = [as_gdf(p).to_crs(crs) for p in paths if p is not None]
    frames = [f[["geometry"]] for f in frames if len(f)]
    if not frames:
        return None
    return gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True), geometry="geometry", crs=crs
    )


def _hits(cells, ring, infra):
    """segment_ids whose zone slice intersects any infrastructure geometry."""
    sl = segment_zones.cut(cells, ring)
    if infra is None or sl.empty:
        return set()
    idx = sl.sindex.query(infra.geometry, predicate="intersects")[1]
    return set(sl.iloc[idx]["segment_id"])


def extract(
    segments, zones_near, zones_wide, flomsikring=(), erosjonssikring=()
):
    """Flomvern and Forbygning columns, keyed on segment_id.

    `zones_near` is (left, right, ring) at 0-10 m, `zones_wide` the same at
    0-20 m.
    """
    seg = as_gdf(segments)
    ids = seg["segment_id"]
    crs = zones_near[0].crs
    flood = _load(*flomsikring, crs=crs)
    erosion = _load(*erosjonssikring, crs=crs)

    out = {}
    for i, side in enumerate(("L", "R")):
        near = _hits(zones_near[i], zones_near[2], flood)
        wide = _hits(zones_wide[i], zones_wide[2], flood)
        out[f"E_{side}_Flomvern"] = pd.Series(
            [
                JA if s in near else (TILBAKETRUKKET if s in wide else NEI)
                for s in ids
            ],
            index=ids.index,
        )
        eros = _hits(zones_wide[i], zones_wide[2], erosion)
        out[f"E_{side}_Forbygning"] = pd.Series(
            [JA if s in eros else NEI for s in ids], index=ids.index
        )

    # Each bank is worth half, so a segment protected on both scores 1.0.
    # A bank abroad has no zone to test and would always read "Nei",
    # capping every border segment at 0.5 however heavily engineered the
    # Norwegian side is. Its half is redistributed instead: the score stays
    # on 0..1 and means the same thing on one bank as on two.
    scale = 2.0 / sides.n_present(seg)
    for name, score in (
        ("Flomvern", {JA: 0.5, TILBAKETRUKKET: 0.25, NEI: 0.0}),
        ("Forbygning", {JA: 0.5, NEI: 0.0}),
    ):
        total = sum(
            sides.mask(seg, s, out[f"E_{s}_{name}"].map(score)).fillna(0.0)
            for s in sides.SIDES
        )
        out[f"E_S_{name}"] = (total * scale).round(4)
    return out


def barrier(segments, claims, dam_line=None, dam_punkt=None):
    """E_S_Barrier: Ja on the one segment each NVE dam sits on.

    One dam, one segment: the aggregated levels count "Ja" segments as
    *Antall* barriers, so a dam line buffered across three consecutive 10 m
    claims — and present again in the point layer — would count three or
    four times and grade every dammed river 1. Dams are collapsed on damNr
    (the point kept, being the more precise location) and each is assigned
    to its nearest channel line within reach.

    Matched against the claims; a segment without one (its stretch never
    routed) falls back to its ELVIS line. Only dams in operation count.
    """
    seg = as_gdf(segments)
    ids = seg["segment_id"]
    crs = seg.crs
    frames = []
    # Points first: with duplicates collapsed on damNr, the point wins.
    for path, reach in (
        (dam_punkt, DAM_POINT_REACH_M),
        (dam_line, DAM_LINE_REACH_M),
    ):
        if path is None:
            continue
        dams = as_gdf(path).to_crs(crs)
        if "status" in dams.columns:
            dams = dams[dams["status"].astype(str).isin(DAM_STATUS_ACTIVE)]
        if len(dams):
            key = dams["damNr"] if "damNr" in dams.columns else dams.index
            frames.append(
                gpd.GeoDataFrame(
                    {"damNr": key.astype(str).values, "_reach": reach},
                    geometry=dams.geometry.values,
                    crs=crs,
                )
            )
    hit: set = set()
    if frames:
        dams = gpd.GeoDataFrame(
            pd.concat(frames, ignore_index=True), geometry="geometry", crs=crs
        ).drop_duplicates("damNr", keep="first")
        cl = as_gdf(claims)
        have = set(cl["segment_id"]) if len(cl) else set()
        rest = seg[~ids.isin(have)]
        lines = gpd.GeoDataFrame(
            pd.concat(
                [
                    cl[["segment_id", "geometry"]]
                    if len(cl)
                    else rest.iloc[0:0][["segment_id", "geometry"]],
                    rest[["segment_id", "geometry"]],
                ],
                ignore_index=True,
            ),
            geometry="geometry",
            crs=crs,
        )
        tree = shapely.STRtree(lines.geometry.values)
        for reach, group in dams.groupby("_reach"):
            near = tree.query_nearest(
                group.geometry.values,
                max_distance=float(reach),
                all_matches=False,
            )
            hit |= set(lines["segment_id"].values[near[1]])
    return {
        "E_S_Barrier": pd.Series(
            [JA if s in hit else NEI for s in ids], index=ids.index
        )
    }
