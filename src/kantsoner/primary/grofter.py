"""Ditches draining into the riparian zone.

A ditch that runs *across* the bank carries field runoff straight past the
vegetation and into the river, so the buffer never gets to filter it. One
running *along* the bank does not. That is why the count is conditioned on
the angle between the ditch and the river rather than being a plain density:
`Grofte_tetthet` is meant to count the ditches that bypass the buffer.

The workbook writes that condition as
``abs(GK_B0-10_orient_groft - K_S_orient_fakt) >= 45``, comparing two
*orientation factors*. Those are `0.7 + 0.3 sin(azimuth)`, so they live in
[0.4, 1.0] and their difference can never reach 45. Implemented on the
azimuths themselves, in degrees, which is what the threshold plainly means.
"""

import numpy as np
import pandas as pd

from ..geometry import segment_zones
from ..io import as_gdf
from ..spec import constants as C
from ..spec.attributes import band_token


def _azimuth(geom):
    """Bearing of a ditch, degrees clockwise from north, folded onto 0-180."""
    coords = list(geom.coords) if geom.geom_type == "LineString" else None
    if not coords or len(coords) < 2:
        return np.nan
    (x0, y0), (x1, y1) = coords[0], coords[-1]
    return float(np.degrees(np.arctan2(x1 - x0, y1 - y0)) % 180.0)


def _angular_gap(a, b):
    """Smallest angle between two undirected bearings, 0-90 degrees."""
    d = (a - b).abs() % 180.0
    return np.minimum(d, 180.0 - d)


def extract(segments, zones, band, ring, grofter, bekk_azimut):
    """Ditch orientation per segment and bypassing-ditch count per side."""
    ids = as_gdf(segments)["segment_id"]
    token = band_token(band)
    empty = {
        f"GK_{token}_groft_azimut": pd.Series(np.nan, index=ids.index),
        f"GK_{token}_orient_groft": pd.Series(np.nan, index=ids.index),
        f"K_L_{token}_Grofte_tetthet": pd.Series(0.0, index=ids.index),
        f"K_R_{token}_Grofte_tetthet": pd.Series(0.0, index=ids.index),
    }
    gr = as_gdf(grofter)
    if gr is None or gr.empty:
        return empty
    gr = gr.to_crs(zones[0].crs)
    gr = gr[gr.geometry.type == "LineString"]
    if gr.empty:
        return empty
    gr = gr.assign(_az=[_azimuth(g) for g in gr.geometry])
    gr = gr[gr["_az"].notna()]
    if gr.empty:
        return empty

    river_az = ids.index.map(pd.Series(bekk_azimut).reindex(ids.index))
    river_az = pd.Series(np.asarray(river_az, dtype="float64"), index=ids.index)

    per_side, az_parts = {}, []
    for side, cells in (("L", zones[0]), ("R", zones[1])):
        sl = segment_zones.cut(cells, ring)
        count = pd.Series(0.0, index=ids.index)
        if not sl.empty:
            hit = gr.sjoin(
                sl[["segment_id", "geometry"]],
                how="inner",
                predicate="intersects",
            )
            if not hit.empty:
                hit = hit.assign(segment_id=hit["segment_id"])
                az = hit.groupby("segment_id")["_az"].mean()
                az_parts.append(az)
                gap = _angular_gap(
                    hit["_az"],
                    hit["segment_id"].map(
                        pd.Series(river_az.values, index=ids.values)
                    ),
                )
                bypass = hit[gap >= C.GROFT_BYPASS_ANGLE_DEG]
                if not bypass.empty:
                    count = (
                        ids.map(bypass.groupby("segment_id").size())
                        .fillna(0.0)
                        .astype(float)
                    )
        per_side[side] = count

    az_all = (
        pd.concat(az_parts).groupby(level=0).mean()
        if az_parts
        else pd.Series(dtype="float64")
    )
    azimut = ids.map(az_all)
    azimut = pd.Series(np.asarray(azimut, dtype="float64"), index=ids.index)
    orient = C.ORIENT_BASE + C.ORIENT_AMPLITUDE * np.sin(np.radians(azimut))

    return {
        f"GK_{token}_groft_azimut": azimut.round(2),
        f"GK_{token}_orient_groft": orient.round(4),
        f"K_L_{token}_Grofte_tetthet": per_side["L"].round(2),
        f"K_R_{token}_Grofte_tetthet": per_side["R"].round(2),
    }
