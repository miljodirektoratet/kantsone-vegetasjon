"""Connectivity and restoration priority.

System B. ``K_lateral`` is local; ``K_oppstroems`` and ``K_nettverk`` need the
directed graph and a spatial index respectively.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from loguru import logger

from ..columns import num
from ..spec import constants as C


def k_lateral(df: pd.DataFrame) -> pd.Series:
    """Share of the 0-20 m zone that is semi-natural, both banks pooled."""
    natural, total = 0.0, 0.0
    for side in ("L", "R"):
        b = f"K_{side}_B0-20"
        for kind in ("skog", "myr", "Aapen"):
            natural = natural + num(df, f"{b}_Areal_{kind}", 0.0)
        total = total + num(df, f"{b}_Areal", 0.0)
    return (natural / total.replace(0, np.nan)).clip(0, 1)


def k_oppstroems(df: pd.DataFrame, graph) -> pd.Series:
    """Distance-decayed mean of Score_total over the upstream network."""
    score = num(df, "E_S_Score_total", 0.0)
    by_id = dict(zip(df["segment_id"], score, strict=False))
    out = {}
    horizon = C.CONNECTIVITY_DECAY_M * 4  # beyond this the weight is < 2 %
    for sid in df["segment_id"]:
        acc, n = 0.0, 0
        for other, dist in graph.walk_upstream(sid, max_distance=horizon):
            if other == sid:
                continue
            acc += by_id.get(other, 0.0) * np.exp(
                -dist / C.CONNECTIVITY_DECAY_M
            )
            n += 1
        out[sid] = acc / n if n else 0.0
    return df["segment_id"].map(out)


def k_nettverk(df: pd.DataFrame, segments) -> pd.Series:
    """Distance-decayed mean of Score_total within a 5 km radius, ignoring flow."""
    score = num(df, "E_S_Score_total", 0.0)
    pts = segments.geometry.centroid
    tree = pts.sindex
    xy = np.column_stack([pts.x.to_numpy(), pts.y.to_numpy()])
    values = score.to_numpy()

    out = np.zeros(len(segments))
    for i, geom in enumerate(pts):
        idx = list(tree.query(geom.buffer(C.NETWORK_RADIUS_M)))
        idx = [j for j in idx if j != i]
        if not idx:
            continue
        d = np.hypot(xy[idx, 0] - xy[i, 0], xy[idx, 1] - xy[i, 1])
        out[i] = np.mean(values[idx] * np.exp(-d / C.CONNECTIVITY_DECAY_M))
    return pd.Series(out, index=df.index)


def compute(df: pd.DataFrame, graph, segments) -> pd.DataFrame:
    out = df.copy()
    out["K_S_K_lateral"] = k_lateral(out).round(4)
    out["K_S_K_oppstroems"] = k_oppstroems(out, graph).round(4)
    out["K_S_K_nettverk"] = k_nettverk(out, segments).round(4)

    out["E_S_K_score"] = (
        out[["K_S_K_lateral", "K_S_K_oppstroems", "K_S_K_nettverk"]]
        .mean(axis=1, skipna=True)
        .fillna(0.0)
        .round(4)
    )
    logger.info("Connectivity computed")
    return out
