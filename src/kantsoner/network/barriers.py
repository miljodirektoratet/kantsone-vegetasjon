"""Barrier-driven continuity: Oppstuving and Fragmentering.

System A. Both need the directed graph from :mod:`kantsoner.geometry.graph`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from loguru import logger

from ..columns import num
from ..spec import constants as C


def _barrier_ids(df: pd.DataFrame) -> set[int]:
    col = df.get("E_S_Barrier")
    if col is None:
        return set()
    hit = col.astype(str).str.strip().str.lower().eq("ja")
    return set(df.loc[hit, "segment_id"])


def oppstuving(df: pd.DataFrame, graph) -> pd.Series:
    """Flag backwater reaches upstream of each barrier.

    Reference width is the mean channel width over 200 m *downstream* of the
    barrier; walking upstream, a segment is impounded while its width is at
    least 1.5x that reference. Gaps of up to two segments below the ratio are
    bridged, so a locally narrow spot does not truncate the pool.
    """
    width = num(df, "K_S_bekk_bredde", None)
    width_by_id = dict(zip(df["segment_id"], width, strict=False))
    flagged: set[int] = set()

    for barrier in _barrier_ids(df):
        downstream = [
            width_by_id.get(sid, np.nan)
            for sid, _ in graph.walk_downstream(
                barrier, max_distance=C.OPPSTUVING_REFERENCE_M
            )
        ]
        reference = np.nanmean(downstream) if downstream else np.nan
        if not np.isfinite(reference) or reference <= 0:
            continue

        gap = 0
        pending: list[int] = []
        for sid, _ in graph.walk_upstream(barrier):
            w = width_by_id.get(sid, np.nan)
            if np.isfinite(w) and w / reference >= C.OPPSTUVING_WIDTH_RATIO:
                flagged.update(pending)
                pending = []
                flagged.add(sid)
                gap = 0
            else:
                gap += 1
                if gap > C.OPPSTUVING_MAX_GAP_SEGMENTS:
                    break
                pending.append(sid)

    logger.info(f"Oppstuving: {len(flagged)} segments flagged")
    return pd.Series(
        np.where(df["segment_id"].isin(flagged), "Ja", "Nei"),
        index=df.index,
    )


def fragmentering(df: pd.DataFrame, graph) -> pd.Series:
    """Classify each segment by how far upstream the nearest barrier is.

    As specified: "for each segment, check distance walking upstream to a
    barrier". A segment below a dam is the fragmented one — cut off from the
    sediment, wood and organisms the reach above would supply. Computed as
    one downstream walk per barrier rather than one upstream walk per
    segment; the minimum is the same.
    """
    barriers = _barrier_ids(df)
    best: dict = {}
    for barrier in barriers:
        for sid, dist in graph.walk_downstream(barrier):
            if dist < best.get(sid, np.inf):
                best[sid] = dist

    def classify(sid) -> str:
        d = best.get(sid)
        if d is None:
            return C.FRAGMENTERING_NO_BARRIER
        for upper, label in C.FRAGMENTERING_CLASSES:
            if d < upper:
                return label
        return C.FRAGMENTERING_CLASSES[-1][1]

    return df["segment_id"].map(classify)


def compute(df: pd.DataFrame, graph) -> pd.DataFrame:
    out = df.copy()
    out["E_S_Oppstuving"] = oppstuving(out, graph)
    out["E_S_Fragmentering"] = fragmentering(out, graph)
    return out
