"""Restoration priority.

System B, prioritet block. Compares the riparian width a segment *has*
(``w_h``) against the width it *ought* to have (``w_opt``), weighted by how
much the gap matters.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import sides
from ..columns import arealbruk as arealbruk_of
from ..columns import is_lake, num
from ..spec import constants as C


def _dominant_arealbruk(df: pd.DataFrame) -> pd.Series:
    left = arealbruk_of(df, "L")
    right = arealbruk_of(df, "R")
    # Take the bank whose land use carries the higher restoration weight.
    # A bank abroad has no land use at all -- not an unrecognised one -- so
    # it must not win the comparison on F_AREALBRUK_DEFAULT, which outranks
    # Skog. Nor may a bank with no zone, which was not measured.
    lw = left.map(C.F_AREALBRUK).fillna(C.F_AREALBRUK_DEFAULT)
    rw = right.map(C.F_AREALBRUK).fillna(C.F_AREALBRUK_DEFAULT)
    lw = lw.mask(sides.unmeasured(df, "L"), -np.inf)
    rw = rw.mask(sides.unmeasured(df, "R"), -np.inf)
    out = pd.Series(np.where(lw >= rw, left, right), index=df.index)
    # No bank with a measurement at all: the report's conservative default.
    nothing = np.isneginf(lw) & np.isneginf(rw)
    return out.mask(nothing, C.AREALBRUK_UNKNOWN)


def fremmedart(df: pd.DataFrame, graph, observations=None):
    """Invasive-species pressure, and the distance it comes from.

    Returns ``(risiko, d_m, arter)``: 0.8 in a segment with an occurrence,
    ``0.8 * exp(-d / 5 km)`` decaying from the nearest one upstream, else 0;
    the network distance in m to that nearest upstream occurrence (0
    in-segment, NaN when there is none within reach); and the species
    recorded in the segment itself.

    Walked *downstream from each observation*, not upstream from each
    segment: a few hundred observations against a million segments is the
    cheap direction, and it gives the same minimum. The walk stops where the
    decay drops below 1 % (FREMMEDART_REACH_KM).
    """
    risiko = pd.Series(0.0, index=df.index)
    d_m = pd.Series(np.nan, index=df.index)
    arter = pd.Series("", index=df.index, dtype="object")
    if observations is None or not len(observations):
        return risiko, d_m, arter

    obs = observations.groupby("segment_id")["art"].agg(
        lambda s: ", ".join(sorted(set(s)))
    )
    best: dict = {}
    reach = C.FREMMEDART_REACH_KM * 1000.0
    for sid in obs.index:
        for other, dist in graph.walk_downstream(sid, max_distance=reach):
            if dist < best.get(other, np.inf):
                best[other] = dist

    ids = df["segment_id"]
    dist = ids.map(best)
    d_m = dist.astype(float)
    risiko = pd.Series(
        np.where(
            dist.notna(),
            C.FREMMEDART_IN_SEGMENT
            * np.exp(-dist.fillna(0.0) / 1000.0 / C.FREMMEDART_DECAY_KM),
            0.0,
        ),
        index=df.index,
    )
    arter = ids.map(obs).fillna("").astype("object")
    return risiko.round(4), d_m.round(1), arter


def compute(df: pd.DataFrame, graph=None, observations=None) -> pd.DataFrame:
    out = df.copy()

    k_score = num(out, "E_S_K_score", 0.0)
    out["E_S_f_konnektivitet"] = (
        C.F_KONNEKTIVITET_BASE - C.F_KONNEKTIVITET_GAIN * k_score
    ).round(4)

    arealbruk = _dominant_arealbruk(out)
    out["E_S_f_arealbruk"] = (
        arealbruk.map(C.F_AREALBRUK).fillna(C.F_AREALBRUK_DEFAULT).round(4)
    )

    slope = (
        pd.concat(
            [
                sides.mask(out, s, num(out, f"K_{s}_B0-10_s", None), measured=True)
                for s in ("L", "R")
            ],
            axis=1,
        )
        .mean(axis=1)
        .fillna(0.0)
    )
    capped = np.minimum(slope, C.F_HELNING_MAX_PCT)
    out["E_S_f_helning"] = (1.0 / np.cos(np.arctan(capped / 100.0))).round(4)

    width = num(out, "K_S_bekk_bredde", 0.0)
    w_basis = np.clip(
        C.W_BASIS_INTERCEPT + C.W_BASIS_WIDTH_GAIN * width,
        C.W_BASIS_MIN_M,
        C.W_BASIS_MAX_M,
    )
    # A lake has no channel width. Left to the line above it would reach
    # the same 10 m only because bekk_bredde is NaN and num() substitutes
    # zero -- the right answer resting on a missing value, which the next
    # change to either would break silently. Pinned instead.
    w_basis = np.where(is_lake(out), C.LAKE_W_BASIS_M, w_basis)
    out["K_S_w_basis"] = pd.Series(w_basis, index=out.index).round(4)

    w_opt = (
        w_basis
        * out["E_S_f_helning"]
        * out["E_S_f_arealbruk"]
        * out["E_S_f_konnektivitet"]
    )
    out["K_S_w_opt"] = w_opt.round(4)

    # w_h is per side; the narrower bank sets the shortfall.
    w_h = (
        pd.concat(
            [num(out, f"K_{s}_w_h", None) for s in ("L", "R")],
            axis=1,
        )
        .min(axis=1)
        .fillna(0.0)
    )
    gap_rel = np.maximum(0.0, (w_opt - w_h) / w_opt.replace(0, np.nan))
    out["K_S_gap_rel"] = (
        pd.Series(gap_rel, index=out.index).fillna(0.0).round(4)
    )

    if graph is not None:
        risiko, d_m, arter = fremmedart(out, graph, observations)
        out["K_S_fremmedart_risiko"] = risiko
        out["K_S_d"] = d_m
        out["K_S_fremmedart_art"] = arter

    score_total = num(out, "E_S_Score_total", 0.0)
    # The last factor is the fremmedart term. A confirmed occurrence in the
    # segment itself (risiko ~ 0.8) nearly doubles the priority, because an
    # invasive stand is a restoration need in its own right, independent of
    # how the kantsone is otherwise doing. Its absence went unnoticed
    # because it is a no-op wherever Artskart has no records -- true of
    # every subarctic catchment, and of Tana in particular, where it was
    # first looked for. It is not true generally: Sandvikselva has 4,840
    # segments with risiko > 0 and Orreaana 3,007, both carrying
    # kjempespringfroe and parkslirekne.
    out["K_S_prioritet"] = (
        out["K_S_gap_rel"]
        * out["E_S_f_arealbruk"]
        * (C.PRIORITET_BASE + C.PRIORITET_GAIN * (1.0 - k_score))
        * (1.0 - score_total)
        * (1.0 + num(out, "K_S_fremmedart_risiko", 0.0))
    ).round(4)
    return out
