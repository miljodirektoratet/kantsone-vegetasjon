"""Shading and stream temperature.

System B, Temperatur block, following the specification formula::

    rekk          = Vhmed * tan(37 deg) * orient_fakt
    Skygge_effekt = min(1, (rekk_L + rekk_R) / min(bekk_bredde, 20))
    Temp_score    = min(1, max(Vp95_L, Vp95_R) / 15)
                    * mean(Lauv_faktor_L, Lauv_faktor_R)
                    * Skygge_effekt

The retired ``calculate_shadow_factor.py`` computed an equivalent
``skygge_effekt`` per nDSM cell, projecting each canopy pixel's shadow and
subtracting its set-back distance. That is better physics but was specified
away in favour of the closed form above; it is recoverable from git history if
the simpler version proves too coarse.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import sides
from ..columns import num
from ..spec import constants as C


def _lauv_faktor(df: pd.DataFrame, side: str) -> pd.Series:
    """1.0 where tall canopy dominates, 0.7 where low growth does, else 0."""
    b = f"K_{side}_B0-10"
    vtot = num(df, f"{b}_Vtot", 0.0)
    tall = num(df, f"{b}_V5-10", 0.0) + num(df, f"{b}_V10pluss", 0.0)
    low = num(df, f"{b}_V1-5", 0.0) + num(df, f"{b}_V0p5-1", 0.0)
    has_cover = vtot >= C.KANTVEG_THRESHOLD_PCT
    return pd.Series(
        np.where(
            has_cover & (tall >= low),
            C.LAUV_FAKTOR_TALL,
            np.where(has_cover, C.LAUV_FAKTOR_LOW, C.LAUV_FAKTOR_NONE),
        ),
        index=df.index,
    )


def compute(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    orient = num(out, "K_S_orient_fakt", 1.0)

    rekk = {}
    for side in ("L", "R"):
        b = f"K_{side}_B0-10"
        vhmed = num(out, f"{b}_Vhmed", 0.0)
        rekk[side] = vhmed * C.TAN_SHADE_ANGLE * orient
        out[f"{b}_rekk"] = rekk[side].round(4)
        out[f"{b}_Lauv_faktor_temp"] = _lauv_faktor(out, side).round(4)

    # Wide channels cannot be shaded across, so the denominator is capped:
    # beyond 20 m extra width no longer dilutes the effect.
    width = num(out, "K_S_bekk_bredde", None)
    denom = np.minimum(width.fillna(C.WIDTH_CAP_M), C.WIDTH_CAP_M)
    denom = denom.where(denom > 0, np.nan)

    # Both banks shade the same channel, so their reach adds. Where one is
    # abroad only the Norwegian bank's reach is known; scaled up to the two
    # banks the formula is written for, rather than counting the unknown one
    # as casting no shadow. n_present is 2 everywhere but a border river.
    scale = 2.0 / sides.n_measured(out)
    reach = sum(
        sides.mask(out, s, rekk[s], measured=True).fillna(0.0)
        for s in sides.SIDES
    )
    skygge = np.minimum(1.0, reach * scale / denom)
    out["K_S_B0-10_Skygge_effekt"] = (
        pd.Series(skygge, index=out.index).fillna(0.0).round(4)
    )

    p95 = (
        pd.concat(
            [
                sides.mask(
                    out, s, num(out, f"K_{s}_B0-10_Vp95", None), measured=True
                )
                for s in sides.SIDES
            ],
            axis=1,
        )
        .max(axis=1)
        .fillna(0.0)
    )
    height_term = np.minimum(1.0, p95 / C.TEMP_REFERENCE_HEIGHT_M)
    lauv_mean = (
        pd.concat(
            [
                sides.mask(
                    out, s, out[f"K_{s}_B0-10_Lauv_faktor_temp"], measured=True
                )
                for s in sides.SIDES
            ],
            axis=1,
        )
        .mean(axis=1, skipna=True)
        .fillna(0.0)
    )

    out["K_S_B0-10_Temp_score"] = (
        height_term * lauv_mean * out["K_S_B0-10_Skygge_effekt"]
    ).round(4)
    return out
