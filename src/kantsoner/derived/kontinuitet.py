"""Lateral continuity of the riparian strip.

System B, Kontinuitet block. Pure tabular.

The workbook's Kontinuitet rows all name ``K_R_B0-10_Loev_kvalitet`` as the
right-side column — a copy-paste slip, since the formulas reference
``Kontin_tre``/``Kontin_busk``. The right side is computed symmetrically here.
"""

from __future__ import annotations

import pandas as pd

from ..columns import num
from ..spec import constants as C


def compute(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for side in ("L", "R"):
        b = f"K_{side}_B0-10"
        tre = (
            num(out, f"{b}_V5-10", 0.0) + num(out, f"{b}_V10pluss", 0.0)
        ) / 100.0
        busk = (
            num(out, f"{b}_V0p5-1", 0.0) + num(out, f"{b}_V1-5", 0.0)
        ) / 100.0
        out[f"{b}_Kontin_tre"] = tre.round(4)
        out[f"{b}_Kontin_busk"] = busk.round(4)
        out[f"{b}_Kontin"] = (
            C.KONTIN_TRE_WEIGHT * tre + C.KONTIN_BUSK_WEIGHT * busk
        ).round(4)
    return out
