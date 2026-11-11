"""Leaf-litter input to the stream (næring).

System B, Næring block. Pure tabular.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..columns import num
from ..spec import constants as C
from .filtrering import _vegetation_situation


def compute(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    width = num(out, "K_S_bekk_bredde", None)
    # Narrow streams receive proportionally more litter per unit length; the
    # ratio saturates at the 20 m cap.
    reach = np.minimum(
        1.0,
        C.WIDTH_CAP_M / np.minimum(width.fillna(C.WIDTH_CAP_M), C.WIDTH_CAP_M),
    )

    for side in ("L", "R"):
        b = f"K_{side}_B0-10"
        situation = _vegetation_situation(out, side)
        kvalitet = situation.map(C.LOEV_KVALITET).astype(float).fillna(0.0)
        out[f"{b}_Loev_kvalitet"] = kvalitet.round(4)

        vtot = num(out, f"{b}_Vtot", 0.0) / 100.0
        score = vtot * kvalitet
        out[f"{b}_Naering_score"] = score.round(4)
        out[f"{b}_Naering_effekt"] = (score * reach).round(4)

    return out
