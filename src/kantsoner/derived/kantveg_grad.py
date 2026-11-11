"""Kantvegetasjon graded: how much of each bank's expected cover is there.

Agreed with the client 2026-10-23, replacing the yes/no reading of a single
Vtot threshold. Per bank, on the 0-10 m zone:

    forventet  = K_{side}_B0-10_Forventet_dekning   (primary/landcover.py)
    vurdert    = Vegdyktig >= C.VEGDYKTIG_MIN_PCT, or Vtot >= forventet
    fjernet    = max(0, (forventet - Vtot) / forventet) x 100

and per segment the two banks weighted half each -- one bank counts at most
half, as the veileder says. A bank that is not assessed contributes 0; a
bank that does not exist (lake shore, the far side of a border river)
leaves the other with full weight. The segment's class is fjernet put
through veilederens tabell 4. Bekkelukking and segments with no zone are
overridden at the aggregation stage, where Beregnet and the culvert flag
are final.

The aggregated value per unit (delivery/aggregate.py) is the length-weighted
mean of fjernet over assessed segments -- not the share of segments in each
class -- classified only where data coverage reaches C.DATADEKNING_MIN.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import sides
from ..columns import num
from ..spec import constants as C
from ..spec import scales

GRADE_NAMES = {
    5: "Svært god",
    4: "God",
    3: "Moderat",
    2: "Dårlig",
    1: "Svært dårlig",
}


def grade_name(pct: pd.Series) -> pd.Series:
    g = scales.pct_to_grade(pct)
    return pd.Series(
        [GRADE_NAMES.get(int(x)) if pd.notna(x) else None for x in g],
        index=pct.index,
    )


def compute(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    fj, vur, absent = {}, {}, {}
    for side in ("L", "R"):
        p = f"K_{side}_B0-10_"
        vtot = num(out, p + "Vtot", None)
        forv = num(out, p + "Forventet_dekning", None)
        vegd = num(out, p + "Vegdyktig", None)
        absent[side] = sides.absent(out, side)
        v = (
            ((vegd >= C.VEGDYKTIG_MIN_PCT) | ((forv > 0) & (vtot >= forv)))
            & vtot.notna()
            & ~absent[side]
        )
        f = ((forv - vtot) / forv.where(forv > 0)).clip(lower=0) * 100.0
        fj[side] = f.where(v).round(1)
        vur[side] = v
        out[f"E_{side}_Kantvegetasjon_fjernet"] = fj[side]
    # Weights: half each; a bank that does not exist hands its half over.
    wL = np.where(absent["R"], 1.0, np.where(absent["L"], 0.0, 0.5))
    wR = 1.0 - wL
    seg = wL * fj["L"].fillna(0.0) + wR * fj["R"].fillna(0.0)
    any_v = vur["L"] | vur["R"]
    seg = pd.Series(seg, index=out.index).where(any_v).round(1)
    out["E_S_Kantvegetasjon_fjernet"] = seg
    out["E_S_Kantvegetasjon_grad"] = grade_name(seg).where(
        any_v, C.KANTVEG_GRAD_IKKE_VURDERT
    )
    return out
