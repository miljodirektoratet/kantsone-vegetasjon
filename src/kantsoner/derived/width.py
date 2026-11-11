"""Area-equivalent riparian width and channel orientation.

Pure tabular. Inputs are the 0-20 m Vtot percentage plus the 0-10 m slope,
both produced by L2.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import sides as side_flags
from ..columns import num
from ..spec import constants as C


def _side_cols(side: str) -> dict[str, str]:
    return {
        "v20": f"K_{side}_B0-20_Vtot",
        "slope": f"K_{side}_B0-10_s",
    }


def _bredde(w_h: pd.Series) -> pd.Series:
    """The delivered class of a width (C.BREDDE_CLASSES)."""
    limits = [lim for lim, _ in C.BREDDE_CLASSES]
    labels = [lbl for _, lbl in C.BREDDE_CLASSES]
    idx = np.searchsorted(limits, w_h.fillna(0.0).to_numpy(), side="right")
    return pd.Series(np.array(labels, dtype=object)[idx], index=w_h.index)


def compute(df: pd.DataFrame) -> pd.DataFrame:
    """Add Bredde, w_h, w_eff (per side) and orient_fakt (combined).

    w_h is the area-equivalent width: WIDTH_BAND_M (20 m) times the share of
    the 0-20 m band that carries crowns. Crowns 0-3 m from the water and
    crowns 10-20 m from it count by their area alone, so a strip set back
    from the water still has a width, and a patchy strip has a fraction of
    one. This replaces a four-rung ladder (0, 3, 10, 20 m) that read 0
    whenever the first 3 m had under 50 % cover.
    """
    out = df.copy()

    for side in ("L", "R"):
        c = _side_cols(side)
        v20 = pd.to_numeric(out.get(c["v20"]), errors="coerce")

        # NaN is "no zone" and reads as no vegetation, as before, rather
        # than propagating NaN into every score.
        share = v20.fillna(0.0).clip(0.0, 100.0) / 100.0

        # A bank clipped away at the national border falls through the same
        # NaN path, and w_h = 0 is exactly the wrong answer for it: the
        # shortfall in priority.py is set by the *narrower* bank, so a
        # foreign bank with no data would decide the restoration priority of
        # every border segment. Blanked instead, which the minimum then
        # skips. A bank with no zone is the same case: nothing was measured.
        gone = side_flags.unmeasured(out, side)
        w_h = (share * C.WIDTH_BAND_M).astype("float64").mask(gone)
        out[f"K_{side}_w_h"] = w_h.round(4)
        out[f"K_{side}_Bredde"] = _bredde(w_h).mask(gone, "")

        # w_eff = w_h / cos(alpha), alpha = arctan(slope_pct / 100): the true
        # along-ground path length a flow crosses on a sloping bank.
        slope = pd.to_numeric(out.get(c["slope"]), errors="coerce").fillna(0.0)
        alpha = np.arctan(slope / 100.0)
        out[f"K_{side}_B0-10_w_eff"] = (w_h / np.cos(alpha)).round(4)

    azimut = num(out, "K_S_bekk_azimut", None)
    out["K_S_orient_fakt"] = (
        C.ORIENT_BASE + C.ORIENT_AMPLITUDE * np.sin(np.radians(azimut))
    ).round(4)
    return out
