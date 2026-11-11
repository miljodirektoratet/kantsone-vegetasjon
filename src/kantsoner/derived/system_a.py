"""System A — the Vannportalen riparian-vegetation classification.

PotVeg's own cutoff (potensiell) is a flat 50%; has_veg's Vtot cutoff is
per segment/side instead, read from K_{side}_{band}_VtotThreshold
(landcover.veg_threshold) -- falls back to the flat 50% wherever that
column is absent or NaN. NaN compares below threshold either way, which
means a side with no riparian zone resolves to "Naturlig fravær av
kantskog" rather than propagating nulls -- Uklassifisert (the actual "we
could not compute this" case) is applied afterward, once Matched is known;
see stages.py:aggregation.

Computed at both B0-10 and B0-20 (0-20's has_veg threshold scaled down per
KANTVEG_WIDE_THRESHOLD_SCALE -- a wider ring dilutes the same real cover
into a lower percentage), then selected per segment into a band-free
E_{side}_Kantvegetasjon_A / Potensiell_Kantvegetasjon pair by
K_S_beregningsbredde (primary/channel.py): whichever band a stretch's
measured channel width says is the right one to measure from. Everything
downstream -- delivery/aggregate.py, the Uklassifisert override in
stages.py:aggregation -- reads the band-free pair, not either raw band.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import sides
from ..columns import num
from ..spec import constants as C
from . import kantveg

THRESH = C.KANTVEG_THRESHOLD_PCT
_BANDS = (("B0-10", 1.0), ("B0-20", C.KANTVEG_WIDE_THRESHOLD_SCALE))


def _compute_band(out: pd.DataFrame, b: str, thresh_scale: float) -> dict:
    """Potensiell_Kantvegetasjon, Kantvegetasjon_A and its arsak/flag, per
    side and combined, for one band."""
    cols = {}
    has_veg_by_side, potensiell_by_side, infra_by_side = {}, {}, {}
    for side in ("L", "R"):
        vtot = num(out, f"K_{side}_{b}_Vtot", None)
        # PotVeg (forest + agriculture + soil-covered Snaumark incl. the
        # frisk classes) supersedes the narrower SkogJordbruk as the
        # potential-vegetation measure; the fallback keeps --allow-missing
        # runs on old attribute files working.
        potveg = num(out, f"K_{side}_{b}_PotVeg", None)
        if potveg.isna().all():
            potveg = num(out, f"K_{side}_{b}_SkogJordbruk", None)

        thresh = num(out, f"K_{side}_{b}_VtotThreshold", THRESH) * thresh_scale
        has_veg = vtot >= thresh
        # Built ground counts as potential too. Asphalt is not forest,
        # agriculture or soil-covered Snaumark, so PotVeg reads it as
        # incapable and the bank fell through to natural absence -- see
        # C.INFRA_POTENSIELL_MIN_PCT. Decided here rather than in the
        # classifier so Potensiell_Kantvegetasjon agrees with the class it
        # produces, and so the aggregated re-derivation can apply the same
        # rule to summed areas.
        built = sum(
            num(out, f"K_{side}_{b}_{c}", 0.0) for c in kantveg.ARSAK_BUILT
        )
        infra_rescued = (
            (built >= C.INFRA_POTENSIELL_MIN_PCT)
            & ~(potveg >= THRESH)
            & ~has_veg
        )
        potensiell = (potveg >= THRESH) | has_veg | infra_rescued
        infra_by_side[side] = infra_rescued
        has_veg_by_side[side] = has_veg
        potensiell_by_side[side] = potensiell
        cols[f"E_{side}_{b}_Potensiell_Kantvegetasjon"] = sides.mask(
            out,
            side,
            pd.Series(np.where(potensiell, "Ja", "Nei"), index=out.index),
        )

    # Status, arsak and the human-intervention flag come out of the same
    # decision per side: which group of competing land uses -- Jordbruk/
    # Skogbruk, or Bygg/Samferdsel/Konstruert -- covers more of a
    # Manglende-equivalent side decides Kantskog_fjernet vs Nedbygd_kantsone;
    # Berg/Myr/Vann/Snoisbre only compete for the Ikke_mulig label. See
    # kantveg.classify_sides.
    for side in ("L", "R"):
        pcts = np.column_stack(
            [
                num(out, f"K_{side}_{b}_{cand}", 0.0).to_numpy()
                for cand in kantveg.ARSAK_LABELS
            ]
        )
        status, arsak, flag = kantveg.classify_sides(
            has_veg_by_side[side], potensiell_by_side[side], pcts
        )
        # A bank clipped away at the national border is blanked rather than
        # classified: its Vtot and PotVeg are 0 because nothing Norwegian
        # maps that ground, and classify_sides would read that as a natural
        # absence of riparian forest. The combiners below then reduce to the
        # Norwegian bank on their own.
        cols[f"E_{side}_{b}_Kantvegetasjon_A"] = sides.mask(out, side, status)
        cols[f"E_{side}_{b}_Kantvegetasjon_A_arsak"] = sides.mask(
            out, side, arsak
        )
        cols[f"E_{side}_{b}_Menneskelig_inngrep"] = sides.mask(out, side, flag)

    left = cols[f"E_L_{b}_Kantvegetasjon_A"]
    right = cols[f"E_R_{b}_Kantvegetasjon_A"]
    cols[f"E_S_{b}_Kantvegetasjon_A"] = kantveg.combine_status(left, right)
    cols[f"E_S_{b}_Kantvegetasjon_A_arsak"] = kantveg.combine_causes(
        cols[f"E_L_{b}_Kantvegetasjon_A_arsak"],
        cols[f"E_R_{b}_Kantvegetasjon_A_arsak"],
    )
    cols[f"E_S_{b}_Menneskelig_inngrep"] = kantveg.combine_flags(
        cols[f"E_L_{b}_Menneskelig_inngrep"],
        cols[f"E_R_{b}_Menneskelig_inngrep"],
        cols[f"E_S_{b}_Kantvegetasjon_A"],
    )
    return cols


def compute(df: pd.DataFrame) -> pd.DataFrame:
    """Add Potensiell_Kantvegetasjon, Kantvegetasjon_A and its arsak/flag,
    per side and combined, at B0-10, at B0-20, and band-free (selected)."""
    out = df.copy()
    for b, scale in _BANDS:
        for name, val in _compute_band(out, b, scale).items():
            out[name] = val

    wide = (
        num(out, "K_S_beregningsbredde", C.PRIMARY_BAND[1]) == C.AREAL_BAND[1]
    )
    for side in ("L", "R"):
        out[f"E_{side}_Potensiell_Kantvegetasjon"] = np.where(
            wide,
            out[f"E_{side}_B0-20_Potensiell_Kantvegetasjon"],
            out[f"E_{side}_B0-10_Potensiell_Kantvegetasjon"],
        )
    for side in ("L", "R", "S"):
        out[f"E_{side}_Kantvegetasjon_A"] = np.where(
            wide,
            out[f"E_{side}_B0-20_Kantvegetasjon_A"],
            out[f"E_{side}_B0-10_Kantvegetasjon_A"],
        )
        out[f"E_{side}_Kantvegetasjon_A_arsak"] = np.where(
            wide,
            out[f"E_{side}_B0-20_Kantvegetasjon_A_arsak"],
            out[f"E_{side}_B0-10_Kantvegetasjon_A_arsak"],
        )
        out[f"E_{side}_Menneskelig_inngrep"] = np.where(
            wide,
            out[f"E_{side}_B0-20_Menneskelig_inngrep"],
            out[f"E_{side}_B0-10_Menneskelig_inngrep"],
        )
    return out
