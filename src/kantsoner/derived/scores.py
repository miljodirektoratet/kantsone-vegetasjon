"""System-B composite scores, per 10 m segment.

Filter_score        = weighted mean of Y_SS_final, Y_P_final, Y_N_total_final and
                      Y_X / 100 (N carries N_sub, the report's Y_N_total)
Filter_score_naering= Filter_score * (0.8 + 0.2 * mean Naering_effekt)
Filter_score_final  = Filter_score_naering * mean Grofte_bypass
Struktur_score      = 0.5 * Temp_score + 0.5 * mean Naering_effekt
Kontinuitet         = 0.7 * mean Kontin_tre + 0.3 * mean Kontin_busk
Score_total         = 0.5 * Filter + 0.3 * Struktur + 0.2 * Kontinuitet
"""

from __future__ import annotations

import pandas as pd

from .. import sides
from ..columns import is_lake, num
from ..spec import constants as C


def _mean_sides(df: pd.DataFrame, suffix: str) -> pd.Series:
    """Mean over the banks that have data.

    The single place every system-B score folds two banks into one, so
    excluding a bank abroad here, or a bank with no zone (not measured, and
    not 0), covers Filter_score, Naering, the grøfte
    bypass and Kontinuitet at once. Without it a border river would be
    scored as having half a riparian zone -- the Norwegian bank's value
    averaged against a zero that means "no data", not "nothing there".
    """
    cols = [
        sides.mask(
            df, side, num(df, f"K_{side}_B0-10_{suffix}", None), measured=True
        )
        for side in sides.SIDES
    ]
    return pd.concat(cols, axis=1).mean(axis=1, skipna=True).fillna(0.0)


def compute(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    w = C.FILTER_WEIGHTS

    # The workbook's Filter_score row mixes L and R columns (SS and N from the
    # left, P from the right); that is a transcription slip. Each efficiency is
    # averaged across both banks here, consistently with every other score.
    filter_score = (
        _mean_sides(out, "Y_SS_final") * w["SS"]
        + _mean_sides(out, "Y_P_final") * w["P"]
        + _mean_sides(out, "Y_N_total_final") * w["N"]
        + _mean_sides(out, "Y_X") * w["X"]
    ) / 100.0
    out["E_S_Filter_score"] = filter_score.round(4)

    naering = _mean_sides(out, "Naering_effekt")
    with_naering = filter_score * (
        C.FILTER_NAERING_BASE + C.FILTER_NAERING_GAIN * naering
    )
    out["E_S_Filter_score_naering"] = with_naering.round(4)

    bypass = _mean_sides(out, "Grofte_bypass")
    # Absent a grøftekart every bypass is 1.0 and this collapses to the
    # previous column, as the spec notes.
    lake = is_lake(out)
    # Naering_effekt scales litter input by channel width, so on a lake it
    # is not small, it is undefined. The naering step is skipped there and
    # the filter score carries through unmodified.
    base = filter_score.where(lake, with_naering)
    out["E_S_Filter_score_final"] = (base * bypass).round(4)

    temp = num(out, "K_S_B0-10_Temp_score", 0.0)
    out["E_S_Struktur_score"] = (
        C.STRUKTUR_TEMP_WEIGHT * temp + C.STRUKTUR_NAERING_WEIGHT * naering
    ).round(4)

    out["E_S_Kontinuitet"] = (
        C.KONTIN_TRE_WEIGHT * _mean_sides(out, "Kontin_tre")
        + C.KONTIN_BUSK_WEIGHT * _mean_sides(out, "Kontin_busk")
    ).round(4)

    sw = C.SCORE_TOTAL_WEIGHTS
    river_total = (
        out["E_S_Filter_score_final"] * sw["filter"]
        + out["E_S_Struktur_score"] * sw["struktur"]
        + out["E_S_Kontinuitet"] * sw["kontinuitet"]
    )
    # Struktur is temperature and naering, both width-derived, so on a lake
    # it is not zero but absent -- and a zero here reads as "worst
    # possible", which is the opposite of what it means. The leg is dropped
    # and the other two restated (see C.LAKE_SCORE_WEIGHTS); Struktur_score
    # itself is blanked so nothing downstream averages a placeholder.
    lw = C.LAKE_SCORE_WEIGHTS
    lake_total = (
        out["E_S_Filter_score_final"] * lw["filter"]
        + out["E_S_Kontinuitet"] * lw["kontinuitet"]
    )
    out["E_S_Score_total"] = river_total.where(~lake, lake_total).round(4)
    out["E_S_Struktur_score"] = out["E_S_Struktur_score"].where(~lake)
    return out
