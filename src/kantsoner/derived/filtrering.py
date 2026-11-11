"""Filtration efficiency — sediment, phosphorus, nitrogen, pesticides.

System B, Filtrering block. Pure tabular.

Note on the source workbook: the ``v_faktor_NP`` row of the Attributes tab is
largely a copy-paste of ``v_faktor_SS`` (its formula assigns to
``v_faktor_SS`` throughout), and the ``H_SS``/``H_PN`` rows carry a stray
``IFELSE`` and an unbalanced bracket. The structure implied by both — a
vegetation-type modifier with the same branches but different values, and a
piecewise slope modifier hinged at 10 % — is what is implemented here, with
the coefficients read off the intended column. Worth confirming against
Bendik's methodology PDF before delivery.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..columns import arealbruk as arealbruk_of
from ..columns import num, text
from ..spec import constants as C


def _cover_terms(df: pd.DataFrame, side: str):
    """(has_cover, tall_dominant, treslag as stored) of one bank.

    A bank has cover when crowns stand on any of the 0-10 m zone. The earlier
    gate was Vtot >= 50 %; it gave "none" (v_faktor 0) to every bank with less
    cover, although the width now counts the area of crowns at any share (see
    width.py). Where the 0-10 m zone has no crowns but the 0-20 m band has,
    the bank is a strip set back from the water, and the height classes of
    the 0-20 m band say whether it is tall or low. The species is the 0-10 m
    one in both cases (the only zone it is measured on).
    """
    b = f"K_{side}_B0-10"
    wide = f"K_{side}_B0-20"
    near_cover = num(df, f"{b}_Vtot", 0.0) > 0
    vtot = num(df, f"{b}_Vtot", 0.0).where(near_cover, num(df, f"{wide}_Vtot", 0.0))
    tall = num(df, f"{b}_V5-10", 0.0) + num(df, f"{b}_V10pluss", 0.0)
    low = num(df, f"{b}_V1-5", 0.0) + num(df, f"{b}_V0p5-1", 0.0)
    tall_wide = num(df, f"{wide}_V5-10", 0.0) + num(df, f"{wide}_V10pluss", 0.0)
    low_wide = num(df, f"{wide}_V1-5", 0.0) + num(df, f"{wide}_V0p5-1", 0.0)
    tall = tall.where(near_cover, tall_wide)
    low = low.where(near_cover, low_wide)
    return vtot > 0, tall >= low, text(df, f"{b}_Treslag")


def species_defaulted(df: pd.DataFrame, side: str) -> pd.Series:
    """True where a tall canopy has no species and løv was assumed."""
    has_cover, tall_dominant, treslag = _cover_terms(df, side)
    return has_cover & tall_dominant & treslag.isin(UNKNOWN_TEXT)


UNKNOWN_TEXT = ("", "<NA>", "nan", "None")


def _vegetation_situation(df: pd.DataFrame, side: str) -> pd.Series:
    """Classify each segment side into the branches the v_faktor tables use.

    Returns one of: tall_lauv, tall_blanding, tall_other, low, <arealbruk>, none.
    A tall canopy with no species is løv (report Tabell C.1).
    """
    has_cover, tall_dominant, treslag = _cover_terms(df, side)
    treslag = treslag.mask(treslag.isin(UNKNOWN_TEXT), C.TRESLAG_UNKNOWN)
    arealbruk = arealbruk_of(df, side)

    out = pd.Series("none", index=df.index, dtype=object)
    # Fall back to arealbruk only where there is no canopy at all.
    out = out.mask(arealbruk.isin(("Grassmark", "Hei og buskmark")), arealbruk)
    out = out.mask(has_cover & ~tall_dominant, "low")
    out = out.mask(has_cover & tall_dominant & (treslag == "Lauv"), "tall_lauv")
    out = out.mask(
        has_cover & tall_dominant & (treslag == "Blandingsskog"),
        "tall_blanding",
    )
    out = out.mask(
        has_cover & tall_dominant & ~treslag.isin(("Lauv", "Blandingsskog")),
        "tall_other",
    )
    return out


def _h_ss(slope: pd.Series) -> pd.Series:
    """Slope modulator for sediment, hinged at 10 % (Sállir A.2 and tabell 9).

    Floored at 0 on both branches. The report lists H_SS as 0-1, but its low
    branch has no limit of its own, and s here is signed: it is negative where
    the ground falls away from the water, and the branch passes 0 at -28.5 %.
    """
    low = (C.H_SS_LOW_INTERCEPT + C.H_SS_LOW_SLOPE * slope) / 100.0
    high = (C.H_SS_HIGH_INTERCEPT + C.H_SS_HIGH_SLOPE * slope) / 100.0
    return pd.Series(
        np.maximum(0.0, np.where(slope <= 10.0, low, high)), index=slope.index
    )


def _h_pn(slope: pd.Series) -> pd.Series:
    """Slope modulator for N/P: flat below 10 %, decaying above."""
    decayed = 1.0 - C.H_PN_DECAY * ((slope - C.H_PN_BREAK_PCT) / 10.0)
    return pd.Series(
        np.where(slope <= C.H_PN_BREAK_PCT, 1.0, np.maximum(0.0, decayed)),
        index=slope.index,
    )


def compute(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    for side in ("L", "R"):
        b = f"K_{side}_B0-10"
        situation = _vegetation_situation(out, side)
        arealbruk = arealbruk_of(out, side)

        v_ss = situation.map(C.V_FAKTOR_SS).astype(float)
        # tall_blanding is a næring-only branch; it filters like tall_other.
        v_ss = v_ss.fillna(
            situation.map({"tall_blanding": C.V_FAKTOR_SS["tall_other"]})
        )
        v_np = situation.map(C.V_FAKTOR_NP).astype(float)
        v_np = v_np.fillna(
            situation.map({"tall_blanding": C.V_FAKTOR_NP["tall_other"]})
        )
        out[f"{b}_v_faktor_SS"] = v_ss.fillna(0.0).round(4)
        out[f"{b}_v_faktor_NP"] = v_np.fillna(0.0).round(4)

        n_sub = arealbruk.map(C.N_SUB).astype(float).fillna(C.N_SUB_DEFAULT)
        out[f"{b}_N_sub"] = n_sub.round(4)

        slope = num(out, f"{b}_s", 0.0)
        out[f"{b}_H_SS"] = _h_ss(slope).round(4)
        out[f"{b}_H_PN"] = _h_pn(slope).round(4)

        w_eff = num(out, f"{b}_w_eff", 0.0)

        y_ss = (
            C.Y_SS_ASYMPTOTE
            * v_ss.fillna(0.0)
            * (1 - np.exp(-C.Y_SS_RATE * w_eff))
        )
        y_x = C.Y_X_ASYMPTOTE * (1 - np.exp(-C.Y_X_RATE * w_eff))
        y_n = (
            C.Y_N_ASYMPTOTE
            * v_np.fillna(0.0)
            * (1 - np.exp(-C.Y_N_RATE * w_eff))
        )
        y_p = (
            C.Y_P_ASYMPTOTE
            * v_np.fillna(0.0)
            * (1 - np.exp(-C.Y_P_RATE * w_eff))
        )

        out[f"{b}_Y_SS"] = y_ss.round(3)
        out[f"{b}_Y_X"] = y_x.round(3)
        out[f"{b}_Y_N"] = y_n.round(3)
        out[f"{b}_Y_P"] = y_p.round(3)
        out[f"{b}_Y_N_total"] = (
            C.Y_N_SURFACE_WEIGHT * y_n + C.Y_N_SUBSURFACE_WEIGHT * n_sub * 100.0
        ).round(3)

        out[f"{b}_Y_SS_final"] = (y_ss * out[f"{b}_H_SS"]).round(3)
        out[f"{b}_Y_P_final"] = (y_p * out[f"{b}_H_PN"]).round(3)
        out[f"{b}_Y_N_final"] = (y_n * out[f"{b}_H_PN"]).round(3)
        # Nitrogen leaves by two ways (report A.2): over the surface, which
        # the width and the slope govern, and through the ground, which
        # N_sub of the land use stands for. This is Y_N_total with the slope
        # modulator on the surface term, and it is what Filter_score reads.
        out[f"{b}_Y_N_total_final"] = (
            C.Y_N_SURFACE_WEIGHT * out[f"{b}_Y_N_final"]
            + C.Y_N_SUBSURFACE_WEIGHT * n_sub * 100.0
        ).round(3)

        # Grøfter short-circuit the filter strip; absent a grøftekart the
        # density is 0 and the bypass factor is a no-op 1.0.
        tetthet = num(out, f"{b}_Grofte_tetthet", 0.0)
        bypass = np.maximum(
            0.0, 1.0 - np.log1p(tetthet) / np.log(C.GROFT_LN_BASE)
        )
        out[f"{b}_Grofte_bypass"] = pd.Series(bypass, index=out.index).round(4)

    return out
