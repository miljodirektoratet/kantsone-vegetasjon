"""Potential-vegetation cover of the riparian zones.

`SkogJordbruk` is the share of a segment's zone slice where vegetation could
grow at all — forest, agriculture, or unspecified open ground on soil. It is
the reference state system A scores against: a bank with no trees is only
"missing" riparian vegetation if vegetation was possible there.
"""

import geopandas as gpd
import numpy as np
import pandas as pd

from ..geometry import segment_zones
from ..io import as_gdf
from ..spec import constants as C
from ..spec.attributes import band_token


def extract(segments, zones, band, potveg):
    """SkogJordbruk columns for one band, both sides, keyed on segment_id."""
    ids = as_gdf(segments)["segment_id"]
    pv = as_gdf(potveg).to_crs(zones[0].crs)

    out = {}
    for side, cells in (("L", zones[0]), ("R", zones[1])):
        prefix = f"K_{side}_{band_token(band)}"
        sl = segment_zones.cut(cells, zones[2])
        atot = ids.map(
            sl.assign(_a=sl.geometry.area).groupby("segment_id")["_a"].sum()
        ).fillna(0.0)

        ov = gpd.overlay(
            sl[["segment_id", "geometry"]],
            pv[["geometry"]],
            how="intersection",
            keep_geom_type=True,
        )
        area = (
            ov.assign(_a=ov.geometry.area).groupby("segment_id")["_a"].sum()
            if not ov.empty
            else pd.Series(dtype=float)
        )
        covered = ids.map(area).fillna(0.0)
        out[f"{prefix}_SkogJordbruk_m2"] = covered.round(2)
        out[f"{prefix}_SkogJordbruk"] = (
            (covered / atot.replace(0, np.nan) * 100)
            .clip(upper=100)
            .fillna(0.0)
            .round(2)
        )
    return out


def _slice_overlay(segments, cells, ring, grunnkart, cols):
    """Zone slice x grunnkart, carrying `cols` and the intersection area."""
    sl = segment_zones.cut(cells, ring)
    if sl.empty:
        return sl.assign(_a=[]), sl
    gk = as_gdf(grunnkart).to_crs(sl.crs)
    keep = [c for c in cols if c in gk.columns]
    ov = gpd.overlay(
        sl[["segment_id", "geometry"]],
        gk[keep + ["geometry"]],
        how="intersection",
        keep_geom_type=True,
    )
    return (ov.assign(_a=ov.geometry.area) if not ov.empty else ov), sl


def classes(segments, zones, band, ring, grunnkart):
    """Arealbruk and Treslag: the class covering most of the zone slice.

    The workbook says "max of okosystemtypeNiva1", which is area-dominance
    rather than a count -- a segment whose slice is mostly forest reads Skog
    even where a dozen small settlement polygons touch it. Grunnkart names the
    ecosystem type in English; the workbook's tables are Norwegian, so the
    class is translated on the way out.
    """
    ids = as_gdf(segments)["segment_id"]
    out = {}
    for side, cells in (("L", zones[0]), ("R", zones[1])):
        prefix = f"K_{side}_{band_token(band)}"
        ov, _ = _slice_overlay(
            segments, cells, ring, grunnkart, ["okosystemtypeniva1", "arealdekkeniva2"]
        )
        if not ov.empty and "arealdekkeniva2" in ov.columns:
            # NaN on ground with no species, which the area sum below skips, so
            # a bank takes the dominant species wherever it has any.
            ov = ov.assign(
                treslag=ov["arealdekkeniva2"].map(C.TRESLAG_FROM_AREALDEKKE)
            )
        for col, name, translate in (
            ("okosystemtypeniva1", "Arealbruk", True),
            ("treslag", "Treslag", False),
        ):
            val = pd.Series(pd.NA, index=ids.index, dtype="object")
            if not ov.empty and col in ov.columns:
                area = ov.groupby(["segment_id", col])["_a"].sum()
                top = (
                    area.reset_index()
                    .sort_values("_a")
                    .groupby("segment_id")
                    .last()
                )
                m = ids.map(top[col])
                if translate:
                    m = m.map(lambda v: C.AREALBRUK_FROM_GRUNNKART.get(v, v))
                val = m
            out[f"{prefix}_{name}"] = val
    return out


def areas(segments, zones, band, ring, grunnkart):
    """Areal and the semi-natural splits, in m2 of zone slice.

    Areal is the slice's own area, not the sum of the classes: grunnkart does
    not tile the zone completely, and using the covered area as denominator
    would score a half-mapped bank as fully natural.
    """
    ids = as_gdf(segments)["segment_id"]
    out = {}
    for side, cells in (("L", zones[0]), ("R", zones[1])):
        prefix = f"K_{side}_{band_token(band)}"
        ov, sl = _slice_overlay(
            segments, cells, ring, grunnkart, ["okosystemtypeniva1"]
        )
        total = (
            ids.map(
                sl.assign(_t=sl.geometry.area).groupby("segment_id")["_t"].sum()
            ).fillna(0.0)
            if not sl.empty
            else pd.Series(0.0, index=ids.index)
        )
        out[f"{prefix}_Areal"] = total.round(2)
        for key, members in C.AREAL_CLASS_FROM_GRUNNKART.items():
            a = pd.Series(0.0, index=ids.index)
            if not ov.empty and "okosystemtypeniva1" in ov.columns:
                hit = ov[ov["okosystemtypeniva1"].isin(members)]
                if not hit.empty:
                    a = ids.map(hit.groupby("segment_id")["_a"].sum()).fillna(
                        0.0
                    )
            out[f"{prefix}_Areal_{key}"] = a.round(2)
    return out


# ── cause and potential-vegetation coverage ────────────────────────────────
# Category masks over the grunnkart rows, all compared lower-cased against
# arealdekkeniva2 (arealdekkeniva1 only for the skog check -- niva2 splits
# Skog into 5 species, so there is no bare "skog" value at that level).
# Settled on the minor_fixes branch, against the old flat arealdekke schema;
# the subtleties carry over unchanged onto niva1/niva2:
#   * Berg excludes Pot_veg ground, or forest on shallow soil (Skog on
#     Grunnlendt / Fjell_i_dagen / Blokkmark) would count as both PotVeg and
#     Berg. It also excludes the human classes so a cause is never both.
#   * The soil-covered Snaumark classes are Pot_veg (vegetatable); the same
#     classes on other ground are neither Pot_veg nor a natural cause.
_BERG_SUBSTRATE = {"fjellidagen", "blokkmark", "grunnlendt"}
_BERG_PROXY_AD = {"snaumarkskrinn", "snaumarkimpediment"}
_HUMAN_AD = {"samferdsel", "snaumarkkonstruert", "bebygd"}
_MYR = {
    "myrapen",
    "myrtresattbarblanding",
    "myrtresattblanding",
    "myrtresattfuru",
    "myrtresattgran",
    "myrtresattlauv",
}
_VANN = {
    "ferskvanninnsjotjern",
    "ferskvanninnsjotjernregulert",
    "ferskvannelvbekk",
    "ferskvannkanal",
    "hav",
}
_SNAUMARK_POTVEG = {
    "snaumarkuspesifisert",
    "snaumarkfrisk",
    "snaumarkmiddelsfrisk",
}
# PotVeg (the delivered vegetation-capable-ground reference state): grunnkart
# okosystemtypeniva1 in these classes, minus skog on impediment-bonitet
# ground -- an impediment stand cannot grow productive forest, so grunnkart
# calling it "skog" should not read as vegetation-capable ground. Supersedes
# the arealdekke-based reference state below (kept only for Berg's own
# exclusion, which predates this and is not being redefined here).
_POTVEG_OKOSYSTEM = {"skog", "dyrketmark", "grasmark"}
_POTVEG_IMPEDIMENT_BONITET = "impediment"


def _cause_masks(gk):
    """{category: boolean mask over grunnkart rows}, mutually consistent.

    Every mask already delivered under the spec workbook's own names (Berg,
    Myr, Vann, Skog, Konstruert, Jordbruk, Snoisbre) keeps its exact existing
    definition -- these are the specification's own terms and are not
    redefined here. Skogbruk/SnaumarkFastfjell/Vatmark are additional,
    arsak-only candidates: Skogbruk aliases Skog (grunnkart still says
    forest, whatever the measured canopy says -- see system_a.py),
    SnaumarkFastfjell is Berg widened to every natural (non-Pot_veg,
    non-Konstruert) Snaumark reading, and Vatmark is Myr | Vann -- the
    workbook attributes they are built from are unaffected.
    """
    ad1 = gk["arealdekkeniva1"].astype("string").str.strip().str.lower()
    ad = gk["arealdekkeniva2"].astype("string").str.strip().str.lower()
    gf = gk["grunnforhold"].astype("string").str.strip().str.lower()
    jordbruk = ad.str.startswith("jordbruk").fillna(False)
    skog = (ad1 == "skog").fillna(False)
    snaumark_jorddekt = (ad.isin(_SNAUMARK_POTVEG) & (gf == "jorddekt")).fillna(
        False
    )
    # arealdekke-based reference state -- kept only for Berg's own exclusion
    # below, not the delivered PotVeg (see _POTVEG_OKOSYSTEM above).
    potveg_arealdekke = skog | jordbruk | snaumark_jorddekt
    human = ad.isin(_HUMAN_AD).fillna(False) | jordbruk
    konstruert = (ad == "snaumarkkonstruert").fillna(False)
    myr = ad.isin(_MYR).fillna(False)
    vann = ad.isin(_VANN).fillna(False)
    berg = (
        (
            gf.isin(_BERG_SUBSTRATE)
            # "ikkeregistrert" is the new schema's rename of the old
            # "ikke_kartlagt" ("not mapped") -- confirmed present under the
            # new spelling, but not confirmed as the same *meaning* ("not
            # registered" could differ); re-check against the source's own
            # documentation if Berg coverage looks off after migration.
            | (ad.isin(_BERG_PROXY_AD) & (gf == "ikkeregistrert"))
        )
        & ~human
        & ~potveg_arealdekke
    ).fillna(False)
    # Every other Snaumark reading: not soil-covered (so not Pot_veg), not
    # Konstruert (so not Nedbygd_kantsone) -- the "rest of snaumark" the
    # Berg substrate test alone does not catch (e.g. Skrinn/Impediment on
    # ground grunnforhold does classify, just not as one of the Berg
    # proxies).
    snaumark_fastfjell = (
        berg
        | (ad.str.startswith("snaumark") & ~snaumark_jorddekt & ~konstruert)
    ).fillna(False)
    ok1 = gk["okosystemtypeniva1"].astype("string").str.strip()
    bonitet = gk["skogbonitet"].astype("string").str.strip()
    impediment_skog = (
        (ok1 == "skog") & (bonitet == _POTVEG_IMPEDIMENT_BONITET)
    ).fillna(False)
    potveg = (ok1.isin(_POTVEG_OKOSYSTEM) & ~impediment_skog).fillna(False)
    return {
        "Konstruert": konstruert,
        "Jordbruk": jordbruk,
        "Berg": berg,
        "Myr": myr,
        "Snoisbre": (ad == "snoisbre").fillna(False),
        "Vann": vann,
        "Skog": skog,
        "SnaumarkJorddekt": snaumark_jorddekt,
        "PotVeg": potveg,
        # arsak-only, not their own workbook attribute -- see docstring.
        "Skogbruk": skog,
        "SnaumarkFastfjell": snaumark_fastfjell,
        "Vatmark": myr | vann,
    }


def causes(segments, zones, band, ring, grunnkart):
    """Cause and Pot_veg coverage per side: why vegetation is absent, or
    could grow. One overlay per side; every category is a mask over its rows.

    Columns per side: Konstruert, Jordbruk, Berg, Myr, Snoisbre, Vann (the
    arsak candidates beyond Bygg/Samferdsel, which vegetation.extract already
    measures), plus the Pot_veg breakdown Skog / SnaumarkJorddekt / PotVeg,
    plus the arsak-only Skogbruk / SnaumarkFastfjell / Vatmark (see
    _cause_masks) — each as a percentage of the zone and as m2.
    """
    ids = as_gdf(segments)["segment_id"]
    out = {}
    for side, cells in (("L", zones[0]), ("R", zones[1])):
        prefix = f"K_{side}_{band_token(band)}"
        ov, sl = _slice_overlay(
            segments,
            cells,
            ring,
            grunnkart,
            [
                "arealdekkeniva1",
                "arealdekkeniva2",
                "grunnforhold",
                "okosystemtypeniva1",
                "skogbonitet",
            ],
        )
        atot = ids.map(
            sl.assign(_a=sl.geometry.area).groupby("segment_id")["_a"].sum()
        ).fillna(0.0)
        masks = _cause_masks(ov) if len(ov) else {}
        for cat in (
            "Konstruert",
            "Jordbruk",
            "Berg",
            "Myr",
            "Snoisbre",
            "Vann",
            "Skog",
            "SnaumarkJorddekt",
            "PotVeg",
            "Skogbruk",
            "SnaumarkFastfjell",
            "Vatmark",
        ):
            if len(ov):
                area = ov[masks[cat]].groupby("segment_id")["_a"].sum()
                covered = ids.map(area).fillna(0.0)
            else:
                covered = pd.Series(0.0, index=ids.index)
            out[f"{prefix}_{cat}_m2"] = covered.round(2)
            out[f"{prefix}_{cat}"] = (
                (covered / atot.replace(0, np.nan) * 100)
                .clip(upper=100)
                .fillna(0.0)
                .round(2)
            )
    return out


# ── dynamic has_veg threshold ───────────────────────────────────────────────
# system_a.py's has_veg check (Vtot >= threshold) uses a per-segment-side
# threshold based on the zone slice's dominant okosystemtypeniva1 class --
# and, where that class is skog, its dominant skogbonitet:
#   skog, bonitet middels/høy   -> 50% (a productive site can support the
#                                   same canopy fraction a lowland stand
#                                   could)
#   skog, any other bonitet     -> 30% (a low-productivity/impediment/
#                                   særs høy-but-otherwise-uncommon stand
#                                   cannot be expected to close canopy the
#                                   same way)
#   heiBuskmark/liteVegetertMark/vatmark -> 30% (naturally sparser cover)
#   everything else              -> 50% (the flat default)
# See derived/system_a.py for how the resulting K_{side}_B0-10_VtotThreshold
# column gets read.
_DYNTHRESH_HIGH_BONITET = {"middelsBonitet", "hoyBonitet"}
_DYNTHRESH_LOW_CLASSES = {"heiBuskmark", "liteVegetertMark", "vatmark"}


def _dominant(ov, ids, col):
    """Per segment_id, the value of `col` covering the most area in `ov`."""
    if ov.empty or col not in ov.columns:
        return pd.Series(pd.NA, index=ids.index, dtype="object")
    top = (
        ov.groupby(["segment_id", col])["_a"]
        .sum()
        .reset_index()
        .sort_values("_a")
        .groupby("segment_id")
        .last()
    )
    return ids.map(top[col])


def veg_threshold(segments, zones, band, ring, grunnkart):
    """K_{side}_{band}_VtotThreshold: the has_veg (Vtot) cutoff for this
    segment/side, per the rule above."""
    ids = as_gdf(segments)["segment_id"]
    out = {}
    for side, cells in (("L", zones[0]), ("R", zones[1])):
        prefix = f"K_{side}_{band_token(band)}"
        ov, _ = _slice_overlay(
            segments,
            cells,
            ring,
            grunnkart,
            ["okosystemtypeniva1", "skogbonitet"],
        )
        thresh = pd.Series(
            C.KANTVEG_THRESHOLD_PCT, index=ids.index, dtype="float64"
        )
        dom_class = _dominant(ov, ids, "okosystemtypeniva1")
        is_skog = dom_class == "skog"
        if is_skog.any() and not ov.empty:
            skog_rows = ov[ov["okosystemtypeniva1"] == "skog"]
            dom_bonitet = _dominant(skog_rows, ids, "skogbonitet")
            high_bonitet = dom_bonitet.isin(_DYNTHRESH_HIGH_BONITET)
            thresh[is_skog & ~high_bonitet] = 30.0
        low_class = dom_class.isin(_DYNTHRESH_LOW_CLASSES)
        thresh[~is_skog & low_class] = 30.0
        out[f"{prefix}_VtotThreshold"] = thresh
    return out


# ── expected canopy cover (graded Kantvegetasjon) ──────────────────────────
_BONITET_KLASSE = {
    "sersHoyBonitet": "svhoy",
    "hoyBonitet": "hoy",
    "middelsBonitet": "middels",
    "lavBonitet": "lav",
    "impediment": "imp",
}


def _forventet_klasse(okosystem, bonitet, arealbruk):
    """Grunnkart polygon -> key of C.FORVENTET_DEKNING_PCT, or None."""
    if okosystem == "skog":
        # skog with bonitet ikkeRelevant (a few %) is read as middels
        return _BONITET_KLASSE.get(bonitet, "middels")
    if okosystem == "dyrketmark":
        return "dyrket"
    if okosystem == "bebygdOpparbeidetAreal":
        return "veg" if arealbruk == "TransportTelek" else "bygg"
    if okosystem == "grasmark":
        return "gras"
    return None


def forventet_dekning(segments, zones, band, ring, grunnkart):
    """K_{side}_{band}_Forventet_dekning and _Vegdyktig, per segment/side.

    Forventet_dekning: the canopy cover the bank would have in its reference
    state -- the area-weighted mean of C.FORVENTET_DEKNING_PCT over the zone
    slice, ground outside those classes counting 0. Vegdyktig: the share of
    the slice in C.VEGDYKTIG_KLASSER. Both are percentages of the same slice
    area Vtot is measured over, so the two compare directly.
    """
    ids = as_gdf(segments)["segment_id"]
    out = {}
    for side, cells in (("L", zones[0]), ("R", zones[1])):
        prefix = f"K_{side}_{band_token(band)}"
        ov, sl = _slice_overlay(
            segments,
            cells,
            ring,
            grunnkart,
            ["okosystemtypeniva1", "skogbonitet", "arealbruklandhovedklasse"],
        )
        forv = pd.Series(np.nan, index=ids.index, dtype="float64")
        vegd = pd.Series(np.nan, index=ids.index, dtype="float64")
        if not sl.empty:
            atot = (
                sl.assign(_a=sl.geometry.area).groupby("segment_id")["_a"].sum()
            )
            if not ov.empty:
                k = [
                    _forventet_klasse(o, b, a)
                    for o, b, a in zip(
                        ov.get("okosystemtypeniva1", pd.Series(index=ov.index)),
                        ov.get("skogbonitet", pd.Series(index=ov.index)),
                        ov.get(
                            "arealbruklandhovedklasse",
                            pd.Series(index=ov.index),
                        ),
                        strict=True,
                    )
                ]
                ov = ov.assign(_k=k)
                ov["_t"] = (
                    ov["_k"].map(C.FORVENTET_DEKNING_PCT).fillna(0.0) * ov["_a"]
                )
                ov["_v"] = ov["_a"].where(
                    ov["_k"].isin(C.VEGDYKTIG_KLASSER), 0.0
                )
                agg = ov.groupby("segment_id")[["_t", "_v"]].sum()
            else:
                agg = pd.DataFrame(columns=["_t", "_v"])
            a = ids.map(atot)
            forv = (ids.map(agg["_t"]).fillna(0.0) / a.replace(0, np.nan)).clip(
                0, 100
            )
            vegd = (
                100 * ids.map(agg["_v"]).fillna(0.0) / a.replace(0, np.nan)
            ).clip(0, 100)
        out[f"{prefix}_Forventet_dekning"] = forv.round(1)
        out[f"{prefix}_Vegdyktig"] = vegd.round(1)
    return out
