"""Roll the official 10 m centerline up to every delivery level.

One parameterised implementation replaces the four near-identical scripts
(100 m / 500 m / 1000 m / river+region) that differed only in block size and
column suffix.

The chain nests by construction — each level is built by summing the previous
level's blocks, never by re-reading the 10 m base::

    10 m -> 100 m (x10) -> 500 m (x5) -> 1000 m (x2)

``river`` and ``region`` collapse straight from the 10 m base.

Only **additive counts** are stored and summed; every ratio is recomputed at
each level from those sums. Averaging ratios of ratios would silently weight
short blocks equally with long ones.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from loguru import logger
from shapely.ops import linemerge

from .. import sides
from ..derived import kantveg
from ..io import save_gpkg
from ..spec import attributes as A
from ..spec import constants as C
from ..spec import scales

RIVER_ID = "river_id_vann-nett"
SEGMENT_ID = "segment_id_vann-nett"
LENGTH = "segment_length_vann-nett"
NAME = "vannforekomstnavn"
# Matched is pass-1 only (a channel could be reconstructed at all) and is
# not what these ratios should exclude on -- Beregnet is pass 2, whether
# this segment actually ended up with a real zone on both sides. See
# stages.py:aggregation for both.
BEREGNET = "Beregnet"
NETWORK_TYPE = "network_type"
CLASS_COL = "E_S_Kantvegetasjon_A"

CLASS_COUNTS = {
    "n_ok": C.OK_KANTVEG,
    "n_kantskog_fjernet": C.KANTSKOG_FJERNET,
    "n_nedbygd_kantsone": C.NEDBYGD_KANTSONE,
    "n_ikke": C.IKKE_MULIG_KANTVEG,
    "n_uklassifisert": C.UKLASSIFISERT_KANTVEG,
}

# Per-segment columns summed on the way up, taken from the registry: adding an
# attribute to the delivery means setting deliver=True there and nothing here.
# All three are summed; they differ in how the sum is reported per level.
SHARE_COLUMNS = tuple(A.delivered("share_ja"))  # -> share of matched length
COUNT_COLUMNS = tuple(A.delivered("count_ja"))  # -> raw count ("Antall")
MEAN_COLUMNS = tuple(A.delivered("mean"))  # -> mean over matched segments
# -> the worst class present in the block. Carried as one count per class
# so it survives the block sums like everything else: the worst is then
# the lowest-graded class with a non-zero count.
WORST_COLUMNS = tuple(A.delivered("worst"))
# 1 is the worst end of the grad-av-endring scale, 5 the best.
WORST_CLASSES = ("Svært dårlig", "Dårlig", "Moderat", "God", "Svært god")

# The Kantvegetasjon cause columns are re-derived per block from summed
# areas rather than averaged: a percentage is sum(area)/sum(Atot), so the
# blocks classify by exactly the same rules as the 10 m segments
# (derived/kantveg.py holds them). These are the per-side categories summed
# up the levels for that; Bygg and Samferdsel ship as percentages only, so
# their areas are reconstructed as pct x Atot / 100. Each category's *source*
# column is picked per segment from B0-10 or B0-20 by K_S_beregningsbredde
# (primary/channel.py) before summing -- a block spanning both narrow and
# wide stretches must sum each segment's own official-width area, not
# whichever band happens to carry the "B0-10" name.
CAUSE_AREA_CATEGORIES = (*kantveg.ARSAK_LABELS[2:], "PotVeg", "Vtot", "Atot")
CAUSE_PCT_CATEGORIES = ("Bygg", "Samferdsel")


@dataclass(frozen=True)
class Level:
    name: str
    parent: (
        str | None
    )  # level rolled up from; None = straight off the 10 m base
    factor: int | None  # how many parent blocks per block
    group_by: tuple[str, ...]  # identity columns defining a block


# Block levels group by chain, not just river: a chain is a run of stretches
# that flow uniquely one into the next (geometry.graph.chains), so a block
# is always consecutive water. Grouping by river alone and sorting the ids
# spliced sibling stretches together.
CHAIN_ID = "chain_id"
CHAIN_POS = "chain_pos"
LEVELS = (
    Level("100m", None, C.ROLLUP_FACTORS["100m"], (RIVER_ID, CHAIN_ID)),
    Level("500m", "100m", C.ROLLUP_FACTORS["500m"], (RIVER_ID, CHAIN_ID)),
    Level("1000m", "500m", C.ROLLUP_FACTORS["1000m"], (RIVER_ID, CHAIN_ID)),
    Level("river", None, None, (RIVER_ID,)),
    Level("region", None, None, ("region_id",)),
)


def _base_counts(gdf: gpd.GeoDataFrame) -> pd.DataFrame:
    """Turn each 10 m segment into its additive contributions.

    Segments are not a fixed 10 m: `_segments` divides each ELVIS stretch
    into round(stretch_length / 10) *equal* pieces, so length is uniform
    within a stretch but varies stretch to stretch. Every share below is
    therefore length-weighted -- the summed length of segments in a state,
    divided by the relevant length total in `_ratios` -- rather than a share
    of segment *count*, so a block's percentage reflects how much of its
    physical network is in each state, not how many pieces it was cut into.
    `n_10m_segments`/`n_matched`/`n_ok`/`n_mangl`/`n_ikke` are kept alongside
    as plain counts for audit, but no ratio is computed from them any more.
    """
    out = pd.DataFrame(index=gdf.index)
    length = pd.to_numeric(gdf.get(LENGTH), errors="coerce").fillna(0.0)
    out["length"] = length
    out["n_10m_segments"] = 1
    matched = (
        gdf[CLASS_COL].notna()
        if CLASS_COL in gdf
        else pd.Series(False, index=gdf.index)
    )
    if BEREGNET in gdf:
        matched = gdf[BEREGNET].astype(str).str.lower().eq("yes")
    out["n_matched"] = matched.astype(int)
    out["length_matched"] = length.where(matched, 0.0)
    # Graded Kantvegetasjon (derived/kantveg_grad.py): the length-weighted
    # sum of each assessed segment's fjernet, and the lengths the unit's
    # data coverage is measured on -- bekkelukking left out of both sides.
    fjs = (
        pd.to_numeric(gdf.get("E_S_Kantvegetasjon_fjernet"), errors="coerce")
        if "E_S_Kantvegetasjon_fjernet" in gdf
        else pd.Series(np.nan, index=gdf.index)
    )
    lukket = (
        gdf["K_S_Rorlagt_mistanke"].astype(str).str.strip().eq("Ja")
        if "K_S_Rorlagt_mistanke" in gdf
        else pd.Series(False, index=gdf.index)
    )
    out["length_grad_vurdert"] = length.where(fjs.notna(), 0.0)
    out["sum_grad_fjernet"] = (fjs * length).where(fjs.notna(), 0.0)
    out["length_ikke_lukket"] = length.where(~lukket, 0.0)
    out["length_data"] = length.where(matched & ~lukket, 0.0)
    for key, label in CLASS_COUNTS.items():
        is_class = (
            (gdf.get(CLASS_COL) == label)
            if CLASS_COL in gdf
            else pd.Series(False, index=gdf.index)
        )
        out[key] = is_class.astype(int)
        out[f"length_{key[2:]}"] = length.where(is_class, 0.0)

    for col in SHARE_COLUMNS:
        if col in gdf:
            is_ja = gdf[col].astype(str).str.strip().str.lower().eq("ja")
            out[f"sum_{col}"] = length.where(is_ja, 0.0)
    for col in WORST_COLUMNS:
        if col in gdf:
            vals = gdf[col].astype("string").str.strip()
            for i, label in enumerate(WORST_CLASSES, start=1):
                # A missing class is "not this class", not an error: a lake
                # shore whose lake has no midline in the network inherits
                # NA here, and the nullable eq() would carry it into the
                # int cast (Hurdalsvassdraget: "cannot convert NA to
                # integer" in aggregate_lakes).
                out[f"n_{col}_{i}"] = vals.eq(label).fillna(False).astype(int)

    for col in COUNT_COLUMNS:
        # A raw occurrence count ("Antall"), not a share of the network --
        # weighting it by length would turn "how many barriers" into some
        # length-scaled quantity that answers no real question.
        if col in gdf:
            out[f"sum_{col}"] = (
                gdf[col]
                .astype(str)
                .str.strip()
                .str.lower()
                .eq("ja")
                .astype(float)
            )
    for col in MEAN_COLUMNS:
        if col in gdf:
            val = pd.to_numeric(gdf[col], errors="coerce").fillna(0.0)
            # Restricted to matched segments, because _ratios divides this
            # by length_matched. Summed over every segment instead, a block
            # where few segments were computed returns a "mean" larger than
            # any value that went into it: the worst in an earlier run was a 99.7 m
            # block with 9.97 m matched, giving E_S_Score_total = 2.48 on a
            # 0-1 scale. An unmatched segment has no score to average in --
            # it reads 0 only because that is what the column defaults to.
            out[f"sum_{col}"] = (val * length).where(matched, 0.0)
    # Each segment's own official band (K_S_beregningsbredde: 10 or 20, see
    # primary/channel.py), not always B0-10 -- a block mixing narrow and wide
    # stretches must sum each segment's correct-width area.
    wide = (
        pd.to_numeric(gdf.get("K_S_beregningsbredde"), errors="coerce")
        == C.AREAL_BAND[1]
    )

    def band_col(side, suffix):
        v10 = gdf.get(f"K_{side}_B0-10_{suffix}")
        v20 = gdf.get(f"K_{side}_B0-20_{suffix}")
        if v10 is None and v20 is None:
            return None
        v10 = (
            pd.to_numeric(v10, errors="coerce")
            if v10 is not None
            else pd.Series(np.nan, index=gdf.index)
        )
        v20 = (
            pd.to_numeric(v20, errors="coerce")
            if v20 is not None
            else pd.Series(np.nan, index=gdf.index)
        )
        return v20.where(wide, v10)

    for side in ("L", "R"):
        # How much of the block's length has this bank in Norway at all. A
        # bank abroad contributes 0 to every area sum, which by itself is
        # indistinguishable from a bank that really is bare -- and the block
        # would then be classified Ikke_mulig on the strength of ground
        # nobody surveyed. _block_causes drops a side whose total is 0.
        out[f"length_norsk_{side}"] = length.where(
            sides.present(gdf, side), 0.0
        )
        for cat in CAUSE_AREA_CATEGORIES:
            val = band_col(side, f"{cat}_m2")
            if val is not None:
                out[f"sum_K_{side}_{cat}_m2"] = val.fillna(0.0)
        atot = band_col(side, "Atot_m2")
        for cat in CAUSE_PCT_CATEGORIES:
            pct = band_col(side, cat)
            if pct is not None and atot is not None:
                out[f"sum_K_{side}_{cat}_m2"] = (
                    pct.fillna(0.0) * atot.fillna(0.0) / 100.0
                )
    return out


def _block_causes(agg: pd.DataFrame) -> pd.DataFrame:
    """The six Kantvegetasjon cause columns, reclassified per block.

    Because a percentage is sum(area) / sum(Atot), a segment counted in
    several blocks does not bias the ratio — numerator and denominator scale
    together. Only the categorical columns are produced; the aggregated
    levels never carry per-side coverage percentages.
    """

    def pct(col, atot):
        if f"sum_{col}" not in agg:
            return pd.Series(0.0, index=agg.index)
        return (agg[f"sum_{col}"] / atot.replace(0, np.nan) * 100).fillna(0.0)

    out = pd.DataFrame(index=agg.index)
    status_by_side = {}
    for side in ("L", "R"):
        atot = agg.get(f"sum_K_{side}_Atot_m2", pd.Series(0.0, index=agg.index))
        vtot = pct(f"K_{side}_Vtot_m2", atot)
        potveg = pct(f"K_{side}_PotVeg_m2", atot)
        has_veg = vtot >= C.KANTVEG_THRESHOLD_PCT
        pcts = np.column_stack(
            [
                pct(f"K_{side}_{cand}_m2", atot).to_numpy()
                for cand in kantveg.ARSAK_LABELS
            ]
        )
        # Same rule as the 10 m layer (derived/system_a.py): built ground
        # is potential ground. Applied to the block's summed areas, so the
        # two cannot classify a bank differently.
        infra_rescued = (
            pd.Series(kantveg.infra_potensiell(pcts), index=has_veg.index)
            & ~(potveg >= C.KANTVEG_THRESHOLD_PCT)
            & ~has_veg
        )
        potensiell = (
            (potveg >= C.KANTVEG_THRESHOLD_PCT) | has_veg | infra_rescued
        )
        status, arsak, flag = kantveg.classify_sides(has_veg, potensiell, pcts)
        # No Norwegian ground on this bank anywhere in the block: blank it,
        # and combine_status/combine_causes/combine_flags reduce to the
        # other bank on their own.
        gone = agg.get(
            f"length_norsk_{side}", pd.Series(1.0, index=agg.index)
        ).le(0)
        status_by_side[side] = status.mask(gone, "")
        out[f"E_{side}_Kantvegetasjon_A_arsak"] = arsak.mask(gone, "")
        out[f"E_{side}_Menneskelig_inngrep"] = flag.mask(gone, "")

    left, right = status_by_side["L"], status_by_side["R"]
    summary = kantveg.combine_status(left, right)
    out["E_S_Kantvegetasjon_A_arsak"] = kantveg.combine_causes(
        out["E_L_Kantvegetasjon_A_arsak"],
        out["E_R_Kantvegetasjon_A_arsak"],
    )
    out["E_S_Menneskelig_inngrep"] = kantveg.combine_flags(
        out["E_L_Menneskelig_inngrep"],
        out["E_R_Menneskelig_inngrep"],
        summary,
    )
    return out


def _sum_columns(frame: pd.DataFrame) -> list[str]:
    return [c for c in frame.columns if c.startswith(("n_", "sum_", "length"))]


def _ratios(block: pd.DataFrame, level: str) -> pd.DataFrame:
    """Recompute every ratio from the summed lengths.

    Every share is length-weighted: the summed length of segments in a
    state, over the relevant length total, not a share of segment count --
    see `_base_counts`. `n_matched`/`n_ok`/`n_kantskog_fjernet`/
    `n_nedbygd_kantsone`/`n_ikke`/`n_uklassifisert` still ride along as
    plain counts (kept for audit), just not used here.
    """
    out = pd.DataFrame(index=block.index)
    length_matched = block["length_matched"]
    length_kantskog_fjernet = block["length_kantskog_fjernet"]
    length_nedbygd_kantsone = block["length_nedbygd_kantsone"]
    length_mangl = length_kantskog_fjernet + length_nedbygd_kantsone
    length_relevant = block["length_ok"] + length_mangl

    out[f"frac_matched_{level}"] = np.where(
        block["length"] > 0,
        length_matched / block["length"],
        np.nan,
    )
    out[f"{level}_Uklassifisert"] = np.where(
        block["length"] > 0,
        block["length_uklassifisert"] / block["length"],
        np.nan,
    )
    # Share of the vegetatable length that is MISSING vegetation, so that,
    # like Forbygning and Flomvern, 0 is the undisturbed state and the 1-5
    # grade reads the same way for all three. Kantskog_fjernet/
    # Nedbygd_kantsone are the same "missing" bucket split by cause -- the
    # grade below is computed from their sum, unaffected by the split.
    # The yes/no reading (share of vegetation-capable length missing kantskog)
    # is kept for reference; the official value is the graded one below.
    out[f"{level}_Kantvegetasjon_janei"] = np.where(
        length_relevant > 0, length_mangl / length_relevant, np.nan
    )
    # Graded: length-weighted mean fjernet over assessed segments, 0-1 like
    # the other pressures, and the unit's data coverage.
    if "length_grad_vurdert" in block:
        lv = block["length_grad_vurdert"]
        out[f"{level}_Kantvegetasjon"] = np.where(
            lv > 0, block["sum_grad_fjernet"] / lv.where(lv > 0) / 100.0, np.nan
        )
        li = block["length_ikke_lukket"]
        out[f"{level}_datadekning"] = np.where(
            li > 0, block["length_data"] / li.where(li > 0), np.nan
        )
    else:
        out[f"{level}_Kantvegetasjon"] = out[f"{level}_Kantvegetasjon_janei"]
    out[f"{level}_Kantskog_fjernet"] = np.where(
        length_relevant > 0, length_kantskog_fjernet / length_relevant, np.nan
    )
    out[f"{level}_Nedbygd_kantsone"] = np.where(
        length_relevant > 0, length_nedbygd_kantsone / length_relevant, np.nan
    )
    out[f"{level}_ikke_mulig_Kantvegetasjon"] = np.where(
        length_matched > 0, block["length_ikke"] / length_matched, np.nan
    )
    for col in (*SHARE_COLUMNS, *MEAN_COLUMNS):
        key = f"sum_{col}"
        if key in block:
            out[f"{level}_{col.replace('E_S_', '')}"] = np.where(
                length_matched > 0, block[key] / length_matched, np.nan
            )
    # Barriers aggregate as a count ("Antall"), not a share — dividing by
    # length_matched would report "barriers per metre", which is meaningless.
    for col in COUNT_COLUMNS:
        key = f"sum_{col}"
        if key in block:
            out[f"{level}_{col.replace('E_S_', '')}"] = block[key]

    for col in WORST_COLUMNS:
        keys = [f"n_{col}_{i}" for i in range(1, len(WORST_CLASSES) + 1)]
        if not all(k in block for k in keys):
            continue
        name = col.replace("E_S_", "")
        counts = block[keys].to_numpy()
        present = counts > 0
        # The lowest-numbered class with a member; NaN where the block has
        # no classified segment at all.
        idx = np.where(present.any(axis=1), present.argmax(axis=1), -1)
        out[f"{level}_{name}"] = [
            WORST_CLASSES[i] if i >= 0 else None for i in idx
        ]
        out[f"{level}_{name}_final"] = scales.fragmentering_to_grade(
            out[f"{level}_{name}"]
        ).to_numpy()

    # Restate on the 1-5 grad-av-endring scale. The scale grades *deviation
    # from the reference state*, so every input must have that polarity.
    #
    # All four measure the presence of a pressure, with 0 % the undisturbed
    # state: Kantvegetasjon is the share of vegetatable length *missing*
    # vegetation (n_mangl / (n_ok + n_mangl)), so one set of bins grades them
    # all the same way. Settled 2026-08-27.
    for short in ("Flomvern", "Forbygning", "Oppstuving"):
        col = f"{level}_{short}"
        if col in out:
            out[f"{col}_final"] = scales.pct_to_grade(out[col] * 100.0)
    veg = f"{level}_Kantvegetasjon"
    if veg in out:
        out[f"{veg}_final"] = scales.pct_to_grade(out[veg] * 100.0)
        # Not classified where too little of the unit was measured.
        dek = out.get(f"{level}_datadekning")
        if dek is not None:
            out[f"{veg}_final"] = out[f"{veg}_final"].where(
                pd.Series(dek, index=out.index).fillna(0) >= C.DATADEKNING_MIN
            )
    if f"{level}_Barrier" in out:
        out[f"{level}_Barrier_final"] = scales.barrier_count_to_grade(
            out[f"{level}_Barrier"]
        )
    if f"{level}_Score_total" in out:
        out[f"{level}_Score_total_final"] = scales.score_b_to_grade(
            out[f"{level}_Score_total"] * 100.0
        )
        # The same rule as Kantvegetasjon: no class on too little data. A
        # unit with no coverage value is wholly piped, so nothing was measured.
        dek = out.get(f"{level}_datadekning")
        if dek is not None:
            out[f"{level}_Score_total_final"] = out[
                f"{level}_Score_total_final"
            ].where(
                pd.Series(dek, index=out.index).fillna(0) >= C.DATADEKNING_MIN
            )
    return out.round(4)


def _merge_lines(geoms):
    """linemerge over a block, flattening multi-part inputs first.

    Each level rolls up the previous one, so a block at 500 m is built from
    100 m geometries that are already MultiLineStrings wherever the river had
    a gap. shapely's linemerge takes a sequence of *single* lines — handed a
    MultiLineString element it raises rather than flattening.
    """
    parts = []
    for g in geoms:
        if g is None or g.is_empty:
            continue
        parts.extend(g.geoms if hasattr(g, "geoms") else [g])
    parts = [p for p in parts if p.geom_type == "LineString"]
    if not parts:
        return None
    return parts[0] if len(parts) == 1 else linemerge(parts)


def _blocks(df: pd.DataFrame, group_by, factor: int | None, order_col: str):
    """Assign a block id: consecutive runs of ``factor`` within each group."""
    df = df.sort_values([*group_by, order_col]).copy()
    if factor is None:
        df["_block"] = df.groupby(list(group_by), dropna=False).ngroup()
        return df
    rank = df.groupby(list(group_by), dropna=False).cumcount()
    df["_block"] = (
        df.groupby(list(group_by), dropna=False).ngroup() * 10**9
        + rank // factor
    )
    return df


def _aggregate_one(
    source: gpd.GeoDataFrame, level: Level, order_col: str, is_base: bool
) -> gpd.GeoDataFrame:
    if is_base:
        counts = _base_counts(source)
        # Only what the roll-up reads: block identity and order, the name,
        # and what gets summed. Copying all ~330 columns of the 10 m table
        # here, and _blocks sorting that copy again, held 3x the table --
        # 24 GB on Tana, ~115 GB projected for Glomma.
        need = {*level.group_by, NAME, order_col, RIVER_ID, "region_id"}
        frame = pd.concat(
            [
                source[
                    [
                        c
                        for c in source.columns
                        if c in need or c.startswith(("n_", "sum_", "length"))
                    ]
                ],
                counts,
            ],
            axis=1,
        )
    else:
        frame = source.drop(columns=["geometry"], errors="ignore").copy()
        # Strip the parent level's suffix so the sum columns line up.
        frame = frame.rename(
            columns={
                c: c.rsplit("_", 1)[0]
                for c in frame.columns
                if c.endswith(f"_{level.parent}")
                and c.startswith(("n_", "sum_", "length_"))
            }
        )
        frame = frame.rename(
            columns={f"aggr_segment_length_vann-nett_{level.parent}": "length"}
        )

    frame["geometry"] = source.geometry.values
    if "region_id" in level.group_by and "region_id" not in frame:
        frame["region_id"] = frame[RIVER_ID].astype(str).str.split("-").str[0]

    frame = _blocks(frame, level.group_by, level.factor, order_col)
    sum_cols = _sum_columns(frame)

    grouped = frame.groupby("_block", sort=True)
    agg = grouped[sum_cols].sum()

    # Identity columns are constant within a block by construction.
    for col in (*level.group_by, NAME):
        if col in frame.columns:
            agg[col] = grouped[col].first()
    if level.name == "region":
        agg["n_rivers"] = grouped[RIVER_ID].nunique()
    if level.parent:
        agg[f"n_{level.parent}_segments"] = grouped.size()

    geoms = grouped["geometry"].apply(_merge_lines)

    out = pd.concat([agg, _ratios(agg, level.name), _block_causes(agg)], axis=1)
    out = out.rename(
        columns={
            **{c: f"{c}_{level.name}" for c in sum_cols if c != "length"},
            "length": f"aggr_segment_length_vann-nett_{level.name}",
        }
    )
    out.insert(
        0, f"aggr_segment_id_vann-nett_{level.name}", range(1, len(out) + 1)
    )
    gdf = gpd.GeoDataFrame(out, geometry=geoms.values, crs=source.crs)
    return gdf.reset_index(drop=True)


def aggregate_all(official_10m: gpd.GeoDataFrame, paths) -> dict[str, Path]:
    """Produce every delivery level from the official 10 m layer."""
    written: dict[str, Path] = {}
    produced: dict[str, gpd.GeoDataFrame] = {}

    for level in LEVELS:
        if level.parent is None:
            source, is_base = official_10m, True
            order_col = CHAIN_POS if CHAIN_POS in source else SEGMENT_ID
        else:
            source = produced[level.parent]
            order_col = f"aggr_segment_id_vann-nett_{level.parent}"
            is_base = False
        gdf = _aggregate_one(source, level, order_col, is_base)
        produced[level.name] = gdf
        path = paths.aggregated(level.name)
        save_gpkg(gdf, path, layer=f"official_centerline_{level.name}")
        written[level.name] = path
        logger.info(f"{level.name}: {len(gdf)} blocks")
    return written
