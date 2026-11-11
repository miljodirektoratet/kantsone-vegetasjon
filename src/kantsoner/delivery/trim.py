"""The delivered column set, written by the pipeline rather than after it.

The trim belongs where the files are made: a run produces `delivery/`
complete, from one keep list.

Two rules decide what stays:

* identity, so a row can be traced back and joined to Vann-Nett;
* every delivered measurement, and nothing it was derived from. The full
  tables stay in `aggregation/` for anyone who needs the workings.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pyogrio
from loguru import logger

from ..io import read_gpkg, save_gpkg

LEVELS = ("10m", "100m", "500m", "1000m", "river", "region")

FIELDS_10M = [
    "segment_id",
    "segment_id_vann-nett",
    "river_id_vann-nett",
    "vannforekomstid",
    "vannforekomstnavn",
    "segment_length_vann-nett",
    "elvenavn",
    "vassdragsnummer",
    "K_S_Kilde",
    "K_S_Usikkerhet",
    "K_S_avstand_beregning",
    "K_S_beregningsbredde",
    "K_S_Laser_aar",
    "K_S_Grenseelv",
    "K_S_Rorlagt_mistanke",
    "K_S_bekk_bredde",
    "K_L_Bredde",
    "K_R_Bredde",
    "E_L_Potensiell_Kantvegetasjon",
    "E_R_Potensiell_Kantvegetasjon",
    "E_S_Kantvegetasjon_A",
    "E_S_Kantvegetasjon_A_arsak",
    "K_L_B0-10_Forventet_dekning",
    "K_R_B0-10_Forventet_dekning",
    "E_L_Kantvegetasjon_fjernet",
    "E_R_Kantvegetasjon_fjernet",
    "E_S_Kantvegetasjon_fjernet",
    "E_S_Kantvegetasjon_grad",
    "E_S_Menneskelig_inngrep",
    "E_L_Flomvern",
    "E_R_Flomvern",
    "E_S_Flomvern",
    "E_S_Forbygning",
    "E_S_Barrier",
    "E_S_Fragmentering",
    "E_S_Filter_score_final",
    "E_S_Struktur_score",
    "E_S_Kontinuitet",
    "E_S_Score_total",
    "E_S_usikkerhet_flagg",
    "E_S_K_score",
    "K_S_prioritet",
    "K_S_fremmedart_risiko",
    "K_S_fremmedart_art",
    "K_L_B0-10_Vtot",
    "K_R_B0-10_Vtot",
    "K_L_B0-20_Vtot",
    "K_R_B0-20_Vtot",
]

FIELDS_LEVEL = [
    "aggr_segment_id_vann-nett_{L}",
    "aggr_segment_length_vann-nett_{L}",
    "river_id_vann-nett",
    "vannforekomstnavn",
    "frac_matched_{L}",
    "{L}_Uklassifisert",
    "{L}_Kantvegetasjon",
    "{L}_Kantvegetasjon_final",
    "{L}_datadekning",
    "{L}_Kantvegetasjon_janei",
    "{L}_Kantskog_fjernet",
    "{L}_Nedbygd_kantsone",
    "{L}_ikke_mulig_Kantvegetasjon",
    "{L}_Flomvern",
    "{L}_Flomvern_final",
    "{L}_Forbygning",
    "{L}_Forbygning_final",
    "{L}_Barrier",
    "{L}_Barrier_final",
    "{L}_Filter_score_final",
    "{L}_Struktur_score",
    "{L}_Kontinuitet",
    "{L}_Score_total",
    "{L}_Score_total_final",
    "E_L_Kantvegetasjon_A_arsak",
    "E_R_Kantvegetasjon_A_arsak",
    "E_S_Kantvegetasjon_A_arsak",
    "E_L_Menneskelig_inngrep",
    "E_R_Menneskelig_inngrep",
    "E_S_Menneskelig_inngrep",
    "length_uklassifisert_{L}",
]

# The lake polygon: who it is, how much shore was measured, and the same
# aggregated values the river levels carry under their own prefix.
FIELDS_LAKE_POLYGON = [
    "lake_id",
    "lake_navn",
    "lake_eu",
    "areal_vannflate_m2",
    "areal_vannforekomst_m2",
    "innsjo_strandlinje_m",
    "innsjo_segmenter",
    "frac_matched_innsjo",
    "innsjo_Uklassifisert",
    "innsjo_Kantvegetasjon",
    "innsjo_Kantvegetasjon_final",
    "innsjo_datadekning",
    "innsjo_Kantvegetasjon_janei",
    "innsjo_Kantskog_fjernet",
    "innsjo_Nedbygd_kantsone",
    "innsjo_ikke_mulig_Kantvegetasjon",
    "innsjo_Flomvern",
    "innsjo_Flomvern_final",
    "innsjo_Forbygning",
    "innsjo_Forbygning_final",
    "innsjo_Barrier",
    "innsjo_Barrier_final",
    "innsjo_Filter_score_final",
    "innsjo_Kontinuitet",
    "innsjo_Score_total",
    "innsjo_Score_total_final",
]

# The shore buffer: the 10 m list, minus what a lake has no answer for --
# channel width and everything keyed on it -- plus the lake it belongs to
# and the area the numbers were measured over.
_LAKE_DROP = {
    "K_S_bekk_bredde",
    "K_S_beregningsbredde",
    "K_S_avstand_beregning",
    "E_S_Struktur_score",
    "K_R_Bredde",
    "E_R_Potensiell_Kantvegetasjon",
    "K_R_B0-10_Forventet_dekning",
    "E_R_Kantvegetasjon_fjernet",
    "E_R_Flomvern",
    "K_R_B0-10_Vtot",
    "K_R_B0-20_Vtot",
    "K_S_Grenseelv",
    "K_S_Rorlagt_mistanke",
}
FIELDS_LAKE_SHORE = ["lake_id", "lake_navn", "area_m2"] + [
    f for f in FIELDS_10M if f not in _LAKE_DROP
]


def fields_for(level: str) -> list[str]:
    if level == "10m":
        return list(FIELDS_10M)
    out = [f.replace("{L}", level) for f in FIELDS_LEVEL]
    if level == "region":
        i = out.index("river_id_vann-nett")
        out[i : i + 1] = ["region_id", "n_rivers"]
    return out


def _write(src: Path, dst: Path, keep: list[str]) -> int:
    # Read only what is kept: the aggregation stage still holds its own copy
    # of the 10 m table while this runs, and all of its ~330 columns again
    # would not fit beside it at Glomma's size.
    fields = set(pyogrio.read_info(src)["fields"])
    have = [c for c in keep if c in fields]
    missing = [c for c in keep if c not in fields]
    gdf = read_gpkg(src, columns=have)
    out = gpd.GeoDataFrame(gdf[have], geometry=gdf.geometry.values, crs=gdf.crs)
    save_gpkg(out, dst)
    if missing:
        logger.debug(f"{dst.name}: not in source, left out: {missing}")
    return len(have)


def run(paths) -> None:
    """Write `delivery/` from `aggregation/`, trimmed."""
    out_dir = paths.root / "delivery"
    out_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for level in LEVELS:
        src = paths.aggregated(level)
        if src.exists():
            _write(src, out_dir / src.name, fields_for(level))
            n += 1
    for name, keep in (
        ("lake_polygon.gpkg", FIELDS_LAKE_POLYGON),
        ("lake_shore.gpkg", FIELDS_LAKE_SHORE),
    ):
        src = paths.aggregation_dir / name
        if src.exists():
            _write(src, out_dir / name, keep)
            n += 1
    logger.info(f"Delivery: {n} trimmed tables -> {out_dir}")
