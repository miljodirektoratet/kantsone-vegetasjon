"""Columns that are not delivered are dropped from the merged files."""

import sqlite3
import sys
from pathlib import Path

import geopandas as gpd
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import merge_delivery as M  # noqa: E402


def _cols(path, table):
    con = sqlite3.connect(path)
    cols = [r[1] for r in con.execute(f'PRAGMA table_info("{table}")')]
    con.close()
    return cols


def test_level_drop_names():
    assert M._level_drop("river") == [
        "frac_matched_river", "river_Kantvegetasjon_janei", "river_Kantskog_fjernet",
        "river_Nedbygd_kantsone", "river_ikke_mulig_Kantvegetasjon",
        "river_Uklassifisert", "length_uklassifisert_river",
    ]


def test_the_condition_and_the_cause_are_never_dropped():
    keep = {"E_S_Kantvegetasjon_fjernet", "E_S_Kantvegetasjon_grad", "E_S_Kantvegetasjon_A",
            "E_S_Kantvegetasjon_A_arsak", "segment_id_vann-nett", "stretch_id",
            "river_Kantvegetasjon", "river_Kantvegetasjon_final", "river_datadekning"}
    for cols in M.DELIVERY_DROP.values():
        assert not keep & set(cols)


def test_drop_on_an_indexed_buffer_keeps_rows_and_other_columns(tmp_path):
    path = tmp_path / "buffer10.gpkg"
    gpd.GeoDataFrame(
        {"segment_id": ["a", "b"], "segment_id_vann-nett": ["a", "b"],
         "E_Potensiell_Kantvegetasjon": ["Ja", "Nei"], "area_m2": [1.0, 2.0]},
        geometry=[box(0, 0, 1, 1), box(1, 0, 2, 1)], crs="EPSG:25833",
    ).to_file(path, layer="buffer10", driver="GPKG")  # with an R-tree and its triggers
    assert M.drop_delivery_columns(path, "buffer10") == ["segment_id", "E_Potensiell_Kantvegetasjon"]
    assert _cols(path, "buffer10") == ["fid", "geom", "segment_id_vann-nett", "area_m2"]
    back = gpd.read_file(path, layer="buffer10")
    assert len(back) == 2 and list(back["area_m2"]) == [1.0, 2.0]


def test_a_product_without_drops_is_untouched(tmp_path):
    assert M.drop_delivery_columns(tmp_path / "nope.gpkg", "single_trees") == []
