"""A merged buffer gets the segment's fjernet and grad from the merged segment layer."""

import sqlite3
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Point, box

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import merge_delivery as M  # noqa: E402

KEY = "segment_id_vann-nett"


def _segments(tmp_path, ids):
    gdf = gpd.GeoDataFrame(
        {
            KEY: ids,
            "E_S_Kantvegetasjon_fjernet": [4.4, np.nan, 43.4][: len(ids)],
            "E_S_Kantvegetasjon_grad": ["Svært god", "Ikke vurdert", "Moderat"][: len(ids)],
        },
        geometry=[Point(i, 0) for i in range(len(ids))],
        crs="EPSG:25833",
    )
    path = tmp_path / "segments.gpkg"
    gdf.to_file(path, layer="official_centerline_10m", driver="GPKG")
    return path


def _buffer(tmp_path, ids):
    gdf = gpd.GeoDataFrame(
        {KEY: ids},
        geometry=[box(i, 0, i + 1, 1) for i in range(len(ids))],
        crs="EPSG:25833",
    )
    path = tmp_path / "buffer10.gpkg"
    gdf.to_file(path, layer="buffer10", driver="GPKG")
    return path


def _rows(path):
    con = sqlite3.connect(path)
    rows = con.execute(
        f'SELECT "{KEY}", E_S_Kantvegetasjon_fjernet, E_S_Kantvegetasjon_grad FROM buffer10 ORDER BY fid'
    ).fetchall()
    con.close()
    return rows


def test_both_banks_of_a_segment_get_its_values_and_a_stranger_gets_null(tmp_path):
    source = _segments(tmp_path, ["a", "b", "c"])
    target = _buffer(tmp_path, ["a", "a", "b", "c", "zzz"])  # two banks of a, one unknown
    got = M.add_segment_condition(target, "buffer10", source, "official_centerline_10m")
    assert got == {"rows": 5, "without_segment": 1}
    rows = _rows(target)
    assert rows[0] == ("a", 4.4, "Svært god") == (rows[1][0], rows[1][1], rows[1][2])
    assert rows[2] == ("b", None, "Ikke vurdert")
    assert rows[3] == ("c", 43.4, "Moderat")
    assert rows[4] == ("zzz", None, None)
    assert not (tmp_path / "buffer10.keys.db").exists()


def test_a_segment_twice_in_the_layer_is_an_error(tmp_path):
    source = _segments(tmp_path, ["a", "a"])
    target = _buffer(tmp_path, ["a"])
    with pytest.raises(RuntimeError, match="twice"):
        M.add_segment_condition(target, "buffer10", source, "official_centerline_10m")


def test_a_buffer_that_has_the_columns_is_not_given_them_twice(tmp_path):
    source = _segments(tmp_path, ["a"])
    target = _buffer(tmp_path, ["a"])
    M.add_segment_condition(target, "buffer10", source, "official_centerline_10m")
    with pytest.raises(RuntimeError, match="already has"):
        M.add_segment_condition(target, "buffer10", source, "official_centerline_10m")


def test_a_missing_segment_layer_says_what_to_merge_first(tmp_path):
    target = _buffer(tmp_path, ["a"])
    with pytest.raises(RuntimeError, match="must be merged first"):
        M.add_segment_condition(target, "buffer10", tmp_path / "nope.gpkg", "official_centerline_10m")
