"""Lake shore segments join the river segments without corrupting their ids.

The geometry stage concatenates the lake track's shore segments onto the
river segments. A river segment with no Vann-Nett match carries a missing
vannforekomstid; if the two frames disagree on the column's dtype, pandas
keeps that as a pd.NA object, and GeoPackage stores it as the text "<NA>",
which then passes every notna() test downstream.
"""

import geopandas as gpd
import pandas as pd
import pyogrio
from shapely.geometry import LineString, box

from kantsoner.geometry import lakes
from kantsoner.io import save_gpkg

CRS = "EPSG:25833"


def test_missing_river_ids_stay_null_next_to_lake_shore(tmp_path):
    water = gpd.GeoDataFrame(geometry=[box(0, 0, 200, 100)], crs=CRS)
    vann_nett = gpd.GeoDataFrame(
        {
            lakes.ID_FIELD: ["006-5-L"],
            lakes.NAME_FIELD: ["Testvatnet"],
            lakes.EU_FIELD: [None],
        },
        geometry=[box(-5, -5, 205, 105)],
        crs=CRS,
    )
    shore, *_ = lakes.build(water, vann_nett, CRS)
    assert len(shore)
    # as stages._with_vannforekomst leaves them: one matched, one not
    rivers = gpd.GeoDataFrame(
        {
            "segment_id": ["a", "b"],
            "vannforekomstid": pd.Series(["006-1-R", pd.NA], dtype="string"),
        },
        geometry=[LineString([(300, 0), (310, 0)])] * 2,
        crs=CRS,
    )
    path = tmp_path / "segments.gpkg"
    save_gpkg(pd.concat([rivers, shore], ignore_index=True), path)
    ids = pyogrio.read_dataframe(path, read_geometry=False)["vannforekomstid"]
    assert ids.iloc[0] == "006-1-R"
    assert pd.isna(ids.iloc[1])
    assert (ids.iloc[2:] == "006-5-L").all()


def test_surface_closes_the_hairline_between_polygons():
    """Two polygons of one lake with a 30 cm gap are one surface, no slit."""
    a, b = box(0, 0, 100, 100), box(100.3, 0, 200, 100)
    merged = lakes._surface([a, b])
    assert merged.geom_type == "Polygon"
    assert len(merged.interiors) == 0
    # the boundary is the outline of the whole, not the outline plus the seam
    assert abs(merged.exterior.length - 600.0) < 1.0
    # a single polygon is left alone
    assert lakes._surface([a]).equals(a)


def test_surface_keeps_polygons_that_are_far_apart():
    """A gap wider than twice CLOSE_GAP_M is land and stays a gap."""
    merged = lakes._surface([box(0, 0, 100, 100), box(102, 0, 200, 100)])
    assert merged.geom_type == "MultiPolygon"
