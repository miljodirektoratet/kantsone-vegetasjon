"""Bounded reads, and the CRS trap in the national height models."""

from __future__ import annotations

import numpy as np
import pytest
import rasterio as rio
from affine import Affine
from rasterio.crs import CRS as RioCRS
from shapely.geometry import box

from kantsoner import io


def _write(tmp_path, name, *, crs, origin=(0.0, 100.0), size=100):
    """A small 1 m raster with the given CRS declaration."""
    path = tmp_path / name
    transform = Affine(1.0, 0.0, origin[0], 0.0, -1.0, origin[1])
    data = np.arange(size * size, dtype="float32").reshape(size, size)
    with rio.open(
        path,
        "w",
        driver="GTiff",
        height=size,
        width=size,
        count=1,
        dtype="float32",
        crs=crs,
        transform=transform,
        nodata=-9999.0,
    ) as dst:
        dst.write(data, 1)
    return path


LOCAL_CS = (
    'LOCAL_CS["ETRS89 / UTM zone 33N",UNIT["metre",1,AUTHORITY["EPSG","9001"]],'
    'AXIS["Easting",EAST],AXIS["Northing",NORTH]]'
)


def test_local_cs_is_treated_as_the_project_crs(tmp_path):
    """The national CHM/DSM/DTM declare LOCAL_CS, which has no datum.

    GDAL cannot build a transformer from it, so any reprojection attempt
    raises — and `rio.open(crs=...)` is ignored in read mode, so it cannot be
    relabelled on open either. An unusable declaration must therefore be
    treated as the project CRS with no transform attempted.
    """
    path = _write(tmp_path, "local.tif", crs=RioCRS.from_wkt(LOCAL_CS))
    src = io.open_raster(path)
    assert src.crs.to_epsg() is None  # the trap
    assert io.effective_crs(src).to_epsg() == 25833

    # Reading against a real EPSG must not raise, and must not shift.
    arr, transform = io.read_window(path, box(10, 10, 20, 20), "EPSG:25833")
    assert arr.size > 0
    assert transform.c == pytest.approx(10.0)


def test_real_crs_still_reprojects(tmp_path):
    """A properly declared, *different* CRS must still be transformed."""
    path = _write(tmp_path, "utm32.tif", crs="EPSG:25832", origin=(0.0, 100.0))
    src = io.open_raster(path)
    assert io.effective_crs(src).to_epsg() == 25832


def test_window_is_clamped_to_the_dataset(tmp_path):
    """A geometry extending past the edge must not read garbage offsets."""
    path = _write(tmp_path, "clamp.tif", crs="EPSG:25833")
    arr, _ = io.read_window(path, box(-500, -500, 50, 50), "EPSG:25833")
    assert arr.shape[0] <= 100 and arr.shape[1] <= 100


def test_mask_burns_the_sentinel(tmp_path):
    import geopandas as gpd

    path = _write(tmp_path, "mask.tif", crs="EPSG:25833")
    footprint = gpd.GeoDataFrame(
        geometry=[box(20, 20, 40, 40)], crs="EPSG:25833"
    )
    arr, _ = io.read_window(
        path, box(0, 0, 100, 100), "EPSG:25833", mask=footprint
    )
    assert (arr == io.MASK_VALUE).any()


def test_clip_drops_sliver_geometries(tmp_path):
    """Clipping polygons leaves LineString/GeometryCollection slivers along the
    cut, and gpd.overlay refuses a frame of mixed families outright."""
    import geopandas as gpd
    from shapely.geometry import LineString, Point, Polygon

    src = tmp_path / "mixed.gpkg"
    gpd.GeoDataFrame(
        {"v": [1, 2, 3, 4]},
        geometry=[
            Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]),
            Polygon([(20, 20), (30, 20), (30, 30), (20, 30)]),
            LineString([(0, 0), (5, 5)]),
            Point(1, 1),
        ],
        crs="EPSG:25833",
    ).to_file(src, driver="GPKG")

    aoi = gpd.GeoDataFrame(
        geometry=[Polygon([(-5, -5), (40, -5), (40, 40), (-5, 40)])],
        crs="EPSG:25833",
    )
    out = io.read_vector(src, aoi, clip=True)
    assert set(out.geom_type.str.replace("Multi", "", regex=False)) == {
        "Polygon"
    }
    assert len(out) == 2


def test_geometry_collections_keep_their_polygonal_parts():
    """A collection must be reduced to its winning-family parts, not dropped —
    otherwise real area disappears from the coverage denominators."""
    import geopandas as gpd
    from shapely.geometry import GeometryCollection, LineString, Polygon

    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    gdf = gpd.GeoDataFrame(
        geometry=[
            poly,
            GeometryCollection(
                [
                    Polygon([(20, 0), (25, 0), (25, 5), (20, 5)]),
                    LineString([(0, 0), (1, 1)]),
                ]
            ),
        ],
        crs="EPSG:25833",
    )
    out = io._dominant_geometry_family(gdf)
    assert len(out) == 2
    assert out.area.sum() == pytest.approx(poly.area + 25.0)


def test_clip_survives_an_invalid_mask(tmp_path):
    """Catchment outlines are not all valid — Lærdal's is self-intersecting —
    and GEOS raises on the intersection rather than skipping it."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    src = tmp_path / "lines.gpkg"
    gpd.GeoDataFrame(
        geometry=[Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])],
        crs="EPSG:25833",
    ).to_file(src, driver="GPKG")

    bowtie = Polygon(
        [(-5, -5), (20, 20), (20, -5), (-5, 20)]
    )  # self-intersecting
    assert not bowtie.is_valid
    aoi = gpd.GeoDataFrame(geometry=[bowtie], crs="EPSG:25833")
    # make_valid turns a bowtie into a GeometryCollection; clip needs polygons.
    out = io.read_vector(src, aoi, clip=True)
    assert len(out) >= 0  # the point is that it does not raise
