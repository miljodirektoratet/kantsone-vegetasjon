"""Reading and writing.

Bounded reads straight from the national datasets — nothing is clipped or
copied to disk. Which FKB/grunnkart tiles cover an AOI comes from the
``gpkg_contents`` table already inside every GeoPackage. A FileGDB source
(grunnkart's national delivery, undivided) has no such table; its layers are
indexed from OGR's own metadata instead — see ``_gdb_layers``.

Writes go via a local temp file, because SQLite-backed formats are unreliable
on the SMB/NFS shares these outputs land on.
"""

from __future__ import annotations

import sqlite3
import tempfile
from functools import cache, lru_cache
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import rasterio as rio
import shapely
from loguru import logger
from pyproj import CRS as PyCRS
from pyproj import Transformer
from rasterio.features import rasterize
from rasterio.windows import Window, from_bounds
from shapely.geometry import box

from .spec.constants import CRS

# Burned into the canopy model where infrastructure covers it.
MASK_VALUE = -1.0


# ── writing ────────────────────────────────────────────────────────────────
def _place(local: Path, dest: Path) -> None:
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    tmp.write_bytes(local.read_bytes())
    tmp.replace(dest)


def save_gpkg(gdf, path: Path | str, layer: str | None = None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        local = Path(tmp) / path.name
        gdf.to_file(local, driver="GPKG", layer=layer or path.stem)
        _place(local, path)
    logger.info(f"Wrote {len(gdf)} features -> {path}")


def _release_memory() -> None:
    """Hand freed heap back to the OS between batches.

    glibc retains freed arenas for reuse, so a worker's RSS ratchets up to
    its largest batch and never comes down — on a national run, six workers
    each pinned ~25 GB of dead heap and the machine ran out of memory while
    most of it was free in the allocator's own books. malloc_trim returns
    it; a no-op on non-glibc platforms.
    """
    import ctypes
    import gc

    gc.collect()
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except OSError:
        pass


def save_gpkg_incremental(frames, path: Path | str, layer=None) -> None:
    """Write frames one at a time so the full table never sits in memory.

    The batched stages produce their output in pieces; accumulating those
    pieces to concatenate at the end would undo the point of batching.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        local = Path(tmp) / path.name
        n = 0
        for gdf in frames:
            gdf.to_file(
                local,
                driver="GPKG",
                layer=layer or path.stem,
                mode="a" if n else "w",
            )
            n += len(gdf)
            del gdf
            _release_memory()
        _place(local, path)
    logger.info(f"Wrote {n} features -> {path}")


# ── vector reading ─────────────────────────────────────────────────────────
def list_layers(path: Path | str) -> list[str]:
    from pyogrio import list_layers as _list

    return [str(row[0]) for row in _list(str(path))]


def input_layers(source) -> list[str | None]:
    """Layer names to iterate for an input that may be a frame or a path.

    An in-memory GeoDataFrame is a single unnamed layer; `as_gdf` ignores
    `layer=None`, so callers keep one loop for both cases.
    """
    if isinstance(source, (gpd.GeoDataFrame, list, tuple)):
        return [None]
    if source is None:
        raise ValueError(
            "input is None — an upstream step returned nothing where a frame "
            "or a path was expected"
        )
    return list_layers(source)


def read_gpkg(path, layer: str | None = None, **kw) -> gpd.GeoDataFrame:
    return gpd.read_file(path, layer=layer, engine="pyogrio", **kw)


def read_all_layers(path, **kw) -> gpd.GeoDataFrame:
    """Union of every layer in a GeoPackage, at the source's own CRS.

    ``gpd.read_file(path)`` on a multi-layer GeoPackage silently reads only
    one layer (with a warning) -- fine for FKB/grunnkart, where a source is
    deliberately handed one named layer at a time, but wrong for a delivery
    that splits one dataset's features across layers by geometry type (e.g.
    Statens Vegvesen's PlastringErosjonssikri: polygon revetment areas and
    line-mapped protection in separate layers of the same file), where
    reading only the default layer drops real features rather than
    unwanted ones. Reprojecting, if needed, is the caller's job -- same as
    every other path handed to `as_gdf`.
    """
    names = pyogrio.list_layers(path)[:, 0]
    frames = [read_gpkg(path, layer=str(name), **kw) for name in names]
    frames = [f for f in frames if len(f)]
    if not frames:
        return gpd.GeoDataFrame(geometry=[])
    return gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True),
        geometry="geometry",
        crs=frames[0].crs,
    )


def as_gdf(source, layer: str | None = None, crs=None) -> gpd.GeoDataFrame:
    """Accept a GeoDataFrame, a path, or a list of either.

    This is what lets `read_vector` hand its result straight to the analysis
    functions — all written to take paths — instead of round-tripping it
    through a working file.
    """
    if isinstance(source, (list, tuple)):
        frames = [f for f in (as_gdf(s, layer, crs) for s in source) if len(f)]
        if not frames:
            return gpd.GeoDataFrame(geometry=[], crs=crs)
        source = gpd.GeoDataFrame(
            pd.concat(frames, ignore_index=True),
            geometry="geometry",
            crs=crs or frames[0].crs,
        )
    elif not isinstance(source, gpd.GeoDataFrame):
        source = read_gpkg(source, layer=layer)
    if crs and source.crs and str(source.crs) != str(crs):
        source = source.to_crs(crs)
    return source


def _ogr_layers(root: Path) -> tuple[tuple, ...]:
    """(path, layer, bounds, srs_id) for every layer in a non-GeoPackage
    single-file vector source (FileGDB, Shapefile, ...).

    These carry no gpkg_contents-style extent table, but OGR reports each
    layer's extent and CRS cheaply from its own metadata (fast_total_bounds,
    no feature decoded) regardless of driver. Each layer becomes one "tile"
    spanning the whole file -- there is no per-region split to reject tiles
    against, so the real narrowing happens in read_vector's bbox pushdown,
    answered from whatever spatial index the driver keeps (FileGDB's .spx;
    a Shapefile's sidecar .qix/.sbn if present, else a full scan).
    """
    out = []
    for name in pyogrio.list_layers(root)[:, 0]:
        info = pyogrio.read_info(root, layer=name)
        epsg = PyCRS.from_user_input(info["crs"]).to_epsg()
        out.append((root, str(name), tuple(info["total_bounds"]), epsg))
    return tuple(out)


@cache
def tile_index(root: Path) -> tuple[tuple, ...]:
    """(path, layer, bounds, srs_id) for every feature layer under ``root``.

    ``root`` is a .gpkg, a directory of them (only files directly inside a
    directory count, matching the per-kommune / per-fylke delivery layout),
    or a single non-GeoPackage vector source (FileGDB, Shapefile, ...) -- a
    national delivery with no per-tile split, handled by ``_ogr_layers``
    instead: it has no gpkg_contents table, so opening it as if it did
    fails outright ("file is not a database"), not just slowly. A FileGDB
    is itself a directory, so this is decided on suffix, not `is_dir()` --
    the flip side of the same "not just single files" fix that added
    Shapefile support: gating on `is_dir()` alone silently routed .gdb back
    into the directory-of-.gpkg-tiles branch below, which found none.
    Read from gpkg_contents with sqlite, no GDAL open: ~18 s cold and ~1.4 s
    warm for FKB's 357 tiles over CIFS, memoised for the run.
    """
    root = Path(root)
    if root.suffix.lower() == ".gdb" or (
        not root.is_dir() and root.suffix.lower() != ".gpkg"
    ):
        return _ogr_layers(root)
    paths = sorted(root.glob("*.gpkg")) if root.is_dir() else [root]
    if not paths:
        raise FileNotFoundError(f"No .gpkg files in {root}")

    out = []
    for path in paths:
        con = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
        try:
            rows = con.execute(
                "SELECT table_name, min_x, min_y, max_x, max_y, srs_id "
                "FROM gpkg_contents WHERE data_type = 'features'"
            ).fetchall()
        finally:
            con.close()
        for name, *bounds, srs_id in rows:
            # NULL bbox means unknown, not empty — keep it as a candidate.
            known = None if None in bounds else tuple(map(float, bounds))
            out.append((path, str(name), known, srs_id))
    return tuple(out)


@cache
def _transformer(src_epsg: int, dst_epsg: int) -> Transformer:
    return Transformer.from_crs(src_epsg, dst_epsg, always_xy=True)


def _bounds_in(bounds, src_crs, dst_srs_id):
    src_epsg = PyCRS.from_user_input(src_crs).to_epsg()
    if not dst_srs_id or not src_epsg or dst_srs_id == src_epsg:
        return bounds
    return _transformer(src_epsg, dst_srs_id).transform_bounds(*bounds)


def _dominant_geometry_family(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Keep only the majority geometry family (polygon / line / point).

    Clipping a polygon layer leaves slivers along the cut — stray LineStrings
    and GeometryCollections — and `gpd.overlay` refuses a frame of mixed
    families outright. Collections are reduced to their parts of the winning
    family rather than dropped, so real area is not lost.
    """
    family = gdf.geom_type.str.replace("Multi", "", regex=False)
    known = family[family.isin(("Polygon", "Point", "LineString"))]
    if known.empty:
        return gdf
    want = known.mode().iloc[0]

    collections = family == "GeometryCollection"
    if collections.any():
        gdf = gdf.copy()
        gdf.loc[collections, "geometry"] = gdf.loc[
            collections, "geometry"
        ].apply(
            lambda g: shapely.union_all(
                [p for p in g.geoms if p.geom_type.replace("Multi", "") == want]
            )
        )
        family = gdf.geom_type.str.replace("Multi", "", regex=False)

    out = gdf[family == want]
    if len(out) < len(gdf):
        logger.debug(f"Dropped {len(gdf) - len(out)} non-{want} clip slivers")
    return out


def read_vector(
    root,
    aoi,
    *,
    layer: str | None = None,
    where: str | None = None,
    columns: list[str] | None = None,
    crs=CRS,
    clip: bool = False,
    exact_mask: bool = False,
) -> gpd.GeoDataFrame:
    """Every feature intersecting ``aoi``, from a file or a directory of tiles.

    ``where`` is pushed into the driver by pyogrio, so attribute filters never
    materialise the rows they reject. The AOI is transformed into *each* tile's
    own CRS — FKB is EPSG:3044 while grunnkart is 25833, so a single
    project-CRS bbox would silently match nothing.
    """
    aoi, aoi_crs = _as_geom(aoi, crs)
    frames = []
    for path, tile_layer, bounds, srs_id in tile_index(Path(root)):
        if layer is not None and tile_layer != layer:
            continue
        bbox = _bounds_in(aoi.bounds, aoi_crs, srs_id)
        if bounds and not box(*bounds).intersects(box(*bbox)):
            continue
        # exact_mask filters driver-side by the AOI *geometry*, so only
        # intersecting features are ever materialised. A bbox read of a
        # riparian AOI pulls the whole rectangle: on a sprawling northern
        # batch that was millions of land-cover polygons -- 25 GB -- read
        # into memory only to be clipped away. bbox and mask are mutually
        # exclusive in pyogrio.
        spatial = (
            {"mask": _geom_in(aoi, aoi_crs, srs_id)}
            if exact_mask
            else {"bbox": bbox}
        )
        gdf = gpd.read_file(
            path,
            layer=tile_layer,
            where=where,
            columns=columns,
            engine="pyogrio",
            **spatial,
        )
        if len(gdf):
            frames.append(gdf.to_crs(crs))

    if not frames:
        logger.warning(f"No features from {Path(root).name} intersect the AOI")
        return gpd.GeoDataFrame(geometry=[], crs=crs)

    out = gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True), geometry="geometry", crs=crs
    )
    if clip:
        # Grunnkart ships some self-intersecting polygons; GEOS raises a
        # "side location conflict" on those rather than skipping them.
        # Both sides can be self-intersecting: grunnkart ships invalid
        # polygons, and some vassdragomr catchment outlines are invalid too
        # (Lærdal is). GEOS raises "side location conflict" rather than
        # skipping them, so repair both before the clip.
        invalid = ~out.geometry.is_valid
        if invalid.any():
            logger.debug(f"Repairing {invalid.sum()} invalid geometries")
            out.loc[invalid, "geometry"] = out.loc[
                invalid, "geometry"
            ].make_valid()
        mask = gpd.GeoSeries([aoi], crs=aoi_crs).to_crs(crs).iloc[0]
        if not mask.is_valid:
            # make_valid on a self-intersecting outline returns a
            # GeometryCollection (polygons plus the offending edges); clip
            # only accepts polygonal input.
            repaired = shapely.make_valid(mask)
            parts = [
                g
                for g in shapely.get_parts(repaired)
                if g.geom_type in ("Polygon", "MultiPolygon")
            ]
            mask = shapely.union_all(parts) if parts else mask
        out = gpd.clip(out, mask)
        out = _dominant_geometry_family(out)
    logger.debug(f"{Path(root).name}: {len(out)} features, {len(frames)} tiles")
    return out.reset_index(drop=True)


# ── raster reading ─────────────────────────────────────────────────────────
@lru_cache(maxsize=8)
def open_raster(path: str) -> rio.DatasetReader:
    """Cached handle — reopening per window over CIFS is what would make
    bounded reads slow. Do not close it."""
    return rio.open(path)


def effective_crs(src) -> rio.crs.CRS:
    """The CRS a dataset's coordinates are really in.

    The national height models declare LOCAL_CS — an engineering CRS with no
    datum — though their coordinates are plain UTM 33N. GDAL cannot build a
    transformer from it, and ``rio.open(crs=...)`` is ignored in read mode, so
    an unresolvable declaration is treated as the project CRS instead.

    Returns a CRS object, not a string: callers reproject with it and read
    `.to_epsg()` off it.
    """
    try:
        epsg = src.crs.to_epsg()
    except Exception:
        epsg = None
    return src.crs if epsg else rio.crs.CRS.from_string(CRS)


def _geom_in(geom, src_crs, dst_srs_id):
    """`geom` transformed into the tile's own CRS, for driver-side masking."""
    src_epsg = PyCRS.from_user_input(src_crs).to_epsg()
    if not dst_srs_id or not src_epsg or dst_srs_id == src_epsg:
        return geom
    return gpd.GeoSeries([geom], crs=src_crs).to_crs(epsg=dst_srs_id).iloc[0]


def _as_geom(geom, crs):
    if isinstance(geom, (gpd.GeoDataFrame, gpd.GeoSeries)):
        return geom.union_all(), (geom.crs or crs)
    return geom, crs


def read_window(path, geom, crs=CRS, *, pad: int = 0, mask=None):
    """The pixels covering ``geom``, plus ``pad`` px of halo.

    ``pad`` matters for gradient work: a slope computed without a one-pixel
    halo is wrong along every polygon edge. ``mask`` burns MASK_VALUE where
    infrastructure polygons cover the window.
    """
    geom, crs = _as_geom(geom, crs)
    src = open_raster(str(path))
    src_crs = effective_crs(src)
    if crs and str(src_crs) != str(crs):
        geom = gpd.GeoSeries([geom], crs=crs).to_crs(src_crs).iloc[0]

    win = from_bounds(*geom.bounds, transform=src.transform)
    win = win.round_offsets(op="floor").round_lengths(op="ceil")
    if pad:
        win = Window(
            win.col_off - pad,
            win.row_off - pad,
            win.width + 2 * pad,
            win.height + 2 * pad,
        )
    win = win.intersection(Window(0, 0, src.width, src.height))
    if win.width <= 0 or win.height <= 0:
        return np.ma.masked_all((0, 0)), src.transform

    arr = src.read(1, window=win, masked=True)
    transform = src.window_transform(win)

    if mask is not None and len(mask):
        shapes = mask.to_crs(src_crs).geometry
        bounds = rio.windows.bounds(win, src.transform)
        shapes = shapes[shapes.intersects(box(*bounds))]
        if len(shapes):
            burned = rasterize(
                ((g, 1) for g in shapes),
                out_shape=arr.shape,
                transform=transform,
                fill=0,
                dtype="uint8",
            )
            arr = np.ma.where(burned == 1, MASK_VALUE, arr)
    return arr, transform
