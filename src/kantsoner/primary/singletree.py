from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio as rio
from loguru import logger
from rasterio.features import geometry_window
from rasterio.features import rasterize as rio_rasterize
from rasterio.features import shapes as rio_shapes
from rasterio.mask import mask as rio_mask
from rasterio.windows import Window
from rasterio.windows import bounds as window_bounds
from rasterstats import zonal_stats
from scipy.ndimage import (
    gaussian_filter,
    maximum,
    maximum_position,
    mean,
    median_filter,
    minimum,
)
from shapely import STRtree, box
from shapely.geometry import Point
from skimage.segmentation import watershed

from ..io import as_gdf, effective_crs, input_layers, save_gpkg

# Defaults for the tiled (memory-bounded) path. A buffer-crop window larger than
# MAX_PIXELS is processed tile-by-tile instead of as one array; smaller windows
# keep the original whole-array path (byte-identical results). See the tiled
# producer below for why the halo must exceed the filter footprint + crown radius.
# Rivers are linear, so the *bounding box* of a detection mask is mostly empty
# — 210 Mpx of bbox for 11 Mpx of riparian zone on the pilot catchment. The
# whole-array path reads all of it; the tiled path skips every tile that does
# not touch the mask. For this workload tiling is not a memory fallback, it is
# the right default, so the threshold is set low enough that it always wins.
# Raise it only to reproduce the old behaviour for comparison.
DEFAULT_MAX_PIXELS = 16_000_000  # ~16 Mpix (~64 MB float32)
DEFAULT_TILE_SIZE = 8192  # core tile edge (px)
DEFAULT_TILE_HALO = 96  # halo margin (px); > 14 px filter reach + crown radius


def singletree(
    ndsm: np.ndarray,
    min_height: float = 0.5,
    compactness: float = 0.1,
    median_size: int = 3,
    gaussian_sigma: float = 1.0,
    tile: gpd.GeoDataFrame = None,
) -> np.ndarray:
    """Process NDSM array to detect single trees using watershed segmentation.

    Args:
        ndsm (np.ndarray): Input NDSM array
        nodata_value: Value indicating nodata in the input array
        min_height (float): Minimum height to consider as tree
        compactness (float): Watershed compactness parameter
        median_size (int): Size of median filter
        gaussian_sigma (float): Sigma for Gaussian filter

    Returns:
        tuple: (result raster with shape (1, h, w), watershed labels with shape (h, w), original ndsm 2d)
    """
    ndsm = ndsm[0]

    ndsm_smoothed = median_filter(ndsm, size=median_size)
    nodata_mask_median = (ndsm_smoothed < min_height) | np.isnan(ndsm_smoothed)
    ndsm_smoothed = gaussian_filter(ndsm_smoothed, sigma=gaussian_sigma)
    nodata_mask = (
        (ndsm_smoothed < min_height) | np.isnan(ndsm_smoothed) | np.isnan(ndsm)
    )

    # Invert for watershed
    elevation_inv = -ndsm_smoothed
    elevation_inv[nodata_mask] = 0

    # Apply watershed segmentation
    segments = watershed(
        elevation_inv,
        mask=~nodata_mask,
        connectivity=1,
        compactness=compactness,
    )

    # Get maximum height per segment using original ndsm
    ndsm = np.nan_to_num(ndsm, nan=0)
    ndsm[nodata_mask] = 0
    result = maximum(ndsm, segments, index=segments)
    result[nodata_mask_median] = np.nan

    # Add a channel dimension
    result_raster = result[None]

    return result_raster, segments, ndsm


# ── Shared helpers (used by both the whole-array and tiled paths) ────────────


def _load_rivers(rivers_path, crs):
    """Read all layers of the rivers GeoPackage, concatenated, reprojected."""
    layers = input_layers(rivers_path)
    gdfs = [as_gdf(rivers_path, layer=ly).to_crs(crs) for ly in layers]
    if len(gdfs) > 1:
        return gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), crs=crs)
    return gdfs[0]


def _rasterize_inside(geoms, out_shape, transform):
    """Boolean mask: True where a pixel is covered by any of ``geoms``."""
    valid = [(g, 1) for g in geoms if g is not None]
    if not valid:
        return np.zeros(out_shape, dtype=bool)
    return rio_rasterize(
        valid, out_shape=out_shape, transform=transform, fill=0, dtype=np.uint8
    ).astype(bool)


def _water_suppress(
    ndsm_band, river_geoms, transform, water_suppression_height
):
    """Zero nDSM water-surface artifacts inside river polygons, in place.

    Pixels inside a river polygon with 0 < ndsm < water_suppression_height are
    set to 0 (then excluded by the watershed min_height threshold); higher
    pixels (real overhanging canopy) are preserved. Returns (n_water,
    n_suppressed). Mirrors the original inline block.
    """
    water_mask = _rasterize_inside(river_geoms, ndsm_band.shape, transform)
    n_water = int(water_mask.sum())
    suppress = (
        water_mask & (ndsm_band > 0) & (ndsm_band < water_suppression_height)
    )
    n_suppressed = int(suppress.sum())
    ndsm_band[suppress] = 0
    return n_water, n_suppressed


def _segment_stats(ndsm_2d, labels, label_ids):
    """Per-segment reductions over the label array (one entry per label id).

    Returns min/mean/max height and the apex (row, col) of each segment, using
    the scipy.ndimage reductions the original code already used.
    """
    mn = np.asarray(minimum(ndsm_2d, labels, index=label_ids), dtype=float)
    mean_h = np.asarray(mean(ndsm_2d, labels, index=label_ids), dtype=float)
    mx = np.asarray(maximum(ndsm_2d, labels, index=label_ids), dtype=float)
    pos = np.atleast_2d(
        np.asarray(maximum_position(ndsm_2d, labels, index=label_ids))
    )
    return {
        "min_height": mn,
        "mean_height": mean_h,
        "max_height": mx,
        "apex_row": pos[:, 0].astype(int),
        "apex_col": pos[:, 1].astype(int),
    }


def _build_top_points_gdf(
    ids, top_x, top_y, tree_height, min_h, mean_h, max_h, dtm_path, ndsm_crs
):
    """Build the tree-top points GeoDataFrame (geometry + all scalar columns).

    Samples the DTM once for bottom_z and derives top_z / height_class exactly
    as the original code did (lines for bottom_z, abs_z, height_class).
    """
    top_x = np.asarray(top_x)
    top_y = np.asarray(top_y)
    tree_height = np.asarray(tree_height)

    with rio.open(dtm_path) as dtm_src:
        bottom_z = np.array(
            list(dtm_src.sample(zip(top_x, top_y, strict=True)))
        ).flatten()

    abs_z = tree_height + bottom_z
    top_pts = [
        Point(x, y, z)
        for x, y, z in zip(top_x, top_y, tree_height, strict=True)
    ]
    gdf = gpd.GeoDataFrame(geometry=top_pts, crs=ndsm_crs.to_epsg())
    gdf["id"] = list(ids)
    gdf["bottom_z"] = bottom_z.astype(np.float64).round(2)
    gdf["tree_height"] = tree_height.astype(np.float64).round(2)
    gdf["top_x"] = top_x.astype(np.float64).round(2)
    gdf["top_y"] = top_y.astype(np.float64).round(2)
    gdf["top_z"] = np.asarray(abs_z).astype(np.float64).round(2)
    gdf["min_height"] = np.asarray(min_h).round(2)
    gdf["mean_height"] = np.asarray(mean_h).round(2)
    gdf["max_height"] = np.asarray(max_h).round(2)
    gdf["height_class"] = np.select(
        [
            gdf["tree_height"] < 1,
            gdf["tree_height"] < 5,
            gdf["tree_height"] < 10,
        ],
        ["0p5-1", "1-5", "5-10"],
        default="10pluss",
    )
    return gdf


def _add_centerline_distances(top_pts_gdf, centerlines_path, ndsm_crs):
    """Add stretch_id, dist_2D, dist_3D via nearest-route join."""
    cl_layers = input_layers(centerlines_path)
    cl_gdfs = [as_gdf(centerlines_path, layer=ly) for ly in cl_layers]
    cl_gdf = (
        gpd.GeoDataFrame(pd.concat(cl_gdfs, ignore_index=True))
        if len(cl_gdfs) > 1
        else cl_gdfs[0]
    )
    cl_gdf = cl_gdf.to_crs(ndsm_crs)

    joined = gpd.sjoin_nearest(
        top_pts_gdf,
        cl_gdf[["stretch_id", "geometry"]],
        how="left",
        distance_col="dist_2D",
    )
    # sjoin_nearest may create duplicates if equidistant; keep first match
    joined = joined[~joined.index.duplicated(keep="first")]
    top_pts_gdf = top_pts_gdf.copy()
    top_pts_gdf["stretch_id"] = joined["stretch_id"].values
    top_pts_gdf["dist_2D"] = np.array(joined["dist_2D"].values).round(2)
    top_pts_gdf["dist_3D"] = np.sqrt(
        top_pts_gdf["dist_2D"] ** 2 + top_pts_gdf["tree_height"] ** 2
    ).round(2)
    logger.info("  Added stretch_id, dist_2D, dist_3D columns")
    return top_pts_gdf


# ── Whole-array path (unchanged behaviour; used when window <= max_pixels) ───


def _empty_segments(ndsm_crs):
    """Canopy polygons frame with no rows but the columns a full one has."""
    segments_gdf = gpd.GeoDataFrame(
        {"id": pd.Series(dtype=int)}, geometry=[], crs=ndsm_crs.to_epsg()
    )
    segments_gdf["ndsm_min"] = pd.Series(dtype=float)
    segments_gdf["ndsm_max"] = pd.Series(dtype=float)
    return segments_gdf


def _process_whole(
    ndsm_src,
    mask_gdf,
    rivers_gdf,
    water_suppression_height,
    dtm_path,
    centerlines_path,
    output_raster_path,
    write_heights_raster,
    ndsm_crs,
    singletree_kwargs,
):
    """Original whole-array path: load the buffer-cropped NDSM, run the kernel
    once, vectorise. Returns (top_pts_gdf, segments_gdf)."""
    if mask_gdf is not None:
        ndsm, ndsm_transform = rio_mask(
            ndsm_src, mask_gdf.geometry, crop=True, nodata=0
        )
        meta = ndsm_src.meta.copy()
        meta.update(
            {
                "transform": ndsm_transform,
                "width": ndsm.shape[2],
                "height": ndsm.shape[1],
                "nodata": np.nan,
            }
        )
    else:
        ndsm = ndsm_src.read()
        ndsm_transform = ndsm_src.transform
        meta = ndsm_src.meta.copy()
        meta["nodata"] = np.nan

    if rivers_gdf is not None:
        logger.info(
            f"Suppressing water-surface artifacts (nDSM < "
            f"{water_suppression_height} m inside river polygons)..."
        )
        n_water, n_suppressed = _water_suppress(
            ndsm[0],
            list(rivers_gdf.geometry),
            ndsm_transform,
            water_suppression_height,
        )
        logger.info(
            f"  Suppressed {n_suppressed} water-artifact pixels "
            f"({n_water} total water pixels in the AOI)"
        )

    logger.info("Running watershed segmentation...")
    result_raster, ws_labels, ndsm_2d = singletree(ndsm, **singletree_kwargs)

    if write_heights_raster and not Path(output_raster_path).exists():
        Path(output_raster_path).parent.mkdir(parents=True, exist_ok=True)
        meta.update(
            {"dtype": result_raster.dtype, "count": 1, "BIGTIFF": "IF_SAFER"}
        )
        with rio.open(output_raster_path, "w", **meta) as dst:
            dst.write(result_raster)
        logger.info(f"Saved raster: {output_raster_path}")

    logger.info("Extracting top points...")
    unique_labels = np.unique(ws_labels)
    unique_labels = unique_labels[unique_labels > 0]

    if unique_labels.size == 0:
        # Bare ground (Finnmark coast, skerries): scipy hands back no apex
        # positions, and _segment_stats cannot index an empty (1, 0) array.
        logger.info("No canopies in the buffer")
        top_pts_gdf = _build_top_points_gdf(
            [], *([np.array([])] * 6), dtm_path, ndsm_crs
        )
        segments_gdf = _empty_segments(ndsm_crs).merge(
            top_pts_gdf.drop(columns=["geometry"]), on="id", how="left"
        )
        return top_pts_gdf, segments_gdf

    stats = _segment_stats(ndsm_2d, ws_labels, unique_labels)
    top_rows, top_cols = stats["apex_row"], stats["apex_col"]
    top_z_vals = ndsm_2d[top_rows, top_cols]
    top_x, top_y = ndsm_transform * (top_cols + 0.5, top_rows + 0.5)

    top_pts_gdf = _build_top_points_gdf(
        unique_labels.tolist(),
        top_x,
        top_y,
        top_z_vals,
        stats["min_height"],
        stats["mean_height"],
        stats["max_height"],
        dtm_path,
        ndsm_crs,
    )
    if centerlines_path is not None:
        logger.info("Calculating distances to river centerlines...")
        top_pts_gdf = _add_centerline_distances(
            top_pts_gdf, centerlines_path, ndsm_crs
        )

    logger.info("Extracting segment polygons...")
    shapes = rio_shapes(ws_labels, connectivity=4, transform=ndsm_transform)
    records = [
        {"geometry": geometry, "properties": {"id": int(value)}}
        for (geometry, value) in shapes
        if value > 0
    ]
    segments_gdf = gpd.GeoDataFrame.from_features(
        records, crs=ndsm_crs.to_epsg()
    )
    segments_gdf = segments_gdf.sort_values(by="id").reset_index(drop=True)
    logger.info(
        f"{len(segments_gdf)} segments from {len(unique_labels)} labels"
    )

    # join top point attributes, then ndsm zonal stats (original order)
    segments_gdf = segments_gdf.merge(
        top_pts_gdf.drop(columns=["geometry"]), on="id", how="left"
    )
    zonal_stats_result = zonal_stats(
        segments_gdf,
        ndsm_2d,
        stats=["min", "max"],
        affine=ndsm_transform,
        nodata=-1,
        geojson_out=True,
    )
    stats_gdf = gpd.GeoDataFrame.from_features(zonal_stats_result)
    segments_gdf[["min", "max"]] = (
        stats_gdf[["min", "max"]].astype(np.float64).round(2)
    )
    segments_gdf.rename(
        columns={"min": "ndsm_min", "max": "ndsm_max"}, inplace=True
    )

    return top_pts_gdf, segments_gdf


# ── Tiled path (memory-bounded; used when window > max_pixels) ───────────────


def _process_tiled(
    ndsm_src,
    mask_gdf,
    buf_window,
    rivers_gdf,
    water_suppression_height,
    dtm_path,
    centerlines_path,
    output_raster_path,
    write_heights_raster,
    tile_size,
    tile_halo,
    ndsm_crs,
    singletree_kwargs,
):
    """Memory-bounded producer: walk the buffer bbox in core tiles (+halo),
    run the same kernel per tile, keep each tree once (apex in core), and merge.
    Returns (top_pts_gdf, segments_gdf)."""
    buf_transform = ndsm_src.window_transform(buf_window)
    buf_col0, buf_row0 = int(buf_window.col_off), int(buf_window.row_off)
    buf_w, buf_h = int(buf_window.width), int(buf_window.height)

    mask_geoms = (
        np.array(list(mask_gdf.geometry), dtype=object)
        if mask_gdf is not None
        else None
    )
    mask_tree = (
        STRtree(mask_geoms)
        if mask_geoms is not None and len(mask_geoms)
        else None
    )
    river_geoms = (
        np.array(list(rivers_gdf.geometry), dtype=object)
        if rivers_gdf is not None
        else None
    )
    river_tree = (
        STRtree(river_geoms)
        if river_geoms is not None and len(river_geoms)
        else None
    )

    heights_dst = None
    if write_heights_raster:
        meta = ndsm_src.meta.copy()
        meta.update(
            {
                "transform": buf_transform,
                "width": buf_w,
                "height": buf_h,
                "count": 1,
                "dtype": "float32",
                "nodata": np.nan,
                "compress": "lzw",
                "tiled": True,
                "blockxsize": 512,
                "blockysize": 512,
                "BIGTIFF": "YES",
            }
        )
        Path(output_raster_path).parent.mkdir(parents=True, exist_ok=True)
        heights_dst = rio.open(output_raster_path, "w", **meta)

    # accumulators across tiles
    ids_l, topx_l, topy_l, th_l, min_l, mean_l, max_l = ([] for _ in range(7))
    seg_frames = []
    id_offset = 0
    n_tiles = n_kept_tiles = 0

    try:
        for r0 in range(0, buf_h, tile_size):
            core_h = min(tile_size, buf_h - r0)
            for c0 in range(0, buf_w, tile_size):
                core_w = min(tile_size, buf_w - c0)
                n_tiles += 1

                core_win = Window(buf_col0 + c0, buf_row0 + r0, core_w, core_h)
                core_box = box(*window_bounds(core_win, ndsm_src.transform))
                # Skip tiles whose core doesn't touch the buffer (thin ribbon →
                # most tiles skipped, so we never even read them).
                if (
                    mask_tree is not None
                    and len(mask_tree.query(core_box)) == 0
                ):
                    continue
                n_kept_tiles += 1

                read_win = Window(
                    buf_col0 + c0 - tile_halo,
                    buf_row0 + r0 - tile_halo,
                    core_w + 2 * tile_halo,
                    core_h + 2 * tile_halo,
                )
                # boundless read: halo beyond the raster edge becomes NaN (nodata)
                ndsm_tile = ndsm_src.read(
                    1, window=read_win, boundless=True, fill_value=np.nan
                ).astype(np.float32)
                read_transform = ndsm_src.window_transform(read_win)
                read_box = box(*window_bounds(read_win, ndsm_src.transform))

                # Burn the buffer: outside-buffer pixels -> 0, mirroring the
                # whole path's rio_mask(nodata=0) so only in-buffer vegetation
                # is segmented.
                if mask_tree is not None:
                    idx = mask_tree.query(read_box)
                    inside = _rasterize_inside(
                        mask_geoms[idx], ndsm_tile.shape, read_transform
                    )
                    ndsm_tile[~inside] = 0.0

                # Water suppression within the read window (only nearby rivers).
                if river_tree is not None:
                    ridx = river_tree.query(read_box)
                    if len(ridx):
                        _water_suppress(
                            ndsm_tile,
                            river_geoms[ridx],
                            read_transform,
                            water_suppression_height,
                        )

                result_raster, lbl, ndsm2d = singletree(
                    ndsm_tile[np.newaxis, ...], **singletree_kwargs
                )

                labels_all = np.unique(lbl)
                labels_all = labels_all[labels_all > 0]
                if labels_all.size == 0:
                    continue

                stats = _segment_stats(ndsm2d, lbl, labels_all)
                ar, ac = stats["apex_row"], stats["apex_col"]
                # Keep segments whose apex falls in the CORE (half-open), so each
                # tree is emitted by exactly one tile.
                in_core = (
                    (ar >= tile_halo)
                    & (ar < tile_halo + core_h)
                    & (ac >= tile_halo)
                    & (ac < tile_halo + core_w)
                )
                keep = np.nonzero(in_core)[0]
                if keep.size == 0:
                    # still write the heights core block (may contain segments
                    # owned by neighbours) before moving on
                    if heights_dst is not None:
                        _write_heights_block(
                            heights_dst,
                            result_raster,
                            tile_halo,
                            core_h,
                            core_w,
                            c0,
                            r0,
                        )
                    continue

                keep_labels = labels_all[keep]
                new_ids = np.arange(id_offset + 1, id_offset + keep.size + 1)
                id_offset += keep.size

                # apex (tile-local) -> buffer-window pixel -> map coords
                buf_r = r0 + (ar[keep] - tile_halo)
                buf_c = c0 + (ac[keep] - tile_halo)
                top_x, top_y = buf_transform * (buf_c + 0.5, buf_r + 0.5)

                ids_l.append(new_ids)
                topx_l.append(np.asarray(top_x))
                topy_l.append(np.asarray(top_y))
                th_l.append(ndsm2d[ar[keep], ac[keep]])
                min_l.append(stats["min_height"][keep])
                mean_l.append(stats["mean_height"][keep])
                max_l.append(stats["max_height"][keep])

                # polygons for kept labels (already in map CRS via read_transform)
                label_to_id = dict(
                    zip(keep_labels.tolist(), new_ids.tolist(), strict=True)
                )
                recs = [
                    {
                        "geometry": geom,
                        "properties": {"id": label_to_id[int(v)]},
                    }
                    for geom, v in rio_shapes(
                        lbl, connectivity=4, transform=read_transform
                    )
                    if int(v) in label_to_id
                ]
                if recs:
                    seg = gpd.GeoDataFrame.from_features(
                        recs, crs=ndsm_crs.to_epsg()
                    )
                    minmax = pd.DataFrame(
                        {
                            "id": new_ids,
                            "ndsm_min": stats["min_height"][keep].round(2),
                            "ndsm_max": stats["max_height"][keep].round(2),
                        }
                    )
                    seg_frames.append(seg.merge(minmax, on="id", how="left"))

                if heights_dst is not None:
                    _write_heights_block(
                        heights_dst,
                        result_raster,
                        tile_halo,
                        core_h,
                        core_w,
                        c0,
                        r0,
                    )
    finally:
        if heights_dst is not None:
            heights_dst.close()

    logger.info(
        f"Tiled singletree: {n_kept_tiles}/{n_tiles} tiles touched the buffer; "
        f"{id_offset} trees detected"
    )

    ids = np.concatenate(ids_l) if ids_l else np.array([], dtype=int)
    top_pts_gdf = _build_top_points_gdf(
        ids.tolist(),
        np.concatenate(topx_l) if topx_l else np.array([]),
        np.concatenate(topy_l) if topy_l else np.array([]),
        np.concatenate(th_l) if th_l else np.array([]),
        np.concatenate(min_l) if min_l else np.array([]),
        np.concatenate(mean_l) if mean_l else np.array([]),
        np.concatenate(max_l) if max_l else np.array([]),
        dtm_path,
        ndsm_crs,
    )
    if centerlines_path and len(top_pts_gdf):
        logger.info("Calculating distances to river centerlines...")
        top_pts_gdf = _add_centerline_distances(
            top_pts_gdf, centerlines_path, ndsm_crs
        )

    if seg_frames:
        segments_gdf = gpd.GeoDataFrame(
            pd.concat(seg_frames, ignore_index=True), crs=ndsm_crs.to_epsg()
        )
    else:
        segments_gdf = _empty_segments(ndsm_crs)

    segments_gdf = segments_gdf.merge(
        top_pts_gdf.drop(columns=["geometry"]), on="id", how="left"
    )
    return top_pts_gdf, segments_gdf


def _write_heights_block(
    heights_dst, result_raster, halo, core_h, core_w, c0, r0
):
    """Write a tile's CORE block of the per-segment-max raster to the streamed
    full-AOI output (window is relative to the buffer-window output grid)."""
    block = result_raster[0, halo : halo + core_h, halo : halo + core_w].astype(
        "float32"
    )
    heights_dst.write(block, 1, window=Window(c0, r0, core_w, core_h))


def run_singletree(
    ndsm_path,
    dtm_path,
    output_raster_path,
    top_points_path,
    segments_path,
    clip_mask_path=None,
    centerlines_path=None,
    rivers_path=None,
    water_suppression_height=2.0,
    max_pixels=DEFAULT_MAX_PIXELS,
    tile_size=DEFAULT_TILE_SIZE,
    tile_halo=DEFAULT_TILE_HALO,
    write_heights_raster=None,
    **kwargs,
):
    """Run single tree detection and extract top points + segment polygons.

    For buffer-crop windows with at most ``max_pixels`` pixels the original
    whole-array path runs (byte-identical results). Larger windows are processed
    tile-by-tile (memory-bounded); see ``_process_tiled``.

    Args:
        ndsm_path: Path to NDSM raster
        dtm_path: Path to DTM raster (for bottom_z / absolute heights)
        output_raster_path: Path for output heights raster
        top_points_path: Path for output tree top points gpkg
        segments_path: Path for output segment polygons gpkg
        clip_mask_path: Optional path to vector mask (e.g. buffer polygons)
        centerlines_path: Optional path to river centerlines gpkg (for distance columns)
        max_pixels: Whole-array path is used at or below this pixel count.
        tile_size: Core tile edge (px) for the tiled path.
        tile_halo: Halo margin (px) for the tiled path.
        write_heights_raster: None -> write in whole path, skip in tiled path
            (diagnostic only, no consumer); True/False forces it.
        **kwargs: Passed to singletree()
    """
    # output_raster_path is optional (diagnostic only, no consumer).
    wanted = [
        p for p in (output_raster_path, top_points_path, segments_path) if p
    ]
    all_exist = bool(wanted) and all(Path(p).exists() for p in wanted)
    if all_exist:
        print("All singletree output files already exist. Skipping processing.")
        return

    logger.info("Loading canopy height model...")
    # open_uncached applies the LOCAL_CS -> EPSG:25833 fix the national height
    # models need; a plain rio.open would report a CRS mismatch against every
    # vector layer and silently reproject.
    with rio.open(ndsm_path) as ndsm_src:
        ndsm_crs = effective_crs(ndsm_src)

        if clip_mask_path is not None:
            mask_gdf = gpd.read_file(clip_mask_path).to_crs(ndsm_crs)
            buf_window = geometry_window(ndsm_src, mask_gdf.geometry)
            n_pixels = int(buf_window.width) * int(buf_window.height)
        else:
            mask_gdf = None
            buf_window = Window(0, 0, ndsm_src.width, ndsm_src.height)
            n_pixels = int(ndsm_src.width) * int(ndsm_src.height)

        rivers_gdf = (
            _load_rivers(rivers_path, ndsm_crs) if rivers_path else None
        )

        use_tiled = n_pixels > max_pixels
        if write_heights_raster is None:
            write_heights_raster = not use_tiled
        logger.info(
            f"NDSM buffer window: {int(buf_window.width)} x "
            f"{int(buf_window.height)} = {n_pixels:,} px; "
            f"path: {'TILED' if use_tiled else 'WHOLE'}"
        )

        if use_tiled:
            top_pts_gdf, segments_gdf = _process_tiled(
                ndsm_src,
                mask_gdf,
                buf_window,
                rivers_gdf,
                water_suppression_height,
                dtm_path,
                centerlines_path,
                output_raster_path,
                write_heights_raster,
                tile_size,
                tile_halo,
                ndsm_crs,
                kwargs,
            )
        else:
            top_pts_gdf, segments_gdf = _process_whole(
                ndsm_src,
                mask_gdf,
                rivers_gdf,
                water_suppression_height,
                dtm_path,
                centerlines_path,
                output_raster_path,
                write_heights_raster,
                ndsm_crs,
                kwargs,
            )

    # Returned rather than written. The canopy polygons are consumed by the
    # riparian extraction in this same stage — 251 MB per region round-tripped
    # through disk to hand a frame to the next call — and the top points are
    # read by nothing at all. Callers that want them on disk can save them.
    if top_points_path is not None:
        save_gpkg(top_pts_gdf, top_points_path)
    if segments_path is not None:
        save_gpkg(segments_gdf, segments_path)
    logger.info(
        f"Detected {len(segments_gdf)} canopies, {len(top_pts_gdf)} tops"
    )
    return top_pts_gdf, segments_gdf
