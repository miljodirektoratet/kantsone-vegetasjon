r"""
Calculate per-side terrain inclination of riparian zones from the DTM.

For each centerline segment and each side (left/right), computes the
mean and median *signed* terrain slope inside the per-side per-segment
buffer polygon. Slope is defined as the DTM gradient projected onto the
perpendicular-to-river direction pointing AWAY from the river:
    + slope rises away from the river  (typical bank)
    - slope drops away from the river  (depression / floodplain dip)

Inputs:
  - DTM raster (bare-earth elevation, projected, e.g. EPSG:25833)
  - Centerline segments GeoPackage (provides per-segment LineString
    geometry, used to derive each segment's perpendicular-to-river
    direction).
  - Per-side per-segment riparian polygons, produced as verification
    output of extract_riparian_zones_attributes.py:
      riparian_zones_outer_polygons_<min>-<max>m_left.gpkg
      riparian_zones_outer_polygons_<min>-<max>m_right.gpkg

Output (8 columns per buffer range, where <P> = K_L_B<min>-<max> or
K_R_B<min>-<max>):
  <P>_Slope_mean_deg, <P>_Slope_mean_pct,
  <P>_Slope_median_deg, <P>_Slope_median_pct

The columns are written to a dedicated GeoPackage (e.g.
centerline_segments_inclination.gpkg). On first call the file is
created from the centerline segments (carrying river_polygon_id,
river_id, segment_length, geometry) and the slope columns are added.
Subsequent calls for additional buffer ranges append more columns to
the same file; existing columns for the same buffer range are
overwritten.
"""

import numpy as np
import pandas as pd
import rasterio as rio
from rasterio.features import rasterize as rio_rasterize
from rasterio.windows import Window, from_bounds

from ..io import as_gdf

# Whole-raster fast path threshold. At or below this pixel count the DTM and
# its two gradient arrays are precomputed once in memory (float64: 24 bytes/px
# for the three arrays — hence a lower threshold than the float32-based one in
# singletree.py; 200 Mpx ≈ 4.8 GB). Above it, each polygon reads only its own
# DTM window, so memory stays bounded at region scale.
DEFAULT_MAX_PIXELS = 200_000_000


def _fmt_dist(d):
    d_int = int(d)
    if d_int == 0:
        return "0"
    if d_int < 10:
        return f"0{d_int}"
    return str(d_int)


def _compute_world_gradient(dtm_arr, transform):
    """Return (dz_dx, dz_dy) arrays in world units (m elevation per m
    horizontal). dz_dx is positive when elevation increases toward the
    east; dz_dy is positive when elevation increases toward the north.
    Row 0 is the northernmost row, so dz/d(row) needs a sign flip to
    become dz/dy.
    """
    px_x = float(transform.a)
    px_y = float(-transform.e)
    grad_row, grad_col = np.gradient(dtm_arr.astype(np.float64))
    dz_dx = grad_col / px_x
    dz_dy = -grad_row / px_y
    return dz_dx, dz_dy


def _per_side_perp(seg_geom, side):
    """Return the unit vector perpendicular to the segment, pointing
    AWAY from the river on the requested side. Convention matches
    extract_riparian_zones_attributes._split_zone_into_sides:
    perp_left = (-uy, ux) where (ux, uy) is the segment direction
    unit vector (start -> end).
    """
    coords = list(seg_geom.coords)
    start = np.asarray(coords[0], dtype=float)
    end = np.asarray(coords[-1], dtype=float)
    direction = end - start
    length = float(np.linalg.norm(direction))
    if length == 0:
        return None
    unit = direction / length
    perp_left = np.array([-unit[1], unit[0]])
    return perp_left if side == "left" else -perp_left


def _polygon_window(polygon, dtm_src, halo):
    """Integer pixel window covering the polygon's bbox plus ``halo`` px,
    clamped to the raster. Returns None when the polygon misses the raster.
    Clamping (rather than a boundless read) avoids rasterio's per-call VRT
    construction and matches the global gradient's one-sided differences at
    the raster's own edge."""
    win = from_bounds(*polygon.bounds, dtm_src.transform)
    win = win.round_offsets(op="floor").round_lengths(op="ceil")
    col0 = max(int(win.col_off) - halo, 0)
    row0 = max(int(win.row_off) - halo, 0)
    col1 = min(int(win.col_off) + int(win.width) + halo, dtm_src.width)
    row1 = min(int(win.row_off) + int(win.height) + halo, dtm_src.height)
    if col0 >= col1 or row0 >= row1:
        return None
    return Window(col0, row0, col1 - col0, row1 - row0)


def _slope_stats_on_block(dz_dx, dz_dy, block_transform, polygon, perp):
    """Project a gradient block onto ``perp`` over the in-polygon pixels
    (center-inside convention, matching the previous rio_mask default of
    all_touched=False) and return (mean, median) of the finite values."""
    inside = rio_rasterize(
        [(polygon, 1)],
        out_shape=dz_dx.shape,
        transform=block_transform,
        fill=0,
        dtype=np.uint8,
    ).astype(bool)
    if not inside.any():
        return np.nan, np.nan
    signed = dz_dx[inside] * perp[0] + dz_dy[inside] * perp[1]
    valid = signed[np.isfinite(signed)]
    if valid.size == 0:
        return np.nan, np.nan
    return float(valid.mean()), float(np.median(valid))


def _polygon_slope_stats_windowed(polygon, perp, dtm_src, dtm_nodata, halo=2):
    """Windowed signed-slope stats for one polygon (region-scale path).

    Reads only the DTM window covering the polygon (plus a ``halo`` margin so
    the local gradient at every in-polygon pixel uses the same neighbours the
    whole-raster gradient would), computes the world gradient on that block,
    and reduces over the in-polygon pixels. Memory is bounded by polygon size
    instead of the full DTM and gradient arrays.
    """
    if polygon is None or polygon.is_empty:
        return np.nan, np.nan
    win = _polygon_window(polygon, dtm_src, halo)
    if win is None:
        return np.nan, np.nan
    block = dtm_src.read(1, window=win).astype(np.float64)
    if dtm_nodata is not None and not np.isnan(dtm_nodata):
        block[block == dtm_nodata] = np.nan
    block_transform = dtm_src.window_transform(win)
    dz_dx, dz_dy = _compute_world_gradient(block, block_transform)
    return _slope_stats_on_block(dz_dx, dz_dy, block_transform, polygon, perp)


def _polygon_slope_stats_preloaded(polygon, perp, dz_dx, dz_dy, dtm_src):
    """Signed-slope stats from a precomputed global gradient (fast path for
    AOI-sized DTMs): slice the in-memory arrays at the polygon's bbox — no
    per-polygon raster I/O. Numerically identical to the original
    global-gradient implementation."""
    if polygon is None or polygon.is_empty:
        return np.nan, np.nan
    win = _polygon_window(polygon, dtm_src, 0)
    if win is None:
        return np.nan, np.nan
    r0, c0 = int(win.row_off), int(win.col_off)
    rh, cw = int(win.height), int(win.width)
    block_transform = dtm_src.window_transform(win)
    return _slope_stats_on_block(
        dz_dx[r0 : r0 + rh, c0 : c0 + cw],
        dz_dy[r0 : r0 + rh, c0 : c0 + cw],
        block_transform,
        polygon,
        perp,
    )


def extract(
    segments,
    zones,
    band,
    dtm_path,
    max_pixels=DEFAULT_MAX_PIXELS,
):
    """Signed slope columns for one band, both sides, keyed on segment_id."""
    buffer_left_path, buffer_right_path = zones
    buffer_min, buffer_max = band
    segments_path = segments

    range_token = f"B{_fmt_dist(buffer_min)}-{_fmt_dist(buffer_max)}"
    col_prefix_L = f"K_L_{range_token}"
    col_prefix_R = f"K_R_{range_token}"
    # `s` is the specification's name for mean slope in percent; the detailed
    # Slope_* stats are kept for inspection. Every derived formula reads `s`.
    new_cols_L = [
        f"{col_prefix_L}_s",
        f"{col_prefix_L}_Slope_mean_deg",
        f"{col_prefix_L}_Slope_mean_pct",
        f"{col_prefix_L}_Slope_median_deg",
        f"{col_prefix_L}_Slope_median_pct",
    ]
    new_cols_R = [
        f"{col_prefix_R}_s",
        f"{col_prefix_R}_Slope_mean_deg",
        f"{col_prefix_R}_Slope_mean_pct",
        f"{col_prefix_R}_Slope_median_deg",
        f"{col_prefix_R}_Slope_median_pct",
    ]
    new_cols = new_cols_L + new_cols_R

    # ── Load segments + per-side polygons ──────────────────────────────────
    segments_gdf = as_gdf(segments_path).reset_index(drop=True)
    segments_gdf["_seg_idx"] = segments_gdf["segment_id"]
    print(f"  Segments: {len(segments_gdf)}")

    print("Loading left-side polygons ...")
    left_gdf = as_gdf(buffer_left_path).to_crs(segments_gdf.crs)
    print(f"  Polygons: {len(left_gdf)}")
    # A single segment can contribute multiple polygon parts to the
    # verification gpkg (e.g. when its zone-of-influence spans a river
    # boundary). Dissolve so each segment_id maps to one geometry.
    left_gdf = left_gdf.dissolve(by="segment_id").reset_index()
    print(f"  Distinct segments: {len(left_gdf)}")
    print("Loading right-side polygons ...")
    right_gdf = as_gdf(buffer_right_path).to_crs(segments_gdf.crs)
    print(f"  Polygons: {len(right_gdf)}")
    right_gdf = right_gdf.dissolve(by="segment_id").reset_index()
    print(f"  Distinct segments: {len(right_gdf)}")

    # ── Open DTM and pick the gradient strategy ─────────────────────────────
    # AOI-sized DTMs (<= max_pixels): precompute the global gradient once in
    # memory, exactly as the script originally did — per-polygon work is then
    # pure array slicing (fast; thousands of small polygons cost no raster
    # I/O). Region-sized DTMs: fall back to windowed per-polygon reads so
    # memory stays bounded (slower per polygon, but never OOM).
    print(f"Loading DTM from {dtm_path} ...")
    results = {}  # segment_id -> {col_name: value}
    with rio.open(dtm_path) as dtm_src:
        dtm_epsg = dtm_src.crs.to_epsg()
        seg_epsg = segments_gdf.crs.to_epsg()
        if dtm_epsg != seg_epsg:
            raise ValueError(
                f"DTM CRS ({dtm_epsg}) != segments CRS ({seg_epsg}). "
                "Reproject the DTM before running this script."
            )
        dtm_nodata = dtm_src.nodata

        n_pixels = dtm_src.width * dtm_src.height
        if n_pixels <= max_pixels:
            print(
                f"  DTM {dtm_src.width} x {dtm_src.height} = {n_pixels:,} px "
                f"<= {max_pixels:,}: precomputing global gradient (FAST path)"
            )
            dtm_arr = dtm_src.read(1).astype(np.float64)
            if dtm_nodata is not None and not np.isnan(dtm_nodata):
                dtm_arr[dtm_arr == dtm_nodata] = np.nan
            dz_dx, dz_dy = _compute_world_gradient(dtm_arr, dtm_src.transform)
            del dtm_arr  # gradients carry all info the loop needs

            def stats_fn(polygon, perp):
                return _polygon_slope_stats_preloaded(
                    polygon, perp, dz_dx, dz_dy, dtm_src
                )

        else:
            print(
                f"  DTM {dtm_src.width} x {dtm_src.height} = {n_pixels:,} px "
                f"> {max_pixels:,}: windowed per-polygon reads "
                "(memory-bounded path)"
            )

            def stats_fn(polygon, perp):
                return _polygon_slope_stats_windowed(
                    polygon, perp, dtm_src, dtm_nodata
                )

        seg_geom_by_idx = dict(
            zip(segments_gdf["_seg_idx"], segments_gdf.geometry, strict=True)
        )

        print("\nComputing per-segment per-side signed slope stats ...")
        for side, side_gdf, cols in [
            ("left", left_gdf, new_cols_L),
            ("right", right_gdf, new_cols_R),
        ]:
            print(f"  {side}: {len(side_gdf)} polygons")
            n_skipped = 0
            for _, row in side_gdf.iterrows():
                sid = row["segment_id"]
                seg_geom = seg_geom_by_idx.get(sid)
                if seg_geom is None:
                    n_skipped += 1
                    continue
                perp = _per_side_perp(seg_geom, side)
                if perp is None:
                    n_skipped += 1
                    continue
                s_mean, s_med = stats_fn(row.geometry, perp)
                if not np.isfinite(s_mean):
                    n_skipped += 1
                    continue
                # Keyed by name, not position: the column list is shared with
                # the header definition above and a reorder there must not
                # silently write values into the wrong columns.
                prefix = cols[0].rsplit("_", 1)[0]
                mean_pct = round(float(s_mean) * 100.0, 4)
                results.setdefault(sid, {}).update(
                    {
                        f"{prefix}_s": mean_pct,
                        f"{prefix}_Slope_mean_deg": round(
                            float(np.degrees(np.arctan(s_mean))), 4
                        ),
                        f"{prefix}_Slope_mean_pct": mean_pct,
                        f"{prefix}_Slope_median_deg": round(
                            float(np.degrees(np.arctan(s_med))), 4
                        ),
                        f"{prefix}_Slope_median_pct": round(
                            float(s_med) * 100.0, 4
                        ),
                    }
                )
            if n_skipped:
                print(
                    f"    (skipped {n_skipped} polygons with no valid pixels)"
                )

    # ── Build per-row arrays aligned to segments_gdf row order ────────────
    # Map each result (keyed by segment_id) back to the matching row in
    # segments_gdf. Output rows are written in segments_gdf order so the
    # arrays are positionally aligned with the output gdf.
    n_seg = len(segments_gdf)
    arrs = {col: np.full(n_seg, np.nan, dtype=np.float64) for col in new_cols}
    row_by_segment_id = {
        sid: i for i, sid in enumerate(segments_gdf["segment_id"])
    }
    for sid, vals in results.items():
        pos = row_by_segment_id.get(sid)
        if pos is None:
            continue
        for col, v in vals.items():
            arrs[col][pos] = v

    ids = segments_gdf["segment_id"]
    return {c: pd.Series(arrs[c], index=ids.index) for c in new_cols}
