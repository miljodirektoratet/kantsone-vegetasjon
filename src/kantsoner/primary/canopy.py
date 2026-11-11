"""Canopy height statistics over the riparian zone: Vp95 and Vhmed.

These are *pixel* statistics, not crown statistics. `vegetation.py` measures
how much ground the detected crowns cover; this measures how tall the canopy
standing on that ground is. A bank with one big tree and a bank thick with
scrub can share a `Vtot` and differ completely here, which is the layering
("sjiktning") the specification is after.

Only pixels above `CANOPY_STAT_MIN_HEIGHT_M` count, per the workbook's
"pikselbasert (høydeverdier > 1 meter)". A slice with no canopy at all reports
0.0 rather than NaN: `derived/temperatur.py` multiplies Vhmed straight into
`rekk`, and no canopy means no shading, which is a real 0 rather than a gap.
"""

import numpy as np
import pandas as pd
import rasterio as rio
from rasterio.features import geometry_mask
from rasterio.windows import Window, from_bounds

from ..geometry import segment_zones
from ..io import as_gdf
from ..spec import constants as C
from ..spec.attributes import band_token

# Above this, the batch window is read per polygon instead of once. Same
# reasoning as terrain.py: a region-sized canopy model does not fit in memory.
DEFAULT_MAX_PIXELS = 200_000_000


def _window(polygon, src, halo=1):
    """Integer pixel window over the polygon's bbox, clamped to the raster."""
    win = from_bounds(*polygon.bounds, src.transform)
    win = win.round_offsets(op="floor").round_lengths(op="ceil")
    col0 = max(int(win.col_off) - halo, 0)
    row0 = max(int(win.row_off) - halo, 0)
    col1 = min(int(win.col_off) + int(win.width) + halo, src.width)
    row1 = min(int(win.row_off) + int(win.height) + halo, src.height)
    if col0 >= col1 or row0 >= row1:
        return None
    return Window(col0, row0, col1 - col0, row1 - row0)


def _stats(polygon, src, nodata):
    """(p95, median) of canopy height inside the polygon, or None."""
    win = _window(polygon, src)
    if win is None:
        return None
    block = src.read(1, window=win)
    if block.size == 0:
        return None
    inside = ~geometry_mask(
        [polygon],
        out_shape=block.shape,
        transform=src.window_transform(win),
        invert=False,
    )
    vals = block[inside]
    if nodata is not None:
        vals = vals[vals != nodata]
    vals = vals[np.isfinite(vals) & (vals > C.CANOPY_STAT_MIN_HEIGHT_M)]
    if vals.size == 0:
        return None
    return float(np.percentile(vals, 95)), float(np.median(vals))


def extract(segments, zones, band, ring, chm_path):
    """Vp95 and Vhmed for one band, both sides, keyed on segment_id."""
    ids = as_gdf(segments)["segment_id"]
    out = {}
    with rio.open(chm_path) as src:
        nodata = src.nodata
        for side, cells in (("L", zones[0]), ("R", zones[1])):
            prefix = f"K_{side}_{band_token(band)}"
            sl = segment_zones.cut(cells, ring)
            p95 = pd.Series(0.0, index=ids.index)
            med = pd.Series(0.0, index=ids.index)
            if not sl.empty:
                pos = {s: i for i, s in enumerate(ids)}
                for sid, geom in zip(
                    sl["segment_id"], sl.geometry, strict=True
                ):
                    i = pos.get(sid)
                    if i is None or geom.is_empty:
                        continue
                    st = _stats(geom, src, nodata)
                    if st is None:
                        continue
                    # A segment can own several disjoint cells; keep the
                    # tallest reading rather than letting the last one win.
                    p95.iloc[i] = max(p95.iloc[i], st[0])
                    med.iloc[i] = max(med.iloc[i], st[1])
            out[f"{prefix}_Vp95"] = p95.round(2)
            out[f"{prefix}_Vhmed"] = med.round(2)
    return out
