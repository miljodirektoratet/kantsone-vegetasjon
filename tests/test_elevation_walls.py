"""Terrain in the routing cost: DTM channels priced as water, contour walls.

Synthetic in-memory DTMs. The route has to cross a strip of dry ground
between two pools; a gully across the strip, off the direct line, gathers
enough of the window to count as a channel, and the router takes it.
"""

import numpy as np
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from shapely.geometry import LineString, Polygon, box

from kantsoner.geometry import elvis, flowacc


def _dtm(minx, maxy, w, h, surface):
    """An in-memory 1 m float32 DTM; `surface(x, y)` gives the elevation."""
    xs = minx + np.arange(w) + 0.5
    ys = maxy - np.arange(h) - 0.5
    z = surface(*np.meshgrid(xs, ys)).astype("float32")
    mem = MemoryFile()
    with mem.open(
        driver="GTiff",
        width=w,
        height=h,
        count=1,
        dtype="float32",
        nodata=-100.0,
        crs="EPSG:25833",
        transform=from_origin(minx, maxy, 1.0, 1.0),
    ) as dst:
        dst.write(z, 1)
    return mem.open()


def _window(line):
    minx, miny, maxx, maxy = line.buffer(elvis.CORRIDOR_M).bounds
    minx, maxy = np.floor(minx), np.ceil(maxy)
    return minx, maxy, int(np.ceil(maxx) - minx), int(maxy - np.floor(miny))


def test_land_crossing_follows_the_dtm_gully_with_a_dtm():
    # Two pools joined only by a 60 m strip of land (x 170..230). The ground
    # is a broad V-valley with its floor at y = 110 and a gentle fall
    # eastwards, so the window drains into a gully along y = 110 -- off the
    # direct line between the anchors at y = 150, but inside the
    # FLOWACC_NEAR_LINE_M margin where DTM channels count.
    pools = np.array(
        [
            Polygon([(0, 0), (170, 0), (170, 300), (0, 300)]),
            Polygon([(230, 0), (400, 0), (400, 300), (230, 300)]),
        ],
        dtype=object,
    )
    line = LineString([(100, 150), (300, 150)])
    minx, maxy, w, h = _window(line)

    def surface(x, y):
        return 0.05 * np.abs(y - 110) - 0.01 * x + 50.0

    dtm = _dtm(minx, maxy, w, h, surface)
    straight = elvis._route(line, pools)
    gully = elvis._route(line, pools, dtm=dtm)
    assert straight is not None and gully is not None
    direct = box(170, 135, 230, 180)
    assert straight.intersects(direct), (
        "without a DTM the direct crossing is cheapest"
    )
    assert not gully.intersects(direct), "with a DTM the gully is cheaper land"
    assert gully.intersects(box(170, 95, 230, 125)), (
        "the crossing runs down the gully"
    )


def test_flow_channels_pick_out_the_gully_and_nothing_on_a_plane():
    minx, maxy, w, h = 0.0, 400.0, 400, 400
    valley = _dtm(
        minx, maxy, w, h, lambda x, y: 0.05 * np.abs(y - 200) - 0.01 * x + 50.0
    )
    chan = elvis._flow_channels(valley, minx, maxy - h, minx + w, maxy, h, w)
    assert chan.shape == (h, w)
    rows = np.flatnonzero(chan.any(axis=1))
    assert len(rows) and abs(rows.mean() - 200) < 30, (
        "channel pixels sit on the valley floor"
    )
    assert chan.mean() < 0.1, "and nowhere else"
    plane = _dtm(minx, maxy, w, h, lambda x, y: 0.02 * x + 0.01 * y)
    assert not elvis._flow_channels(
        plane, minx, maxy - h, minx + w, maxy, h, w
    ).any()


def test_accumulation_counts_the_whole_window_at_the_outlet():
    # A plane tilted toward one corner: every cell drains to it.
    z = (np.arange(50)[:, None] + np.arange(50)[None, :]).astype("float64")
    acc = flowacc.accumulation(z, res=1.0)
    assert acc.max() == 50 * 50
    assert np.unravel_index(acc.argmax(), acc.shape) == (0, 0)


def test_walls_mark_level_changes_only():
    minx, maxy, w, h = 0.0, 100.0, 100, 100
    dtm = _dtm(minx, maxy, w, h, lambda x, y: x / 10.0)
    zc = elvis._coarse_elevation(dtm, minx, maxy - h, minx + w, maxy, h, w)
    walls = elvis._elevation_walls(zc, h, w)
    col = walls.any(axis=0)
    assert np.array_equal(walls, np.broadcast_to(col, walls.shape))
    assert 0 < col.sum() < w
    flat = _dtm(minx, maxy, w, h, lambda x, y: np.full_like(x, 5.0))
    zf = elvis._coarse_elevation(flat, minx, maxy - h, minx + w, maxy, h, w)
    assert not elvis._elevation_walls(zf, h, w).any()
