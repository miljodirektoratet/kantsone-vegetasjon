"""D8 flow accumulation on a small elevation window.

Enough to say where the DTM thinks water runs inside one stretch's routing
window: depressions filled by priority flood (Barnes et al. 2014), steepest-
descent receivers, upslope area summed in one pass down the filled surface.
Pure numpy plus a heap; a corridor window at 8 m is tens of thousands of
cells, ~0.1 s, so there is nothing to gain from a raster library here.

Window-local by construction: a cell's upslope area only counts what lies
inside the window, so values are a floor on the true accumulation. A
threshold applied to them therefore picks out channels that gather their
water nearby -- creeks, gullies, the reach of a river inside the window --
which is exactly what the router needs where the water polygons stop.
"""

from __future__ import annotations

import heapq

import numpy as np

_NEIGHBOURS = (
    (-1, -1),
    (-1, 0),
    (-1, 1),
    (0, -1),
    (0, 1),
    (1, -1),
    (1, 0),
    (1, 1),
)


def fill_depressions(z: np.ndarray, eps: float = 1e-3) -> np.ndarray:
    """The surface with every pit raised until it drains off the edge,
    each raised cell `eps` above the one it drains to, so no flat is left
    for D8 to stall on."""
    h, w = z.shape
    filled = z.astype("float64", copy=True)
    seen = np.zeros(z.shape, dtype=bool)
    heap: list[tuple[float, int, int]] = []
    for i in range(h):
        for j in (0, w - 1):
            if not seen[i, j]:
                heapq.heappush(heap, (float(filled[i, j]), i, j))
                seen[i, j] = True
    for j in range(w):
        for i in (0, h - 1):
            if not seen[i, j]:
                heapq.heappush(heap, (float(filled[i, j]), i, j))
                seen[i, j] = True
    while heap:
        zc, i, j = heapq.heappop(heap)
        for di, dj in _NEIGHBOURS:
            a, b = i + di, j + dj
            if 0 <= a < h and 0 <= b < w and not seen[a, b]:
                seen[a, b] = True
                v = max(float(filled[a, b]), zc + eps)
                filled[a, b] = v
                heapq.heappush(heap, (v, a, b))
    return filled


def accumulation(z: np.ndarray, res: float) -> np.ndarray:
    """Upslope area per cell in m2, D8 on the filled surface."""
    filled = fill_depressions(z)
    h, w = filled.shape
    best = np.full((h, w), -np.inf)
    ri = np.repeat(np.arange(h)[:, None], w, axis=1)
    rj = np.repeat(np.arange(w)[None, :], h, axis=0)
    rec = ri * w + rj
    for di, dj in _NEIGHBOURS:
        s = np.full((h, w), -np.inf)
        src = filled[max(0, -di) : h - max(0, di), max(0, -dj) : w - max(0, dj)]
        dst = filled[max(0, di) : h - max(0, -di), max(0, dj) : w - max(0, -dj)]
        s[max(0, -di) : h - max(0, di), max(0, -dj) : w - max(0, dj)] = (
            src - dst
        ) / np.hypot(di, dj)
        better = s > best
        best[better] = s[better]
        rec[better] = ((ri + di) * w + (rj + dj))[better]
    acc = np.ones(h * w)
    recf, bestf = rec.ravel(), best.ravel()
    for k in np.argsort(-filled.ravel(), kind="stable"):
        if bestf[k] > 0:
            acc[recf[k]] += acc[k]
    return acc.reshape(h, w) * res * res
