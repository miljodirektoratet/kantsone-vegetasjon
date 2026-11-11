"""River network topology — the prerequisite for every L4 attribute.

Nothing in the previous codebase built this, which is why `Oppstuving`,
`Fragmentering` and the whole konnektivitet / prioritet block were missing
from both scoring systems.

The graph is built over the 10 m segments: nodes are segment endpoints
(snapped to a tolerance), edges are the segments themselves. Flow direction
is the digitising order (see `build`).
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field

import geopandas as gpd
import numpy as np
import pandas as pd
from loguru import logger

SNAP_TOLERANCE_M = 1.0


def _node_key(x: float, y: float, tol: float) -> tuple[int, int]:
    return (int(round(x / tol)), int(round(y / tol)))


@dataclass
class RiverGraph:
    """Directed segment graph. Edges point downstream."""

    segment_ids: np.ndarray
    downstream: dict[int, list[int]] = field(default_factory=dict)
    upstream: dict[int, list[int]] = field(default_factory=dict)
    length: dict[int, float] = field(default_factory=dict)

    def walk_upstream(self, start: int, max_distance: float | None = None):
        """Yield (segment_id, distance_from_start) breadth-first upstream."""
        seen = {start}
        queue = deque([(start, 0.0)])
        while queue:
            sid, dist = queue.popleft()
            yield sid, dist
            for nxt in self.upstream.get(sid, ()):
                if nxt in seen:
                    continue
                d = dist + self.length.get(nxt, 0.0)
                if max_distance is not None and d > max_distance:
                    continue
                seen.add(nxt)
                queue.append((nxt, d))

    def walk_downstream(self, start: int, max_distance: float | None = None):
        seen = {start}
        queue = deque([(start, 0.0)])
        while queue:
            sid, dist = queue.popleft()
            yield sid, dist
            for nxt in self.downstream.get(sid, ()):
                if nxt in seen:
                    continue
                d = dist + self.length.get(nxt, 0.0)
                if max_distance is not None and d > max_distance:
                    continue
                seen.add(nxt)
                queue.append((nxt, d))


def build(
    segments: gpd.GeoDataFrame,
    *,
    tolerance: float = SNAP_TOLERANCE_M,
    id_col: str = "segment_id",
) -> RiverGraph:
    """Build the directed graph from 10 m segments.

    A segment flows in its own coordinate order. That is
    the right default for elvenett: NVE digitises the network downstream —
    measured on Hurumlandet, 96 % of ElvBekk stretches drop from start to
    end, 0.7 % rise, the rest are flat lake midlines — whereas guessing
    direction per 10 m from the DTM broke a third of the chain (9,796
    directed edges for 14,534 segments) and every network walk died within
    a few segments.
    """

    starts: dict[tuple[int, int], list[int]] = defaultdict(list)
    ends: dict[tuple[int, int], list[int]] = defaultdict(list)
    length: dict[int, float] = {}

    for sid, geom in zip(segments[id_col], segments.geometry, strict=False):
        coords = list(geom.coords)
        a, b = coords[0], coords[-1]
        starts[_node_key(*a, tolerance)].append(sid)
        ends[_node_key(*b, tolerance)].append(sid)
        length[sid] = float(geom.length)

    downstream: dict[int, list[int]] = defaultdict(list)
    upstream: dict[int, list[int]] = defaultdict(list)
    for node, incoming in ends.items():
        for outgoing in starts.get(node, ()):
            for sid in incoming:
                if sid == outgoing:
                    continue
                downstream[sid].append(outgoing)
                upstream[outgoing].append(sid)

    graph = RiverGraph(
        segment_ids=segments[id_col].to_numpy(),
        downstream=dict(downstream),
        upstream=dict(upstream),
        length=length,
    )
    logger.info(
        f"Graph: {len(graph.segment_ids)} segments, "
        f"{sum(len(v) for v in downstream.values())} directed edges"
    )
    return graph


def chains(
    segments: gpd.GeoDataFrame, graph: RiverGraph, river_col: str = "elvID"
) -> pd.DataFrame:
    """Order the segments along the flow, for block aggregation.

    Returns a frame indexed like ``segments`` with ``chain_id`` and
    ``chain_pos``. A chain is a run of stretches of one river where each
    flows *uniquely* into the next — the last segment of A has exactly one
    downstream neighbour, it is the first segment of B, and B has no other
    upstream feeder within the river. Chains break at every junction and
    divergence, so a block can only ever contain consecutive water.

    ``river_col`` decides what "one river" means -- delivery/aggregate.py
    groups blocks by the Vann-Nett vannforekomst, which can span several
    elvID reaches or split a single elvID, so stages.py passes it a column
    already resolved to vannforekomstid (falling back to elvID where a
    segment has none) rather than the default.

    Needed because the 10 m ids sort along a stretch but not along the
    river: elvenett strekningLnr is not spatial order, and blocks cut from
    a lexical sort spliced the end of one stretch onto the start of a
    sibling 370 m away.
    """
    sid = segments["segment_id"].to_numpy()
    stretch = segments["stretch_id"].to_numpy()
    river = segments[river_col].astype(str).to_numpy()
    order = np.argsort(sid, kind="stable")
    by_stretch: dict[str, list[int]] = defaultdict(list)
    for i in order:
        by_stretch[stretch[i]].append(i)
    first = {st: rows[0] for st, rows in by_stretch.items()}
    last = {st: rows[-1] for st, rows in by_stretch.items()}
    row_of = {sid[i]: i for i in range(len(sid))}

    # A -> B when A's last segment flows only into B's first, B's first is
    # fed only by A, and both belong to the same river.
    succ: dict[str, str] = {}
    pred: dict[str, str] = {}
    for st, i_last in last.items():
        down = graph.downstream.get(sid[i_last], [])
        if len(down) != 1:
            continue
        j = row_of.get(down[0])
        if j is None:
            continue
        nxt = stretch[j]
        if nxt == st or river[j] != river[i_last] or first[nxt] != j:
            continue
        if len(graph.upstream.get(sid[j], [])) != 1:
            continue
        succ[st] = nxt
        pred[nxt] = st

    chain_id = np.zeros(len(sid), dtype=np.int64)
    chain_pos = np.zeros(len(sid), dtype=np.int64)
    cid = 0
    seen: set[str] = set()
    for st in sorted(by_stretch):
        if st in pred or st in seen:
            continue
        cid += 1
        pos = 0
        cur: str | None = st
        while cur is not None and cur not in seen:
            seen.add(cur)
            for i in by_stretch[cur]:
                chain_id[i] = cid
                chain_pos[i] = pos
                pos += 1
            cur = succ.get(cur)
    # Anything left is a cycle; each becomes its own chain.
    for st in sorted(by_stretch):
        if st in seen:
            continue
        cid += 1
        pos = 0
        cur = st
        while cur is not None and cur not in seen:
            seen.add(cur)
            for i in by_stretch[cur]:
                chain_id[i] = cid
                chain_pos[i] = pos
                pos += 1
            cur = succ.get(cur)
    return pd.DataFrame(
        {"chain_id": chain_id, "chain_pos": chain_pos}, index=segments.index
    )
