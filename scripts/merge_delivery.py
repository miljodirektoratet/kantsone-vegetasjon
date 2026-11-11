#!/usr/bin/env python3
"""Merge a finished run's region deliveries into one GeoPackage per product.

The run is a finished `kantsoner run` (--source, its base_dir: one folder per
region). It is only read, and the aggregation rules are the pipeline's own
code (--code, <repo>/src), called rather than copied.

Output, one file and one layer per product, in -o/--output:

  official_centerline_10m/100m/500m/1000m/region, lake_shore, lake_polygon
      The delivered rows and measurements, plus eUVannomradeId (the vannområde
      holding the largest part of the feature) and region_folder (the region it
      came from; aggregated ids restart in every region, so they are unique only
      together with it). Four things differ from the regional files, below:
      the vannforekomst, Score_total_final, water delivered twice, and the lake
      a shore belongs to.
  official_centerline_river
      The vannforekomst level: one row per vannforekomst in Miljødirektoratet's
      file, from the 10 m segments, by the run's roll-up: additive sums, then
      every ratio and class recomputed from them. Lake shore is left out, as it
      is from the other centerline levels.
  official_centerline_vannomraade
      The same roll-up, one row per eUVannomradeId.
  buffer10
      The 0-10 m zone of the rivers, one polygon per segment and bank that has
      any: zones_left/right cut to the stretch's buffer10 ring, the same cut the
      attributes were measured on. See BUFFER_FIELDS. Each polygon also carries
      the segment's E_S_Kantvegetasjon_fjernet and E_S_Kantvegetasjon_grad, copied
      from its row in the merged official_centerline_10m (SEGMENT_VALUES), so the
      zones can be drawn by condition class without a join.
  lake_buffer10
      The same 0-10 m zone around the lakes, one polygon per shore segment, with
      the same two values from the merged lake_shore. See LAKE_BUFFER_FIELDS.
      lake_shore is the 0-20 m zone these were cut from.
  single_trees
      The crown polygons, see TREE_FIELDS. Not in the default set: run it on
      its own once the export is finished. Mean height is measured here from
      the CHM, because the export did not keep it. A crown two regions found,
      the same apex, is kept once. The outlines are simplified by TREE_SIMPLIFY_M
      to take the pixel stairs out; area and mean height are measured before that.

The vannforekomst comes from the Vann-Nett download of 2026-10-26 (--vannforekomster),
the one the run was made on, so the run's ids are the file's. (The February
export of Miljødirektoratet has the same lines, but water bodies, mostly in the
north, split and renumbered: 316 of its river ids are gone and 903 are new, and
16 lakes are new. Merged against it, 3.1 % of the network got
another id than the run's.) Each segment takes the line it runs on, by the
longest overlap within 1 m; this also settles the segments at a junction, which
the run gave to whichever line came first. vannforekomstid, river_id_vann-nett
and vannforekomstnavn in every layer are the file's. A lake the file has as a
lake keeps its id; one it does not have belongs to the river that runs through it.

Score_total_final follows Sállir tabell 8 and chapter 3.11 of the report: Svært
dårlig below 0.20 is 1, Svært god from 0.80 is 5, and no grade under 50 %
datadekning. A region whose grades run the other way (made by older code) is
turned here.

Water delivered twice: a stretch crossing a region border was run in both
regions, so about 5,000 segments exist twice, and the copies differ. A copy is
the same geometry in two regions, whatever its id, and is kept once, from the
region that owns its vannforekomst (the first part of the id is the region's
number). The other copy goes from every layer, and the rows it was part of are
recomputed from the 10 m segments by the run's own rules: the region rows of the
water involved, and the block rows of that river in that region. region_folder
names the regions holding real length, longest first.

A lake with shore in two regions is not a copy: each region numbered its own
shore from 0, so the same id can be different shore. Both stay, the second
region's numbered on from the first's in lake_shore and lake_buffer10, and the
lake has one row with the sums of both.

Which lake a shore belongs to: the run gave each mapped water surface to one lake,
so a surface holding two lakes of the file had all its shore under the first and
none under the second. Here a surface is shared by the lakes of the file that
really lie in it. One such lake takes all its shore; several share it, each piece
going to the lake nearest through the water. A piece on a narrow arm of another
lake reaching into the surface takes that lake, if the one just chosen is over
50 m away. Shore never leaves its surface, so it cannot cross land. lake_id and
the vannforekomst follow in lake_shore and lake_buffer10, and the lake's inherited
values in lake_shore. The rows of lake_polygon that change are rebuilt from the
run's own sums, and a lake gaining its first shore gets a row with the file's
polygon.
segment_id stays the run's, as on the rivers, so its first part can name another
lake.

Double lines in the lake shore: some 10 m pieces are not shore but a line drawn
through open water where two mapped water surfaces meet, once from each side, the
copies on top of each other. Grunnkart cuts a lake at a kommune border into two
polygons that stop centimetres short of each other (both copies under one run
lake), and the shore rule above brings the copies of two run lakes that met along
a seam under one lake of the file. A piece is a double line when, at 0.2, 0.5 and
0.8 of its length, another piece of the same lake lies within 0.5 m (see
double_lines). Both copies are removed, from lake_shore and lake_buffer10 and,
those with no zone too, from the sums of lake_polygon, which are recomputed from
the shore that remains. They are counted in the log and the README, and listed in
lake_shore_double_lines.csv.

At the end every line of the file is checked against the result, in short pieces:
a piece is present if a delivered segment lies within 1 m of it. A piece that is
not there but lies in a lake is not missing: lakes are delivered by lake_polygon
and lake_shore. Written: vannforekomst_lines_check.csv, one row per vannforekomst,
with missing_km outside lakes and in_lake_km apart; vannforekomst_lines_missing.gpkg,
the pieces that are not there and not in a lake; and lakes_not_delivered.csv, the
lakes of the file with no row in lake_polygon.

Also written: README.md with every file's layer, rows and fields, and how the
vannområde ids were assigned. It is rebuilt from what is in the folder on each
run, so a later single_trees run adds to it.

    uv run python scripts/merge_delivery.py -s <RUN> -o <OUT> --work <WORK> \
        --vannomrader <DATA>/geonorge/natur/vannomraader_miljoedir_25833_20260227.gpkg \
        --vannforekomster <DATA>/vannnett/vannforekomster_20261026.gpkg -w 3
    uv run python scripts/merge_delivery.py ... --products single_trees \
        --chm <DATA>/hoydedata/chm_norway_1m_epsg25833_20260911.tif

Work happens in --work, a local disk, because SQLite and network shares do not
mix; each finished file
is then copied to the output, SHA-256 checked and removed locally. Regions
finished in a stopped run are not redone, so one --work belongs to one set of
inputs and one version of this script, and refuses any other: after a change,
use a new one. A product already in the output is skipped unless --overwrite is
given. -r is for trying regions out, in an output folder of their own.

Memory: each worker is capped by --max-worker-gb, and the whole run stops itself
below --stop-below-gb; rerun to resume. The workers are also the kernel's first
choice if memory runs out, so a merge cannot take down what runs beside it.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import multiprocessing as mp
import os
import pickle
import resource
import shutil
import signal
import sqlite3
import struct
import subprocess
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import quote

# The pipeline's code, which the merge calls for the lake sums and the national
# roll-ups: it has to be the code that made the region files.
CODE = Path(__file__).resolve().parent.parent / "src"

VO_ID = "eUVannomradeId"
FOLDER = "region_folder"
MD_KEY = "md_id"  # the working column for the file's vannforekomst
MD_RIVERS, MD_LAKES = "VannforekomstElv", "VannforekomstInnsjø"
MD_ID, MD_NAME = "identifikasjon_lokalId", "navn"
BLOCKS = {
    "official_centerline_100m": "100m",
    "official_centerline_500m": "500m",
    "official_centerline_1000m": "1000m",
}
DELIVERED = (
    "official_centerline_10m",
    *BLOCKS,
    "official_centerline_region",
    "lake_shore",
    "lake_polygon",
)
GRADED = {**BLOCKS, "official_centerline_region": "region", "lake_polygon": "innsjo"}
RIVER = "official_centerline_river"
VO_PRODUCT = "official_centerline_vannomraade"
ROLLUPS = (RIVER, VO_PRODUCT)
BUFFER = "buffer10"
LAKE_BUFFER = "lake_buffer10"
TREES = "single_trees"
RING = "buffer10"  # the source file and layer of the 0-10 m ring
DEFAULT_PRODUCTS = (*DELIVERED, *ROLLUPS, BUFFER, LAKE_BUFFER)
# Buffer -> the merged layer whose rows hold the segment's condition (SEGMENT_VALUES).
SEGMENT_SOURCE = {BUFFER: "official_centerline_10m", LAKE_BUFFER: "lake_shore"}


def _level_drop(level):
    """Columns of an aggregated level that are not delivered: the yes/no-rule shares
    (the condition is <level>_Kantvegetasjon and _final, from the graded fjernet) and
    the coverage figures that repeat <level>_datadekning."""
    shares = ("Kantvegetasjon_janei", "Kantskog_fjernet", "Nedbygd_kantsone",
              "ikke_mulig_Kantvegetasjon", "Uklassifisert")
    return [f"frac_matched_{level}", *(f"{level}_{s}" for s in shares),
            f"length_uklassifisert_{level}"]


# Columns taken out of each merged product just before it is published (decided with
# the client's delivery, 2026-11-11). The merge itself still uses some of them (the
# copies are removed by segment_id), so they go at the very end. A name a product
# does not have is skipped.
DELIVERY_DROP = {
    # constant (Kilde "Elv", beregningsbredde 10), the old rule's potential, and
    # segment_id, which repeats segment_id_vann-nett
    "official_centerline_10m": ["segment_id", "K_S_Kilde", "K_S_beregningsbredde",
                                "E_L_Potensiell_Kantvegetasjon",
                                "E_R_Potensiell_Kantvegetasjon"],
    **{f"official_centerline_{b}": _level_drop(b) for b in ("100m", "500m", "1000m")},
    "official_centerline_river": _level_drop("river"),
    # a region has no vannforekomstnavn: the column is empty
    "official_centerline_region": [*_level_drop("region"), "vannforekomstnavn"],
    "official_centerline_vannomraade": _level_drop("vannomraade"),
    # lake_eu is empty in Vann-Nett's lakes
    "lake_polygon": [*_level_drop("innsjo"), "lake_eu"],
    BUFFER: ["segment_id", "E_Potensiell_Kantvegetasjon"],
    LAKE_BUFFER: ["segment_id", "E_Potensiell_Kantvegetasjon"],
}
ALL_PRODUCTS = (*DEFAULT_PRODUCTS, TREES)

BANK_VALUES = {
    "area_m2": "area of this polygon",
    "K_B0-10_Vtot": "K_{side}_B0-10_Vtot: vegetation cover in the 0-10 m zone, %",
    "E_Potensiell_Kantvegetasjon": "E_{side}_Potensiell_Kantvegetasjon",
    "E_Kantvegetasjon_A": "E_{side}_Kantvegetasjon_A: class of this bank",
    "E_Kantvegetasjon_A_arsak": "E_{side}_Kantvegetasjon_A_arsak: cause",
}
# The segment's own condition, the same on both of its banks and the same as in its row
# of the merged segment layer: the share of the expected canopy cover that is gone, and
# the class that share gives (kantveg_grad.py). Not made per region: add_segment_condition
# copies them from the merged official_centerline_10m / lake_shore after the buffer is merged.
SEGMENT_VALUES = {
    "E_S_Kantvegetasjon_fjernet": (
        "E_S_Kantvegetasjon_fjernet: share of the expected canopy cover that is gone, "
        "0-100 %, both banks half each; empty when the segment is not assessed"
    ),
    "E_S_Kantvegetasjon_grad": (
        "E_S_Kantvegetasjon_grad: class of that share, Svært god (under 5 %) to "
        "Svært dårlig (80 % or more), or Ikke vurdert, Ikke vurdert (bekkelukking), "
        "Ikke beregnet"
    ),
}
BUFFER_FIELDS = {
    FOLDER: "region the segment came from",
    VO_ID: "the segment's vannområde",
    "segment_id": "the 10 m segment; joins official_centerline_10m",
    "stretch_id": "the ELVIS stretch",
    "side": "L or R bank",
    "segment_id_vann-nett": "from the segment",
    "river_id_vann-nett": "the segment's vannforekomst",
    "vannforekomstid": "the segment's vannforekomst",
    "vannforekomstnavn": "its name in Vann-Nett",
    "elvenavn": "the river's name in ELVIS, from the segment",
    **BANK_VALUES,
}
LAKE_BUFFER_FIELDS = {
    FOLDER: "region the shore segment came from",
    VO_ID: "the shore segment's vannområde",
    "lake_id": "the lake; joins lake_polygon",
    "segment_id": "the shore segment; joins lake_shore",
    "stretch_id": "the run's lake this shore was cut from",
    "side": "L: a shore has one bank",
    "segment_id_vann-nett": "from the shore segment",
    "river_id_vann-nett": "the shore segment's vannforekomst",
    "vannforekomstid": "the shore segment's vannforekomst",
    "vannforekomstnavn": "its name in Vann-Nett",
    "elvenavn": "the lake's name in ELVIS, from the shore segment (lake_shore's elvenavn)",
    **{
        name: text.replace("{side}", "L").replace("this bank", "this shore")
        for name, text in BANK_VALUES.items()
    },
}
TREE_FIELDS = {
    FOLDER: "region the crown came from; tree_id is unique only with it",
    VO_ID: "the vannområde holding the crown's apex",
    "tree_id": "id within the region",
    "stretch_id": "nearest ELVIS stretch of the region the crown came from",
    "max_height": "maks høyde: CHM at the apex, m",
    "mean_height": (
        "middelhøyde: mean CHM over the crown's pixels, m. Detection's own mean "
        "also counted pixels it had zeroed (outside its mask, low in rivers), "
        "so this reads higher on about 1 crown in 10"
    ),
    "crown_area_m2": "kroneareal, m², of the crown before its outline was simplified",
    "height_class": "høydeklasse: 1-5, 5-10 or 10pluss",
}

ADDED = {
    "lake_id": "the lake in Vann-Nett; joins lake_polygon and lake_shore",
    VO_ID: "the vannområde holding the largest part of the feature",
    FOLDER: "the region it came from",
    "vannforekomstid": "the vannforekomst in Vann-Nett",
    "river_id_vann-nett": "the vannforekomst in Vann-Nett",
    "vannforekomstnavn": "its name in Vann-Nett",
}
# Stated, not inherited: regions hold a mix of single and multi parts.
GEOMETRY = {
    "official_centerline_10m": "LINESTRING",
    "official_centerline_100m": "MULTILINESTRING",
    "official_centerline_500m": "MULTILINESTRING",
    "official_centerline_1000m": "MULTILINESTRING",
    "official_centerline_river": "MULTILINESTRING",
    "official_centerline_region": "MULTILINESTRING",
    "lake_shore": "MULTIPOLYGON",
    "lake_polygon": "MULTIPOLYGON",
    BUFFER: "MULTIPOLYGON",
    LAKE_BUFFER: "MULTIPOLYGON",
    TREES: "MULTIPOLYGON",
}

CHUNK = 500_000  # rows read at a time; every sum below is additive
NEW_SHORE_IDS = "new_shore_ids.pkl"  # in --work, next to copies.pkl
LAKE_INDEX = "lakes.pkl"  # in --work: the file's lakes, for the shore rule
SHORE_PLAN = "shore_plan.pkl"  # in --work: what the rule changes, for every worker
TREE_COPIES = "tree_copies.pkl"  # in --work: region -> the apexes it does not keep
# Half the diagonal of a 1 m CHM pixel: the stairs of a vectorised crown go, which
# halves its vertices, and no vertex that is kept moves. The merged file is about
# 270 million crowns.
TREE_SIMPLIFY_M = 0.71
LAKE_KEY = "lake_id"
INHERITED = ("E_S_Barrier", "E_S_Fragmentering")
REACH_M = 20.0  # a lake claims a surface if its overlap holds a circle of this radius
TIP_M = 5.0  # an arm holds a circle of this radius; a piece this near lies on it
FAR_M = 50.0  # only a piece this far from the lake it was given can follow an arm
WATER_CELL_M = 5.0  # raster cell of the walk through the water
WATER_MAX_CELLS = 6_000_000  # very large surfaces get a coarser cell
TREE_CELL_M = 2048  # crowns are measured a CHM window of this size at a time
LINE_TOLERANCE_M = 1.0  # the run's own distance for tagging a segment to a line
# A shore piece is a double line when, at each of these fractions of its length,
# another piece of the same lake lies within DOUBLE_M. See double_lines.
DOUBLE_STATIONS = (0.2, 0.5, 0.8)
DOUBLE_M = 0.5
PIECE_VERTICES = 64  # the file's lines are indexed in pieces of this many
FINE_M = 10.0  # and checked against the result in pieces of at most this length
# A 10 m row the run did not compute (Beregnet = no) reads Uklassifisert in
# E_S_Kantvegetasjon_A (kantsoner.spec.constants.UKLASSIFISERT_KANTVEG; the merge
# does not import the run's code at this level). Aggregation keeps such rows out of
# every mean, so only the 10 m table needs its score and priority emptied.
UNCLASSIFIED = "Uklassifisert"
BLANK_UNCLASSIFIED = ("E_S_Score_total", "K_S_prioritet")
CRS = "EPSG:25833"
VO = None  # the vannområde index, set in each worker
MD = None  # the vannforekomst index, set in each worker
COPIES = {}  # region -> segment ids it does not keep, set in each worker
COPIED_RIVERS = {}  # region -> river ids that lost segments there, likewise
SHORE_COPIES = {}  # region -> shore segment ids it does not keep, likewise
DOUBLES = {}  # region -> shore segment ids that are double lines, likewise
LAKES = None  # the file's lakes, set in each worker of the shore rule
MOVED = {}  # region -> {shore segment id: the lake it moves to}, set in each worker
TARGETS = {}  # lake -> (elvenavn, *INHERITED) for the shore it receives, likewise
LAKE_SUMS = {}  # region -> the lakes whose sums it must give, likewise
TREE_DROPS = {}  # region -> the apexes of crowns it does not keep, likewise


# ── which vannområde ───────────────────────────────────────────────────────
def _polygons(geom):
    import shapely

    parts = shapely.get_parts(geom)
    nested = shapely.get_type_id(parts) >= 6  # multipolygon or collection
    if nested.any():
        parts = shapely.get_parts(
            [p for g in parts for p in shapely.get_parts(g)]
        )
    return [
        p for p in parts if shapely.get_type_id(p) == 3 and not p.is_empty
    ]


def _polygon_part(geom):
    """The polygons of a geometry as one MultiPolygon; None if it holds none."""
    import shapely

    polys = _polygons(geom)
    return shapely.multipolygons(polys) if polys else None


def _subdivide(geom, max_vertices=2000):
    """Halve a polygon until each piece is small enough to test quickly."""
    import shapely

    if shapely.get_num_coordinates(geom) <= max_vertices:
        return [geom]
    x0, y0, x1, y1 = geom.bounds
    if x1 - x0 > y1 - y0:
        xm = (x0 + x1) / 2
        halves = [shapely.box(x0, y0, xm, y1), shapely.box(xm, y0, x1, y1)]
    else:
        ym = (y0 + y1) / 2
        halves = [shapely.box(x0, y0, x1, ym), shapely.box(x0, ym, x1, y1)]
    out = []
    for half in halves:
        for piece in _polygons(shapely.intersection(geom, half)):
            out.extend(_subdivide(piece, max_vertices))
    return out


def build_vo_index(path: Path, cache: Path) -> None:
    import pyogrio
    import shapely

    if cache.exists():
        return
    gdf = pyogrio.read_dataframe(path)
    table = gdf.drop(columns="geometry")
    pieces, owner = [], []
    for i, geom in enumerate(shapely.make_valid(gdf.geometry.values)):
        for poly in _polygons(geom):
            found = _subdivide(poly)
            pieces.extend(found)
            owner.extend([i] * len(found))
    tmp = cache.with_suffix(".tmp")
    with open(tmp, "wb") as f:
        pickle.dump((table, shapely.to_wkb(pieces), owner), f)
    tmp.rename(cache)


class Vannomrader:
    def __init__(self, cache: Path):
        import numpy as np
        import shapely

        with open(cache, "rb") as f:
            self.table, wkb, owner = pickle.load(f)
        self.pieces = shapely.from_wkb(wkb)
        self.owner = np.asarray(owner)
        self.tree = shapely.STRtree(self.pieces)
        self.ids = self.table[VO_ID].to_numpy()

    def assign(self, geoms):
        """eUVannomradeId per geometry, and the counts of how it was decided.

        Inside one vannområde: that one. Across a boundary: the one holding
        the most length (lines) or area (polygons). Touching none, as on the
        outer coast: the nearest.
        """
        import numpy as np
        import pandas as pd
        import shapely

        def size(g):
            area = shapely.area(g)
            return np.where(area > 0, area, shapely.length(g))

        # By part, not by feature: a region row is millions of vertices, and
        # testing it whole against every piece it touches means thousands of
        # copies of it. Only parts that cross a boundary are cut.
        geoms = np.asarray(geoms)
        parts, feature = shapely.get_parts(geoms, return_index=True)
        out = np.full(len(geoms), -1)
        src, piece = self.tree.query(parts, predicate="intersects")
        pairs = pd.DataFrame({"part": src, "vo": self.owner[piece], "p": piece})
        n_vo = pairs.groupby("part")["vo"].nunique()
        crossing = pairs["part"].map(n_vo).to_numpy() > 1

        whole = pairs[~crossing].drop_duplicates("part")
        shares = [whole.assign(size=size(parts[whole["part"].to_numpy()]))]
        if crossing.any():
            cut = pairs[crossing]
            valid = shapely.make_valid(parts[cut["part"].to_numpy()])
            inside = shapely.intersection(
                valid, self.pieces[cut["p"].to_numpy()]
            )
            shares.append(cut.assign(size=size(inside)))
        shares = pd.concat(shares, ignore_index=True)
        shares["feature"] = feature[shares["part"].to_numpy()]
        total = shares.groupby(["feature", "vo"])["size"].sum().reset_index()
        n_found = total.groupby("feature")["vo"].size()
        best = total.sort_values(["feature", "size"]).drop_duplicates(
            "feature", keep="last"
        )
        out[best["feature"].to_numpy()] = best["vo"].to_numpy()
        several = n_found.index[n_found > 1]

        outside = np.nonzero(out < 0)[0]
        max_m = 0.0
        if len(outside):
            i, near = self.tree.query_nearest(geoms[outside], all_matches=False)
            out[outside[i]] = self.owner[near]
            max_m = float(
                shapely.distance(geoms[outside[i]], self.pieces[near]).max()
            )
        counts = {
            "rows": len(geoms),
            "boundary": len(several),
            "outside": len(outside),
            "outside_max_m": round(max_m, 1),
        }
        return self.ids[out], counts


# ── which vannforekomst ────────────────────────────────────────────────────
def build_md_index(path: Path, cache: Path) -> None:
    """Miljødirektoratet's rivers as short pieces of line, ready for a region.

    Whole water bodies run to hundreds of kilometres, so a spatial index on
    them says nothing; pieces of a few hundred metres keep every test local.
    """
    import numpy as np
    import pyogrio
    import shapely

    if cache.exists():
        return
    rivers = pyogrio.read_dataframe(path, layer=MD_RIVERS, columns=[MD_ID, MD_NAME])
    lakes = pyogrio.read_dataframe(
        path, layer=MD_LAKES, columns=[MD_ID, MD_NAME], read_geometry=False
    )
    parts, owner = shapely.get_parts(rivers.geometry.values, return_index=True)
    xy, at = shapely.get_coordinates(parts, return_index=True)
    starts = np.searchsorted(at, np.arange(len(parts)))
    ends = np.append(starts[1:], len(xy))
    pieces, owners = [], []
    for part, (a, b) in enumerate(zip(starts, ends, strict=True)):
        for first in range(a, b - 1, PIECE_VERTICES):
            last = min(first + PIECE_VERTICES + 1, b)
            pieces.append(shapely.linestrings(xy[first:last]))
            owners.append(owner[part])
    pieces = np.array(pieces, dtype=object)
    fine = np.zeros(len(pieces), dtype=np.int64)
    for i in range(0, len(pieces), 100_000):
        cut = shapely.segmentize(pieces[i : i + 100_000], FINE_M)
        fine[i : i + 100_000] = shapely.get_num_coordinates(cut) - 1
    index = {
        "ids": rivers[MD_ID].to_numpy(dtype=object),
        "names": rivers[MD_NAME].to_numpy(dtype=object),
        "lake_names": dict(zip(lakes[MD_ID], lakes[MD_NAME], strict=True)),
        "wkb": shapely.to_wkb(pieces),
        "bounds": shapely.bounds(pieces),
        "owner": np.asarray(owners, dtype=np.int32),
        "offset": np.concatenate([[0], np.cumsum(fine)]),
    }
    tmp = cache.with_suffix(".tmp")
    with open(tmp, "wb") as f:
        pickle.dump(index, f)
    tmp.rename(cache)


class Vannforekomster:
    def __init__(self, cache: Path):
        with open(cache, "rb") as f:
            self.__dict__.update(pickle.load(f))
        self.index = {v: i for i, v in enumerate(self.ids)}
        self.river_names = dict(zip(self.ids, self.names, strict=True))

    def name(self, vannforekomst):
        return self.river_names.get(vannforekomst) or self.lake_names.get(
            vannforekomst
        )

    def positions(self, ids):
        """Where each id sits in the file; -1 for one it does not have."""
        import numpy as np

        return np.array([self.index.get(i, -1) for i in ids])

    def near(self, bounds, margin=5000.0):
        return NearbyLines(self, bounds, margin)


class NearbyLines:
    """The pieces around one region: no worker holds the whole country."""

    def __init__(self, md, bounds, margin):
        import numpy as np
        import shapely

        x0, y0, x1, y1 = bounds
        b = md.bounds
        self.ids = md.ids
        self.keep = np.nonzero(
            (b[:, 0] <= x1 + margin)
            & (b[:, 2] >= x0 - margin)
            & (b[:, 1] <= y1 + margin)
            & (b[:, 3] >= y0 - margin)
        )[0]
        self.geoms = shapely.from_wkb(md.wkb[self.keep])
        self.owner = md.owner[self.keep]
        self.tree = shapely.STRtree(self.geoms)

    def assign(self, geoms, current=None):
        """The vannforekomst each line runs on, and how that was decided.

        Near one: that one. Near several, as at a junction: the one it runs
        along for longest; a tie, where Vann-Nett draws two bodies on the same
        line, stays with `current`, the position the feature already has. Near
        none: the nearest.
        """
        import numpy as np
        import pandas as pd
        import shapely

        geoms = np.asarray(geoms)
        parts, feature = shapely.get_parts(geoms, return_index=True)
        src, piece = self.tree.query(
            parts, predicate="dwithin", distance=LINE_TOLERANCE_M
        )
        pairs = pd.DataFrame({"part": src, "owner": self.owner[piece], "piece": piece})
        n_owners = pairs.groupby("part")["owner"].nunique()
        several = pairs["part"].map(n_owners).to_numpy() > 1
        whole = pairs[~several].drop_duplicates("part")
        at = whole["part"].to_numpy()
        shares = [
            pd.DataFrame(
                {
                    "part": at,
                    "owner": whole["owner"].to_numpy(),
                    "size": shapely.length(parts[at]),
                }
            )
        ]
        rows = []
        for (part, owner), group in pairs[several].groupby(["part", "owner"]):
            line = shapely.union_all(self.geoms[group["piece"].to_numpy()])
            near = shapely.buffer(line, LINE_TOLERANCE_M)
            rows.append(
                (part, owner, shapely.length(shapely.intersection(parts[part], near)))
            )
        if rows:
            shares.append(pd.DataFrame(rows, columns=["part", "owner", "size"]))
        shares = pd.concat(shares, ignore_index=True)
        shares["feature"] = feature[shares["part"].to_numpy()]
        total = shares.groupby(["feature", "owner"])["size"].sum().reset_index()
        if current is not None:
            same = total["owner"].to_numpy() == current[total["feature"].to_numpy()]
            total["size"] += same * 1e-3
        n_found = total.groupby("feature")["owner"].size()
        best = total.sort_values(["feature", "size"]).drop_duplicates(
            "feature", keep="last"
        )
        out = np.full(len(geoms), -1)
        out[best["feature"].to_numpy()] = best["owner"].to_numpy()
        outside = np.nonzero(out < 0)[0]
        max_m = 0.0
        if len(outside):
            i, near = self.tree.query_nearest(geoms[outside], all_matches=False)
            out[outside[i]] = self.owner[near]
            max_m = float(shapely.distance(geoms[outside[i]], self.geoms[near]).max())
        counts = {
            "rows": len(geoms),
            "boundary": int((n_found > 1).sum()),
            "outside": len(outside),
            "outside_max_m": round(max_m, 1),
        }
        return self.ids[out], counts


def lake_vannforekomst(src):
    """lake_id -> its vannforekomst in the file, for one region.

    A lake the file has as a lake is itself. One it does not have, as the
    later download's new lakes in Finnmark, is part of the river running
    through it: the one with the most line inside it.
    """
    import numpy as np
    import pandas as pd
    import pyogrio
    import shapely

    path = src / "delivery" / "lake_polygon.gpkg"
    if not _layer_rows(path, "lake_polygon"):
        return {}
    lakes = pyogrio.read_dataframe(path, columns=["lake_id"])
    known = lakes["lake_id"].isin(MD.lake_names)
    found = dict(zip(lakes["lake_id"][known], lakes["lake_id"][known], strict=True))
    rest = lakes[~known]
    if rest.empty:
        return found
    geoms = shapely.make_valid(rest.geometry.values)
    lines = MD.near(shapely.total_bounds(geoms))
    for lake, geom in zip(rest["lake_id"], geoms, strict=True):
        hit = lines.tree.query(geom, predicate="intersects")
        if len(hit):
            inside = shapely.length(shapely.intersection(lines.geoms[hit], geom))
            per = pd.Series(inside).groupby(lines.owner[hit]).sum()
            found[lake] = MD.ids[per.idxmax()]
        else:
            near = lines.tree.query_nearest(geom, all_matches=False)
            found[lake] = MD.ids[lines.owner[np.atleast_1d(near)[0]]]
    return found


# ── which lake ─────────────────────────────────────────────────────────────
def build_lake_index(path: Path, cache: Path) -> None:
    """The file's lakes as one geometry per id: it has four ids twice."""
    import pyogrio
    import shapely

    if cache.exists():
        return
    rows = pyogrio.read_dataframe(path, layer=MD_LAKES, columns=[MD_ID])
    rows["geometry"] = shapely.make_valid(rows.geometry.values)
    ids, geoms = [], []
    for lake, group in rows.groupby(MD_ID, sort=True):
        ids.append(lake)
        values = group.geometry.values
        geoms.append(values[0] if len(values) == 1 else shapely.union_all(values))
    tmp = cache.with_suffix(".tmp")
    with open(tmp, "wb") as f:
        pickle.dump((ids, shapely.to_wkb(geoms)), f)
    tmp.rename(cache)


class FileLakes:
    def __init__(self, cache: Path):
        import numpy as np
        import shapely

        with open(cache, "rb") as f:
            ids, wkb = pickle.load(f)
        self.ids = np.array(ids, dtype=object)
        self.geoms = shapely.from_wkb(wkb)
        self.index = {v: i for i, v in enumerate(ids)}
        self.tree = shapely.STRtree(self.geoms)

    def geom(self, lake):
        return self.geoms[self.index[lake]]


def _surfaces(lines):
    """The water of one run lake, from its shore pieces.

    polygonize gives every face inside the rings, islands too. A face inside an
    odd number of rings is water. Undelivered pieces count: one missing leaves a
    ring open.
    """
    import numpy as np
    import shapely

    faces = shapely.get_parts(shapely.polygonize(list(lines)))
    if not len(faces):
        return faces
    shells = shapely.polygons(shapely.get_exterior_ring(faces))
    inside, _ = shapely.STRtree(shells).query(
        shapely.point_on_surface(faces), predicate="within"
    )
    return faces[np.bincount(inside, minlength=len(faces)) % 2 == 1]


def _reach(geom):
    """Radius of the largest circle inside the eight largest polygons of geom."""
    import numpy as np
    import shapely

    best = 0.0
    parts = [p for p in _polygons(geom) if p.area > 0]
    for part in sorted(parts, key=lambda p: -p.area)[:8]:
        if np.sqrt(part.area / np.pi) > best:
            circle = shapely.maximum_inscribed_circle(part, tolerance=0.5)
            best = max(best, shapely.length(circle))
    return best


def _wide_parts(geom, limit):
    """The polygons of an overlap that themselves reach limit, or None."""
    import shapely

    keep = [p for p in _polygons(geom) if p.area > 0 and _reach(p) >= limit]
    return shapely.multipolygons(keep) if keep else None


def _distance_to(poly, points):
    """Metres from each point to a polygon, 0 inside it; through an index on its edges."""
    import numpy as np
    import shapely

    poly = shapely.make_valid(poly)
    rings = shapely.get_rings(np.array(_polygons(poly), dtype=object))
    xy, ring = shapely.get_coordinates(rings, return_index=True)
    same = ring[:-1] == ring[1:]
    edges = shapely.linestrings(np.stack([xy[:-1][same], xy[1:][same]], axis=1))
    (at, _), metres = shapely.STRtree(edges).query_nearest(
        points, return_distance=True, all_matches=False
    )
    out = np.full(len(points), np.inf)
    out[at] = metres
    x, y = shapely.get_coordinates(points).T
    shapely.prepare(poly)
    out[shapely.contains_xy(poly, x, y)] = 0.0
    return out


def _water_walk(surface, claims, points):
    """Metres from each point to each claim, walking on the surface's water.

    And each claim's share of the water, the cells nearest it: no line in the map
    says where one lake stops, so the area is divided as the shore is.
    """
    import numpy as np
    import rasterio.features
    import rasterio.transform
    import shapely
    from skimage.graph import MCP_Geometric

    x0, y0, x1, y1 = surface.bounds
    cell = max(WATER_CELL_M, float(np.sqrt((x1 - x0) * (y1 - y0) / WATER_MAX_CELLS)))
    width = int(np.ceil((x1 - x0) / cell)) + 2
    height = int(np.ceil((y1 - y0) / cell)) + 2
    transform = rasterio.transform.from_origin(x0 - cell, y1 + cell, cell, cell)

    def burn(geom, touched):
        return rasterio.features.rasterize(
            [geom],
            out_shape=(height, width),
            transform=transform,
            all_touched=touched,
            dtype="uint8",
        ).astype(bool)

    water = burn(surface, True)
    costs = np.where(water, 1.0, -1.0)  # below zero is impassable
    x, y = shapely.get_coordinates(points).T
    col = np.clip(((x - (x0 - cell)) / cell).astype(int), 0, width - 1)
    row = np.clip((((y1 + cell) - y) / cell).astype(int), 0, height - 1)
    walk = np.full((len(claims), len(points)), np.inf)
    nearest = np.full(water.shape, np.inf)
    owner = np.full(water.shape, -1, dtype=np.int32)
    for k, claim in enumerate(claims):
        start = burn(claim, False) & water
        if not start.any():
            start = burn(claim, True) & water
        if not start.any():
            continue
        rim = start & ~(
            np.roll(start, 1, 0)
            & np.roll(start, -1, 0)
            & np.roll(start, 1, 1)
            & np.roll(start, -1, 1)
        )
        cost, _ = MCP_Geometric(costs, fully_connected=True).find_costs(
            starts=np.argwhere(rim).tolist()
        )
        cost[start] = 0.0
        metres = cost[row, col] * cell
        walk[k] = np.where(np.isfinite(metres), metres, np.inf)
        closer = water & (cost < nearest)
        nearest[closer] = cost[closer]
        owner[closer] = k
    cells = np.bincount(owner[owner >= 0], minlength=len(claims))
    return walk, cells / max(int(cells.sum()), 1)


def _follow_arms(arms, points, labels):
    """Pieces on a narrow arm of another lake take it, if their lake is far off."""
    import numpy as np
    import pandas as pd
    import shapely

    out = labels.copy()
    best = np.full(len(points), np.inf)
    tree = shapely.STRtree(points)
    for lake, overlap, area in sorted(arms, key=lambda a: -a[2]):
        near = tree.query(overlap, predicate="dwithin", distance=TIP_M)
        near = near[labels[near] != lake]
        if not len(near):
            continue
        to_arm = _distance_to(overlap, points[near])
        to_label = np.full(len(near), -1.0)
        for label in pd.unique(labels[near]):
            if label in LAKES.index:
                same = labels[near] == label
                to_label[same] = _distance_to(LAKES.geom(label), points[near][same])
        take = (to_label > FAR_M) & (to_arm <= TIP_M) & (to_arm < best[near])
        out[near[take]] = lake
        best[near[take]] = to_arm[take]
    return out


def _lake_rows(lake, lines, points):
    """One run lake's shore pieces: the lake of each, how, and the area given away.

    A piece lying on no surface goes to the nearest one. A run lake with no surface
    at all, or one whose surface holds no lake of the file, keeps its id: that is
    the 15 lakes the file does not have.
    """
    import numpy as np
    import shapely

    labels = np.full(len(points), lake, dtype=object)
    basis = np.full(len(points), "no_surface", dtype=object)
    given = []
    water = _surfaces(lines)
    if not len(water):
        return labels, basis, given
    tree = shapely.STRtree(water)
    which = np.full(len(points), -1)
    i, j = tree.query(points, predicate="dwithin", distance=0.01)
    which[i] = j
    rest = np.flatnonzero(which < 0)
    if len(rest):
        i, j = tree.query_nearest(points[rest], all_matches=False)
        which[rest[i]] = j
    own = LAKES.index.get(lake)
    for at, surface in enumerate(water):
        here = np.flatnonzero(which == at)
        if not len(here):
            continue
        hits = LAKES.tree.query(surface, predicate="intersects")
        names, overlaps, wide, arms = [], [], [], []
        if any(h != own for h in hits):
            for h in hits:
                overlap = shapely.intersection(surface, LAKES.geoms[h])
                area = shapely.area(overlap)
                if area <= 0:
                    continue
                reach = _reach(overlap)
                if h == own or reach >= REACH_M:
                    names.append(LAKES.ids[h])
                    overlaps.append(overlap)
                    wide.append(_wide_parts(overlap, REACH_M) or overlap)
                if reach >= TIP_M:
                    arms.append((LAKES.ids[h], overlap, area))
        elif own in hits:
            names = [lake]
        if not names:
            basis[here] = "not_in_file"
            continue
        if len(names) == 1:
            basis[here] = "one_lake"
            labels[here] = names[0]
            shares = [1.0]
        else:
            basis[here] = "split"
            order = np.argsort([-o.area for o in overlaps])  # a tie: the larger overlap
            names = np.array(names, dtype=object)[order]
            claims = [wide[k] for k in order]
            q = points[here]
            straight = np.stack([_distance_to(c, q) for c in claims]).argmin(0)
            walk, shares = _water_walk(surface, claims, q)
            reached = np.isfinite(walk)
            by_water = np.where(reached, walk, 1e12).argmin(0)
            labels[here] = names[np.where(reached.any(0), by_water, straight)]
        given += [
            (lake, name, surface.area * share)
            for name, share in zip(names, shares, strict=True)
            if name != lake and share > 0
        ]
        if arms:
            done = _follow_arms(arms, points[here], labels[here])
            basis[here[done != labels[here]]] = "arm"
            labels[here] = done
    return labels, basis, given


def double_lines(lines, lake, run, skip=None):
    """Which shore pieces are a line drawn through open water, from one side.

    Where two mapped water surfaces meet, each side draws its own shore, so the
    line exists twice, one copy on top of the other. Grunnkart does it where it cuts
    a lake at a kommune border into two polygons that stop centimetres short of each
    other (both copies are then under one run lake), and the shore rule above does
    it where it moves the shore of two run lakes that met along a seam to one lake
    of the file. Neither is shore: a piece is a double line when, at three stations
    along it (0.2, 0.5 and 0.8 of its length), another piece of the same lake
    lies within 0.5 m. All three must: a piece that only touches another at its
    midpoint, a neck or a point, is real shore and stays.

    The pieces are tested against the lake they now belong to (lake), not the run
    lake they were cut from. A copy that another region keeps (skip) is gone from
    the result, so it neither is tested nor counts as the other piece. The other
    piece may differ from station to station, and a piece cannot be its own pair.

    Neighbours along a ring are no pair when they reach a station only through
    the joint. Pieces of 10 m have their stations 2 m or more from both ends, so
    the next piece in the ring cannot reach them; only a piece under 2.5 m, a
    ring's last, would see its neighbour at a station. So a piece that meets
    a piece of its own run lake end to end is not counted as the other piece
    at a station unless it lies nearer that station than the joint does. A
    hairpin, a ring folding back over itself, is a line drawn twice like any
    other: its two legs are different pieces, and when one runs within 0.5 m
    of the other at all three stations both go, even if they meet end to end at
    the turn or where the ring closes.

    Both copies go. The test is the same from either side, but the two copies
    can be of different length, so a piece that a flagged piece stood on at all
    three stations is flagged too, even if its own stations lie beyond the other. A lake whose every piece would go (the
    ring of a lake drawn only by its seam) is left whole: there would be nothing
    to count. Returns the flags, and the lakes left whole.
    """
    import numpy as np
    import pandas as pd
    import shapely

    n = len(lines)
    flag = np.zeros(n, dtype=bool)
    if n < 2:
        return flag, []
    used = np.ones(n, dtype=bool) if skip is None else ~np.asarray(skip, dtype=bool)
    at = np.flatnonzero(used)
    if len(at) < 2:
        return flag, []
    part, owner, mark = lines[at], np.asarray(lake)[at], np.asarray(run)[at]
    xy, which = shapely.get_coordinates(part, return_index=True)
    first = np.searchsorted(which, np.arange(len(part)))
    last = np.searchsorted(which, np.arange(len(part)), side="right") - 1
    ends = np.stack([xy[first], xy[last]], axis=1)  # piece, start or end, x or y
    tree = shapely.STRtree(part)
    size = len(part)
    stations = np.zeros((len(DOUBLE_STATIONS), size), dtype=bool)
    pairs = None
    for k, fraction in enumerate(DOUBLE_STATIONS):
        points = shapely.line_interpolate_point(part, fraction, normalized=True)
        xy_at = shapely.get_coordinates(points)
        i, j = tree.query(points, predicate="dwithin", distance=DOUBLE_M)
        keep = (i != j) & (owner[i] == owner[j])
        i, j = i[keep], j[keep]
        # A neighbour in the ring (end to end, one run lake) is the other piece only
        # if it reaches the station other than through the joint: a station
        # near a joint is near the neighbour, and that is no line drawn twice.
        joint = np.full(len(i), np.inf)
        for a in (0, 1):
            for b in (0, 1):
                same = np.abs(ends[i, a] - ends[j, b]).max(axis=1) < 1e-6
                reach = np.hypot(*(xy_at[i] - ends[i, a]).T)
                joint = np.where(same, np.minimum(joint, reach), joint)
        along = shapely.distance(points[i], part[j])
        through = (mark[i] == mark[j]) & (along >= joint - 1e-6)
        i, j = i[~through], j[~through]
        stations[k, i] = True
        both = i.astype(np.int64) * size + j  # the pairs of this station
        pairs = both if pairs is None else np.intersect1d(pairs, both)
    double = stations.all(axis=0)
    # The other copy: a piece that a flagged piece stood on at all three stations.
    twin = np.zeros(size, dtype=bool)
    twin[(pairs % size)[double[pairs // size]]] = True
    double |= twin
    if double.any():
        every = pd.Series(double).groupby(owner).transform("all").to_numpy()
        whole = sorted(set(owner[every]))
        double &= ~every
    else:
        whole = []
    flag[at[double]] = True
    return flag, whole


def part_shore(src, target, min_free):
    """Which lake of the file each shore segment of a region belongs to.

    Every shore segment, those with no zone too: they are in the lake's sums.
    Kept: the rows that move, the rows per lake before and after, and the area.
    """
    import numpy as np
    import pandas as pd
    import pyogrio
    import shapely

    wait_for_memory(min_free)
    t0 = time.time()
    path = src / "geometry" / "segments_10m.gpkg"
    found = pd.DataFrame(columns=["segment_id", LAKE_KEY, "lake", "basis", "m"])
    doubles = pd.DataFrame(columns=["segment_id", "lake", "m"])
    whole = []
    given = []
    if LAKE_KEY in pyogrio.read_info(path)["fields"]:  # not in a region with no lake
        rows = pyogrio.read_dataframe(
            path, where="objektType = 'Innsjøstrand'", columns=["segment_id", LAKE_KEY]
        )
        lines = rows.geometry.values.to_numpy()
        points = shapely.line_interpolate_point(lines, 0.5, normalized=True)
        lake = rows[LAKE_KEY].to_numpy(dtype=object).copy()
        basis = np.empty(len(rows), dtype=object)
        for name, at in rows.groupby(LAKE_KEY, sort=False).indices.items():
            lake[at], basis[at], parts = _lake_rows(name, lines[at], points[at])
            given += parts
        found = pd.DataFrame(
            {
                "segment_id": rows["segment_id"].to_numpy(),
                LAKE_KEY: rows[LAKE_KEY].to_numpy(),
                "lake": lake,
                "basis": basis,
                "m": shapely.length(lines),
            }
        )
        # A copy another region keeps is not here for the test, nor in the sums.
        copies = found["segment_id"].isin(SHORE_COPIES.get(src.name, ())).to_numpy()
        double, whole = double_lines(lines, lake, rows[LAKE_KEY].to_numpy(), copies)
        doubles = found.loc[double, ["segment_id", "lake", "m"]]
    tmp = target.with_suffix(".tmp")
    pd.to_pickle(
        {
            "moved": found[found["lake"] != found[LAKE_KEY]],
            "before": found[LAKE_KEY].value_counts(),
            # Without the double lines: what the lake keeps in this region.
            "after": found.loc[~found["segment_id"].isin(doubles["segment_id"]), "lake"]
            .value_counts(),
            "basis": found["basis"].value_counts(),
            "given": pd.DataFrame(given, columns=[LAKE_KEY, "lake", "m2"]),
            "doubles": doubles,
            "whole": whole,
        },
        tmp,
    )
    tmp.rename(target)
    moved = int((found["lake"] != found[LAKE_KEY]).sum())
    return (
        f"{target.stem}: {len(found):,} shore, {moved:,} moved, {len(doubles):,} double, "
        f"{time.time() - t0:.0f} s"
    )


def _midline_values(lake, geom, regions, source):
    """The values of a lake's own midline, for its shore.

    For a lake the run did not select, so nothing tagged the midline: tag_inside's
    test, the midpoint within the polygon, then the run's own inheritance.
    """
    import pandas as pd
    import pyogrio
    import shapely
    from kantsoner import stages
    from kantsoner.spec import constants as C

    found = []
    for region in sorted(regions):
        rows = pyogrio.read_dataframe(
            source / region / "attributes" / "derived.gpkg",
            bbox=tuple(geom.bounds),
            columns=["segment_id", C.KILDE_COL, "objektType", *INHERITED],
        )
        rows = rows[rows[C.KILDE_COL] != C.KILDE_LAKE]
        if len(rows):
            middle = shapely.line_interpolate_point(
                rows.geometry.values, 0.5, normalized=True
            )
            inside = rows[shapely.within(middle, geom)]
            found.append(pd.DataFrame(inside.drop(columns="geometry")))
    shore = pd.DataFrame(
        {
            "segment_id": ["shore"],
            C.KILDE_COL: [C.KILDE_LAKE],
            "objektType": ["Innsjøstrand"],
            **{c: [None] for c in INHERITED},
        }
    )
    frame = pd.concat([*found, shore], ignore_index=True).assign(**{LAKE_KEY: lake})
    last = stages._lake_network_inheritance(frame).iloc[-1]
    return tuple(None if pd.isna(last[c]) else last[c] for c in INHERITED)


def _shore_values(moved, before, source, md, lakes):
    """lake -> (elvenavn, *INHERITED) for the shore each lake receives.

    Barrier and Fragmentering belong to the lake, not to the segment,
    so a moved segment takes them from the lake's own shore, or its midline.
    """
    targets = set(moved["lake"])
    values = {}
    wanted = '"lake_id", "elvenavn", ' + ", ".join(f'"{c}"' for c in INHERITED)
    held = before[before["lake"].isin(targets)]
    for region, group in held.groupby("region"):
        file = source / region / "delivery" / "lake_shore.gpkg"
        if not file.exists():
            continue
        names = ", ".join(f"'{lake}'" for lake in group["lake"])
        where = f"WHERE lake_id IN ({names})"
        for lake, *row in _rows(file, "lake_shore", f"DISTINCT {wanted}", where):
            if values.setdefault(lake, tuple(row)) != tuple(row):
                raise RuntimeError(f"{lake}: its shore has two sets of lake values")
    where = moved.groupby("lake")["region"].agg(set)
    for lake in sorted(targets - set(values)):
        found = _midline_values(lake, lakes.geom(lake), where[lake], source)
        values[lake] = (md.lake_names[lake], *found)
    return values


def plan_shore(regions, source, work, md_cache) -> None:
    """What the shore rule changes, for the workers and the national rebuild."""
    import pandas as pd

    path = work / SHORE_PLAN
    if path.exists():
        return

    def counts(found, name):
        return pd.concat(
            [
                f[name].rename("rows").rename_axis("lake").reset_index().assign(region=r)
                for r, f in found.items()
            ],
            ignore_index=True,
        )

    found = {r: pd.read_pickle(work / "shore" / f"{r}.pkl") for r in regions}
    moved = pd.concat(
        [f["moved"].assign(region=r) for r, f in found.items()], ignore_index=True
    )
    before, after = counts(found, "before"), counts(found, "after")
    gone = set(before["lake"]) - set(after["lake"])
    if gone:
        raise RuntimeError(f"lakes that lose all their shore: {sorted(gone)}")
    md, lakes = Vannforekomster(md_cache), FileLakes(work / LAKE_INDEX)
    unknown = set(moved["lake"]) - set(lakes.index)
    if unknown:
        raise RuntimeError(f"moved to lakes the file does not have: {sorted(unknown)}")
    given = pd.concat([f["given"] for f in found.values()], ignore_index=True)
    area = given.groupby("lake")["m2"].sum().sub(
        given.groupby(LAKE_KEY)["m2"].sum(), fill_value=0.0
    )
    doubles = pd.concat(
        [f["doubles"].assign(region=r) for r, f in found.items()], ignore_index=True
    )
    regions_of = after.groupby("lake")["region"].nunique()
    touched = {
        *moved[LAKE_KEY],
        *moved["lake"],
        *regions_of.index[regions_of > 1],
        *doubles["lake"],  # a lake that loses a line has new sums
    }
    wanted = after[after["lake"].isin(touched)].groupby("region")["lake"].agg(set)
    plan = {
        "moved": {
            r: dict(zip(g["segment_id"], g["lake"], strict=True))
            for r, g in moved.groupby("region")
        },
        "values": _shore_values(moved, before, source, md, lakes),
        "area": area.to_dict(),
        "wanted": wanted.to_dict(),
        "doubles": {
            r: list(g["segment_id"]) for r, g in doubles.groupby("region")
        },
    }
    doubles.rename(columns={"lake": LAKE_KEY, "region": FOLDER})[
        [FOLDER, "segment_id", LAKE_KEY, "m"]
    ].round({"m": 3}).to_csv(work / "double_lines.csv", index=False)
    basis = pd.concat([f["basis"] for f in found.values()]).groupby(level=0).sum()
    (work / "shore_moved.json").write_text(
        json.dumps(
            {
                "rows": len(moved),
                "km": float(moved["m"].sum() / 1000),
                "arm_rows": int((moved["basis"] == "arm").sum()),
                "lakes_gaining": int(moved["lake"].nunique()),
                "lakes_losing": int(moved[LAKE_KEY].nunique()),
                "new_lakes": len(set(moved["lake"]) - set(before["lake"])),
                "basis": {k: int(v) for k, v in basis.items()},
                "doubles": {
                    "rows": len(doubles),
                    "km": float(doubles["m"].sum() / 1000),
                    "lakes": int(doubles["lake"].nunique()),
                    "regions": int(doubles["region"].nunique()),
                    "left_whole": sorted({w for f in found.values() for w in f["whole"]}),
                },
            },
            indent=1,
        )
    )
    tmp = path.with_suffix(".tmp")
    with open(tmp, "wb") as f:
        pickle.dump(plan, f)
    tmp.rename(path)


def _add_counts(total, counts):
    for k in ("rows", "boundary", "outside"):
        total[k] = total.get(k, 0) + counts[k]
    total["outside_max_m"] = max(
        total.get("outside_max_m", 0.0), counts["outside_max_m"]
    )
    return total


def _chunks(path, layer, **kw):
    """Read a layer CHUNK rows at a time, with the fid as the index."""
    import pyogrio

    n = pyogrio.read_info(path, layer=layer)["features"]
    for skip in range(0, n, CHUNK):
        yield pyogrio.read_dataframe(
            path,
            layer=layer,
            skip_features=skip,
            max_features=CHUNK,
            fid_as_index=True,
            **kw,
        )


def _layer_rows(path: Path, layer: str) -> int:
    import pyogrio

    if not path.exists():
        return 0
    return pyogrio.read_info(path, layer=layer)["features"]


def _save_stats(part: Path, stats: dict) -> None:
    part.with_suffix(".json").write_text(json.dumps(stats))


def _rows(path, table, columns, where=""):
    """Rows straight from SQLite, for the few columns that are all we need."""
    return [r for chunk in _row_chunks(path, table, columns, where) for r in chunk]


def _row_chunks(path, table, columns, where=""):
    uri = "file:" + quote(str(path), safe="/:") + "?mode=ro&immutable=1"
    con = sqlite3.connect(uri, uri=True)
    try:
        rows = con.execute(f'SELECT {columns} FROM "{table}" {where}')
        while chunk := rows.fetchmany(CHUNK):
            yield chunk
    finally:
        con.close()


def _hash(values):
    import numpy as np
    from pandas.util import hash_array

    return hash_array(np.array(values, dtype=object), categorize=False)


def _in_two_regions(hashes, places):
    """The hashes seen in more than one region."""
    import numpy as np

    order = np.lexsort((places, hashes))
    h, p = hashes[order], places[order]
    return np.unique(h[1:][(h[1:] == h[:-1]) & (p[1:] != p[:-1])])


def _shared_rows(source, regions, name, water):
    """The rows of a table whose id or geometry is also in another region.

    Hashes find them, in memory a national list of strings would not fit; the
    rows themselves are then compared by the caller.
    """
    import numpy as np
    import pandas as pd

    files = {i: source / r / "delivery" / f"{name}.gpkg" for i, r in enumerate(regions)}
    files = {i: f for i, f in files.items() if f.exists()}
    ids, shapes, places = [], [], []
    for i, file in files.items():
        for chunk in _row_chunks(file, name, '"segment_id", geom'):
            segment, shape = zip(*chunk, strict=True)
            ids.append(_hash(segment))
            shapes.append(_hash(shape))
            places.append(np.full(len(chunk), i, dtype=np.int32))
    if not ids:
        return None
    ids, shapes, places = (np.concatenate(x) for x in (ids, shapes, places))
    same_id = _in_two_regions(ids, places)
    same_shape = _in_two_regions(shapes, places)
    found = []
    columns = ["segment_id", "water", "shape"]
    for i in np.unique(places[np.isin(ids, same_id) | np.isin(shapes, same_shape)]):
        for chunk in _row_chunks(files[i], name, f'"segment_id", "{water}", geom'):
            rows = pd.DataFrame(chunk, columns=columns)
            hit = np.isin(_hash(rows["segment_id"]), same_id) | np.isin(
                _hash(rows["shape"]), same_shape
            )
            found.append(rows[hit].assign(table=name, region=regions[i]))
    return pd.concat(found, ignore_index=True) if found else None


def _new_shore_ids(source, clash):
    """(region, id, new id) for a lake's shore outside the region owning it.

    Each region numbered its shore from 0, so the owner's ids stay and the others
    go on from its last number.
    """
    new = []

    def ids(region, lake):
        file = source / region / "delivery" / "lake_shore.gpkg"
        where = f"WHERE lake_id = '{lake}'"
        return sorted(r[0] for r in _rows(file, "lake_shore", '"segment_id"', where))

    for lake, group in clash.groupby("water"):
        owner = group.loc[group["own"], "region"].unique()
        if len(owner) != 1:
            raise RuntimeError(f"{lake}: shore in two regions, and no one owner")
        top = int(ids(owner[0], lake)[-1].rsplit("_", 1)[1])
        for region in sorted(set(group.loc[~group["own"], "region"])):
            for old in ids(region, lake):
                top += 1
                new.append((region, old, f"{lake}_{top:06d}"))
    return new


def find_copies(source, regions, path) -> None:
    """Water two regions both delivered, and which copy of each goes.

    A copy is the same geometry in two regions. The id alone proves nothing
    either way: shore is numbered from 0 in each region, so one id can be two
    different shores, and a stretch can be called _p1 in one region and not in
    the next. What only shares its id is left, and gets an id of its own.
    """
    import pandas as pd

    if path.exists():
        return
    columns = ["segment_id", "water", "table", "region", "own"]
    found = [
        _shared_rows(source, regions, name, water)
        for name, water in (
            ("official_centerline_10m", "vannforekomstid"),
            ("lake_shore", "lake_id"),
        )
    ]
    found = [f for f in found if f is not None]
    copies, new_ids = pd.DataFrame(columns=columns), []
    if found:
        rows = pd.concat(found, ignore_index=True)
        rows["own"] = rows["region"].str[:3] == rows["water"].str.split("-").str[0]
        twice = rows.groupby(["table", "shape"])["region"].transform("nunique") > 1
        copies = rows[twice]
        owners = (
            copies[copies["own"]].groupby(["table", "shape"])["region"].nunique()
        ).reindex(copies.groupby(["table", "shape"]).size().index, fill_value=0)
        if (owners != 1).any():
            raise RuntimeError(f"{(owners != 1).sum()} copies have no one owner")
        others = rows[~twice]
        clash = others[
            others.groupby(["table", "segment_id"])["region"].transform("nunique") > 1
        ]
        if (clash["table"] != "lake_shore").any():
            raise RuntimeError("river segments share an id but not their geometry")
        new_ids = _new_shore_ids(source, clash)
        copies = copies[~copies["own"]][columns]
    with open(path.with_name(NEW_SHORE_IDS), "wb") as f:
        pickle.dump(new_ids, f)
    part = path.with_suffix(".tmp")
    with open(part, "wb") as f:
        pickle.dump(copies, f)
    part.rename(path)


def segment_ids(src: Path, cache: Path, lakes: dict, moved: dict):
    """Per delivered segment and shore: its vannområde, vannforekomst and lake.

    One lookup, so a segment, its buffer polygons and its place in the
    roll-ups can never disagree. Returns it, and the counts per product.
    """
    import pandas as pd
    import pyogrio

    if cache.exists():
        return pd.read_pickle(cache)
    frames, stats = [], {}
    name = "official_centerline_10m"
    path = src / "delivery" / f"{name}.gpkg"
    if _layer_rows(path, name):
        bounds = pyogrio.read_info(path, layer=name, force_total_bounds=True)
        lines = MD.near(bounds["total_bounds"])
        total, on_lines = {}, {}
        columns = ["segment_id", "vannforekomstid"]
        for chunk in _chunks(path, name, columns=columns):
            geoms = chunk.geometry.values.to_numpy()
            vo, counts = VO.assign(geoms)
            _add_counts(total, counts)
            current = MD.positions(chunk["vannforekomstid"])
            md, counts = lines.assign(geoms, current)
            _add_counts(on_lines, counts)
            frames.append(
                pd.DataFrame(
                    {"vo": vo, "md": md}, index=chunk["segment_id"].to_numpy()
                )
            )
        stats[name] = {**total, "md": on_lines}
    name = "lake_shore"
    path = src / "delivery" / f"{name}.gpkg"
    if _layer_rows(path, name):
        total = {}
        for chunk in _chunks(path, name, columns=["segment_id", "lake_id"]):
            vo, counts = VO.assign(chunk.geometry.values)
            _add_counts(total, counts)
            lake = chunk["segment_id"].map(moved).astype(object)
            frames.append(
                pd.DataFrame(
                    {
                        "vo": vo,
                        "md": lake.fillna(chunk["lake_id"].map(lakes)).to_numpy(),
                        "lake": lake.fillna(chunk["lake_id"]).to_numpy(),
                    },
                    index=chunk["segment_id"].to_numpy(),
                )
            )
        stats[name] = total
    columns = ["vo", "md", "lake"]
    seg = pd.concat(frames) if frames else pd.DataFrame(columns=columns)
    seg = seg.reindex(columns=columns).astype(dict.fromkeys(columns, "category"))
    tmp = cache.with_suffix(".tmp")
    pd.to_pickle((seg, stats), tmp)
    tmp.rename(cache)
    return seg, stats


# ── one product of one region ──────────────────────────────────────────────
def correct_grades(con, table: str, level: str) -> None:
    """Score_total_final as Sállir tabell 8 and chapter 3.11 of the report have it.

    A region made by older code has the grade reversed and is turned here; the
    sign of the grade against the score tells, so that a region made by the
    current code stays as it is. 6 minus the grade, not a new grade from the
    score: the score is rounded and the grade is not, and a block of 0.19998
    reads 0.2000.
    A Kantvegetasjon value with no class marks too little data exactly, where
    the rounded datadekning could not.
    """
    from kantsoner.spec import constants as C

    score, grade = f"{level}_Score_total", f"{level}_Score_total_final"
    cover, veg = f"{level}_datadekning", f"{level}_Kantvegetasjon"
    sign = con.execute(
        f'SELECT sum(("{grade}" - 3) * ("{score}" - 0.5)) FROM "{table}"'
    ).fetchone()[0]
    if sign is not None and sign < 0:
        con.execute(
            f'UPDATE "{table}" SET "{grade}" = 6 - "{grade}" '
            f'WHERE "{grade}" IS NOT NULL'
        )
    con.execute(
        f'UPDATE "{table}" SET "{grade}" = NULL WHERE "{cover}" IS NULL '
        f'OR "{cover}" < ? OR ("{veg}" IS NOT NULL AND "{veg}_final" IS NULL)',
        (C.DATADEKNING_MIN,),
    )


def grade_rule(frame, level: str) -> None:
    """The same for rows computed here. The scale itself is set in tabell_8."""
    import pandas as pd
    from kantsoner.spec import constants as C

    cover = pd.to_numeric(frame[f"{level}_datadekning"], errors="coerce").fillna(0)
    grade, veg = f"{level}_Score_total_final", f"{level}_Kantvegetasjon"
    refused = frame[veg].notna() & frame[f"{veg}_final"].isna()
    frame[grade] = frame[grade].where((cover >= C.DATADEKNING_MIN) & ~refused)


def tabell_8() -> None:
    """The scale of Sállir tabell 8 for rows computed here, whatever --code has
    (older code had it reversed)."""
    from kantsoner.spec import scales

    scales._SCORE_B_SCORES = [1, 2, 3, 4, 5]


def _gpkg_blob(geom, srs_id):
    """A geometry as GeoPackage stores it: header, envelope, little-endian WKB."""
    import shapely

    x0, y0, x1, y1 = shapely.bounds(geom)
    head = b"GP\x00\x03" + struct.pack("<i4d", srs_id, x0, x1, y0, y1)
    return head + shapely.to_wkb(geom, byte_order=1)


def _only_polygons(con, path, name):
    """A few shore rows are polygons plus a stray line: keep the polygons.

    The layer is MultiPolygon, and a strict reader skips anything else. The
    line has no area, so every measurement on the row stays true.
    """
    import numpy as np
    import shapely

    column, srs_id = con.execute(
        "SELECT column_name, srs_id FROM gpkg_geometry_columns WHERE table_name = ?",
        (name,),
    ).fetchone()
    for chunk in _chunks(path, name, columns=[]):
        geoms = chunk.geometry.values.to_numpy()
        for i in np.nonzero(shapely.get_type_id(geoms) == 7)[0]:
            geom = _polygon_part(geoms[i])
            if not geom.is_valid:
                geom = _polygon_part(shapely.make_valid(geom))
            con.execute(
                f'UPDATE "{name}" SET "{column}" = ? WHERE fid = ?',
                (_gpkg_blob(geom, srs_id), int(chunk.index[i])),
            )


def not_computed(aggregation, table) -> list:
    """The segment_id of every row the run did not compute (Beregnet = no).

    From the region's aggregation table: the delivery has no Beregnet. Those rows
    read Uklassifisert, or Ikke vurdert (bekkelukking) where the stream is a
    culvert and was not computed either. A culvert that was computed is not here.
    """
    return [r[0] for r in _rows(aggregation, table, '"segment_id"', "WHERE \"Beregnet\" = 'no'")]


def blank_unclassified(con, table, not_computed=()) -> int:
    """Empty the score and the priority of the rows that were not computed.

    A row is not computed when it reads Uklassifisert, or when its segment_id is
    in `not_computed` (Beregnet = no in the aggregation, which also catches the
    culverts that were not computed). Returns the number of rows. A table
    without these columns raises: the delivery has them, and a silent skip would
    leave the values in.
    """
    con.execute("CREATE TEMP TABLE IF NOT EXISTS _not_computed (segment_id)")
    con.execute("DELETE FROM _not_computed")
    con.executemany(
        "INSERT INTO _not_computed VALUES (?)", ((i,) for i in not_computed)
    )
    sets = ", ".join(f'"{c}" = NULL' for c in BLANK_UNCLASSIFIED)
    return con.execute(
        f'UPDATE "{table}" SET {sets} WHERE "E_S_Kantvegetasjon_A" = ? '
        "OR segment_id IN (SELECT segment_id FROM _not_computed)",
        (UNCLASSIFIED,),
    ).rowcount


def part_delivered(name, src, part, region, seg, stats, lakes):
    """Copy the file and set what differs from it.

    Copied, not rewritten, so the rows, types, geometry and every measurement
    stay exactly what was delivered, except that the 10 m rows that were not
    computed lose their score and priority (blank_unclassified).
    """
    import pyogrio

    path = src / "delivery" / f"{name}.gpkg"
    if not _layer_rows(path, name):  # no file, or 0 rows: nothing to merge
        return
    tmp = part.with_name(part.stem + ".tmp.gpkg")
    shutil.copyfile(path, tmp)
    con = sqlite3.connect(tmp)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    # The R-tree triggers call SpatiaLite functions sqlite3 does not have, and
    # a part needs no index: it is only read through once.
    triggers = con.execute(
        "SELECT name FROM sqlite_master WHERE type='trigger' AND name LIKE ?",
        (f"rtree_{name}_%",),
    ).fetchall()
    for (trigger,) in triggers:
        con.execute(f'DROP TRIGGER "{trigger}"')
    con.execute(f'ALTER TABLE "{name}" ADD COLUMN "{VO_ID}" TEXT')
    con.execute(f'ALTER TABLE "{name}" ADD COLUMN "{FOLDER}" TEXT')

    def update(columns, rows):
        sets = ", ".join(f'"{c}" = ?' for c in columns)
        con.executemany(f'UPDATE "{name}" SET {sets} WHERE fid = ?', rows)

    out = {}
    if name in ("official_centerline_10m", "lake_shore"):
        out = stats[name]
        columns = [VO_ID, FOLDER, "vannforekomstid", "river_id_vann-nett",
                   "vannforekomstnavn"]
        for chunk in _chunks(path, name, columns=["segment_id"], read_geometry=False):
            look = seg.reindex(chunk["segment_id"].to_numpy())
            vo, md = _text(look["vo"]), _text(look["md"])
            update(
                columns,
                (
                    (v, region, m, m, MD.name(m), fid)
                    for v, m, fid in zip(vo, md, chunk.index, strict=True)
                ),
            )
        if name == "official_centerline_10m":
            skipped = not_computed(src / "aggregation" / f"{name}.gpkg", name)
            out = {**out, "blanked_unclassified": blank_unclassified(con, name, skipped)}
        if name == "lake_shore":
            _only_polygons(con, path, name)
            if region in MOVED:
                for chunk in _chunks(
                    path, name, columns=["segment_id"], read_geometry=False
                ):
                    lake = chunk["segment_id"].map(MOVED[region]).dropna()
                    update(
                        [LAKE_KEY, "elvenavn", *INHERITED],
                        ((new, *TARGETS[new], fid) for fid, new in lake.items()),
                    )
    elif name in BLOCKS:
        # Blocks of a river that lost segments here are built again at the end.
        gone = sorted(COPIED_RIVERS.get(region, ()))
        marks = ",".join("?" * len(gone))
        removed = con.execute(
            f'DELETE FROM "{name}" WHERE "river_id_vann-nett" IN ({marks})', gone
        ).rowcount
        bounds = pyogrio.read_info(path, layer=name, force_total_bounds=True)
        lines = MD.near(bounds["total_bounds"])
        columns = [VO_ID, FOLDER, "river_id_vann-nett", "vannforekomstnavn"]
        for chunk in _chunks(path, name, columns=["river_id_vann-nett"]):
            geoms = chunk.geometry.values.to_numpy()
            vo, counts = VO.assign(geoms)
            _add_counts(out, counts)
            current = MD.positions(chunk["river_id_vann-nett"])
            md, _ = lines.assign(geoms, current)
            update(
                columns,
                (
                    (v, region, m, MD.name(m), fid)
                    for v, m, fid in zip(vo, md, chunk.index, strict=True)
                ),
            )
        out["rows"] -= removed
        out["removed"] = removed
    elif name == "lake_polygon":
        con.execute(f'ALTER TABLE "{name}" ADD COLUMN "vannforekomstid" TEXT')
        con.execute(f'ALTER TABLE "{name}" ADD COLUMN "vannforekomstnavn" TEXT')
        columns = [VO_ID, FOLDER, "vannforekomstid", "vannforekomstnavn", "lake_navn"]
        for chunk in _chunks(path, name, columns=["lake_id", "lake_navn"]):
            vo, counts = VO.assign(chunk.geometry.values)
            _add_counts(out, counts)
            update(
                columns,
                (
                    (
                        v,
                        region,
                        lakes.get(lake),
                        MD.name(lakes.get(lake)),
                        MD.lake_names.get(lake, navn),
                        fid,
                    )
                    for v, lake, navn, fid in zip(
                        vo, chunk["lake_id"], chunk["lake_navn"], chunk.index,
                        strict=True,
                    )
                ),
            )
    else:  # a region: its first river's name says nothing about it
        columns = [VO_ID, FOLDER, "vannforekomstnavn"]
        for chunk in _chunks(path, name, columns=[]):
            vo, counts = VO.assign(chunk.geometry.values)
            _add_counts(out, counts)
            update(
                columns,
                ((v, region, None, fid) for v, fid in zip(vo, chunk.index, strict=True)),
            )
    if name in GRADED:
        correct_grades(con, name, GRADED[name])
    con.commit()
    con.close()
    _save_stats(part, out)
    tmp.rename(part)


def _base_columns(path, aggregate, sides, extra=()):
    """What _base_counts reads, so ~130 of the 330 columns come off disk."""
    import pyogrio

    want = {
        aggregate.SEGMENT_ID,
        aggregate.RIVER_ID,
        aggregate.LENGTH,
        aggregate.NAME,
        aggregate.BEREGNET,
        aggregate.CLASS_COL,
        "E_S_Kantvegetasjon_fjernet",
        "K_S_Rorlagt_mistanke",
        "K_S_beregningsbredde",
        *sides.OUTSIDE_COLS.values(),
        *sides.NO_BANK_COLS.values(),
        *aggregate.SHARE_COLUMNS,
        *aggregate.COUNT_COLUMNS,
        *aggregate.MEAN_COLUMNS,
        *aggregate.WORST_COLUMNS,
        aggregate.CHAIN_ID,
        aggregate.CHAIN_POS,
        *extra,
    }
    for side in sides.SIDES:
        for band in ("B0-10", "B0-20"):
            for cat in aggregate.CAUSE_AREA_CATEGORIES:
                want.add(f"K_{side}_{band}_{cat}_m2")
            for cat in aggregate.CAUSE_PCT_CATEGORIES:
                want.add(f"K_{side}_{band}_{cat}")
    fields = pyogrio.read_info(path)["fields"]
    return [c for c in fields if c in want]


def part_rollups(src, part, region, seg, work):
    """The additive sums of this region, per vannområde and vannforekomst.

    Without its copies. Sums over regions and over chunks of a region are the
    same sums, so the national roll-ups get exactly what one pass over all
    segments would. Stray water, whose id belongs to another region, is also
    kept aside for the rows that are recomputed at the end, under the run's id,
    which is what those rows were built on.
    """
    import pandas as pd
    from kantsoner import sides
    from kantsoner.delivery import aggregate

    name = "official_centerline_10m"
    path = src / "aggregation" / f"{name}.gpkg"
    (work / "side_rows").mkdir(exist_ok=True)
    (work / "id_map").mkdir(exist_ok=True)
    if not _layer_rows(path, name):
        (work / "side_rows" / f"{region}.done").touch()
        return
    columns = _base_columns(path, aggregate, sides)
    pieces, aside, ids = [], [], []
    gone = COPIES.get(region, set())
    level = aggregate.Level("part", None, None, (VO_ID, MD_KEY))
    for chunk in _chunks(path, name, columns=columns):
        chunk = chunk[~chunk[aggregate.SEGMENT_ID].isin(gone)]
        if chunk.empty:
            continue
        look = seg.reindex(chunk[aggregate.SEGMENT_ID].to_numpy())
        chunk[VO_ID] = look["vo"].astype(object).to_numpy()
        chunk[MD_KEY] = look["md"].astype(object).to_numpy()
        if chunk[VO_ID].isna().any() or chunk[MD_KEY].isna().any():
            raise RuntimeError(f"{region}: segment without an id")
        river = chunk[aggregate.RIVER_ID].astype(str)
        stray = river.str.split("-").str[0] != region[:3]
        if stray.any():
            aside.append(chunk[stray].assign(**{FOLDER: region}))
        ids.append(
            chunk.groupby([aggregate.RIVER_ID, MD_KEY])[aggregate.LENGTH].sum()
        )
        pieces.append(
            aggregate._aggregate_one(chunk, level, aggregate.SEGMENT_ID, True)
        )
    if aside:
        pd.concat(aside).to_pickle(work / "side_rows" / f"{region}.pkl")
    (work / "side_rows" / f"{region}.done").touch()
    if not pieces:
        return
    pd.concat(ids).groupby(level=[0, 1]).sum().to_pickle(
        work / "id_map" / f"{region}.pkl"
    )
    tmp = part.with_suffix(".tmp")
    pd.concat(pieces, ignore_index=True).assign(**{FOLDER: region}).to_pickle(tmp)
    tmp.rename(part)


def _run_shore_table(src):
    """Every shore segment as the run rolled it up, those with no zone too.

    Only stages.aggregation has that table: it fills in what a segment with no
    zone is, and neither its input nor its output can stand in. So it is run here on
    the shore rows of derived.gpkg, and stopped where it hands them to the roll-up.
    """
    import pyogrio
    from kantsoner import sides, stages
    from kantsoner.delivery import aggregate
    from kantsoner.delivery import lakes as run
    from kantsoner.spec import constants as C

    derived = src / "attributes" / "derived.gpkg"
    columns = _base_columns(
        derived,
        aggregate,
        sides,
        ("segment_id", LAKE_KEY, "routed", "fallback", "elvID", "elvenavn", "objektType",
         "vannforekomstid", "stretch_id", C.KILDE_COL, "E_S_Kantvegetasjon_grad"),
    )

    class Reached(Exception):
        pass

    def read(path, **kw):
        if path != derived:
            return real_read(path, **kw)
        where = f"\"{C.KILDE_COL}\" = '{C.KILDE_LAKE}'"
        return pyogrio.read_dataframe(path, where=where, columns=columns)

    def stop(gdf, polys):
        raise Reached(gdf)

    real_read, real_rollup = stages.read_gpkg, run.aggregate_lakes
    stages.read_gpkg, run.aggregate_lakes = read, stop
    try:
        stages.aggregation(
            None, SimpleNamespace(derived=derived, geometry_dir=src / "geometry")
        )
    except Reached as found:
        return found.args[0]
    finally:
        stages.read_gpkg, run.aggregate_lakes = real_read, real_rollup
    raise RuntimeError(f"{src}: the run never reached its lakes")


def part_lake_sums(src, region, work):
    """The sums of the lakes that change, from this region's shore as it now stands.

    A lake's sums are its shore's, so nothing but the shore's lake changes, and the
    inherited values with it. The table must be the run's before anything is done
    with it, so its sums are held against the run's first.
    """
    import pandas as pd
    import pyogrio
    from kantsoner.delivery import aggregate

    table = _run_shore_table(src)
    run_rows = pyogrio.read_dataframe(
        src / "aggregation" / "lake_polygon.gpkg", read_geometry=False
    ).drop_duplicates(LAKE_KEY).set_index(LAKE_KEY)

    def sums(frame):
        base = aggregate._base_counts(frame)
        base[LAKE_KEY] = frame[LAKE_KEY].to_numpy()
        return base.groupby(LAKE_KEY).sum(numeric_only=True)

    columns = aggregate._sum_columns(run_rows)
    before = sums(table)[columns].astype(float)
    worst = (before - run_rows[columns].astype(float).reindex(before.index)).abs()
    if set(before.index) != set(run_rows.index) or not worst.max().max() < 1e-6:
        raise RuntimeError(f"{region}: the shore does not give the run's sums")
    # The rows the merge does not keep are not in the sums either: the copies another
    # region holds, and the double lines, with or without a zone.
    out_of = set(SHORE_COPIES.get(region, ())) | DOUBLES.get(region, set())
    table = table[~table["segment_id"].isin(out_of)].copy()
    moved = set(MOVED.get(region, {})) - out_of
    if not moved <= set(table["segment_id"]):
        raise RuntimeError(f"{region}: moved shore that the table does not have")
    lake = table["segment_id"].map(MOVED.get(region, {}))
    changed = lake.notna().to_numpy()
    table.loc[changed, LAKE_KEY] = lake[changed]
    for at, col in enumerate(INHERITED, start=1):
        table[col] = table[col].astype(object)
        table.loc[changed, col] = pd.Series(
            [TARGETS[new][at] for new in lake[changed]],
            index=table.index[changed],
            dtype=object,
        )
    out = sums(table[table[LAKE_KEY].isin(LAKE_SUMS[region])])
    out["areal_vannflate_m2"] = (
        run_rows["areal_vannflate_m2"].astype(float).reindex(out.index).fillna(0.0)
    )
    target = work / "lake_sums" / f"{region}.pkl"
    tmp = target.with_suffix(".tmp")
    out.to_pickle(tmp)
    tmp.rename(target)


def part_coverage(src, region, work):
    """Which pieces of the file's lines have a delivered segment of this region.

    Also what this region delivered off the lines, measured at the segment's
    middle: a segment only touching a line at its end is not on it.
    """
    import numpy as np
    import pyogrio
    import shapely

    name = "official_centerline_10m"
    path = src / "delivery" / f"{name}.gpkg"
    target = work / "coverage" / f"{region}.npy"
    if not _layer_rows(path, name):
        np.save(target, np.array([], dtype=np.int64))
        return
    bounds = pyogrio.read_info(path, layer=name, force_total_bounds=True)
    lines = MD.near(bounds["total_bounds"])
    geoms, touched, off_n, off_m, away = [], [], 0, 0.0, []
    for chunk in _chunks(path, name, columns=["segment_id"]):
        g = chunk.geometry.values.to_numpy()
        geoms.append(g)
        touched.append(
            np.unique(
                lines.tree.query(
                    g, predicate="dwithin", distance=LINE_TOLERANCE_M
                )[1]
            )
        )
        middle = shapely.line_interpolate_point(g, 0.5, normalized=True)
        on = lines.tree.query(
            middle, predicate="dwithin", distance=LINE_TOLERANCE_M
        )[0]
        off = np.ones(len(g), dtype=bool)
        off[on] = False
        off_n += int(off.sum())
        off_m += float(shapely.length(g[off]).sum())
        if off.any():
            _, distance = lines.tree.query_nearest(
                middle[off], return_distance=True, all_matches=False
            )
            away.append(
                {"segment_id": chunk["segment_id"].to_numpy()[off], "m": distance}
            )
    tree = shapely.STRtree(np.concatenate(geoms))
    touched = np.unique(np.concatenate(touched))
    present = []
    for block in np.array_split(touched, max(1, len(touched) // 20_000)):
        if not len(block):
            continue
        fine, piece, k = _fine_pieces(lines.geoms[block])
        hit = np.unique(
            tree.query(fine, predicate="dwithin", distance=LINE_TOLERANCE_M)[0]
        )
        where = MD.offset[lines.keep[block[piece]]] + k
        present.append(where[hit])
    tmp = target.with_suffix(".tmp.npy")
    np.save(tmp, np.concatenate(present) if present else np.array([], dtype=np.int64))
    target.with_suffix(".json").write_text(
        json.dumps({"off_line_segments": off_n, "off_line_m": off_m})
    )
    if away:
        with open(target.with_suffix(".off.pkl"), "wb") as f:
            pickle.dump(away, f)
    tmp.rename(target)


def _fine_pieces(lines):
    """Midpoints of the short pieces of each line, with their line and number."""
    import numpy as np
    import shapely

    xy, at = shapely.get_coordinates(shapely.segmentize(lines, FINE_M), return_index=True)
    same = at[1:] == at[:-1]
    first = np.searchsorted(at, np.arange(len(lines)))
    piece = at[:-1][same]
    k = np.arange(len(xy) - 1)[same] - first[piece]
    return shapely.points((xy[:-1][same] + xy[1:][same]) / 2), piece, k


def _cut_to_ring(zones, rings):
    """zone ∩ ring per row; keeps the polygon part of what comes back."""
    import numpy as np
    import shapely

    try:
        out = shapely.intersection(zones, rings)
    except shapely.errors.GEOSException:
        zones, rings = shapely.make_valid(zones), shapely.make_valid(rings)
        try:
            out = shapely.intersection(zones, rings)
        except shapely.errors.GEOSException:
            out = shapely.intersection(zones, rings, grid_size=0.01)
    # A zone that only touches the ring comes back as a line or a point.
    kind = shapely.get_type_id(out)
    out[~np.isin(kind, (3, 6, 7))] = None
    for i in np.nonzero(kind == 7)[0]:
        out[i] = _polygon_part(out[i])
    # Now and then one polygon comes back as two, the second inside the first.
    for i in np.nonzero(~shapely.is_valid(out) & ~shapely.is_missing(out))[0]:
        out[i] = _polygon_part(shapely.make_valid(out[i]))
    return out


def _bank_table(src, name, ident, banks):
    """The columns a bank inherits, one row per segment.

    Text is held as categories: on Glomma, 5 million rows of plain strings in
    a dozen columns took more memory than the rest of the job together.
    """
    import numpy as np
    import pandas as pd
    import pyogrio
    from pandas.api.types import union_categoricals

    path = src / "aggregation" / f"{name}.gpkg"
    if not _layer_rows(path, name):
        return None
    per_side = [
        "K_{s}_B0-10_Vtot",
        "K_{s}_B0-10_Atot_m2",
        "E_{s}_Potensiell_Kantvegetasjon",
        "E_{s}_Kantvegetasjon_A",
        "E_{s}_Kantvegetasjon_A_arsak",
    ]
    numbers = {c.format(s=s) for s in banks for c in per_side[:2]}
    cols = ["segment_id", *ident]
    cols += [c.format(s=s) for s in banks for c in per_side]
    missing = set(cols) - set(pyogrio.read_info(path, layer=name)["fields"])
    if missing:
        raise RuntimeError(f"{path}: no {sorted(missing)}")
    frames = []
    for chunk in _chunks(path, name, columns=cols, read_geometry=False):
        for col in cols:
            if col in numbers:
                chunk[col] = pd.to_numeric(chunk[col], errors="coerce")
            elif col != "segment_id":
                chunk[col] = chunk[col].astype("string").astype("category")
        frames.append(chunk.reset_index(drop=True))
    table = {}
    for col in cols:
        pieces = [f[col] for f in frames]
        if col in numbers or col == "segment_id":
            table[col] = np.concatenate([p.to_numpy() for p in pieces])
        else:
            table[col] = union_categoricals(pieces)
    table = pd.DataFrame(table)
    if table["segment_id"].duplicated().any():
        raise RuntimeError(f"{path}: a segment twice")
    return table.set_index("segment_id")


def _text(column):
    """A category column as written: plain text, or null."""
    return column.astype(object).where(column.notna(), None).to_numpy()


def part_bank_zones(src, part, region, seg, product):
    """buffer10 or lake_buffer10: the 0-10 m cut of each segment's zones."""
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    import pyogrio
    import shapely

    if product == LAKE_BUFFER:
        name, fields = "lake_shore", LAKE_BUFFER_FIELDS
        sides = (("L", "left"),)  # a shore has one bank
        ident = ["lake_id", "elvenavn"]
    else:
        name, fields = "official_centerline_10m", BUFFER_FIELDS
        sides = (("L", "left"), ("R", "right"))
        ident = ["elvenavn"]
    banks = _bank_table(src, name, ident, "".join(c for c, _ in sides))
    if banks is None:
        return
    # Rings stay as WKB and only those a chunk needs become geometries.
    path = src / "geometry" / f"{RING}.gpkg"
    crs = pyogrio.read_info(path)["crs"]
    ids, wkb = [], []
    for chunk in _chunks(path, RING, columns=["id"]):
        ids.append(chunk["id"].to_numpy())
        wkb.append(shapely.to_wkb(chunk.geometry.values.to_numpy()))
    ids, wkb = np.concatenate(ids), np.concatenate(wkb)
    position = pd.Series(np.arange(len(ids)), index=ids)
    tmp = part.with_name(part.stem + ".tmp.gpkg")
    tmp.unlink(missing_ok=True)
    stats = {"rows": 0, "no_ring": 0, "empty": 0, "area_m2": 0.0, "atot_m2": 0.0}
    for code, side in sides:
        layer = f"zones_{side}"
        for chunk in _chunks(
            src / "geometry" / f"{layer}.gpkg",
            layer,
            columns=["segment_id", "stretch_id"],
        ):
            chunk = chunk[chunk["segment_id"].isin(banks.index)]
            at = chunk["stretch_id"].map(position)
            stats["no_ring"] += int(at.isna().sum())
            chunk, at = chunk[at.notna()], at.dropna().astype(int).to_numpy()
            if chunk.empty:
                continue
            needed, inverse = np.unique(at, return_inverse=True)
            geoms = _cut_to_ring(
                chunk.geometry.values.to_numpy(),
                shapely.from_wkb(wkb[needed])[inverse],
            )
            keep = ~(shapely.is_missing(geoms) | shapely.is_empty(geoms))
            stats["empty"] += int((~keep).sum())
            chunk, geoms = chunk[keep], geoms[keep]
            if chunk.empty:
                continue
            row = banks.reindex(chunk["segment_id"].to_numpy())
            area = shapely.area(geoms)
            segment = chunk["segment_id"].to_numpy()
            look = seg.reindex(segment)
            water = _text(look["md"])
            values = {
                FOLDER: region,
                VO_ID: _text(look["vo"]),
                "river_id_vann-nett": water,
                "vannforekomstid": water,
                "vannforekomstnavn": [MD.name(w) for w in water],
                "segment_id": segment,
                "stretch_id": chunk["stretch_id"].to_numpy(),
                "side": code,
                # The delivery's name for segment_id: aggregation sets one to the other.
                "segment_id_vann-nett": segment,
                "area_m2": np.round(area, 2),
                "K_B0-10_Vtot": row[f"K_{code}_B0-10_Vtot"].to_numpy(),
                "E_Potensiell_Kantvegetasjon": _text(
                    row[f"E_{code}_Potensiell_Kantvegetasjon"]
                ),
                "E_Kantvegetasjon_A": _text(row[f"E_{code}_Kantvegetasjon_A"]),
                "E_Kantvegetasjon_A_arsak": _text(
                    row[f"E_{code}_Kantvegetasjon_A_arsak"]
                ),
                **{col: _text(row[col]) for col in ident},
            }
            if product == LAKE_BUFFER:
                values[LAKE_KEY] = _text(look["lake"])
                moved = MOVED.get(region)
                if moved:
                    # A shore segment moved to another lake has that lake's elvenavn,
                    # as in lake_shore (part_delivered).
                    new = pd.Series(segment).map(moved)
                    if new.notna().any():
                        values["elvenavn"] = [
                            TARGETS[n][0] if pd.notna(n) else e
                            for n, e in zip(new, values["elvenavn"], strict=True)
                        ]
            stats["rows"] += len(segment)
            stats["area_m2"] += float(area.sum())
            stats["atot_m2"] += float(
                row[f"K_{code}_B0-10_Atot_m2"].fillna(0).sum()
            )
            gdf = gpd.GeoDataFrame(
                pd.DataFrame(values)[list(fields)], geometry=geoms, crs=crs
            )
            pyogrio.write_dataframe(
                gdf,
                tmp,
                layer=product,
                append=tmp.exists(),
                promote_to_multi=True,
                geometry_type="MultiPolygon",
                layer_options={"SPATIAL_INDEX": "NO"},
            )
    if tmp.exists():
        _save_stats(part, stats)
        tmp.rename(part)


def crown_means(chm, geoms, top_x, top_y):
    """Mean CHM per crown, over the pixels the crown polygon covers.

    Crowns are vectorised from CHM pixels, so a pixel's centre being inside
    the polygon is exactly membership of the crown. Crowns are burned a window
    at a time; one partly hidden under another, as at a batch edge, is measured
    again on its own.
    """
    import numpy as np
    import pandas as pd
    import shapely
    from affine import Affine
    from rasterio import features, windows

    out = np.full(len(geoms), np.nan)
    cells = pd.Series(
        np.floor(top_x / TREE_CELL_M) * 1e7 + np.floor(top_y / TREE_CELL_M)
    )
    for _, rows in cells.groupby(cells).groups.items():
        rows = np.asarray(rows)
        x0, y0, x1, y1 = shapely.total_bounds(geoms[rows])
        col0, row0 = np.floor(~chm.transform * (x0, y1)).astype(int)
        col1, row1 = np.ceil(~chm.transform * (x1, y0)).astype(int)
        win = windows.Window(col0, row0, col1 - col0, row1 - row0)
        height = chm.read(1, window=win, boundless=True, fill_value=np.nan)
        height = height.astype(np.float64)
        valid = np.isfinite(height) & (height != chm.nodata)
        transform = chm.window_transform(win)

        labels = features.rasterize(
            zip(geoms[rows], range(1, len(rows) + 1)),
            out_shape=height.shape,
            transform=transform,
            fill=0,
            dtype="int32",
        )
        burned = np.bincount(labels.ravel(), minlength=len(rows) + 1)[1:]
        ok = valid & (labels > 0)
        n = np.bincount(labels[ok], minlength=len(rows) + 1)[1:]
        total = np.bincount(
            labels[ok], weights=height[ok], minlength=len(rows) + 1
        )[1:]
        with np.errstate(invalid="ignore", divide="ignore"):
            out[rows] = total / n
        # A crown is whole pixels, so its area is its pixel count.
        pixels = np.rint(shapely.area(geoms[rows]))
        for k in np.nonzero(burned < pixels)[0]:
            x0, y0, x1, y1 = geoms[rows[k]].bounds
            c0, r0 = np.floor(~transform * (x0, y1)).astype(int)
            c1, r1 = np.ceil(~transform * (x1, y0)).astype(int)
            mask = features.rasterize(
                [(geoms[rows[k]], 1)],
                out_shape=(r1 - r0, c1 - c0),
                transform=transform * Affine.translation(c0, r0),
                fill=0,
                dtype="uint8",
            ).astype(bool)
            mask &= valid[r0:r1, c0:c1]
            if mask.any():
                out[rows[k]] = height[r0:r1, c0:c1][mask].mean()
    return out


def apex_keys(x, y):
    """A crown's apex is one pixel, whichever region found it."""
    import numpy as np

    return np.rint(x * 100).astype(np.int64) * (1 << 31) + np.rint(y * 100).astype(
        np.int64
    )


def part_apex(src, target):
    """Every crown's apex, how far its stretch is, and whether the region owns it."""
    import numpy as np

    path = src / "delivery" / f"{TREES}.gpkg"
    keys, dist, own = [], [], []
    if path.exists():
        columns = '"top_x", "top_y", "dist_2D", "stretch_id"'
        for chunk in _row_chunks(path, TREES, columns):
            x, y, d, stretch = zip(*chunk, strict=True)
            keys.append(apex_keys(np.array(x, dtype=float), np.array(y, dtype=float)))
            dist.append(np.array(d, dtype=np.float32))
            own.append(np.array([s[:3] == src.name[:3] for s in stretch], dtype=bool))
    tmp = target.with_suffix(".tmp.npz")
    np.savez(
        tmp,
        key=np.concatenate(keys) if keys else np.empty(0, np.int64),
        dist=np.concatenate(dist) if dist else np.empty(0, np.float32),
        own=np.concatenate(own) if own else np.empty(0, bool),
    )
    tmp.rename(target)
    return f"{src.name}: apexes read"


def find_tree_copies(regions, work) -> None:
    """Crowns two regions found, and which copy of each goes.

    A crown near a border is detected by both regions, at the same pixel and the
    same height. Which stays: the copy nearest its stretch, so that stretch_id
    holds; at equal distance, which is a stretch in both regions, the region owning
    it; then the lower number.
    """
    import numpy as np
    import pandas as pd

    files = [work / "apex" / f"{r}.npz" for r in regions]
    keys = [np.load(f)["key"] for f in files]
    places = [np.full(len(k), i, dtype=np.int16) for i, k in enumerate(keys)]
    shared = _in_two_regions(np.concatenate(keys), np.concatenate(places))
    del places
    found = []
    for i, (file, key) in enumerate(zip(files, keys, strict=True)):
        if not len(shared) or not len(key):
            continue
        at = np.minimum(np.searchsorted(shared, key), len(shared) - 1)
        hit = shared[at] == key
        if hit.any():
            with np.load(file) as z:
                dist, own = z["dist"][hit], z["own"][hit]
            found.append(
                pd.DataFrame({"key": key[hit], "dist": dist, "own": own, "place": i})
            )
    drops = {}
    if found:
        rows = pd.concat(found, ignore_index=True).sort_values(
            ["key", "dist", "own", "place"], ascending=[True, True, False, True]
        )
        gone = rows[rows.duplicated("key")]
        drops = {
            regions[i]: np.sort(g["key"].to_numpy()) for i, g in gone.groupby("place")
        }
    tmp = work / f"{TREE_COPIES}.tmp"
    with open(tmp, "wb") as f:
        pickle.dump(drops, f)
    tmp.rename(work / TREE_COPIES)


def part_trees(src, part, region, seg, chm_path):
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    import pyogrio
    import rasterio
    import shapely

    path = src / "delivery" / f"{TREES}.gpkg"
    if not _layer_rows(path, TREES):
        return
    tmp = part.with_name(part.stem + ".tmp.gpkg")
    tmp.unlink(missing_ok=True)
    stats = {}
    no_mean = 0
    cols = ["tree_id", "tree_height", "height_class", "top_x", "top_y", "stretch_id"]
    with rasterio.open(chm_path) as chm:
        for chunk in _chunks(path, TREES, columns=cols):
            x, y = chunk["top_x"].to_numpy(), chunk["top_y"].to_numpy()
            if region in TREE_DROPS:
                chunk = chunk[~np.isin(apex_keys(x, y), TREE_DROPS[region])]
                if chunk.empty:
                    continue
                x, y = chunk["top_x"].to_numpy(), chunk["top_y"].to_numpy()
            ids, counts = VO.assign(shapely.points(x, y))
            _add_counts(stats, counts)
            geoms = chunk.geometry.values.to_numpy()
            mean = crown_means(chm, geoms, x, y)
            no_mean += int(np.isnan(mean).sum())
            out = pd.DataFrame(
                {
                    FOLDER: region,
                    VO_ID: ids,
                    "tree_id": chunk["tree_id"].to_numpy(),
                    "stretch_id": chunk["stretch_id"].to_numpy(),
                    "max_height": chunk["tree_height"].to_numpy(),
                    "mean_height": np.round(mean, 2),
                    "crown_area_m2": np.round(shapely.area(geoms), 1),
                    "height_class": chunk["height_class"].to_numpy(),
                }
            )
            outline = shapely.simplify(geoms, TREE_SIMPLIFY_M, preserve_topology=True)
            pyogrio.write_dataframe(
                gpd.GeoDataFrame(out, geometry=outline, crs=chunk.crs),
                tmp,
                layer=TREES,
                append=tmp.exists(),
                layer_options={"SPATIAL_INDEX": "NO"},
            )
    if not tmp.exists():  # every crown of the region belongs to another
        return
    stats["no_mean"] = no_mean
    _save_stats(part, stats)
    tmp.rename(part)


# ── one region ─────────────────────────────────────────────────────────────
def available_gb() -> float:
    with open("/proc/meminfo") as f:
        info = dict(line.split(":", 1) for line in f)
    return int(info["MemAvailable"].split()[0]) / 2**20


def wait_for_memory(min_gb: float) -> None:
    while available_gb() < min_gb:
        time.sleep(30)


def stop_before_the_machine_does(stop_gb: float) -> None:
    """Stop the whole run while memory is still there.

    The machine may be running something that cannot be restarted so cheaply,
    and a stopped merge only redoes the regions it was in the middle of.
    """
    while True:
        if available_gb() < stop_gb:
            try:
                print(f"only {available_gb():.0f} GB free: stopping", flush=True)
            finally:  # a broken stdout must not cost the machine
                os.killpg(0, signal.SIGKILL)
        time.sleep(2)


def die_with_parent() -> None:
    """A worker must not outlive main: a relaunch would then write its parts
    at the same time, and nothing would show that rows had gone."""
    parent = os.getppid()

    def watch():
        while os.getppid() == parent:
            time.sleep(2)
        os._exit(1)

    threading.Thread(target=watch, daemon=True).start()


def limit_worker(code, max_gb) -> None:
    os.environ.setdefault("GDAL_NUM_THREADS", "1")
    os.environ.setdefault("GDAL_CACHEMAX", "256")
    # If memory does run out, the kernel takes these first; and one worker
    # cannot grow past its share however large the region.
    die_with_parent()
    Path("/proc/self/oom_score_adj").write_text("1000")
    limit = int(max_gb * 2**30)
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    sys.path.insert(0, code)
    from loguru import logger

    logger.remove()


def init_shore_worker(code, lake_cache, copies, max_gb) -> None:
    global LAKES, SHORE_COPIES
    limit_worker(code, max_gb)
    LAKES = FileLakes(Path(lake_cache))
    with open(copies, "rb") as f:
        found = pickle.load(f)
    shores = found[found["table"] == "lake_shore"]
    SHORE_COPIES = {r: set(g["segment_id"]) for r, g in shores.groupby("region")}


def init_worker(code, cache, md_cache, copies, plan, max_gb) -> None:
    global VO, MD, COPIES, COPIED_RIVERS, SHORE_COPIES, MOVED, TARGETS, LAKE_SUMS
    global TREE_DROPS, DOUBLES
    limit_worker(code, max_gb)
    tabell_8()
    VO = Vannomrader(Path(cache))
    MD = Vannforekomster(Path(md_cache))
    with open(copies, "rb") as f:
        found = pickle.load(f)
    rivers = found[found["table"] == "official_centerline_10m"]
    shores = found[found["table"] == "lake_shore"]
    COPIES = {r: set(g["segment_id"]) for r, g in rivers.groupby("region")}
    COPIED_RIVERS = {r: set(g["water"]) for r, g in rivers.groupby("region")}
    SHORE_COPIES = {r: set(g["segment_id"]) for r, g in shores.groupby("region")}
    if Path(plan).exists():  # not there when only the trees are merged
        with open(plan, "rb") as f:
            found = pickle.load(f)
        MOVED, TARGETS, LAKE_SUMS = found["moved"], found["values"], found["wanted"]
        DOUBLES = {r: set(ids) for r, ids in found["doubles"].items()}
    tree_copies = Path(copies).with_name(TREE_COPIES)
    if tree_copies.exists():
        with open(tree_copies, "rb") as f:
            TREE_DROPS = pickle.load(f)


def part_path(work: Path, product: str, region: str) -> Path:
    if product in ROLLUPS:  # one pass makes both
        return work / "parts" / "rollups" / f"{region}.pkl"
    return work / "parts" / product / f"{region}.gpkg"


def do_region(region, source, work, products, chm, min_free):
    t0 = time.time()
    src = source / region
    done = []
    seg = stats = lakes = None
    for product in products:
        part = part_path(work, product, region)
        if part.exists():
            continue
        wait_for_memory(min_free)
        if lakes is None and product != TREES:
            lakes = lake_vannforekomst(src)
        if seg is None and product != TREES and product != "lake_polygon":
            seg, stats = segment_ids(
                src, work / "segments" / f"{region}.pkl", lakes, MOVED.get(region, {})
            )
        if product in ROLLUPS:
            part_rollups(src, part, region, seg, work)
        elif product in (BUFFER, LAKE_BUFFER):
            part_bank_zones(src, part, region, seg, product)
        elif product == TREES:
            part_trees(src, part, region, seg, chm)
        else:
            part_delivered(product, src, part, region, seg, stats, lakes)
        done.append(product)
    covered = work / "coverage" / f"{region}.npy"
    if "official_centerline_10m" in products and not covered.exists():
        wait_for_memory(min_free)
        part_coverage(src, region, work)
    sums = work / "lake_sums" / f"{region}.pkl"
    if "lake_polygon" in products and region in LAKE_SUMS and not sums.exists():
        wait_for_memory(min_free)
        part_lake_sums(src, region, work)
    return f"{region}: {len(done)} products, {(time.time() - t0) / 60:.1f} min"


# ── national ───────────────────────────────────────────────────────────────
def finish_rollup(product, work, regions, out, cache, md_cache) -> dict:
    """The national roll-up: the run's own second step, on the region sums."""
    import geopandas as gpd
    import pandas as pd
    import pyogrio
    from kantsoner.delivery import aggregate, trim

    frames = [
        pd.read_pickle(part_path(work, product, r))
        for r in sorted(regions)
        if part_path(work, product, r).exists()  # a region with no segments
    ]
    source = gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True), geometry="geometry", crs=CRS
    )
    name, key = ("river", MD_KEY) if product == RIVER else ("vannomraade", VO_ID)
    level = aggregate.Level(name, "part", None, (key,))
    gdf = aggregate._aggregate_one(
        source, level, "aggr_segment_id_vann-nett_part", False
    )
    grade_rule(gdf, name)
    keep = trim.fields_for(name)
    if product == RIVER:
        md = Vannforekomster(md_cache)
        length = "aggr_segment_length_vann-nett_part"
        gdf[aggregate.RIVER_ID] = gdf[MD_KEY]
        gdf[aggregate.NAME] = [md.name(i) for i in gdf[MD_KEY]]
        gdf[VO_ID] = Vannomrader(cache).assign(gdf.geometry.values)[0]
        gdf[FOLDER] = gdf[MD_KEY].map(_folders(source, MD_KEY, length))
        keep += [VO_ID, FOLDER]
    else:
        with open(cache, "rb") as f:
            table = pickle.load(f)[0].set_index(VO_ID)
        gdf["vannomradeId"] = gdf[VO_ID].map(table["vannomradeId"])
        gdf["vannomradenavn"] = gdf[VO_ID].map(table["navn"])
        gdf["n_rivers"] = gdf[VO_ID].map(source.groupby(VO_ID)[MD_KEY].nunique())
        # The region list, with the vannområde in place of the region.
        at = keep.index(aggregate.RIVER_ID)
        keep[at : at + 2] = [VO_ID, "vannomradeId", "vannomradenavn", "n_rivers"]
    gdf = gpd.GeoDataFrame(gdf[keep], geometry=gdf.geometry.values, crs=CRS)
    pyogrio.write_dataframe(
        gdf,
        out,
        layer=product,
        promote_to_multi=True,
        geometry_type="MultiLineString",
    )
    return {"rows": len(gdf)}


def append_parts(product: str, parts: list[Path], out: Path, consume=False) -> dict:
    """Native append: nothing is re-encoded, and the first part sets the types.

    The largest part goes first. Its columns have both values and gaps, so it
    has the widest types: K_S_Laser_aar is REAL there but INTEGER in a region
    with no gaps, and vannforekomstnavn has no type in one with no names.

    With consume a part goes as soon as it is in, so the disk holds the parts or
    the file, not both. A stopped merge then starts over, after do_region has made
    the missing parts again.
    """
    for part in parts:  # a short part would otherwise pass unnoticed
        stats = json.loads(part.with_suffix(".json").read_text())
        if stats["rows"] != _layer_rows(part, product):
            raise RuntimeError(f"{part}: rows differ from what was written")
    parts = sorted(parts, key=lambda p: (-p.stat().st_size, p.name))
    # Not the Arrow path: it refuses a part whose types differ from the first.
    env = {
        **os.environ,
        "OGR_SQLITE_SYNCHRONOUS": "OFF",
        "OGR2OGR_USE_ARROW_API": "NO",
    }
    out.unlink(missing_ok=True)
    rows = 0
    for i, part in enumerate(parts):
        rows += _layer_rows(part, product)
        subprocess.run(
            [
                "ogr2ogr",
                "-f",
                "GPKG",
                *([] if i == 0 else ["-update", "-append"]),
                "-nln",
                product,
                "-nlt",
                GEOMETRY[product],
                "-gt",
                "100000",
                "-unsetFid",
                "-lco",
                "SPATIAL_INDEX=NO",
                str(out),
                str(part),
                product,
            ],
            check=True,
            env=env,
            capture_output=True,
            text=True,
        )
        if consume:
            part.unlink()
    got = _layer_rows(out, product)
    if got != rows:
        raise RuntimeError(f"{product}: {got} rows merged of {rows}")
    return {"rows": rows}


def spatial_index(out: Path, product: str) -> None:
    subprocess.run(
        [
            "ogrinfo",
            "-q",
            str(out),
            "-sql",
            f"SELECT CreateSpatialIndex('{product}', 'geom')",
        ],
        check=True,
        env={**os.environ, "OGR_SQLITE_SYNCHRONOUS": "OFF"},
        capture_output=True,
        text=True,
    )


def delete_rows(path, product, key, rows) -> int:
    """Delete the rows named by (region_folder, key) pairs. Returns how many."""
    con = sqlite3.connect(path)
    con.execute("PRAGMA synchronous=OFF")
    con.execute("CREATE TEMP TABLE gone (folder TEXT, id TEXT, PRIMARY KEY (folder, id))")
    con.executemany("INSERT OR IGNORE INTO gone VALUES (?, ?)", rows)
    n = con.execute(
        f'DELETE FROM "{product}" WHERE EXISTS (SELECT 1 FROM gone WHERE '
        f'gone.folder = "{product}"."{FOLDER}" AND gone.id = "{product}"."{key}")'
    ).rowcount
    con.commit()
    con.close()
    return n


def _append(path, product, rows, cache) -> None:
    """Add recomputed rows, in the layer's columns, with their vannområde."""
    import geopandas as gpd
    import pandas as pd
    import pyogrio

    if rows.empty:
        return
    geoms = rows.geometry.values
    rows = pd.DataFrame(rows.drop(columns="geometry"))
    rows[VO_ID] = Vannomrader(cache).assign(geoms)[0]
    fields = pyogrio.read_info(path, layer=product)["fields"]
    pyogrio.write_dataframe(
        gpd.GeoDataFrame(rows[list(fields)], geometry=geoms, crs=CRS),
        path,
        layer=product,
        append=True,
        promote_to_multi=True,
    )


def _side_rows(work, regions, copies):
    """The segments kept aside, already without their copies."""
    import geopandas as gpd
    import pandas as pd

    folder = work / "side_rows"
    missing = [
        r
        for r in set(copies["region"])
        if r in regions and not (folder / f"{r}.pkl").exists()
    ]
    frames = [
        pd.read_pickle(folder / f"{r}.pkl")
        for r in sorted(regions)
        if (folder / f"{r}.pkl").exists()
    ]
    if missing or (len(copies) and not frames):
        # Not every copy leaves something behind, but its region's part ran.
        have = [r for r in missing if part_ran(work, r)]
        if len(have) != len(missing):
            raise RuntimeError(
                "side rows are missing: run official_centerline_vannomraade too"
            )
    if not frames:
        return None
    return gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True), geometry="geometry", crs=CRS
    )


def part_ran(work, region) -> bool:
    return (work / "side_rows" / f"{region}.done").exists()


def _folders(rows, key, length):
    """Per id, the regions holding its segments, longest first."""
    per = rows.groupby([key, FOLDER])[length].sum().reset_index()
    per = per.sort_values([key, length], ascending=[True, False])
    return per.groupby(key)[FOLDER].agg("; ".join)


def _rivers_per_region(work, regions):
    """Vannforekomster delivered, counted by the region their id belongs to."""
    import pandas as pd

    ids = set()
    for region in regions:
        path = work / "id_map" / f"{region}.pkl"
        if path.exists():
            ids.update(pd.read_pickle(path).index.get_level_values(1))
    return pd.Series(sorted(ids), dtype=object).str.split("-").str[0].value_counts()


def _rebuild_regions(path, product, aside, source, per_region, cache) -> int:
    """A region's stray rows join its own: the sums of what remains of them."""
    import pandas as pd
    from kantsoner.delivery import aggregate

    key = "region_id"
    table = pd.DataFrame(
        _rows(path, product, f'"{key}", "{FOLDER}"'), columns=["id", "folder"]
    )
    is_own = table["folder"].str[:3] == table["id"]
    stray = table[~is_own]
    delete_rows(path, product, key, stray[["folder", "id"]].values.tolist())
    kept = None
    if aside is not None:
        prefix = aside[aggregate.RIVER_ID].astype(str).str.split("-").str[0]
        kept = aside[prefix != aside[FOLDER].str[:3]].assign(**{key: prefix})
    if kept is not None and len(kept):
        _join_stray_rows(path, product, table[is_own], kept, source, per_region, cache)
    con = sqlite3.connect(path)
    con.executemany(
        f'UPDATE "{product}" SET n_rivers = ? WHERE "{key}" = ?',
        [(int(n), i) for i, n in per_region.items()],
    )
    con.commit()
    con.close()
    return len(stray)


def _join_stray_rows(path, product, owners, kept, source, per_region, cache):
    import geopandas as gpd
    import pandas as pd
    import pyogrio
    from kantsoner.delivery import aggregate

    key, id_col = "region_id", "aggr_segment_id_vann-nett_region"
    length = "aggr_segment_length_vann-nett_region"
    own = owners.set_index("id")["folder"]
    level = aggregate.Level("region", None, None, (key,))
    extra = aggregate._aggregate_one(kept, level, aggregate.SEGMENT_ID, True)
    frames, numbers, names = [extra], {}, {}
    for region_id in extra[key]:
        if region_id not in own.index:  # its own folder is not in this -r run
            continue
        file = source / own[region_id] / "aggregation" / f"{product}.gpkg"
        row = pyogrio.read_dataframe(file, where=f"{key} = '{region_id}'")
        frames.append(row)
        numbers[region_id] = int(row[id_col].iloc[0])
        names[region_id] = own[region_id]
    both = gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True), geometry="geometry", crs=CRS
    ).drop(columns="n_rivers")
    both["_order"] = -both[length]
    both[aggregate.RIVER_ID] = ""  # the roll-up counts these; n_rivers is below
    level = aggregate.Level("region", "region", None, (key,))
    out = aggregate._aggregate_one(both, level, "_order", False)
    grade_rule(out, "region")
    out[id_col] = [numbers.get(i, 0) for i in out[key]]
    out["n_rivers"] = out[key].map(per_region)
    out[aggregate.NAME] = None  # as on every region row
    others = _folders(kept, key, aggregate.LENGTH)
    out[FOLDER] = [
        "; ".join(filter(None, (names.get(i), others[i]))) for i in out[key]
    ]
    delete_rows(path, product, key, [[f, i] for i, f in names.items()])
    _append(path, product, out, cache)


def _rebuild_blocks(path, product, aside, rivers, work, regions, cache, md_cache):
    """The blocks of a river in the region that lost its copies.

    Their old rows went when the part was made. The chain is the run's own:
    100 m from the segments, then 500 and 1000 m from those. The owner's
    blocks are not touched. Returns how many rows fewer the product has.
    """
    import pandas as pd
    import shapely
    from kantsoner.delivery import aggregate

    name = BLOCKS[product]
    id_col = f"aggr_segment_id_vann-nett_{name}"
    key = aggregate.RIVER_ID
    removed = sum(
        json.loads(p.read_text()).get("removed", 0)
        for p in (work / "parts" / product).glob("*.json")
    )
    if aside is None:
        return removed
    con = sqlite3.connect(path)
    top = dict(
        con.execute(
            f'SELECT "{FOLDER}", max("{id_col}") FROM "{product}" '
            f'GROUP BY "{FOLDER}"'
        ).fetchall()
    )
    con.close()
    md = Vannforekomster(md_cache)
    levels = {level.name: level for level in aggregate.LEVELS}
    pairs = rivers[["region", "water"]].drop_duplicates()
    rebuilt = []
    for region, group in pairs.groupby("region"):
        rows = aside[
            (aside[FOLDER] == region) & aside[key].isin(group["water"])
        ]
        if rows.empty:
            continue
        out, order, base = rows, aggregate.CHAIN_POS, True
        for step in ("100m", "500m", "1000m"):
            out = aggregate._aggregate_one(out, levels[step], order, base)
            order, base = f"aggr_segment_id_vann-nett_{step}", False
            if step == name:
                break
        grade_rule(out, name)
        geoms = out.geometry.values.to_numpy()
        lines = md.near(shapely.total_bounds(geoms))
        ids, _ = lines.assign(geoms, md.positions(out[key]))
        out[key] = ids
        out[aggregate.NAME] = [md.name(i) for i in ids]
        out[id_col] = range(top.get(region, 0) + 1, top.get(region, 0) + 1 + len(out))
        rebuilt.append(out.assign(**{FOLDER: region}))
    if rebuilt:
        _append(path, product, pd.concat(rebuilt, ignore_index=True), cache)
    return removed - sum(len(r) for r in rebuilt)


def _one_row_per_lake(path, product) -> int:
    """One row per lake.

    A lake twice in one region is the same numbers on two near-equal polygons, and
    keeps the larger. In two regions, each holds only its own shore: the owner's
    row stays, and rebuild_lakes gives it the sums of both.
    """
    import pyogrio

    rows = pyogrio.read_dataframe(
        path, layer=product, columns=[LAKE_KEY, FOLDER], fid_as_index=True
    )
    rows["own"] = rows[FOLDER].str[:3] == rows[LAKE_KEY].str.split("-").str[0]
    rows["area"] = rows.geometry.area
    rows = rows[rows[LAKE_KEY].duplicated(keep=False)]
    keep = rows.sort_values([LAKE_KEY, "own", "area"]).drop_duplicates(
        LAKE_KEY, keep="last"
    )
    gone = rows.drop(index=keep.index)
    con = sqlite3.connect(path)
    con.executemany(
        f'DELETE FROM "{product}" WHERE fid = ?', [(int(i),) for i in gone.index]
    )
    con.commit()
    con.close()
    return len(gone)


def rebuild_lakes(path, product, work, regions, cache, md_cache) -> None:
    """The rows of the lakes whose shore moved or lies in two regions.

    The roll-up is the run's: additive sums, then the ratios. A lake's polygon stays;
    one that gains its first shore gets a row with the file's. region_folder names
    the regions holding its shore, longest first.
    """
    import geopandas as gpd
    import pandas as pd
    import pyogrio
    from kantsoner.delivery import aggregate
    from kantsoner.delivery import lakes as run

    frames = [
        pd.read_pickle(work / "lake_sums" / f"{r}.pkl").assign(**{FOLDER: r})
        for r in sorted(regions)
        if (work / "lake_sums" / f"{r}.pkl").exists()
    ]
    if not frames:
        return
    with open(work / SHORE_PLAN, "rb") as f:
        plan = pickle.load(f)
    area = pd.Series(plan["area"], dtype=float)
    per = pd.concat(frames).rename_axis(LAKE_KEY).reset_index()
    block = per.groupby(LAKE_KEY)[aggregate._sum_columns(per)].sum()
    wanted = {lake for r in regions for lake in plan["wanted"].get(r, ())}
    if set(block.index) != wanted:
        raise RuntimeError(f"{len(set(block.index) ^ wanted)} lakes differ from the plan")
    out = aggregate._ratios(block, run.LEVEL)
    grade_rule(out, run.LEVEL)
    out[f"{run.LEVEL}_strandlinje_m"] = block["length"].round(1)
    out[f"{run.LEVEL}_segmenter"] = block["n_10m_segments"].astype("int64")
    surface = per.groupby(LAKE_KEY)["areal_vannflate_m2"].sum()
    out["areal_vannflate_m2"] = (
        surface.add(area, fill_value=0.0).reindex(out.index).round(1)
    )
    longest = per.sort_values("length", ascending=False)
    out[FOLDER] = longest.groupby(LAKE_KEY)[FOLDER].agg("; ".join)

    def plain(value):
        if pd.isna(value):
            return None
        return value.item() if hasattr(value, "item") else value

    have = dict(_rows(path, product, f'"{LAKE_KEY}", fid'))
    fields = set(pyogrio.read_info(path, layer=product)["fields"])
    columns = [c for c in out.columns if c in fields]
    sets = ", ".join(f'"{c}" = ?' for c in columns)
    con = sqlite3.connect(path)
    con.executemany(
        f'UPDATE "{product}" SET {sets} WHERE fid = ?',
        (
            [*(plain(v) for v in out.loc[lake, columns]), have[lake]]
            for lake in out.index
            if lake in have
        ),
    )
    con.commit()
    con.close()
    new = out[~out.index.isin(list(have))].rename_axis(LAKE_KEY).reset_index()
    if new.empty:
        return
    md, lakes = Vannforekomster(md_cache), FileLakes(work / LAKE_INDEX)
    geoms = [_polygon_part(lakes.geom(lake)) for lake in new[LAKE_KEY]]
    new["lake_navn"] = new["vannforekomstnavn"] = [
        md.lake_names[lake] for lake in new[LAKE_KEY]
    ]
    new["lake_eu"] = None  # empty on every lake of the file
    new["areal_vannforekomst_m2"] = [round(g.area, 1) for g in geoms]
    new["vannforekomstid"] = new[LAKE_KEY]
    _append(path, product, gpd.GeoDataFrame(new, geometry=geoms, crs=CRS), cache)


def renumber_shore(path, product, work) -> int:
    """Shore numbered alike in two regions gets ids of its own. Returns how many."""
    with open(work / NEW_SHORE_IDS, "rb") as f:
        rows = pickle.load(f)
    con = sqlite3.connect(path)
    con.execute("PRAGMA synchronous=OFF")
    con.execute(
        "CREATE TEMP TABLE new_id (folder TEXT, id TEXT, new TEXT, "
        "PRIMARY KEY (folder, id))"
    )
    con.executemany("INSERT INTO new_id VALUES (?, ?, ?)", rows)
    match = f'folder = "{product}"."{FOLDER}" AND id = "{product}"."segment_id"'
    n = con.execute(
        f'UPDATE "{product}" SET ("segment_id", "segment_id_vann-nett") = '
        f"(SELECT new, new FROM new_id WHERE {match}) "
        f"WHERE EXISTS (SELECT 1 FROM new_id WHERE {match})"
    ).rowcount
    con.commit()
    con.close()
    return n


def remove_copies(product, path, work, regions, source, cache, md_cache) -> int:
    """Take the second copy of everything out of one merged product."""
    with open(work / "copies.pkl", "rb") as f:
        copies = pickle.load(f)
    rivers = copies[copies["table"] == "official_centerline_10m"]
    shores = copies[copies["table"] == "lake_shore"]
    if product in ("official_centerline_10m", BUFFER):
        pairs = rivers[["region", "segment_id"]].values.tolist()
        return delete_rows(path, product, "segment_id", pairs)
    if product in ("lake_shore", LAKE_BUFFER):
        pairs = shores[["region", "segment_id"]].values.tolist()
        removed = delete_rows(path, product, "segment_id", pairs)
        # The double lines go before the renumbering, by the run's ids. Only the
        # ones with a zone are rows here; the rest are in the sums alone. They are
        # counted apart from the copies, which the README's count of those is for.
        with open(work / SHORE_PLAN, "rb") as f:
            plan = pickle.load(f)
        lines = [(r, i) for r, ids in plan["doubles"].items() for i in ids]
        gone = delete_rows(path, product, "segment_id", lines)
        (work / f"doubles_{product}.txt").write_text(str(gone))
        renumber_shore(path, product, work)
        return removed
    if product == "lake_polygon":
        removed = _one_row_per_lake(path, product)
        rebuild_lakes(path, product, work, regions, cache, md_cache)
        return removed
    if product == TREES:  # already left out of the parts
        with open(work / TREE_COPIES, "rb") as f:
            return sum(len(keys) for keys in pickle.load(f).values())
    if product == "official_centerline_region":
        aside = _side_rows(work, regions, rivers)
        per_region = _rivers_per_region(work, regions)
        return _rebuild_regions(path, product, aside, source, per_region, cache)
    if product in BLOCKS:
        aside = _side_rows(work, regions, rivers)
        return _rebuild_blocks(
            path, product, aside, rivers, work, regions, cache, md_cache
        )
    return 0


def pieces_in_lakes(a, b, lakes):
    """Which pieces (from a[i] to b[i]) have their middle inside one of the lakes."""
    import numpy as np
    import shapely

    tree = shapely.STRtree(np.array(lakes, dtype=object))
    middle = shapely.points((a + b) / 2)
    inside = np.zeros(len(middle), dtype=bool)
    inside[tree.query(middle, predicate="intersects")[0]] = True
    return inside


def outside_lakes(table):
    """Split a check table's gap into what is missing and what a lake delivers.

    A piece of a vannforekomst's line that no delivered segment lies on, but that
    lies in a lake, is not missing: the lake is delivered by lake_polygon and
    lake_shore, and the line in it carries the id of the lake or crosses another.
    On entry missing_km is every piece not there and in_lake_km those of them in a
    lake; on return missing_km is the pieces outside lakes alone, and present_pct
    is the share of the line outside lakes that is there. In-lake km stays in its
    own column, so nothing is hidden.
    """
    import numpy as np

    table["missing_km"] = (table["missing_km"] - table["in_lake_km"]).clip(lower=0.0)
    outside = table["line_km"] - table["in_lake_km"]
    share = table["present_km"] / outside.where(outside > 1e-6) * 100
    table["present_pct"] = share.where(outside > 1e-6, np.nan).round(3)
    return table


def present_line(c):
    """The comparison in one line, for the log and the README."""
    outside = c["line_km"] - c["in_lake_km"]
    # Only a partial run can have no line outside lakes at all.
    pct = f"{c['present_km'] / outside * 100:.3f} %" if outside > 1e-6 else "n/a"
    return (
        f"Present: {pct} of {outside:,.1f} km; "
        f"missing outside lakes {c['missing_km']:,.2f} km; "
        f"inside lakes {c['in_lake_km']:,.2f} km (delivered as lakes)"
    )


def compare_with_file(work, regions, output, md_path, md_cache, partial=False):
    """Is every line of the file in the result, and is anything else?

    A short piece of line is present if a delivered segment lies within 1 m of
    it. Writes one row per vannforekomst, and the pieces that are not there and
    not in a lake. A piece in a lake is not missing: a lake that is its own
    vannforekomst is delivered by its shore, so its midline is left out on
    purpose, and the line in any lake is delivered by the lake. Its length is
    kept in in_lake_km. A partial run answers only for the water bodies of its
    own regions.
    """
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    import pyogrio
    import shapely

    md = Vannforekomster(md_cache)
    numbers = {r[:3] for r in regions}
    wanted = np.array([not partial or i[:3] in numbers for i in md.ids])
    present = np.zeros(int(md.offset[-1]), dtype=bool)
    off_n, off_m = 0, 0.0
    for region in regions:
        file = work / "coverage" / f"{region}.npy"
        present[np.load(file)] = True
        note = file.with_suffix(".json")
        if note.exists():
            found = json.loads(note.read_text())
            off_n += found["off_line_segments"]
            off_m += found["off_line_m"]
    n = len(md.ids)
    line, there, gaps = np.zeros(n), np.zeros(n), []
    for start in range(0, len(md.wkb), 100_000):
        geoms = shapely.from_wkb(md.wkb[start : start + 100_000])
        xy, at = shapely.get_coordinates(
            shapely.segmentize(geoms, FINE_M), return_index=True
        )
        same = at[1:] == at[:-1]
        a, b, piece = xy[:-1][same], xy[1:][same], at[:-1][same]
        k = np.arange(len(xy) - 1)[same] - np.searchsorted(
            at, np.arange(len(geoms))
        )[piece]
        owner = md.owner[start + piece]
        keep = wanted[owner]
        a, b, piece, k, owner = a[keep], b[keep], piece[keep], k[keep], owner[keep]
        hit = present[md.offset[start + piece] + k]
        length = np.hypot(*(b - a).T)
        line += np.bincount(owner, weights=length, minlength=n)
        there += np.bincount(owner[hit], weights=length[hit], minlength=n)
        if (~hit).any():
            gaps.append((a[~hit], b[~hit], owner[~hit], length[~hit]))
    table = pd.DataFrame(
        {
            "vannforekomstid": md.ids,
            "vannforekomstnavn": md.names,
            "line_km": line / 1000,
            "present_km": there / 1000,
        }
    )
    table["missing_km"] = table["line_km"] - table["present_km"]
    table["in_lake_km"] = 0.0
    table["wanted"] = wanted
    if gaps:
        a, b, owner, length = (np.concatenate(x) for x in zip(*gaps, strict=True))
        lakes = pyogrio.read_dataframe(md_path, layer=MD_LAKES, columns=[])
        polys = list(shapely.make_valid(lakes.geometry.values))
        delivered = output / "lake_polygon.gpkg"
        if delivered.exists():
            own = pyogrio.read_dataframe(delivered, columns=[])
            polys += list(shapely.make_valid(own.geometry.values))
        in_lake = pieces_in_lakes(a, b, polys)
        table["in_lake_km"] = (
            np.bincount(owner[in_lake], weights=length[in_lake], minlength=n) / 1000
        )
        # Only the pieces outside lakes are missing, and only they are written.
        out = ~in_lake
        pieces = shapely.linestrings(np.stack([a[out], b[out]], axis=1))
        rows = pd.DataFrame({"owner": owner[out], "m": length[out]})
        merged = [
            (md.ids[o][:3], md.ids[o], md.names[o], group["m"].sum() / 1000,
             shapely.line_merge(shapely.multilinestrings(pieces[group.index])))
            for o, group in rows.groupby("owner")
        ]
        if merged:
            gdf = gpd.GeoDataFrame(
                pd.DataFrame(
                    merged,
                    columns=["region_nr", "vannforekomstid", "vannforekomstnavn",
                             "missing_km", "geometry"],
                ),
                geometry="geometry",
                crs=CRS,
            )
            local = work / "merged" / "vannforekomst_lines_missing.gpkg"
            local.unlink(missing_ok=True)
            pyogrio.write_dataframe(
                gdf,
                local,
                layer="vannforekomst_lines_missing",
                promote_to_multi=True,
                geometry_type="MultiLineString",
            )
            publish(local, output / local.name)
        else:
            # Nothing missing outside lakes: no file, and none left from before.
            (output / "vannforekomst_lines_missing.gpkg").unlink(missing_ok=True)
    table = outside_lakes(table[table.pop("wanted")].copy())
    table = table.sort_values("missing_km", ascending=False)
    table.round(4).to_csv(
        output / "vannforekomst_lines_check.csv", index=False, encoding="utf-8-sig"
    )
    no_line = table["present_km"] < 0.001
    all_in_lake = table["in_lake_km"] >= table["line_km"] - 0.001
    summary = {
        "vannforekomster": len(table),
        "line_km": float(table["line_km"].sum()),
        "present_km": float(table["present_km"].sum()),
        "missing_km": float(table["missing_km"].sum()),
        "in_lake_km": float(table["in_lake_km"].sum()),
        "complete": int((table["missing_km"] < 0.001).sum()),
        "with_gaps": int((table["missing_km"] >= 0.001).sum()),
        "only_by_shore": int((no_line & all_in_lake).sum()),
        "nothing_delivered": int((no_line & ~all_in_lake).sum()),
        "off_line_segments": off_n,
        "off_line_km": off_m / 1000,
    }
    if (output / "lake_polygon.gpkg").exists():
        delivered = {
            r[0]
            for r in _rows(output / "lake_polygon.gpkg", "lake_polygon", "lake_id")
        }
        lakes = {
            lake
            for lake in md.lake_names
            if not partial or lake[:3] in numbers
        }
        file_lakes = FileLakes(work / LAKE_INDEX)
        absent = sorted(lakes - delivered)
        pd.DataFrame(
            {
                "vannforekomstid": absent,
                "vannforekomstnavn": [md.lake_names[i] for i in absent],
                "km2": [round(file_lakes.geom(i).area / 1e6, 4) for i in absent],
            }
        ).sort_values("km2", ascending=False).to_csv(
            output / "lakes_not_delivered.csv", index=False, encoding="utf-8-sig"
        )
        summary["lakes_in_file"] = len(lakes)
        summary["lakes_not_delivered"] = len(absent)
        summary["lakes_not_in_file"] = len(delivered - set(md.lake_names))
    (work / "comparison.json").write_text(json.dumps(summary, indent=1))
    return summary


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def publish(local: Path, target: Path) -> None:
    part = target.with_suffix(".part")
    shutil.copyfile(local, part)
    if sha256(local) != sha256(part):
        part.unlink()
        raise RuntimeError(f"{local.name}: checksum differs after the copy")
    part.rename(target)
    local.unlink()


def add_segment_condition(target: Path, product: str, source: Path, layer: str) -> dict:
    """Give a merged buffer the segment's E_S_Kantvegetasjon_fjernet and _grad.

    Joined on segment_id_vann-nett from `layer` in the merged `source` file (the merged
    official_centerline_10m for buffer10, lake_shore for lake_buffer10), not made per
    region: the values are the ones in the delivered segment rows. A polygon whose
    segment is not in the layer gets null. Returns the row counts.
    """
    key = '"segment_id_vann-nett"'
    if not source.exists():
        raise RuntimeError(f"{product}: {source.name} must be merged first, it is not in {source.parent}")
    con = sqlite3.connect(target, uri=True)
    con.isolation_level = None
    # The R-tree triggers of an indexed file call functions sqlite3 does not have. They
    # fire on geometry changes only; the names must exist for the UPDATE to compile.
    for name in ("ST_IsEmpty", "ST_MinX", "ST_MaxX", "ST_MinY", "ST_MaxY"):
        con.create_function(name, 1, lambda g: None)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA cache_size=-4000000")
    con.execute(f"ATTACH DATABASE '{source.resolve().as_uri()}?mode=ro&immutable=1' AS src")
    have = {r[1] for r in con.execute(f'PRAGMA main.table_info("{product}")')}
    wanted = {r[1] for r in con.execute(f'PRAGMA src.table_info("{layer}")')}
    if "segment_id_vann-nett" not in have:
        raise RuntimeError(f"{target.name}: no segment_id_vann-nett")
    if not {"segment_id_vann-nett", *SEGMENT_VALUES} <= wanted:
        raise RuntimeError(f"{source.name}: needs segment_id_vann-nett and {sorted(SEGMENT_VALUES)}")
    if set(SEGMENT_VALUES) & have:
        raise RuntimeError(f"{target.name}: already has {sorted(set(SEGMENT_VALUES) & have)}")
    # The keys (48 million for the rivers) go in a file next to the target, on the work
    # disk, not in the temp directory.
    keys = target.with_name(target.stem + ".keys.db")
    keys.unlink(missing_ok=True)
    con.execute(f"ATTACH DATABASE '{keys}' AS kd")
    con.execute("PRAGMA kd.journal_mode=OFF")
    con.execute("PRAGMA kd.synchronous=OFF")
    con.execute("CREATE TABLE kd.segment_condition (sid TEXT PRIMARY KEY, fj REAL, gr TEXT) WITHOUT ROWID")
    try:
        con.execute(
            f"INSERT INTO kd.segment_condition SELECT {key}, E_S_Kantvegetasjon_fjernet, "
            f'E_S_Kantvegetasjon_grad FROM src."{layer}" WHERE {key} IS NOT NULL'
        )
    except sqlite3.IntegrityError as e:
        raise RuntimeError(f"{source.name}: a segment_id_vann-nett twice ({e})") from e
    con.execute(f'ALTER TABLE "{product}" ADD COLUMN E_S_Kantvegetasjon_fjernet REAL')
    con.execute(f'ALTER TABLE "{product}" ADD COLUMN E_S_Kantvegetasjon_grad TEXT')
    first, last = con.execute(f'SELECT min(fid), max(fid) FROM "{product}"').fetchone()
    for start in range(first or 1, (last or 0) + 1, 5_000_000):  # a transaction per 5 million
        con.execute("BEGIN")
        con.execute(
            f'UPDATE "{product}" SET (E_S_Kantvegetasjon_fjernet, E_S_Kantvegetasjon_grad) = '
            f"(SELECT fj, gr FROM kd.segment_condition WHERE sid = \"{product}\".{key}) "
            "WHERE fid >= ? AND fid < ?",
            (start, start + 5_000_000),
        )
        con.execute("COMMIT")
    rows, with_grad = con.execute(
        f'SELECT count(*), count(E_S_Kantvegetasjon_grad) FROM "{product}"'
    ).fetchone()
    con.close()
    keys.unlink(missing_ok=True)
    return {"rows": rows, "without_segment": rows - with_grad}


def drop_delivery_columns(path: Path, product: str) -> list[str]:
    """Drop DELIVERY_DROP[product] from a merged file. Returns the names dropped."""
    wanted = DELIVERY_DROP.get(product, [])
    if not wanted:
        return []
    con = sqlite3.connect(path)
    con.isolation_level = None
    # ALTER TABLE re-reads the schema, whose R-tree triggers name ST_* functions.
    for name in ("ST_IsEmpty", "ST_MinX", "ST_MaxX", "ST_MinY", "ST_MaxY"):
        con.create_function(name, 1, lambda g: None)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    have = {r[1] for r in con.execute(f'PRAGMA table_info("{product}")')}
    dropped = [c for c in wanted if c in have]
    for col in dropped:
        con.execute(f'ALTER TABLE "{product}" DROP COLUMN "{col}"')
    con.close()
    return dropped


def finish_product(product, work, regions, output, cache, md_cache):
    t0 = time.time()
    local = work / "merged" / f"{product}.gpkg"
    local.parent.mkdir(exist_ok=True)
    if product in ROLLUPS:
        result = finish_rollup(product, work, regions, local, cache, md_cache)
    else:
        parts = [part_path(work, product, r) for r in sorted(regions)]
        parts = [p for p in parts if p.exists()]
        if not parts:
            raise RuntimeError(f"{product}: no regions to merge")
        result = append_parts(product, parts, local, consume=product == TREES)
        source = parts_source(work)
        removed = remove_copies(
            product, local, work, regions, source, cache, md_cache
        )
        result["rows"] = _layer_rows(local, product)
        (work / f"removed_{product}.txt").write_text(str(removed))
        if product in SEGMENT_SOURCE:
            layer = SEGMENT_SOURCE[product]
            added = add_segment_condition(local, product, output / f"{layer}.gpkg", layer)
            (work / f"segment_condition_{product}.txt").write_text(str(added))
        spatial_index(local, product)
    dropped = drop_delivery_columns(local, product)
    (work / f"dropped_{product}.txt").write_text(", ".join(dropped))
    size = local.stat().st_size / 1e9
    publish(local, output / f"{product}.gpkg")
    # Sidecars stay: README.md reads them on later runs too.
    if product not in ROLLUPS:
        for region in regions:
            part_path(work, product, region).unlink(missing_ok=True)
    minutes = (time.time() - t0) / 60
    return f"{product}: {result['rows']:,} rows, {size:.1f} GB, {minutes:.1f} min"


def parts_source(work: Path) -> Path:
    return Path(json.loads((work / "inputs.json").read_text())["source"])


def double_lines_readme(c, work, output):
    """The README's section on the double lines, from what plan_shore counted."""
    left = {}
    for name in ("lake_shore", LAKE_BUFFER):
        count = work / f"doubles_{name}.txt"
        if count.exists():
            left[name] = int(count.read_text())
    out = [
        "",
        "## Double lines in the lake shore",
        "",
        "Where two mapped water surfaces meet, each draws its own shore, so a "
        "line through open water exists twice, one copy on top of the other. "
        "Grunnkart does it where it cuts a lake at a kommune border into two "
        "polygons that stop centimetres short of each other; the rule above "
        "does it where it brings the shore of two run lakes that met along a "
        "seam to one lake of the file. Such a piece is not shore. A 10 m "
        "piece is a double line when, at 0.2, 0.5 and 0.8 of its length, "
        "another piece of the same lake lies within 0.5 m; all three must, so "
        "a piece that only touches another at one point stays. A piece is not "
        "its own pair, and the next piece in its ring counts only where it lies "
        "nearer the station than the joint between them, as where the ring "
        "folds back on itself. Pieces are compared within a region. Both copies "
        "are removed, "
        "those with no zone too, which are in the lake's sums but not in "
        "lake_shore. innsjo_strandlinje_m, innsjo_segmenter, "
        "innsjo_datadekning, the shares, innsjo_Score_total and its grade are "
        "recomputed from the shore that remains, by the same rules as before. "
        "The rule is geometric only: a piece on a kommune border stays unless "
        "it passes it. The pieces are listed in lake_shore_double_lines.csv "
        "when that file is in the folder.",
        "",
        f"- Double lines removed: {c['rows']:,} rows, {c['km']:,.2f} km, in "
        f"{c['lakes']:,} lakes in {c['regions']:,} regions",
    ]
    if left:
        out.append(
            "- Of these, the ones with a zone, which were rows in "
            + " and ".join(f"{n} ({k:,} rows)" for n, k in left.items())
        )
    if c.get("left_whole"):
        out.append(
            f"- Lakes left whole, because every piece of them was a double line: "
            f"{', '.join(c['left_whole'])}"
        )
    return out


def write_readme(output: Path, work: Path) -> None:
    import pyogrio

    lines = [
        "# Merged delivery",
        "",
        f"Written {time.strftime('%F %T')} by merge_delivery.py.",
        "",
        "## Files",
        "",
        "| File | Layer | Geometry | Rows | Size GB |",
        "|---|---|---|---:|---:|",
    ]
    fields = []
    for path in sorted(output.glob("*.gpkg")):
        layer = pyogrio.list_layers(path)[0][0]
        info = pyogrio.read_info(path, layer=layer)
        size = path.stat().st_size / 1e9
        lines.append(
            f"| {path.name} | {layer} | {info['geometry_type']} | "
            f"{info['features']:,} | {size:.2f} |"
        )
        fields.append(f"### {path.name}\n")
        for name, dtype in zip(info["fields"], info["dtypes"], strict=True):
            about = {
                BUFFER: {**BUFFER_FIELDS, **SEGMENT_VALUES},
                LAKE_BUFFER: {**LAKE_BUFFER_FIELDS, **SEGMENT_VALUES},
                TREES: TREE_FIELDS,
            }.get(layer, ADDED)
            about = about.get(name)
            fields.append(f"- `{name}` ({dtype})" + (f": {about}" if about else ""))
        fields.append("")
    lines += [
        "",
        "## How eUVannomradeId was assigned",
        "",
        "Inside one vannområde: that one. Across a boundary: the one with the "
        "largest part of the feature. Touching none: the nearest. Trees use "
        "the apex; the two buffer layers and the vannområde roll-up use the "
        "segment's id.",
        "",
        "| Product | Rows | Across a boundary | Outside all | Furthest m |",
        "|---|---:|---:|---:|---:|",
    ]
    for product in (*DELIVERED, TREES):
        total = {}
        for sidecar in (work / "parts" / product).glob("*.json"):
            _add_counts(total, json.loads(sidecar.read_text()))
        if total:
            lines.append(
                f"| {product} | {total['rows']:,} | {total['boundary']:,} | "
                f"{total['outside']:,} | {total['outside_max_m']:,} |"
            )
    blanked = sum(
        json.loads(p.read_text()).get("blanked_unclassified", 0)
        for p in (work / "parts" / "official_centerline_10m").glob("*.json")
    )
    if blanked:
        lines += [
            "",
            f"official_centerline_10m: E_S_Score_total and K_S_prioritet are empty on "
            f"{blanked:,} rows that were not computed (Uklassifisert in "
            "E_S_Kantvegetasjon_A, and culverts the run did not compute). The "
            "roll-ups already leave these rows out, so no other table changes.",
        ]
    # The file --vannforekomster named (stamped into --work before any README).
    file = Path(
        json.loads((work / "inputs.json").read_text())["vannforekomster"]
    ).name
    lines += [
        "",
        "## Vannforekomster",
        "",
        "vannforekomstid, river_id_vann-nett and vannforekomstnavn come from "
        f"{file}, which is Vann-Nett. Each "
        "segment takes the line it runs on, by the longest overlap within 1 m. "
        "official_centerline_river has one row per vannforekomst; the blocks take "
        "the line holding most of them, and region rows carry no name. A lake the "
        "file has as a lake keeps its id; one it does not have takes the river "
        "running through it, and keeps its lake_id. region_folder and "
        "the blocks still follow the regions of the run.",
    ]
    done = work / "comparison.json"
    if done.exists():
        c = json.loads(done.read_text())
        lines += [
            "",
            "## Comparison with the file",
            "",
            "A short piece of line is present if a delivered segment lies within "
            "1 m of it. A piece that is not there but lies in a lake is not "
            "missing: lakes are delivered by lake_polygon and lake_shore, and the "
            "line in a lake carries the lake's own id or crosses another lake. "
            "Per vannforekomst in vannforekomst_lines_check.csv: missing_km and "
            "present_pct leave those pieces out, and in_lake_km holds their "
            "length. What is missing outside lakes is in "
            "vannforekomst_lines_missing.gpkg, with region_nr, the region "
            "owning the vannforekomst (the first part of its id).",
            "",
            f"- Lines in the file: {c['line_km']:,.1f} km in "
            f"{c['vannforekomster']:,} vannforekomster",
            f"- {present_line(c)}",
            f"- Vannforekomster complete outside lakes: {c['complete']:,}, of which "
            f"{c.get('only_by_shore', 0):,} are lakes whose whole line lies in the "
            "lake and which are delivered by their shore; with a gap outside "
            f"lakes: {c['with_gaps']:,}, of which with nothing delivered: "
            f"{c['nothing_delivered']:,}",
            f"- Delivered off every line: {c['off_line_segments']:,} segments, "
            f"{c['off_line_km']:,.2f} km",
        ]
        if "lakes_in_file" in c:
            lines.append(
                f"- Lakes in the file: {c['lakes_in_file']:,}; not delivered: "
                f"{c['lakes_not_delivered']:,}. Delivered lakes not in the file: "
                f"{c['lakes_not_in_file']:,}"
            )
    moved = work / "shore_moved.json"
    if moved.exists():
        c = json.loads(moved.read_text())
        zoned = ""
        if (output / "lake_shore.gpkg").exists():
            other_lake = "WHERE substr(segment_id, 1, length(lake_id)) <> lake_id"
            count = _rows(output / "lake_shore.gpkg", "lake_shore", "count(*)", other_lake)
            zoned = (
                f". Of the moved segments, {count[0][0]:,} have a zone, and so a "
                "row in lake_shore"
            )
        lines += [
            "",
            "## Which lake a shore belongs to",
            "",
            "The run gave each mapped water surface to one lake, so a surface "
            "holding two lakes of the file had all its shore under the first. "
            "Here a surface is shared by the lakes of the file that really lie in "
            "it: the run's own lake, and any lake whose overlap with the surface "
            "holds a circle of radius 20 m. One such lake takes all the surface's "
            "shore. Several share it, each 10 m piece going to the lake nearest "
            "through the water. A piece on a narrow arm of another lake takes that "
            "lake, if the one just chosen is over 50 m away. Shore never leaves "
            "its surface, so it cannot cross land.",
            "",
            f"- Moved: {c['rows']:,} shore segments, {c['km']:,.1f} km, to "
            f"{c['lakes_gaining']:,} lakes from {c['lakes_losing']:,}; "
            f"{c['arm_rows']:,} of them by the arm{zoned}",
            f"- Lakes that got their first shore, and so a row in lake_polygon: "
            f"{c['new_lakes']:,}",
            "",
            "Where the file draws no line between two lakes in a surface, the "
            "front lies halfway through the water between them, and is arbitrary "
            "by the length of that water: nationally under 100 m at 72 of the 100 "
            "pairs of neighbouring lakes, up to 1 km at 24, and more at 4, at most "
            "5.9 km on the Øyeren delta. Shore in water no outline holds, such as "
            "channels on a delta, goes by distance through the water alone.",
            "",
            "In lake_shore and lake_buffer10 a moved segment has the new lake's "
            "lake_id, vannforekomstid, vannforekomstnavn, river_id_vann-nett and "
            "elvenavn; in lake_shore also E_S_Barrier and "
            "E_S_Fragmentering, which belong to the lake. segment_id stays the "
            "run's, as on the rivers, so its first part can name another lake. "
            "stretch_id in lake_buffer10 stays the run's lake too.",
            "",
            "The rows of lake_polygon that change are rebuilt from the run's sums, "
            "with the segments that have no zone. areal_vannflate_m2 is the "
            "lake's part of a shared surface, divided the same way through the "
            "water, along a line that is not mapped. Shore moved by the arm brings "
            "no water with it.",
            "",
        ]
        if "doubles" in c:
            lines += double_lines_readme(c["doubles"], work, output)
        if (output / "lakes_not_delivered.csv").exists():
            lines += [
                "",
                "A lake of the file listed in lakes_not_delivered.csv has no shore "
                "in the run, and no moving of shore can give it any: its water "
                "polygon is under the run's limit of 5,000 m², it is under 30 % of "
                "a larger polygon that no lake reaches 30 % of, or Grunnkart maps "
                "it as sea.",
            ]
    removed = [
        f"- {product}: {(work / f'removed_{product}.txt').read_text()} rows fewer"
        for product in (*DELIVERED, BUFFER, LAKE_BUFFER, TREES)
        if (work / f"removed_{product}.txt").exists()
    ]
    if removed:
        lines += [
            "",
            "## Water delivered by two regions",
            "",
            "A segment with the same geometry in two regions is kept once, from "
            "the region that owns its vannforekomst in the run, and the rows it "
            "was part of are recomputed from the 10 m segments. A lake with two "
            "polygons in one region keeps the larger. A lake with shore in two "
            "regions has one row with the sums of both, and the second region's "
            "shore is numbered on from the first's in lake_shore and "
            "lake_buffer10. region_folder names the regions holding real length, "
            "longest first. Apart from the shore's lake above, all other rows are "
            "as delivered, including water drawn twice under two ids inside one "
            "region, and about 100 m in four places that two regions cut "
            "differently. A crown with the same apex in two regions is kept once: "
            "the copy nearest its stretch, at equal distance the region owning the "
            "stretch, then the lower number. It keeps that copy's outline, which "
            "the other region can have drawn differently. Crowns found with "
            "different apexes stay, so a few overlap. The vannområde table above "
            "counts rows before all this, and before the rows of the new lakes, "
            "except that the old block rows of a rebuilt river and the crowns not "
            "kept are already out of it.",
            "",
            *removed,
        ]
    lines += [
        "",
        "## Reading the layers together",
        "",
        "- buffer10 and lake_buffer10 have a polygon only where a bank has area "
        "in the 0-10 m ring. About 5 % of the river banks have none, most of "
        "them on segments that were not computed, and about 0.2 % of the whole "
        "shore, nearly all of it with no zone and so no row in lake_shore either.",
        "- innsjo_segmenter counts a lake's whole shore; lake_shore holds the "
        "segments that have a zone, so it can have fewer.",
        "- Join lake_shore and lake_buffer10 to lake_polygon on lake_id, not on "
        "the start of segment_id.",
        "- A row of official_centerline_region has no vannforekomstnavn: a region "
        "has no one name.",
        f"- single_trees: outlines are simplified by {TREE_SIMPLIFY_M} m, half a CHM "
        "pixel's diagonal, which takes the pixel stairs out. crown_area_m2 and "
        "mean_height are from the whole crown, before that. About 1 crown in 7,000 "
        "reads over 50 m, up to 1,100 m: steep terrain, where the CHM is wrong. They "
        "are kept as detected.",
        "",
        "## Fields",
        "",
        *fields,
    ]
    (output / "README.md").write_text("\n".join(lines), encoding="utf-8")


def code_id(root: Path) -> str:
    """A SHA-256 of the path and bytes of every file of kantsoner/ under --code."""
    files = sorted(
        p
        for p in (root / "kantsoner").rglob("*")
        if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"
    )
    h = hashlib.sha256()
    for p in files:
        h.update(str(p.relative_to(root)).encode() + b"\0" + p.read_bytes() + b"\0")
    return h.hexdigest()[:16]


def check_source(source: Path, regions: list, config) -> list[str]:
    """What keeps --source from being merged, if anything: with `config` (a
    run's YAML, --regions-from) every region of it must have a folder."""
    import re

    import yaml

    problems = []
    if config is not None:
        listed = yaml.safe_load(Path(config).read_text(encoding="utf-8"))["regions"]
        names = [r if isinstance(r, str) else r["name"] for r in listed]
        have = {f.split("_", 1)[-1].replace(":", "") for f in regions}
        lost = [
            n
            for n in names
            if re.sub(r"[\s/]+", "_", n.strip()).replace(":", "") not in have
        ]
        if lost:
            problems.append(
                f"{len(lost)} region(s) of {Path(config).name} have no folder in "
                f"{source}: {'; '.join(lost[:10])}"
            )
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        epilog=f"Products: {', '.join(ALL_PRODUCTS)}",
    )
    ap.add_argument("-o", "--output", type=Path, required=True)
    ap.add_argument(
        "--products",
        nargs="+",
        choices=ALL_PRODUCTS,
        default=list(DEFAULT_PRODUCTS),
        help="default: everything except single_trees",
    )
    ap.add_argument("-w", "--workers", type=int, default=3)
    ap.add_argument("-r", "--region", action="append", help="region folder")
    ap.add_argument(
        "-s",
        "--source",
        type=Path,
        required=True,
        help="a finished run: its base_dir (<RUN>), one folder per region",
    )
    ap.add_argument(
        "--regions-from",
        type=Path,
        help="a run config: refuse to start unless --source holds all its regions",
    )
    ap.add_argument(
        "--work", type=Path, required=True, help="a working folder on a local disk"
    )
    ap.add_argument(
        "--code",
        default=str(CODE),
        help="the folder holding the kantsoner package that made the run "
        "(<repo>/src); default: this repository's",
    )
    ap.add_argument(
        "--vannomrader",
        type=Path,
        required=True,
        help="<DATA>/geonorge/natur/vannomraader_miljoedir_25833_20260227.gpkg",
    )
    ap.add_argument(
        "--vannforekomster",
        type=Path,
        required=True,
        help="the Vann-Nett download the run was made on, with the layers and "
        "fields of Miljødirektoratet's file: "
        "<DATA>/vannnett/vannforekomster_20261026.gpkg",
    )
    ap.add_argument(
        "--chm",
        type=Path,
        help="for single_trees: <DATA>/hoydedata/chm_norway_1m_epsg25833_20260911.tif",
    )
    ap.add_argument(
        "--min-free-gb",
        type=float,
        default=25.0,
        help="a region's next step waits until this much memory is free",
    )
    ap.add_argument(
        "--stop-below-gb",
        type=float,
        default=16.0,
        help="the whole run stops if free memory falls below this",
    )
    ap.add_argument(
        "--max-worker-gb",
        type=float,
        default=12.0,
        help="hard limit on one worker's address space",
    )
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument(
        "--allow-missing-trees",
        action="store_true",
        help="merge single_trees though some regions have no file",
    )
    a = ap.parse_args()
    if TREES in a.products and a.chm is None:
        ap.error("single_trees needs --chm")

    sys.path.insert(0, a.code)
    if os.getpgrp() != os.getpid():  # a session leader already leads its group
        os.setpgrp()
    Path("/proc/self/oom_score_adj").write_text("1000")
    threading.Thread(
        target=stop_before_the_machine_does,
        args=(a.stop_below_gb,),
        daemon=True,
    ).start()
    a.work.mkdir(parents=True, exist_ok=True)
    lock = open(a.work / "lock", "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        print(f"another run is using {a.work}", flush=True)
        return 2
    # Before anything is stamped: a refused run must leave --work as it was.
    listed = a.region or sorted(
        p.name for p in a.source.glob("[0-9][0-9][0-9]_*") if p.is_dir()
    )
    problems = check_source(a.source, listed, None if a.region else a.regions_from)
    if problems:
        for line in problems:
            print(line, flush=True)
        return 2
    # Parts and caches are only good for the inputs that made them.
    inputs = json.dumps(
        {
            "source": str(a.source),
            "code": a.code,
            "vannomrader": str(a.vannomrader),
            "modified": a.vannomrader.stat().st_mtime_ns,
            "vannforekomster": str(a.vannforekomster),
            "file_modified": a.vannforekomster.stat().st_mtime_ns,
            "regions": sorted(a.region or []),
            "script": sha256(Path(__file__)),
            # What the path holds: changed code must not reuse the parts of
            # the code before.
            "code_id": code_id(Path(a.code))
            if (Path(a.code) / "kantsoner").is_dir()
            else None,
        }
    )
    stamp = a.work / "inputs.json"
    if stamp.exists() and stamp.read_text() != inputs:
        print(
            f"{a.work} was made by other inputs or another version of this "
            "script: use a new --work",
            flush=True,
        )
        return 2
    stamp.write_text(inputs)
    if a.region and any((a.output / f"{p}.gpkg").exists() for p in a.products):
        print("-r makes a partial product: use an output folder of its own")
        return 2
    for p in a.products:
        part_path(a.work, p, "x").parent.mkdir(parents=True, exist_ok=True)
    for sub in ("segments", "coverage", "merged", "lake_sums"):
        (a.work / sub).mkdir(parents=True, exist_ok=True)
    a.output.mkdir(parents=True, exist_ok=True)

    products = [
        p
        for p in a.products
        if a.overwrite or not (a.output / f"{p}.gpkg").exists()
    ]
    skipped = sorted(set(a.products) - set(products))
    if skipped:
        print(f"already in the output, skipped: {skipped}", flush=True)
    for p in products:
        need = SEGMENT_SOURCE.get(p)
        if need and need not in products and not (a.output / f"{need}.gpkg").exists():
            print(f"{p} takes the segment's condition from {need}: merge it too, or have it in {a.output}")
            return 2
    regions = a.region or sorted(
        p.name for p in a.source.glob("[0-9][0-9][0-9]_*") if p.is_dir()
    )
    done = a.work / "comparison.json"
    stale = not done.exists() or (
        (a.output / "lake_polygon.gpkg").exists()
        and "lakes_in_file" not in json.loads(done.read_text())
    )
    compare = "official_centerline_10m" in products or (
        "official_centerline_10m" in a.products
        and stale
        and all((a.work / "coverage" / f"{r}.npy").exists() for r in regions)
    )
    if not products and not compare:
        write_readme(a.output, a.work)
        return 0

    if TREES in products and not a.allow_missing_trees:
        missing = [
            r
            for r in regions
            if not (a.source / r / "delivery" / f"{TREES}.gpkg").exists()
        ]
        if missing:
            print(
                f"{len(missing)} regions have no {TREES}.gpkg, e.g. "
                f"{missing[:3]}. If the export has finished and they have no "
                f"crowns, use --allow-missing-trees.",
                flush=True,
            )
            return 1

    cache = a.work / "vannomrader.pkl"
    build_vo_index(a.vannomrader, cache)
    md_cache = a.work / "vannforekomster.pkl"
    build_md_index(a.vannforekomster, md_cache)
    tabell_8()

    def size(region):
        path = a.source / region / "geometry" / "zones_left.gpkg"
        return path.stat().st_size if path.exists() else 0

    regions = sorted(regions, key=size, reverse=True)  # largest first
    print(f"{len(regions)} regions, products: {products}", flush=True)
    copies = a.work / "copies.pkl"
    find_copies(a.source, sorted(regions), copies)

    failed = []
    ctx = mp.get_context("spawn")
    if TREES in products:
        # Simplified, parts or file come to about 3/4 of the source.
        (a.work / "merged" / f"{TREES}.gpkg").unlink(missing_ok=True)
        files = [a.source / r / "delivery" / f"{TREES}.gpkg" for r in regions]
        held = sum(p.stat().st_size for p in (a.work / "parts" / TREES).glob("*.gpkg"))
        need = sum(f.stat().st_size for f in files if f.exists()) - held
        free = shutil.disk_usage(a.work).free
        if free < need:
            print(f"{need / 1e9:.0f} GB needed in {a.work}, {free / 1e9:.0f} free")
            return 1
    if TREES in products and not (a.work / TREE_COPIES).exists():
        (a.work / "apex").mkdir(exist_ok=True)
        with ProcessPoolExecutor(
            a.workers,
            mp_context=ctx,
            initializer=limit_worker,
            initargs=(a.code, a.max_worker_gb),
        ) as pool:
            jobs = {
                pool.submit(
                    part_apex, a.source / r, a.work / "apex" / f"{r}.npz"
                ): r
                for r in regions
                if not (a.work / "apex" / f"{r}.npz").exists()
            }
            for i, job in enumerate(as_completed(jobs), 1):
                try:
                    print(time.strftime("%F %T"), f"{i}/{len(jobs)}", job.result(), flush=True)
                except Exception as e:  # noqa: BLE001 - report, keep the others going
                    failed.append(jobs[job])
                    print(f"FAILED {jobs[job]}: {e!r}", flush=True)
        if failed:
            print(f"{len(failed)} regions failed, nothing merged: {failed}", flush=True)
            return 1
        find_tree_copies(sorted(regions), a.work)
        print(time.strftime("%F %T"), "crowns in two regions found", flush=True)
    plan = a.work / SHORE_PLAN
    lake_cache = a.work / LAKE_INDEX
    build_lake_index(a.vannforekomster, lake_cache)
    if products != [TREES] and not plan.exists():
        (a.work / "shore").mkdir(exist_ok=True)
        with ProcessPoolExecutor(
            a.workers,
            mp_context=ctx,
            initializer=init_shore_worker,
            initargs=(a.code, str(lake_cache), str(copies), a.max_worker_gb),
        ) as pool:
            jobs = {
                pool.submit(
                    part_shore,
                    a.source / r,
                    a.work / "shore" / f"{r}.pkl",
                    a.min_free_gb,
                ): r
                for r in regions
                if not (a.work / "shore" / f"{r}.pkl").exists()
            }
            for i, job in enumerate(as_completed(jobs), 1):
                try:
                    print(time.strftime("%F %T"), f"{i}/{len(jobs)}", job.result(), flush=True)
                except Exception as e:  # noqa: BLE001 - report, keep the others going
                    failed.append(jobs[job])
                    print(f"FAILED {jobs[job]}: {e!r}", flush=True)
        if failed:
            print(f"{len(failed)} regions failed, nothing merged: {failed}", flush=True)
            return 1
        plan_shore(regions, a.source, a.work, md_cache)
        print(time.strftime("%F %T"), "shore planned", flush=True)
    counted = a.work / "shore_moved.json"
    if products != [TREES] and counted.exists():
        d = json.loads(counted.read_text())["doubles"]
        print(
            time.strftime("%F %T"),
            f"Double lines removed: {d['rows']:,} rows, {d['km']:,.2f} km, "
            f"in {d['lakes']:,} lakes",
            flush=True,
        )
    with ProcessPoolExecutor(
        a.workers,
        mp_context=ctx,
        initializer=init_worker,
        initargs=(
            a.code, str(cache), str(md_cache), str(copies), str(plan), a.max_worker_gb
        ),
    ) as pool:
        jobs = {
            pool.submit(
                do_region,
                r,
                a.source,
                a.work,
                products,
                a.chm,
                a.min_free_gb,
            ): r
            for r in regions
        }
        for i, job in enumerate(as_completed(jobs), 1):
            try:
                print(time.strftime("%F %T"), f"{i}/{len(jobs)}", job.result(), flush=True)
            except Exception as e:  # noqa: BLE001 - report, keep the others going
                failed.append(jobs[job])
                print(f"FAILED {jobs[job]}: {e!r}", flush=True)
    if failed:
        print(f"{len(failed)} regions failed, nothing merged: {failed}", flush=True)
        return 1

    from loguru import logger

    logger.remove()
    errors = []
    # The buffers take the segment's condition from the merged segment layer in the
    # output, so they are finished after it.
    for batch in (
        [p for p in products if p not in SEGMENT_SOURCE],
        [p for p in products if p in SEGMENT_SOURCE],
    ):
        with ThreadPoolExecutor(3) as pool:
            jobs = {
                pool.submit(
                    finish_product,
                    p,
                    a.work,
                    regions,
                    a.output,
                    cache,
                    md_cache,
                ): p
                for p in batch
            }
            for job in as_completed(jobs):
                try:
                    print(time.strftime("%F %T"), job.result(), flush=True)
                except Exception as e:  # noqa: BLE001
                    errors.append(jobs[job])
                    detail = getattr(e, "stderr", "") or ""
                    print(f"FAILED {jobs[job]}: {e!r} {detail[-500:]}", flush=True)
    if not errors and compare:
        summary = compare_with_file(
            a.work, regions, a.output, a.vannforekomster, md_cache, bool(a.region)
        )
        print(time.strftime("%F %T"), "compared with the file:", summary, flush=True)
        print(time.strftime("%F %T"), present_line(summary), flush=True)
    listed = a.work / "double_lines.csv"
    if not errors and "lake_shore" in products and listed.exists():
        local = a.work / "merged" / "lake_shore_double_lines.csv"
        shutil.copyfile(listed, local)
        publish(local, a.output / local.name)
    write_readme(a.output, a.work)
    print(f"finished, {len(errors)} failed: {errors}", flush=True)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
