#!/usr/bin/env python
"""Local GeoPackages of a geonorge Vann-Nett download, in the layout the pipeline reads.

    uv run python scripts/vannnett_from_geonorge.py \
        vannforekomster_0000_norge_25833_FILEGDB.gdb <DATA>/vannnett 20261026

Miljødirektoratet publishes the water bodies as one FileGDB with four layers
(rivers, lakes, coast, groundwater), each keyed on ``identifikasjon_lokalId``.
The pipeline wants two files from it:

* ``elvevannforekomster_<date>.gpkg``: the river lines (VannforekomstElv) as
  the only layer, with the id renamed to ``vannforekomstid`` and the name to
  ``vannforekomstnavn``. io.read_vector reads every layer of a file, and the
  segment tagging (stages._with_vannforekomst) and the ELVIS patch
  (geometry/vannnett_patch.py) read those column names.
* ``innsjovannforekomster_<date>.gpkg``: the lake polygons
  (VannforekomstInnsjø) as published; geometry/lakes.py reads that layer and
  its own field names.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pyogrio
from loguru import logger

RIVERS, LAKES = "VannforekomstElv", "VannforekomstInnsjø"
RENAME = {
    "identifikasjon_lokalId": "vannforekomstid",
    "navn": "vannforekomstnavn",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "gdb", type=Path, help="the geonorge vannforekomster FileGDB"
    )
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("date", help="download date, yyyymmdd, for the file names")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rivers = pyogrio.read_dataframe(args.gdb, layer=RIVERS).rename(
        columns=RENAME
    )
    if rivers["vannforekomstid"].duplicated().any():
        raise ValueError(f"{RIVERS}: vannforekomstid is not unique")
    out = args.out_dir / f"elvevannforekomster_{args.date}.gpkg"
    out.unlink(missing_ok=True)
    rivers.to_file(out, layer="elvevannforekomster", engine="pyogrio")
    logger.info(f"{out}: {len(rivers):,} river water bodies")

    lakes = pyogrio.read_dataframe(args.gdb, layer=LAKES)
    out = args.out_dir / f"innsjovannforekomster_{args.date}.gpkg"
    out.unlink(missing_ok=True)
    lakes.to_file(out, layer=LAKES, engine="pyogrio")
    logger.info(f"{out}: {len(lakes):,} lake water bodies")


if __name__ == "__main__":
    main()
