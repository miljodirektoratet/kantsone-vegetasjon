#!/usr/bin/env python
"""Write ELVIS patched with Vann-Nett (kantsoner/geometry/vannnett_patch.py),
completed with every Vann-Nett line the patch left out
(kantsoner/geometry/vannnett_complete.py).

    uv run python scripts/patch_elvis.py --elvis <DATA>/elvenett/elvenett_nve.gpkg \
        --vannnett <DATA>/vannnett/elvevannforekomster_20261026.gpkg \
        --regine <DATA>/nve_felt/regineenhet_nve.gpkg --out elvenett_vn_patched.gpkg
    uv run python scripts/patch_elvis.py ... --bbox 320000 7010000 330000 7020000 --out test.gpkg
    uv run python scripts/patch_elvis.py ... --complete elvenett_vn_patched_20261027.gpkg --out new.gpkg

``--complete`` takes an existing patched network and only completes it, so
every row it does not touch stays byte-identical to that file.

The output has one layer, ``elvenett_nve``, with ELVIS's own columns plus
``patch`` and ``patch_vn``, so ``paths.elvenett`` can point straight at it
(io.read_vector reads every layer of a file, so the QA layers go to a
separate file, <out>_qa.gpkg).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pyogrio
from loguru import logger
from shapely.geometry import box

from kantsoner.geometry import vannnett_complete as vc
from kantsoner.geometry import vannnett_patch as vp


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--elvis",
        type=Path,
        help="NVE's ELVIS river network, <DATA>/elvenett/elvenett_nve.gpkg "
        "(not read with --complete)",
    )
    ap.add_argument(
        "--vannnett",
        type=Path,
        required=True,
        help="<DATA>/vannnett/elvevannforekomster_20261026.gpkg "
        "(scripts/vannnett_from_geonorge.py)",
    )
    ap.add_argument(
        "--regine",
        type=Path,
        required=True,
        help="NVE's REGINE units, for the river names: "
        "<DATA>/nve_felt/regineenhet_nve.gpkg",
    )
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--complete",
        type=Path,
        help="an existing patched network: only complete it with Vann-Nett",
    )
    ap.add_argument(
        "--bbox", type=float, nargs=4, help="xmin ymin xmax ymax, for testing"
    )
    args = ap.parse_args()
    if not (args.complete or args.elvis):
        ap.error("give --elvis, or --complete")

    t0 = time.time()
    bbox = tuple(args.bbox) if args.bbox else None
    elvis = pyogrio.read_dataframe(args.complete or args.elvis, bbox=bbox)
    if args.complete:
        # an integer field with NULLs (OBJECTID of the additions) reads as
        # float; give it back its type so the file is written as it was
        info = pyogrio.read_info(args.complete)
        for f, dt in zip(info["fields"], info["dtypes"], strict=True):
            if (
                dt.startswith("int")
                and f in elvis
                and elvis[f].dtype.kind == "f"
            ):
                elvis[f] = elvis[f].round().astype("Int64")
    vn = pyogrio.read_dataframe(args.vannnett, bbox=bbox)
    if bbox:
        # a water body is one multipart feature reaching far outside a test
        # box; keep only its parts that touch the box, as for ELVIS
        vn = vn.explode(index_parts=False)
        vn = vn[vn.intersects(box(*bbox))]
    regine = pyogrio.read_dataframe(
        args.regine,
        bbox=bbox,
        columns=["vassdragsnummer", "elvNavnHierarki"],
    )
    logger.info(
        f"read {len(elvis):,} ELVIS, {len(vn):,} Vann-Nett features in {time.time() - t0:.0f} s"
    )
    if args.complete:
        res = vc.complete(elvis, vn, regine)
    else:
        res = vp.patch(elvis, vn, regine)
        done = vc.complete(res.network, vn, regine)
        res = vp.Result(
            done.network,
            {**res.qa, **{f"completion_{k}": v for k, v in done.qa.items()}},
            {
                **res.stats,
                **{f"completion_{k}": v for k, v in done.stats.items()},
            },
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    qa_path = args.out.with_name(args.out.stem + "_qa.gpkg")
    for p in (args.out, qa_path):
        p.unlink(missing_ok=True)
    res.network.to_file(args.out, layer="elvenett_nve", engine="pyogrio")
    for name, gdf in res.qa.items():
        if len(gdf):
            gdf.to_file(qa_path, layer=name, engine="pyogrio")
    res.stats["seconds"] = round(time.time() - t0)
    args.out.with_name(args.out.stem + "_stats.json").write_text(
        json.dumps(res.stats, indent=1)
    )
    logger.info(f"wrote {args.out} ({len(res.network):,} rows): {res.stats}")


if __name__ == "__main__":
    main()
