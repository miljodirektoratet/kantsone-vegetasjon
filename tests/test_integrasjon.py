"""The whole chain on the small test area in data/testomraade/.

Runs exactly what README section 6 runs, as separate commands, into a
temporary folder: `kantsoner run` (all four stages, crown export on), then
scripts/merge_delivery.py for the national products and for single_trees.
Checks that every delivered file is there, has rows, and carries computed
scores, not only blanks. About half a minute.
"""

import os
import subprocess
import sys
from pathlib import Path

import pyogrio
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "testomraade"
CONFIG = ROOT / "config" / "testomraade.yaml"
REGION = "012_Testområde_Korsrudtjernet"
BIN = Path(sys.executable).parent

PRODUCTS = {
    "official_centerline_10m": "Score_total",
    "official_centerline_100m": "Score_total",
    "official_centerline_500m": "Score_total",
    "official_centerline_1000m": "Score_total",
    "official_centerline_river": "Score_total",
    "official_centerline_vannomraade": "Score_total",
    "official_centerline_region": "Score_total",
    "lake_shore": None,
    "lake_polygon": None,
    "buffer10": None,
    "lake_buffer10": None,
    "single_trees": "max_height",
}  # product: a column (suffix) that must carry values

pytestmark = pytest.mark.skipif(
    not DATA.exists(), reason="data/testomraade/ is not present"
)


def _run(args, tmp_path):
    env = {**os.environ, "KANTSONER_EXPORT_SINGLETREE": "1"}
    done = subprocess.run(  # noqa: S603 -- our own commands, no shell
        args, cwd=ROOT, env=env, capture_output=True, text=True, check=False
    )
    (tmp_path / "log.txt").open("a").write(done.stdout + done.stderr)
    assert done.returncode == 0, done.stdout[-3000:] + done.stderr[-3000:]


@pytest.fixture(scope="module")
def delivery(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("integrasjon")
    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    cfg["paths"]["base_dir"] = str(tmp / "run")
    config = tmp / "config.yaml"
    config.write_text(yaml.safe_dump(cfg, allow_unicode=True), encoding="utf-8")

    _run([BIN / "kantsoner", "run", "-c", config], tmp)
    merge = [
        sys.executable,
        ROOT / "scripts" / "merge_delivery.py",
        "-s", tmp / "run",
        "-o", tmp / "leveranse",
        "--work", tmp / "merge",
        "-w", "1",
        "--min-free-gb", "0",
        "--stop-below-gb", "0",
        "--vannomrader",
        DATA / "vannomraader" / "vannomraader_miljoedir_25833_20260227.gpkg",
        "--vannforekomster", DATA / "vannnett" / "vannforekomster_20261026.gpkg",
    ]  # fmt: skip
    _run(merge, tmp)
    chm = DATA / "hoydedata" / "chm_norway_1m_epsg25833_20260911.tif"
    _run([*merge, "--products", "single_trees", "--chm", chm], tmp)
    return tmp


def test_region_run_writes_every_stage(delivery):
    region = delivery / "run" / REGION
    for f in (
        "geometry/segments_10m.gpkg",
        "attributes/attributes.gpkg",
        "attributes/derived.gpkg",
        "aggregation/official_centerline_10m.gpkg",
        "delivery/official_centerline_10m.gpkg",
        "delivery/single_trees.gpkg",
    ):
        assert (region / f).exists(), f


@pytest.mark.parametrize("product", PRODUCTS)
def test_delivered_product(delivery, product):
    path = delivery / "leveranse" / f"{product}.gpkg"
    assert path.exists()
    gdf = pyogrio.read_dataframe(path)
    assert len(gdf) > 0
    assert gdf.crs is not None and gdf.crs.to_epsg() == 25833
    assert gdf.geometry.notna().all()
    suffix = PRODUCTS[product]
    if suffix:
        # Levels above 10 m prefix their columns: 100m_Score_total, ...
        cols = [c for c in gdf.columns if c.endswith(suffix)]
        assert cols, f"{product}: no *{suffix} column"
        assert gdf[cols[0]].notna().any(), f"{product}: {cols[0]} is all blank"


def test_10m_has_measured_values(delivery):
    gdf = pyogrio.read_dataframe(
        delivery / "leveranse" / "official_centerline_10m.gpkg",
        read_geometry=False,
    )
    for col in ("K_L_B0-10_Vtot", "K_S_Laser_aar", "K_S_prioritet"):
        assert gdf[col].notna().mean() > 0.9, col
