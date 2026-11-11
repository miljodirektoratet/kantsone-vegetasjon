"""Typed configuration and per-region path resolution.

Only two kinds of value live in ``config.yaml``: site paths (where the
national datasets are mounted) and run settings (which regions, which bands).
Method constants — coefficients, thresholds, class tables — are code, in
:mod:`kantsoner.spec.constants`.

Derived paths are code too. The per-region output tree is built here from
``base_dir`` + the region name rather than being restated in YAML.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Literal

import geopandas as gpd
import yaml
from loguru import logger
from pydantic import BaseModel, Field, field_validator

from .spec import constants as C


class FkbPaths(BaseModel):
    bygning_dir: Path
    bygnanlegg_dir: Path
    vann_dir: Path


class NvePaths(BaseModel):
    """NVE sikringstiltak, delivered split by geometry type."""

    flomsikring_line: Path | None = None
    flomsikring_punkt: Path | None = None
    erosjonssikring_line: Path | None = None
    erosjonssikring_punkt: Path | None = None
    dam_line: Path | None = None
    dam_punkt: Path | None = None


class Paths(BaseModel):
    base_dir: Path
    vassdragomr: Path
    # Two more layers a region can be cut from -- see RegionSpec. Optional:
    # a config that names only vassdragsområder never opens them.
    vannomraader: Path | None = None  # Miljødirektoratet vannområder
    regineenhet: Path | None = None  # NVE REGINE units, elvNavnHierarki
    chm: Path
    dtm: Path
    grunnkart_dir: Path
    elvenett: Path
    fkb: FkbPaths
    nve: NvePaths = Field(default_factory=NvePaths)
    vegvesen: Path | None = None  # SVV plastring/erosjonssikring
    grofter: Path | None = None
    fremmedart: Path | None = None
    elvevannforekomster: Path | None = None  # Vann-Nett water-body lines
    # Vann-Nett lake water bodies -- the reporting unit the lake track
    # aggregates to, and the filter deciding which water surfaces count as
    # a lake at all. Identity only: the polygons are generalised (perimeter
    # differs from the mapped shoreline by up to 3.4x), so the shoreline
    # buffered is grunnkart/FKB's, matched to these by largest overlap.
    # Layer VannforekomstInnsjø; absent means the lake track is skipped.
    vannforekomster: Path | None = None
    nhm_prosjektdekning: Path | None = None  # laser survey year -> Laser_aar
    # Norway's administrative outline. Only border catchments need it: it
    # keeps a buffer on a border river off the far bank, which every
    # Norwegian source stops short of. See geometry/landegrense.py.
    landegrense: Path | None = None


class RunSettings(BaseModel):
    """What varies between runs and machines.

    Nothing here changes the method — bands and detection thresholds are the
    specification and live in `spec/constants.py`, where a typo cannot
    silently desynchronise them from the attribute registry.
    """

    include_creeks: bool = True
    workers: int = 1
    # Segments per attributes batch. Bounds peak memory per region; whole
    # stretches are kept together, so a catchment smaller than this runs as
    # one. 50k keeps a giant's batch under ~6 GB.
    batch_segments: int = 50_000
    # Do not start another region while less than this is free. A coarse
    # brake against the giants stacking up; the interleaved order does the
    # real work.
    min_free_gb: float = 25.0
    # GB the pool's stages may hold together, by the per-stage estimates in
    # cli.STAGE_GB_PER_1000. A region waits for room before each stage
    # instead of only at its start, so the long, light attributes stages run
    # side by side while geometry and derived queue for memory. 0 = off: the
    # region-level brake above and the interleaved order only.
    mem_budget_gb: float = 0.0


RegionSource = Literal["vassdragomr", "vannomraader", "regine"]


class RegionSpec(BaseModel):
    """One region: a name and the national layer it is looked up in.

    ``vassdragomr`` (the default, so a bare string in the YAML still means a
    vassdragsområde) matches ``vassomr``; ``vannomraader`` matches the
    Miljødirektoratet ``navn``; ``regine`` takes every REGINE unit whose
    ``elvNavnHierarki`` has the river as one of its levels -- the river's
    whole catchment, tributaries included -- and dissolves them.
    """

    name: str
    source: RegionSource = "vassdragomr"

    model_config = {"frozen": True}


class Config(BaseModel):
    crs: str = C.CRS
    paths: Paths
    regions: list[RegionSpec]
    run: RunSettings = Field(default_factory=RunSettings)

    @field_validator("regions", mode="before")
    @classmethod
    def _non_empty(cls, v):
        if not v:
            raise ValueError("config.regions is empty — nothing to do")
        # A bare string is a vassdragsområde; a mapping names its source.
        return [{"name": x} if isinstance(x, str) else x for x in v]

    @field_validator("regions")
    @classmethod
    def _unique_names(cls, v):
        names = [r.name for r in v]
        dupes = sorted({n for n in names if names.count(n) > 1})
        if dupes:
            raise ValueError(f"config.regions repeats: {', '.join(dupes)}")
        return v

    @property
    def region_names(self) -> list[str]:
        return [r.name for r in self.regions]

    def region_spec(self, name: str) -> RegionSpec:
        for r in self.regions:
            if r.name == name:
                return r
        raise KeyError(f"{name!r} is not in config.regions")

    @classmethod
    def load(cls, path: Path | str) -> Config:
        with Path(path).open(encoding="utf-8") as fh:
            return cls.model_validate(yaml.safe_load(fh))

    def region_paths(self, region: str) -> RegionPaths:
        return RegionPaths.resolve(self, region)


@lru_cache(maxsize=4)
def _layer(path: str) -> gpd.GeoDataFrame:
    """A national region layer, read once per process."""
    return gpd.read_file(path)


def locate_region(
    cfg: Config, spec: RegionSpec
) -> tuple[str, gpd.GeoDataFrame]:
    """The region's number and its one-row outline, from its source layer.

    Reads only; `RegionPaths.resolve` is what writes the AOI. The number is
    the source's own id (vassomrnr, vannomradeId, or the REGINE
    nedbfeltHavVassdragNr of the river's catchment) so folders sort with the
    national numbering rather than an invented one.
    """
    if spec.source == "vassdragomr":
        gdf = _layer(str(cfg.paths.vassdragomr))
        match = gdf[gdf["vassomr"] == spec.name]
        nr_col = "vassomrnr"
    elif spec.source == "vannomraader":
        if cfg.paths.vannomraader is None:
            raise ValueError("paths.vannomraader is not configured")
        gdf = _layer(str(cfg.paths.vannomraader))
        match = gdf[gdf["navn"] == spec.name]
        nr_col = "vannomradeId"
    else:
        if cfg.paths.regineenhet is None:
            raise ValueError("paths.regineenhet is not configured")
        gdf = _layer(str(cfg.paths.regineenhet))
        levels = gdf["elvNavnHierarki"].fillna("").str.split("/")
        units = gdf[levels.apply(lambda p: spec.name in p)]
        if units.empty:
            raise ValueError(
                f"Region {spec.name!r}: no REGINE unit has it in elvNavnHierarki"
            )
        nrs = sorted(
            units["nedbfeltHavVassdragNr"].dropna().astype(str).unique()
        )
        if len(nrs) != 1:
            # The same river name in two catchments: nothing to dissolve
            # them by without a qualifier.
            raise ValueError(
                f"Region {spec.name!r} spans {len(nrs)} REGINE catchments "
                f"({', '.join(nrs)}); name is ambiguous"
            )
        match = gpd.GeoDataFrame(
            {"nedbfeltHavVassdragNr": nrs, "n_units": [len(units)]},
            geometry=[units.geometry.union_all()],
            crs=units.crs,
        )
        nr_col = "nedbfeltHavVassdragNr"
    if len(match) != 1:
        raise ValueError(
            f"Region {spec.name!r} matched {len(match)} rows in the "
            f"{spec.source} layer (expected exactly 1)"
        )
    nr = str(match.iloc[0][nr_col])
    out = gpd.GeoDataFrame(
        {"region": [spec.name], "region_nr": [nr], "source": [spec.source]},
        geometry=match.geometry.values,
        crs=match.crs,
    )
    return nr, out


def sanitise(name: str) -> str:
    """`Lærdalsvassdraget/Lærdalsfjorden` -> `Lærdalsvassdraget_Lærdalsfjorden`."""
    return re.sub(r"[\s/]+", "_", name.strip())


class RegionPaths(BaseModel):
    """Every output path for one region. Built, never configured."""

    region: str
    folder: str
    root: Path
    aoi: Path
    area_km2: float = 0.0  # of the outline; orders a parallel run

    model_config = {"arbitrary_types_allowed": True}

    @classmethod
    def resolve(cls, cfg: Config, region: str) -> RegionPaths:
        """Locate the region in its national layer and build its tree.

        Writes the AOI GeoPackage on first use. This was
        ``resolve_river_region.py``; it is not a pipeline step because it
        produces no analysis output, only the region's own extent.
        """
        nr, match = locate_region(cfg, cfg.region_spec(region))
        match = match.to_crs(cfg.crs)
        folder = f"{nr}_{sanitise(region)}"
        root = Path(cfg.paths.base_dir) / folder
        aoi = root / f"{folder}.gpkg"
        if not aoi.exists():
            root.mkdir(parents=True, exist_ok=True)
            from .io import save_gpkg

            save_gpkg(match, aoi)
            logger.info(f"Wrote AOI for {region} -> {aoi}")
        return cls(
            region=region,
            folder=folder,
            root=root,
            aoi=aoi,
            area_km2=float(match.geometry.area.sum() / 1e6),
        )

    # ── stage output trees ────────────────────────────────────────────────
    @property
    def geometry_dir(self) -> Path:
        return self.root / "geometry"

    @property
    def attributes_dir(self) -> Path:
        return self.root / "attributes"

    @property
    def aggregation_dir(self) -> Path:
        return self.root / "aggregation"

    # ── L1 geometry ───────────────────────────────────────────────────────
    @property
    def rivers(self) -> Path:
        return self.geometry_dir / "rivers.gpkg"

    @property
    def segments(self) -> Path:
        """The 10 m ELVIS segments — the delivery keys, cut in geometry."""
        return self.geometry_dir / "segments_10m.gpkg"

    @property
    def routes(self) -> Path:
        """Reconstructed channel line per ELVIS stretch."""
        return self.geometry_dir / "routes.gpkg"

    @property
    def claims(self) -> Path:
        """Each segment's piece of route — the segment/channel link."""
        return self.geometry_dir / "claims.gpkg"

    def zone_side(self, side: str) -> Path:
        """Per-segment slice of the riparian zone on one bank.

        Cut at the widest band; each extractor intersects it with its own
        ring. Written so the partition can be inspected — left/right
        assignment and bank ownership at river splits are the most
        error-prone geometry in the pipeline.
        """
        return self.geometry_dir / f"zones_{side}.gpkg"

    def zones(self, band: tuple[int, int]) -> Path:
        """The riparian band as a buffer outwards from the river edge.

        Named for its outer distance: a band always starts at the water, so
        `buffer20` says what it is more plainly than `zones_0-20m` did.
        """
        return self.geometry_dir / f"buffer{band[1]}.gpkg"

    # ── L2/L3 attributes ──────────────────────────────────────────────────
    # One file per stage, named for the stage that wrote it. Each is a
    # complete, self-contained layer of the 10 m segments with everything
    # computed up to that point — so a stage boundary is visible on disk and
    # any of them opens in QGIS on its own. The duplicated geometry (~7 MB a
    # copy) is the price of that explicitness.
    @property
    def attributes_file(self) -> Path:
        """L2: measured attributes."""
        return self.attributes_dir / "attributes.gpkg"

    @property
    def derived(self) -> Path:
        """L3: + system A and B formulas."""
        return self.attributes_dir / "derived.gpkg"

    # ── delivery ───────────────────────────────────────────────────────
    @property
    def official_10m(self) -> Path:
        return self.aggregation_dir / "official_centerline_10m.gpkg"

    def aggregated(self, level: str) -> Path:
        return self.aggregation_dir / f"official_centerline_{level}.gpkg"

    def mkdirs(self) -> None:
        for d in (self.geometry_dir, self.attributes_dir, self.aggregation_dir):
            d.mkdir(parents=True, exist_ok=True)
