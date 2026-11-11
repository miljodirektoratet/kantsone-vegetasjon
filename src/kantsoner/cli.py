"""The single entry point.

One command per layer, named after its module, plus ``run`` to chain them.
There is no L0 command: bounded reads live in :mod:`kantsoner.io`, a library
not a stage — nothing is clipped or materialised, so there is nothing to run.
"""

from __future__ import annotations

import ctypes
import gc
import json
import multiprocessing as mp
import os
import sys
import threading
import time
from collections import deque
from concurrent.futures import (
    FIRST_COMPLETED,
    ProcessPoolExecutor,
    wait,
)
from multiprocessing.managers import BaseManager
from pathlib import Path

import click
from loguru import logger

from . import __version__
from .config import Config, RegionPaths
from .spec import attributes as A

STAGES = ("geometry", "attributes", "derived", "aggregation")


def _setup_logging(verbose: bool) -> None:
    logger.remove()
    logger.add(
        sys.stderr,
        level="DEBUG" if verbose else "INFO",
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <7}</level> | {message}",
    )


def _regions(cfg: Config, selected: tuple[str, ...]) -> list[RegionPaths]:
    names = list(selected) if selected else cfg.region_names
    unknown = [n for n in names if n not in cfg.region_names]
    if unknown:
        raise click.BadParameter(f"Not in config.regions: {', '.join(unknown)}")
    out, bad = [], []
    for n in names:
        try:
            out.append(cfg.region_paths(n))
        except Exception as exc:
            # A national run must not die on one name the catchment layer
            # cannot resolve; report it in the final summary instead.
            logger.warning(f"cannot resolve region {n!r}: {exc}")
            bad.append(n)
    if bad:
        logger.warning(
            f"skipping {len(bad)} unresolvable region(s): {', '.join(bad)}"
        )
    return out


_config_option = click.option(
    "-c",
    "--config",
    "config_path",
    type=click.Path(exists=True, path_type=Path),
    default="config/config.yaml",
    show_default=True,
    help="Path to config.yaml.",
)
_region_option = click.option(
    "-r",
    "--region",
    "regions",
    multiple=True,
    help="Region to process. Repeatable; defaults to every region in the config.",
)
_dry_run_option = click.option(
    "--dry-run", is_flag=True, help="Print the resolved plan and exit."
)


_allow_missing_option = click.option(
    "--allow-missing",
    is_flag=True,
    help=(
        "Compute derived scores even when required L2 inputs are absent. "
        "Those scores will be WRONG (missing inputs read as zero); for "
        "debugging only."
    ),
)


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__)
@click.option("-v", "--verbose", is_flag=True, help="Debug-level logging.")
@click.pass_context
def cli(ctx, verbose):
    """Kantsoner — riparian zone morphology from remote sensing."""
    ctx.ensure_object(dict)
    _setup_logging(verbose)


def _stage_command(name: str, help_text: str):
    """Build one stage subcommand; they share every option and the region loop."""

    @cli.command(name=name, help=help_text)
    @_config_option
    @_region_option
    @_dry_run_option
    @_allow_missing_option
    def _cmd(config_path, regions, dry_run, allow_missing):
        cfg = Config.load(config_path)
        for paths in _regions(cfg, regions):
            if dry_run:
                click.echo(f"{name}: {paths.region} -> {paths.root}")
                continue
            paths.mkdirs()
            logger.info(f"── {name}: {paths.region}")
            _run_stage(name, cfg, paths, allow_missing=allow_missing)

    return _cmd


def _run_stage(name: str, cfg: Config, paths: RegionPaths, **kw) -> None:
    """Call a stage, passing only the options it accepts."""
    import inspect

    from . import stages

    fn = getattr(stages, name)
    accepted = inspect.signature(fn).parameters
    fn(cfg, paths, **{k: v for k, v in kw.items() if k in accepted})


_stage_command(
    "geometry",
    "L1 — routes, kantsoner and buffers on the 10 m ELVIS segments.",
)
_stage_command("attributes", "L2 — every measured attribute, per batch.")
_stage_command(
    "derived", "L3 — system A + B formulas plus the topological attributes."
)
_stage_command("aggregation", "L4 — roll up to every reporting level.")


def _run_region(
    config_path,
    region,
    selected,
    allow_missing,
    budget=None,
    stretches=0,
    resume=False,
):
    """One region, start to finish. Runs in a worker process.

    With a ``budget`` (a _Budget proxy) each stage first takes its memory
    estimate from the shared pool and hands it back when done.
    """
    import os

    # Each worker gets its own GDAL; let the pool provide the parallelism
    # rather than every reader spawning threads and cache of its own.
    os.environ.setdefault("GDAL_NUM_THREADS", "1")
    os.environ.setdefault("GDAL_CACHEMAX", "256")

    cfg = Config.load(config_path)
    paths = cfg.region_paths(region)
    paths.mkdirs()
    if budget is not None:
        from . import membudget

        membudget.bind(budget, region, stretches, cfg.run.min_free_gb)
    for stage in selected:
        # A region can legitimately hold no network at all (Røst: islands
        # with no mapped rivers). Geometry then writes empty layers, and
        # every later stage would only trip over the empty input.
        if stage != "geometry" and paths.segments.exists():
            from pyogrio import read_info

            if not read_info(paths.segments).get("features"):
                logger.warning(
                    f"── {stage}: {region} skipped — 0 segments "
                    "(no ELVIS network in the region)"
                )
                continue
        if resume and stage == "geometry" and _geometry_done(paths):
            logger.info(f"── geometry: {region} done before, kept")
            continue
        if budget is not None:
            need = _stage_gb(stage, stretches, cfg.run.mem_budget_gb)
            budget.acquire(region, need)
            logger.info(
                f"── {stage}: {region} ({need:.0f} GB of the budget, "
                f"{budget.used():.0f}/{cfg.run.mem_budget_gb:.0f} in use)"
            )
        else:
            logger.info(f"── {stage}: {region}")
        try:
            _run_stage(stage, cfg, paths, allow_missing=allow_missing)
        finally:
            if budget is not None:
                budget.release(region)
        _give_back_memory()


@cli.command()
@_config_option
@_region_option
@_dry_run_option
@_allow_missing_option
@click.option(
    "--from", "from_stage", type=click.Choice(STAGES), default=STAGES[0]
)
@click.option("--to", "to_stage", type=click.Choice(STAGES), default=STAGES[-1])
@click.option(
    "--resume",
    is_flag=True,
    help="Keep a region's geometry if a previous run finished it.",
)
@click.option(
    "-w",
    "--workers",
    type=int,
    default=None,
    help="Catchments to run in parallel. Overrides run.workers in the config.",
)
def run(
    config_path,
    regions,
    dry_run,
    allow_missing,
    from_stage,
    to_stage,
    workers,
    resume,
):
    """Run every stage, or a slice of them, for each region."""
    cfg = Config.load(config_path)
    lo, hi = STAGES.index(from_stage), STAGES.index(to_stage)
    if lo > hi:
        raise click.BadParameter(
            f"--from {from_stage} comes after --to {to_stage}"
        )
    selected = STAGES[lo : hi + 1]
    chosen = _regions(cfg, regions)

    if dry_run:
        for paths in chosen:
            click.echo(
                f"{paths.region}: {' -> '.join(selected)} in {paths.root}"
            )
        return

    n = workers if workers is not None else cfg.run.workers
    os.environ["KANTSONER_OUTER_WORKERS"] = str(max(n, 1))
    if n <= 1:
        failed = []
        for paths in chosen:
            try:
                _run_region(
                    config_path,
                    paths.region,
                    selected,
                    allow_missing,
                    resume=resume,
                )
            except Exception as exc:
                logger.exception(f"{paths.region} failed: {exc}")
                failed.append(paths.region)
        _report(failed)
        return

    ctx = mp.get_context("spawn")
    budget, manager, sizes = None, None, {}
    if cfg.run.mem_budget_gb > 0:
        # Largest first by network size: area misleads (Pasvikelva is 19,000
        # km2 with 8,000 stretches, Drammen 17,000 km2 with 67,000). The
        # budget, not the order, keeps the giants' heavy stages apart.
        sizes = _stretch_counts(cfg, chosen)
        order = sorted(chosen, key=lambda p: -sizes[p.region])
        manager = _BudgetManager(ctx=ctx)
        manager.start()
        budget = manager.Budget(cfg.run.mem_budget_gb, cfg.run.min_free_gb)
        logger.info(
            f"Running {len(order)} regions on {n} workers, largest network "
            f"first, stages admitted within {cfg.run.mem_budget_gb:.0f} GB"
        )
    else:
        # Largest first: catchments span 62 to 43,000 km2, so a big one
        # started late leaves the pool idle while it finishes alone.
        order = _interleaved(cfg, chosen)
        logger.info(
            f"Running {len(order)} regions on {n} workers, "
            f"giants interleaved with small ones, "
            f"next region admitted at >{cfg.run.min_free_gb} GB free"
        )
    failed = []
    queue = deque(order)
    running: dict = {}
    # A fresh process per region: memory a finished region's worker could
    # not hand back to the OS must not count against the next one.
    with ProcessPoolExecutor(
        max_workers=n, mp_context=ctx, max_tasks_per_child=1
    ) as pool:
        while queue or running:
            # Admit work lazily. Submitting everything up front let the six
            # largest catchments start together and their geometry state
            # alone (~20 GB each) exceeded RAM; the OOM killer then took one
            # worker and the pool with it. A region is only started when
            # there is room for it to grow -- or, with a budget, started
            # freely and held at each stage until there is.
            while (
                queue
                and len(running) < n
                and (budget is not None or _memory_ok(cfg, running))
            ):
                p = queue.popleft()
                fut = pool.submit(
                    _run_region,
                    str(config_path),
                    p.region,
                    selected,
                    allow_missing,
                    budget,
                    sizes.get(p.region, 0),
                    resume,
                )
                running[fut] = p.region
                logger.info(f"▶ {p.region} ({len(running)} running)")
            if not running:
                time.sleep(30)
                continue
            done, _ = wait(running, timeout=60, return_when=FIRST_COMPLETED)
            for fut in done:
                region = running.pop(fut)
                try:
                    fut.result()
                    logger.info(f"✓ {region}")
                except Exception as exc:
                    logger.error(f"✗ {region}: {exc}")
                    failed.append(region)
                    if budget is not None:
                        # A killed worker never reached its finally.
                        budget.release(region)
    if manager is not None:
        manager.shutdown()
    _report(failed)


def _free_gb() -> float:
    with Path("/proc/meminfo").open() as fh:
        for line in fh:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024**2
    return float("inf")


def _memory_ok(cfg, running) -> bool:
    """Room to start one more region.

    The floor applies once anything is running; the first region always
    starts. A region ramps up over minutes, so this is a coarse brake, not
    a guarantee — the interleaved order is what keeps the giants apart.
    """
    if not running:
        return True
    return _free_gb() > cfg.run.min_free_gb


# Memory a region's stage holds, in GB per 1,000 ELVIS stretches in the
# region, on top of BUDGET_BASE_GB. Measured on a national run: Glomma (104,975
# stretches) peaked at 83.6 GB in geometry (zoning, 8 workers) and held 2-21
# GB through attributes; derived peaked at 93 GB, aggregation at 69 GB.
# The 150 smallest regions (<= 14,600 stretches), four at a time, peaked at
# 18 GB together. Scheduling estimates, not part of the method.
STAGE_GB_PER_1000 = {
    "geometry": 0.45,  # + elvis.ZONING_GB_PER_1000 while zoning
    "attributes": 0.25,
    "derived": 0.9,
    "aggregation": 0.7,
}
BUDGET_BASE_GB = 1.5
# A stage that has waited this long for memory reserves its room: smaller
# ones may still overtake it, but only into what is left beside it, so a
# giant's derived is not starved by a stream of small regions slipping into
# the gaps -- and does not idle the pool while it cannot fit yet.
STARVE_S = 3600


def _stage_gb(stage, stretches, total):
    """A stage's estimate, capped at the budget so it can always run alone."""
    need = BUDGET_BASE_GB + STAGE_GB_PER_1000[stage] * stretches / 1000
    return min(need, total)


class _Budget:
    """Memory tokens shared by a pool's workers, served by a manager.

    Smaller requests may overtake a waiting one that does not fit yet;
    once it has waited STARVE_S they may only use the room it leaves. Nothing is admitted below ``floor`` GB
    actually free unless the pool holds nothing: the estimates are only
    estimates, and other jobs may share the machine.
    """

    def __init__(self, total, floor):
        self.total, self.floor = total, floor
        self.held: dict = {}
        self.waiting: dict = {}
        self.cv = threading.Condition()

    def _fits(self, key, need):
        used = sum(self.held.values())
        if not self.held:
            return True
        if used + need > self.total or _free_gb() < self.floor:
            return False
        now = time.monotonic()
        starved = [
            (since, k, n)
            for k, (since, n) in self.waiting.items()
            if now - since > STARVE_S
        ]
        if not starved:
            return True
        _, oldest, reserved = min(starved)
        return oldest == key or used + need + reserved <= self.total

    def acquire(self, key, need):
        self.try_acquire(key, need, None)

    def try_acquire(self, key, need, timeout):
        """Take ``need`` GB for ``key``; False if it did not fit in time."""
        deadline = None if timeout is None else time.monotonic() + timeout
        with self.cv:
            self.waiting.setdefault(key, (time.monotonic(), need))
            while not self._fits(key, need):
                left = 30 if deadline is None else deadline - time.monotonic()
                if left <= 0:
                    del self.waiting[key]
                    return False
                # Woken by every release; the timeout re-reads free memory.
                self.cv.wait(min(left, 30))
            del self.waiting[key]
            self.held[key] = need
            self.cv.notify_all()
            return True

    def release(self, key):
        with self.cv:
            self.held.pop(key, None)
            self.cv.notify_all()

    def used(self):
        with self.cv:
            return sum(self.held.values())


class _BudgetManager(BaseManager):
    pass


_BudgetManager.register("Budget", _Budget)


def _stretch_counts(cfg, chosen) -> dict:
    """ELVIS stretches per region (midpoint inside the outline), cached.

    The size that drives a region's memory and time; the whole national
    file reads in seconds.
    """
    cache = Path(cfg.paths.base_dir) / "region_sizes.json"
    counts = json.loads(cache.read_text()) if cache.exists() else {}
    todo = [p for p in chosen if p.region not in counts]
    if todo:
        import geopandas as gpd
        import pyogrio
        import shapely

        lines = pyogrio.read_dataframe(str(cfg.paths.elvenett), columns=[])
        mid = shapely.line_interpolate_point(
            lines.geometry.values, 0.5, normalized=True
        )
        tree = shapely.STRtree(mid)
        for p in todo:
            aoi = gpd.read_file(p.aoi).to_crs(lines.crs).geometry.union_all()
            counts[p.region] = int(len(tree.query(aoi, predicate="contains")))
        cache.write_text(json.dumps(counts, ensure_ascii=False, indent=0))
    return counts


def _geometry_done(paths) -> bool:
    """Geometry's last write is the right-bank zones, staged then renamed."""
    last = paths.zone_side("right")
    return (
        last.exists()
        and paths.segments.exists()
        and last.stat().st_mtime >= paths.segments.stat().st_mtime
    )


def _give_back_memory():
    """Return freed heap to the OS between stages (glibc keeps it otherwise)."""
    gc.collect()
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except OSError:
        pass


def _interleaved(cfg, chosen):
    """Largest, smallest, second-largest, second-smallest, ...

    The long poles still start at the front so they never finish alone at
    the end, but no more than about half the pool is ever giants at once —
    strictly largest-first put every giant in the machine together.
    """
    big = deque(sorted(chosen, key=lambda p: -p.area_km2))
    out = []
    while big:
        out.append(big.popleft())
        if big:
            out.append(big.pop())
    return out


def _report(failed):
    if failed:
        raise SystemExit(f"Failed regions: {', '.join(failed)}")


@cli.group()
def spec():
    """Inspect the attribute registry."""


@spec.command("list")
@click.option("--system", type=click.Choice(["A", "B"]), default=None)
@click.option("--layer", type=click.Choice(["primary", "derived", "network"]))
@click.option("--producer", default=None)
def attributes_list(system, layer, producer):
    """Print every declared attribute — diff this against the Excel."""
    rows = A.select(system=system, layer=layer, producer=producer)
    width = max((len(c) for a in rows for c in a.columns()), default=10)
    for attr in rows:
        for col in attr.columns():
            click.echo(
                f"{col:<{width}}  {attr.layer:<8} {attr.producer:<14} "
                f"{attr.system:<2} {attr.aggregation}"
            )
    click.echo(
        f"\n{len(rows)} attributes, {len(A.columns(system=system, layer=layer, producer=producer))} columns"
    )


@cli.group()
def config():
    """Inspect configuration."""


@config.command("check")
@_config_option
@_region_option
def config_check(config_path, regions):
    """Resolve every source and covering tile, and open each one.

    This is what makes a missing kommune tile or an unmounted share fail in
    seconds instead of two hours into `primary`.
    """
    from .io import tile_index

    cfg = Config.load(config_path)
    problems = []

    single = {
        "chm": cfg.paths.chm,
        "dtm": cfg.paths.dtm,
        "elvenett": cfg.paths.elvenett,
        "vassdragomr": cfg.paths.vassdragomr,
        "vannomraader": cfg.paths.vannomraader,
        "regineenhet": cfg.paths.regineenhet,
        "vegvesen": cfg.paths.vegvesen,
        "grofter": cfg.paths.grofter,
        "fremmedart": cfg.paths.fremmedart,
        "elvevannforekomster": cfg.paths.elvevannforekomster,
        "nhm_prosjektdekning": cfg.paths.nhm_prosjektdekning,
        "landegrense": cfg.paths.landegrense,
        **{f"nve.{k}": v for k, v in cfg.paths.nve.model_dump().items()},
    }
    for name, path in single.items():
        if path is None:
            click.echo(f"  - {name:<20} not configured (optional)")
        elif Path(path).exists():
            click.echo(f"  ✓ {name:<20} {path}")
        else:
            click.echo(f"  ✗ {name:<20} MISSING {path}")
            problems.append(name)

    dirs = {
        "grunnkart": cfg.paths.grunnkart_dir,
        "fkb.bygning": cfg.paths.fkb.bygning_dir,
        "fkb.bygnanlegg": cfg.paths.fkb.bygnanlegg_dir,
        "fkb.vann": cfg.paths.fkb.vann_dir,
    }
    for name, root in dirs.items():
        try:
            entries = tile_index(Path(root))
            click.echo(f"  ✓ {name:<20} {len(entries)} layers indexed")
        except Exception as exc:
            click.echo(f"  ✗ {name:<20} {exc}")
            problems.append(name)

    # Resolve regions without writing anything — a check should not create
    # directories or AOI files as a side effect.
    from .config import locate_region

    for region in regions or cfg.region_names:
        try:
            spec = cfg.region_spec(region)
            nr, outline = locate_region(cfg, spec)
        except Exception as exc:
            click.echo(f"  ✗ region {region}: {exc}")
            problems.append(region)
            continue
        km2 = outline.to_crs(cfg.crs).geometry.area.sum() / 1e6
        click.echo(f"  ✓ region {region}  [{spec.source} {nr}, {km2:,.0f} km2]")

    if problems:
        raise SystemExit(f"{len(problems)} problem(s): {', '.join(problems)}")
    click.echo("\nAll sources resolved.")


def main() -> None:
    cli(obj={})


if __name__ == "__main__":
    main()
