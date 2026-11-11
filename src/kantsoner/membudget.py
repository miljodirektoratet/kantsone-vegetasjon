"""A region worker's handle on the parallel run's memory budget.

cli._run_region binds it when the run has a budget (cli._Budget). Code that
knows it is about to peak -- geometry's zoning pool -- takes a surcharge
for the duration, so the stage as a whole can be booked at what it holds
for hours rather than at what it holds for minutes. Unbound (a region run
on its own) every call is a no-op.
"""

from __future__ import annotations

from contextlib import contextmanager

from loguru import logger

_bound = None  # (budget proxy, region, ELVIS stretches, free-GB floor)


def bind(budget, region, stretches, floor):
    global _bound
    _bound = (budget, region, stretches, floor)


@contextmanager
def surcharge(gb_per_1000: float, what: str):
    """Hold extra budget while the block runs.

    Never blocks for good: when the budget has no room, the block goes
    ahead untracked as soon as that much memory is actually free above the
    floor -- the stage's own booking already holds the region's base, and a
    worker that waits forever inside a stage holds that booking forever.
    """
    if _bound is None:
        yield
        return
    from .cli import _free_gb

    budget, region, stretches, floor = _bound
    need = gb_per_1000 * stretches / 1000
    key = f"{region} [{what}]"
    held = False
    while True:
        if budget.try_acquire(key, need, 60):
            held = True
            break
        if _free_gb() > need + floor:
            logger.info(
                f"{what}: no budget for {need:.0f} GB; {_free_gb():.0f} GB free, going ahead"
            )
            break
    try:
        yield
    finally:
        if held:
            budget.release(key)
