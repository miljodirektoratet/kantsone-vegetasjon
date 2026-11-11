"""L3 — attributes derivable from L2 without touching a raster again.

Every module here exposes ``compute(df) -> df`` and is pure tabular: no I/O,
no geometry, no CRS. That is what makes the whole system-B method testable
with a ten-row fixture and no GIS data, and it is enforced by
``tests/test_layering.py``.

Order matters — ``width`` produces ``w_eff`` for ``filtrering``, and
``filtrering``/``temperatur``/``naering``/``kontinuitet`` all feed ``scores``.
"""

from __future__ import annotations

import pandas as pd

from . import (
    filtrering,
    kantveg_grad,
    kontinuitet,
    naering,
    scores,
    system_a,
    temperatur,
    usikkerhet,
    width,
)

PIPELINE = (
    system_a,
    kantveg_grad,
    width,
    filtrering,
    temperatur,
    naering,
    kontinuitet,
    scores,
    usikkerhet,
)


# Columns the derived formulas cannot do without. `_columns.num` defaults a
# missing column to 0.0 — right for genuinely optional sources (no grøftekart
# means no Grofte_tetthet), but it would also silently zero a *required* input
# and hand back plausible scores computed from nothing. These are checked.
REQUIRED = (
    "K_{side}_B0-3_Vtot",
    "K_{side}_B0-10_Vtot",
    "K_{side}_B0-20_Vtot",
    "K_{side}_B0-10_V0p5-1",
    "K_{side}_B0-10_V1-5",
    "K_{side}_B0-10_V5-10",
    "K_{side}_B0-10_V10pluss",
    "K_{side}_B0-10_Vhmed",
    "K_{side}_B0-10_Vp95",
    "K_{side}_B0-10_s",
    "K_{side}_B0-10_SkogJordbruk",
    "K_{side}_B0-10_PotVeg",
    # system_a classifies at 0-20 m too, for stretches beregningsbredde
    # picks that band for; a silent zero default there would misclassify
    # every wide stretch the same way a missing 0-10 PotVeg/SkogJordbruk
    # would misclassify a narrow one.
    "K_{side}_B0-20_SkogJordbruk",
    "K_{side}_B0-20_PotVeg",
)

# Side-less inputs, same rule.
REQUIRED_COMBINED = (
    "K_S_bekk_bredde",
    "K_S_bekk_azimut",
)


def missing_inputs(df: pd.DataFrame) -> list[str]:
    """Required columns absent from ``df``."""
    per_side = [
        col
        for template in REQUIRED
        for col in (template.format(side=s) for s in ("L", "R"))
        if col not in df.columns
    ]
    return per_side + [c for c in REQUIRED_COMBINED if c not in df.columns]


def compute_all(
    df: pd.DataFrame, *, allow_missing: bool = False
) -> pd.DataFrame:
    """Run every derived module in dependency order.

    Raises when a required L2 input is absent, rather than quietly treating it
    as zero. Pass ``allow_missing=True`` to compute anyway — the affected
    scores will be wrong, so this is for debugging only.
    """
    absent = missing_inputs(df)
    if absent and not allow_missing:
        raise ValueError(
            f"{len(absent)} required L2 column(s) missing, so the derived "
            f"scores would be computed from zeros: {', '.join(absent)}. "
            "Fix the primary stage, or pass allow_missing=True to proceed "
            "knowing the output is wrong."
        )
    for module in PIPELINE:
        df = module.compute(df)
    return df


__all__ = [
    "PIPELINE",
    "compute_all",
    *(m.__name__.split(".")[-1] for m in PIPELINE),
]
