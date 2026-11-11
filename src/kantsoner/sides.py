"""Which banks of a segment carry data, for the formulas that combine both.

Almost every column is per bank and needs nothing from here. The handful
that fold two banks into one number do: a minimum (the width shortfall), a
mean (every system-B score), a sum of two halves (Flomvern/Forbygning), a
worst-of (Kantvegetasjon_A). Each of those assumes two banks exist, and
each is wrong on a border river, where one bank is foreign ground the
clip in :mod:`kantsoner.geometry.landegrense` has removed.

The distinction that matters is between a bank measured as bare and a bank
not measured at all. Both arrive as zero or NaN; only the second may be
skipped. ``absent`` is the one place that tells them apart, and it answers
False when the border column is missing -- an old attributes file, or a
catchment nowhere near a border -- so nothing changes for the 99.9 % of
the network this was never about.
"""

from __future__ import annotations

import pandas as pd

SIDES = ("L", "R")
YES, NO = "Ja", "Nei"
# Per-side marker written onto the segments layer by geometry/landegrense.py:
# this bank's buffer was (almost) entirely outside Norway and has been
# removed, so it carries no measurement at all rather than a measurement of
# zero. Declared here rather than there so that derived/ can read it without
# importing anything geometric -- see tests/test_layering.py.
OUTSIDE_COLS = {"L": "K_L_Utenfor_Norge", "R": "K_R_Utenfor_Norge"}
# The second way a bank can carry no measurement: there is no such bank.
# A lake shore is one-sided -- the zone lies outside the water and there is
# nothing on the other side of the line to measure. Kept separate from
# OUTSIDE_COLS rather than reusing it, because that one also drives
# K_S_Grenseelv, and marking every lake in Norway as a border river would
# be a plain untruth in the delivery.
NO_BANK_COLS = {"L": "K_L_Ingen_bredd", "R": "K_R_Ingen_bredd"}
# Set when either bank is. The one column a reader of the delivery needs to
# spot a border reach without knowing which side of it is abroad.
BORDER_COL = "K_S_Grenseelv"


def _flag(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    # Plain bool, not pandas' nullable one: these feed .mask(), np.where and
    # arithmetic, and a pd.NA leaking into any of those is a silent wrong
    # answer rather than an error.
    return (
        df[col].astype("string").str.strip().eq(YES).fillna(False).astype(bool)
    )


def absent(df: pd.DataFrame, side: str) -> pd.Series:
    """True where this bank carries no measurement and must be skipped.

    Either because the ground is foreign and was clipped away, or because
    the feature has only one bank at all. The formulas that fold two banks
    into one do not care which: both mean "do not count this side".
    """
    return _flag(df, OUTSIDE_COLS[side]) | _flag(df, NO_BANK_COLS[side])


def present(df: pd.DataFrame, side: str) -> pd.Series:
    return ~absent(df, side)


def no_zone(df: pd.DataFrame, side: str) -> pd.Series:
    """True where this bank has no zone area at all, so nothing was measured.

    The bank exists (it is not abroad, not the missing side of a lake shore)
    but the geometry stage could not cut a zone for it: the first or last
    segment of a stretch often loses one side to the bank split. Every
    percentage is then 0 of nothing, which is "not measured", not "nothing
    there". The test is the widest band, the one aggregation reads to decide
    whether a segment was computed (stages.aggregation); an attributes file
    without the area column answers False, so an old file reads as before.
    """
    for col in (f"K_{side}_B0-20_Atot_m2", f"K_{side}_B0-10_Atot_m2"):
        if col in df.columns:
            # NaN is unknown, not "no zone": only a recorded 0 counts.
            return pd.to_numeric(df[col], errors="coerce") <= 0
    return pd.Series(False, index=df.index)


def unmeasured(df: pd.DataFrame, side: str) -> pd.Series:
    """A bank that carries no measurement: absent, or a bank with no zone."""
    return absent(df, side) | no_zone(df, side)


def n_measured(df: pd.DataFrame) -> pd.Series:
    """Banks with a measurement, 1 to 2, as n_present but without the empty ones."""
    n = sum((~unmeasured(df, side)).astype(int) for side in SIDES)
    return n.clip(lower=1)


def n_present(df: pd.DataFrame) -> pd.Series:
    """How many banks a segment has data for: 2 normally, 1 on a border.

    Never 0: a segment with both banks abroad is not in a Norwegian
    catchment at all, and dividing by it would be worse than the answer it
    would replace.
    """
    n = sum(present(df, side).astype(int) for side in SIDES)
    return n.clip(lower=1)


def mask(
    df: pd.DataFrame, side: str, values: pd.Series, *, measured: bool = False
) -> pd.Series:
    """``values`` with the absent bank's entries blanked.

    With ``measured`` a bank with no zone is blanked too: the two-bank means
    of the scores skip it instead of counting it as 0.

    NaN for numbers and "" for text -- both of which the combining code
    already skips, so a border segment reduces to its Norwegian bank
    without any of it learning what a border is.
    """
    out = values.copy()
    gone = unmeasured(df, side) if measured else absent(df, side)
    if not gone.any():
        return out
    blank = float("nan") if pd.api.types.is_numeric_dtype(out) else ""
    return out.mask(gone.to_numpy(), blank)
