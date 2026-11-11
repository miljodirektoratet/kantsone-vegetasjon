"""Column accessors that survive a missing input column.

Not every source is always configured — no grøftekart means no
`Grofte_tetthet` — and a segment side may have no riparian zone at all.
`pd.to_numeric(df.get(name))` silently returns a *scalar* when the column is
absent, which then breaks on the next Series method.
"""

from __future__ import annotations

import pandas as pd


def num(df: pd.DataFrame, name: str, default: float | None = 0.0) -> pd.Series:
    """Numeric column. `default=None` keeps NaN where the value is missing."""
    if name not in df.columns:
        return pd.Series(
            default if default is not None else float("nan"),
            index=df.index,
            dtype="float64",
        )
    out = pd.to_numeric(df[name], errors="coerce")
    return out if default is None else out.fillna(default)


def text(df: pd.DataFrame, name: str, default: str = "") -> pd.Series:
    if name not in df.columns:
        return pd.Series(default, index=df.index, dtype=object)
    return df[name].fillna(default).astype(object)


def arealbruk(df: pd.DataFrame, side: str) -> pd.Series:
    """The land use of one bank. A bank that exists but has none is åker.

    Report Tabell C.1: with no arealbruk the default is åker, the most
    conservative class (highest width requirement, lowest N_sub). A bank that
    does not exist (abroad, the missing side of a lake shore) or has no zone
    (nothing was measured) has no land use at all, not an unknown one, and
    stays empty.
    """
    from . import sides
    from .spec import constants as C

    a = text(df, f"K_{side}_B0-10_Arealbruk")
    a = a.mask(a.isin(("<NA>", "nan", "None")), "")
    return a.mask((a == "") & ~sides.unmeasured(df, side), C.AREALBRUK_UNKNOWN)


def is_lake(df) -> pd.Series:
    """True for a lake shore segment.

    Reads the source marker the geometry stage writes. Falls back to
    objektType so an attributes file made before the marker existed still
    classifies correctly rather than silently treating shore as riverbank.
    """
    import pandas as pd

    from .spec import constants as C

    if C.KILDE_COL in df.columns:
        return (
            df[C.KILDE_COL]
            .astype("string")
            .str.strip()
            .eq(C.KILDE_LAKE)
            .fillna(False)
            .astype(bool)
        )
    if "objektType" in df.columns:
        return (
            df["objektType"].astype("string").eq("Innsjøstrand").fillna(False)
        )
    return pd.Series(False, index=df.index)
