"""usikkerhet_flagg: which default values a segment's scores rest on.

Report appendix C: every calculation that uses a default value is flagged, so
the user can tell which segments are well covered by data and which are not.
The flag is an integer bitmask (Tabell C.2); a segment that used no default
is 0 and is tested bit by bit in a GIS (``usikkerhet_flagg & 2 != 0``).

Only the defaults this pipeline really applies are set:

* bit 1 (value 2): arealbruk missing, åker used. A bank with a zone but no
  land use read (columns.arealbruk). A bank with no zone was not measured and
  sets nothing.
* bit 2 (value 4): løv/bar missing, løv used. A bank with a tall canopy and
  no species (filtrering.species_defaulted).

The other bits of Tabell C.2 are never set here, because the pipeline does
not use that default: it computes K_score, CHM height, azimuth and
continuity instead of defaulting them, and reads a missing slope as 0 %, not
the report's 5 % (bit 0), and a missing ditch map as no ditches without
saying so (bit 8; the grøfte column is 0 for both).

Pure tabular, like the rest of L3.
"""

from __future__ import annotations

import pandas as pd

from .. import sides
from ..columns import text
from . import filtrering

BIT_AREALBRUK = 2
BIT_LOEV = 4


def compute(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    flag = pd.Series(0, index=out.index, dtype="int64")
    for side in sides.SIDES:
        exists = ~sides.unmeasured(out, side)
        blank = text(out, f"K_{side}_B0-10_Arealbruk").isin(
            filtrering.UNKNOWN_TEXT
        )
        flag |= (blank & exists).astype("int64") * BIT_AREALBRUK
        flag |= (
            filtrering.species_defaulted(out, side) & exists
        ).astype("int64") * BIT_LOEV
    out["E_S_usikkerhet_flagg"] = flag
    return out
