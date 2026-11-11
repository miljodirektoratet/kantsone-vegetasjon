"""The 1-5 'grad av endring' scale.

Every aggregated parameter is restated on the Vannportalen five-point scale
(5 = best, 1 = worst) as a ``*_final`` column. See the Aggregation tab of the
v03 workbook.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Percent-of-length bins shared by Kantvegetasjon, Forbygning and Flomvern.
# Read as: < 5 % affected -> 5 (best); >= 80 % -> 1 (worst).
_PCT_BINS = [-np.inf, 5, 20, 50, 80, np.inf]
_PCT_SCORES = [5, 4, 3, 2, 1]

# Barrier count per river: none -> 5, 1-2 -> 3, 3 -> 2, >=4 -> 1.
# Bins are left-closed, so the edges are the first count of each class.
_BARRIER_BINS = [-np.inf, 1, 3, 4, np.inf]
_BARRIER_SCORES = [5, 3, 2, 1]

# System-B total score, 0-100. High is good here, where the bins above measure
# damage, so the grades run the other way: Sállir tabell 8, from Svært dårlig
# below 0.20 to Svært god from 0.80, on Vannportalen's 5 = Svært god.
_SCORE_B_BINS = [-np.inf, 20, 40, 60, 80, np.inf]
_SCORE_B_SCORES = [1, 2, 3, 4, 5]


def _bin(values, bins, scores) -> pd.Series:
    s = pd.to_numeric(pd.Series(values), errors="coerce")
    out = pd.cut(s, bins=bins, labels=scores, right=False, ordered=False)
    return pd.Series(out, index=s.index).astype("Float64")


def pct_to_grade(values) -> pd.Series:
    """Percent of length affected (0-100) -> 1-5. Used by Kantvegetasjon,
    Forbygning and Flomvern."""
    return _bin(values, _PCT_BINS, _PCT_SCORES)


def barrier_count_to_grade(values) -> pd.Series:
    """Number of barriers on a watercourse -> 1-5."""
    return _bin(values, _BARRIER_BINS, _BARRIER_SCORES)


def score_b_to_grade(values) -> pd.Series:
    """System-B aggregated score (0-100) -> 1-5."""
    return _bin(values, _SCORE_B_BINS, _SCORE_B_SCORES)


def fragmentering_to_grade(labels) -> pd.Series:
    """The Fragmentering class labels -> 1-5."""
    mapping = {
        "Svært god": 5,
        "God": 4,
        "Moderat": 3,
        "Dårlig": 2,
        "Svært dårlig": 1,
    }
    return pd.Series(labels).map(mapping).astype("Float64")
