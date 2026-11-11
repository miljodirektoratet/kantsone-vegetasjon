"""Shared Kantvegetasjon status + cause rules.

Why riparian vegetation is missing, decided from the competing land-use
coverages. One module holds the rules because two places classify: the
per-segment evaluation in :mod:`derived.system_a`, and the per-block
recompute in :mod:`delivery.aggregate` -- the aggregated levels re-derive the
status and cause from summed areas rather than averaging percentages, and
keeping the rules here guarantees the 10 m and the aggregated levels
classify identically.

Kantvegetasjon_A is five-valued. `has_veg`/`potensiell` (the Vtot/Pot_veg
50 % threshold test) still decide the broad bucket, unchanged from before:

    has_veg                              -> OK_KANTVEG
    not has_veg, potensiell               -> Kantskog_fjernet or Nedbygd_kantsone
    not has_veg, not potensiell            -> Ikke_mulig (natural absence)

What's new is the middle bucket splits in two, by which *group* of
candidates covers more of the side: Jordbruk+Skogbruk ("agri") against
Bygg+Samferdsel+Konstruert ("built"). The group split happens on group
*totals* first, deliberately -- letting every one of the eight candidates
compete freely for "biggest single one" could hand a Manglende-status
segment to Berg or Vann just because a sliver of bedrock shares the slice,
which has no home in a taxonomy where Kantskog_fjernet/Nedbygd_kantsone only
exist under `potensiell`. Once the group is decided, `side_cause` picks the
specific label to report *within* that group only.

Ported from the ``minor_fixes`` branch (``utils/kantveg_core.py``), where the
original three-state rules were settled; only the plumbing differs there.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..spec import constants as C

YES, NO = "Ja", "Nei"
HUMAN_NA = "N/A"

# Candidate coverages, each typed by group, in tie-break priority order
# within their group (first listed wins a tie). Coverages handed to
# side_cause must align with whichever group's own order is passed in.
ARSAK_BUILT = ("Bygg", "Samferdsel", "Konstruert")
ARSAK_AGRI = ("Jordbruk", "Skogbruk")
ARSAK_NATURAL = ("SnaumarkFastfjell", "Vatmark", "Snoisbre")
ARSAK_LABELS = ARSAK_BUILT + ARSAK_AGRI + ARSAK_NATURAL
ARSAK_TYPES = (
    ("built",) * len(ARSAK_BUILT)
    + ("agri",) * len(ARSAK_AGRI)
    + ("natural",) * len(ARSAK_NATURAL)
)
# A same-group cause within this many percentage points of the dominant one
# is reported alongside it ("Bygg/Samferdsel").
ARSAK_WINDOW_PCT = 10.0


def side_cause(
    values, labels=ARSAK_LABELS, types=ARSAK_TYPES
) -> tuple[str, str]:
    """``(label, type)`` for one side's candidate coverages.

    The dominant cause plus any same-type cause within ``ARSAK_WINDOW_PCT``
    of it, joined with '/' highest-first. Ties resolve by candidate order
    (first in ``labels``). ``("", "")`` when every coverage is 0.

    ``labels``/``types`` default to every candidate (used for the Ikke_mulig
    natural-cause pick); pass a narrower pair -- e.g. just ``ARSAK_AGRI`` and
    its types -- to pick a label from within one already-decided group.
    """
    values = [0.0 if v is None or np.isnan(v) else float(v) for v in values]
    highest = max(values)
    if highest <= 0:
        return "", ""
    dom_i = max(range(len(values)), key=lambda i: (values[i], -i))
    dom_type = types[dom_i]
    picked = [
        (values[i], i, labels[i])
        for i in range(len(values))
        if highest - values[i] <= ARSAK_WINDOW_PCT and types[i] == dom_type
    ]
    picked.sort(key=lambda t: (-t[0], t[1]))
    return "/".join(label for _, _, label in picked), dom_type


def _side_group_totals(pcts):
    """(agri_total, built_total) row-wise sums over an (n, 8) pct array."""
    idx = {label: i for i, label in enumerate(ARSAK_LABELS)}
    agri = sum(pcts[:, idx[label]] for label in ARSAK_AGRI)
    built = sum(pcts[:, idx[label]] for label in ARSAK_BUILT)
    return agri, built


def infra_potensiell(cause_pcts) -> np.ndarray:
    """True where built cover alone makes the bank count as potential.

    Ground under asphalt has not stopped being capable of vegetation; it
    has been built on. PotVeg cannot see that -- it counts forest,
    agriculture and soil-covered Snaumark -- so without this a paved bank
    failed the potential test and was filed as natural absence. The cause
    reported stays Bygg / Samferdsel / Konstruert, which say more than a
    single "infrastructure" would. See C.INFRA_POTENSIELL_MIN_PCT.
    """
    vals = np.asarray(cause_pcts, dtype=float)
    _, built = _side_group_totals(vals)
    return built >= C.INFRA_POTENSIELL_MIN_PCT


def classify_sides(
    has_veg: pd.Series, potensiell: pd.Series, cause_pcts
) -> tuple[pd.Series, pd.Series, pd.Series]:
    """(status, arsak, flag) Series for one side of many segments.

    ``cause_pcts`` is an (n, 8) array-like aligned with ``ARSAK_LABELS``.
    """
    vals = np.asarray(cause_pcts, dtype=float)
    idx = {label: i for i, label in enumerate(ARSAK_LABELS)}
    agri_total, built_total = _side_group_totals(vals)

    status = pd.Series(
        C.IKKE_MULIG_KANTVEG, index=has_veg.index, dtype="object"
    )
    status[has_veg] = C.OK_KANTVEG
    manglende = ~has_veg & potensiell
    is_agri = manglende & (agri_total >= built_total)
    is_built = manglende & ~is_agri
    status[is_agri] = C.KANTSKOG_FJERNET
    status[is_built] = C.NEDBYGD_KANTSONE

    arsak = pd.Series("", index=has_veg.index, dtype="object")
    flag = pd.Series("", index=has_veg.index, dtype="object")
    flag[has_veg] = HUMAN_NA

    agri_i = [idx[label] for label in ARSAK_AGRI]
    built_i = [idx[label] for label in ARSAK_BUILT]
    nat_i = [idx[label] for label in ARSAK_NATURAL]
    for pos in np.flatnonzero(is_agri.to_numpy()):
        label, _ = side_cause(vals[pos, agri_i], ARSAK_AGRI, ("agri",) * 2)
        arsak.iat[pos] = label
        flag.iat[pos] = YES if label else ""
    for pos in np.flatnonzero(is_built.to_numpy()):
        label, _ = side_cause(vals[pos, built_i], ARSAK_BUILT, ("built",) * 3)
        arsak.iat[pos] = label
        flag.iat[pos] = YES
    for pos in np.flatnonzero((~has_veg & ~potensiell).to_numpy()):
        label, _ = side_cause(vals[pos, nat_i], ARSAK_NATURAL, ("natural",) * 3)
        arsak.iat[pos] = label
        flag.iat[pos] = NO if label else ""
    return status, arsak, flag


def combine_causes(left: pd.Series, right: pd.Series) -> pd.Series:
    """Both sides' cause labels, deduplicated, 'left, right' when they differ."""

    def one(a, b):
        seen = []
        for cause in (a, b):
            if cause and cause not in seen:
                seen.append(cause)
        return ", ".join(seen)

    return pd.Series(
        [one(a, b) for a, b in zip(left, right, strict=True)],
        index=left.index,
        dtype="object",
    )


def combine_status(left: pd.Series, right: pd.Series) -> pd.Series:
    """Combined Kantvegetasjon_A: OK only if both sides are; Ikke_mulig if
    either side is (it dominates); otherwise whichever of
    Kantskog_fjernet/Nedbygd_kantsone a negative side reports, agreement
    when both report one, 'left, right' when they differ.

    A blank side is one with no data rather than a bad score -- a bank
    abroad, clipped away at the national border -- and is left out of the
    vote entirely. ``classify_sides`` never produces a blank, so this only
    fires where something upstream deliberately cleared a side.
    """

    def one(*vals):
        vals = [v for v in vals if v]
        if not vals:
            return ""
        if all(v == C.OK_KANTVEG for v in vals):
            return C.OK_KANTVEG
        if any(v == C.IKKE_MULIG_KANTVEG for v in vals):
            return C.IKKE_MULIG_KANTVEG
        seen = []
        for s in vals:
            if s != C.OK_KANTVEG and s not in seen:
                seen.append(s)
        return ", ".join(seen)

    return pd.Series(
        [one(a, b) for a, b in zip(left, right, strict=True)],
        index=left.index,
        dtype="object",
    )


def combine_flags(
    left: pd.Series, right: pd.Series, summary_status: pd.Series
) -> pd.Series:
    """Combined human flag: "N/A" when the summary is OK, agreement when the
    sides agree (ignoring N/A), "Ja/Nei" when human and natural mix."""

    def one(a, b, status):
        if status == C.OK_KANTVEG:
            return HUMAN_NA
        vals = {f for f in (a, b) if f and f != HUMAN_NA}
        if not vals:
            return ""
        if vals == {YES}:
            return YES
        if vals == {NO}:
            return NO
        return f"{YES}/{NO}"

    return pd.Series(
        [
            one(a, b, s)
            for a, b, s in zip(left, right, summary_status, strict=True)
        ],
        index=left.index,
        dtype="object",
    )
