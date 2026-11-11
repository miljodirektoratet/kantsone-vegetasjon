"""The attribute registry — a transcription of the specification workbooks.

Every per-segment attribute is declared exactly once here. Extraction
and aggregation both read from this list instead of carrying their own
hardcoded column names.

Column-name grammar::

    <family>_<side>_B<min>-<max>_<name>
       K         L       0-10       Vtot

``side`` and ``band`` are omitted when an attribute has none, e.g.
``E_L_Flomvern`` or ``K_S_bekk_bredde``.

Run ``kantsoner spec list`` to print this for diffing against the Excel.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import constants as C

Band = tuple[int, int]


def band_token(band: Band) -> str:
    """`(0, 3)` -> `B0-3`.

    Matches the workbook. An earlier version zero-padded single digits to
    `B0-03` because that is what the then-current data carried; the
    specification is the contract, so the data was renamed instead.
    """
    return f"B{band[0]}-{band[1]}"


@dataclass(frozen=True)
class Attribute:
    """One declared attribute, expanded over its sides and bands."""

    name: str
    family: str  # "K" measured | "E" evaluated
    layer: str  # "primary" | "derived" | "network"
    producer: str  # module under that layer which computes it
    system: str  # "A" | "B" | "AB"
    sides: tuple[str, ...] = ("L", "R")  # () when side-less
    bands: tuple[Band, ...] = ()  # () when band-less
    dtype: str = "float"
    aggregation: str = "mean"  # see delivery/aggregate.py
    deliver: bool = False  # rolled up to 100 m / 500 m / ... in the delivery
    note: str = ""

    def columns(self) -> list[str]:
        sides = self.sides or ("",)
        bands = self.bands or (None,)
        out = []
        for band in bands:
            for side in sides:
                parts = [self.family]
                if side:
                    parts.append(side)
                if band is not None:
                    parts.append(band_token(band))
                parts.append(self.name)
                out.append("_".join(parts))
        return out


DELIVERED_SCORES = (
    "Filter_score_final",
    "Struktur_score",
    "Kontinuitet",
    "Score_total",
)

P10 = (C.PRIMARY_BAND,)
P20 = (C.AREAL_BAND,)
ALL_BANDS = C.BANDS

REGISTRY: list[Attribute] = [
    # ── L2 primary: vegetation (nDSM) ──────────────────────────────────────
    *[
        Attribute(
            f"V{hc}",
            "K",
            "primary",
            "vegetation",
            "B",
            bands=P10,
            note="% cover of the zone by this height class",
        )
        for hc in C.HEIGHT_CLASSES
    ],
    Attribute(
        "Vtot",
        "K",
        "primary",
        "vegetation",
        "AB",
        bands=ALL_BANDS,
        note="% cover by trees >= 1 m",
    ),
    Attribute("Vtot_m2", "K", "primary", "vegetation", "AB", bands=P10),
    Attribute(
        "Atot_m2",
        "K",
        "primary",
        "vegetation",
        "AB",
        bands=ALL_BANDS,
        note="zone area; denominator for every percentage",
    ),
    Attribute(
        "Vp95",
        "K",
        "primary",
        "canopy",
        "B",
        bands=P10,
        note="95th pct of canopy height > 1 m",
    ),
    Attribute(
        "Vhmed",
        "K",
        "primary",
        "canopy",
        "B",
        bands=P10,
        note="median canopy height > 1 m; drives rekk",
    ),
    Attribute("Vstdev", "K", "primary", "vegetation", "B", bands=P10),
    Attribute("Bygg", "K", "primary", "vegetation", "B", bands=P10),
    Attribute("Samferdsel", "K", "primary", "vegetation", "B", bands=P10),
    # ── L2 primary: terrain (DTM) ──────────────────────────────────────────
    Attribute(
        "s",
        "K",
        "primary",
        "terrain",
        "B",
        bands=P10,
        note="mean slope, percent",
    ),
    # ── L2 primary: land cover (grunnkart) ─────────────────────────────────
    Attribute(
        "SkogJordbruk",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone with vegetation potential (Pot_veg)",
    ),
    Attribute("SkogJordbruk_m2", "K", "primary", "landcover", "A", bands=P10),
    # Why vegetation is absent, or could grow — the arsak candidate
    # coverages (Bygg and Samferdsel come from vegetation) and the Pot_veg
    # breakdown. Ported from the minor_fixes branch.
    Attribute(
        "Konstruert",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone that is Snaumark_konstruert",
    ),
    Attribute("Konstruert_m2", "K", "primary", "landcover", "A", bands=P10),
    Attribute(
        "Jordbruk",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone under agriculture",
    ),
    Attribute("Jordbruk_m2", "K", "primary", "landcover", "A", bands=P10),
    Attribute(
        "Berg",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone on bedrock / blocky / thin soil",
    ),
    Attribute("Berg_m2", "K", "primary", "landcover", "A", bands=P10),
    Attribute(
        "Myr",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone that is mire",
    ),
    Attribute("Myr_m2", "K", "primary", "landcover", "A", bands=P10),
    Attribute(
        "Snoisbre",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone under permanent snow/ice",
    ),
    Attribute("Snoisbre_m2", "K", "primary", "landcover", "A", bands=P10),
    Attribute(
        "Vann",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone that is water",
    ),
    Attribute("Vann_m2", "K", "primary", "landcover", "A", bands=P10),
    Attribute(
        "Skog",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone that is forest",
    ),
    Attribute("Skog_m2", "K", "primary", "landcover", "A", bands=P10),
    Attribute(
        "SnaumarkJorddekt",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone soil-covered treeless ground",
    ),
    Attribute(
        "SnaumarkJorddekt_m2", "K", "primary", "landcover", "A", bands=P10
    ),
    Attribute(
        "PotVeg",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone with vegetation potential (Skog+Jordbruk+SnaumarkJorddekt)",
    ),
    Attribute("PotVeg_m2", "K", "primary", "landcover", "A", bands=P10),
    # arsak-only candidates, not their own workbook term: Skogbruk aliases
    # Skog above; SnaumarkFastfjell widens Berg to every natural Snaumark
    # reading; Vatmark is Myr | Vann. See primary/landcover.py:_cause_masks.
    Attribute(
        "Skogbruk",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone that is forest (alias of Skog, for arsak)",
    ),
    Attribute("Skogbruk_m2", "K", "primary", "landcover", "A", bands=P10),
    Attribute(
        "SnaumarkFastfjell",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone that is bedrock/blocky/thin-soil or other natural Snaumark",
    ),
    Attribute(
        "SnaumarkFastfjell_m2", "K", "primary", "landcover", "A", bands=P10
    ),
    Attribute(
        "Vatmark",
        "K",
        "primary",
        "landcover",
        "A",
        bands=P10,
        note="% of zone that is mire or water (Myr | Vann)",
    ),
    Attribute("Vatmark_m2", "K", "primary", "landcover", "A", bands=P10),
    Attribute(
        "Arealbruk",
        "K",
        "primary",
        "landcover",
        "B",
        bands=P10,
        dtype="str",
        aggregation="mode",
        note="dominant okosystemtypeNiva1",
    ),
    Attribute(
        "Treslag",
        "K",
        "primary",
        "landcover",
        "B",
        bands=P10,
        dtype="str",
        aggregation="mode",
    ),
    *[
        Attribute(f"Areal_{k}", "K", "primary", "landcover", "B", bands=P20)
        for k in ("skog", "myr", "Aapen")
    ],
    Attribute(
        "Areal",
        "K",
        "primary",
        "landcover",
        "B",
        bands=P20,
        note="total zone area at 0-20 m",
    ),
    # ── L2 primary: channel geometry ───────────────────────────────────────
    Attribute("bekk_bredde", "K", "primary", "channel", "B", sides=("S",)),
    Attribute("bekk_azimut", "K", "primary", "channel", "B", sides=("S",)),
    Attribute(
        "bekk_bredde_stretch",
        "K",
        "primary",
        "channel",
        "B",
        sides=("S",),
        note="median bekk_bredde over the stretch; decides beregningsbredde",
    ),
    Attribute(
        "beregningsbredde",
        "K",
        "primary",
        "channel",
        "B",
        sides=("S",),
        dtype="int",
        aggregation="mode",
        note="10 or 20 m: the band this stretch's Kantvegetasjon_A is drawn from",
    ),
    Attribute(
        "avstand_beregning",
        "K",
        "primary",
        "channel",
        "B",
        sides=("S",),
        note="distance from the segment's own line to the claim its bands were built from",
    ),
    Attribute(
        "Usikkerhet",
        "K",
        "primary",
        "channel",
        "B",
        sides=("S",),
        dtype="str",
        aggregation="share_ja",
        note="Ja when avstand_beregning exceeds max(50, 2x bekk_bredde) -- likely the wrong water body",
    ),
    Attribute(
        "Rorlagt_mistanke",
        "K",
        "primary",
        "rorlagt",
        "B",
        sides=("S",),
        dtype="str",
        aggregation="share_ja",
        note="Ja: a run of segments with no mapped water nearby, over jordbruk-dominated land -- suspected culvert, not a detection",
    ),
    # ── L1 geometry: the national border (geometry/landegrense.py) ─────────
    Attribute(
        "Utenfor_Norge",
        "K",
        "primary",
        "landegrense",
        "AB",
        dtype="str",
        note="Ja: this bank's buffer is foreign ground, clipped away -- it "
        "carries no measurement and the both-bank formulas skip it",
    ),
    Attribute(
        "Grenseelv",
        "K",
        "primary",
        "landegrense",
        "AB",
        sides=("S",),
        dtype="str",
        aggregation="share_ja",
        note="Ja where either bank is abroad: the segment is measured from "
        "the Norwegian side only",
    ),
    # ── L2 primary: laser survey metadata (NHM_ProsjektDekning) ────────────
    Attribute(
        "Laser_aar",
        "K",
        "primary",
        "laser_metadata",
        "B",
        sides=("S",),
        dtype="int",
        note="AARSTALL of the NHM_ProsjektDekning polygon covering the "
        "greatest length of the segment",
    ),
    # ── L2 primary: infrastructure (NVE / SVV) ─────────────────────────────
    Attribute(
        "Flomvern",
        "E",
        "primary",
        "infrastructure",
        "A",
        dtype="str",
        aggregation="share_ja",
    ),
    Attribute(
        "Forbygning",
        "E",
        "primary",
        "infrastructure",
        "A",
        dtype="str",
        aggregation="share_ja",
    ),
    Attribute(
        "Flomvern",
        "E",
        "primary",
        "infrastructure",
        "A",
        sides=("S",),
        dtype="float",
        aggregation="mean",
        note="0-1, each bank worth half, as in the veileder",
        deliver=True,
    ),
    Attribute(
        "Forbygning",
        "E",
        "primary",
        "infrastructure",
        "A",
        sides=("S",),
        dtype="float",
        aggregation="mean",
        note="0-1, each bank worth half, as in the veileder",
        deliver=True,
    ),
    Attribute(
        "Barrier",
        "E",
        "primary",
        "infrastructure",
        "A",
        sides=("S",),
        dtype="str",
        aggregation="count_ja",
        deliver=True,
    ),
    # ── L2 primary: grøfter ────────────────────────────────────────────────
    Attribute(
        "groft_azimut",
        "GK",
        "primary",
        "grofter",
        "B",
        sides=(),
        bands=P10,
        note="mean bearing of ditches meeting the zone",
    ),
    Attribute(
        "orient_groft",
        "GK",
        "primary",
        "grofter",
        "B",
        sides=(),
        bands=P10,
    ),
    Attribute(
        "Grofte_tetthet",
        "K",
        "primary",
        "grofter",
        "B",
        bands=P10,
        aggregation="sum",
        note="count of ditches crossing the bank, not running along it",
    ),
    # ── L3 derived: system A ───────────────────────────────────────────────
    # Computed at B0-10 and B0-20 (P10/P20, raw/QA, not delivered on their
    # own), then selected per segment into a band-free pair by
    # K_S_beregningsbredde -- see derived/system_a.py. The band-free entries
    # (no `bands=`) are the ones delivery/aggregate.py actually reads.
    *[
        Attribute(
            "Potensiell_Kantvegetasjon",
            "E",
            "derived",
            "system_a",
            "A",
            bands=bands,
            dtype="str",
            aggregation="share_ja",
        )
        for bands in (P10, P20, ())
    ],
    *[
        Attribute(
            "Kantvegetasjon_A",
            "E",
            "derived",
            "system_a",
            "A",
            bands=bands,
            dtype="str",
            aggregation="kantveg_ratio",
        )
        for bands in (P10, P20, ())
    ],
    *[
        Attribute(
            "Kantvegetasjon_A_arsak",
            "E",
            "derived",
            "system_a",
            "A",
            bands=bands,
            sides=("L", "R", "S"),
            dtype="str",
            aggregation="block_cause",
            note="dominant cause of a negative status; blocks reclassify from summed areas",
        )
        for bands in (P10, P20, ())
    ],
    *[
        Attribute(
            "Menneskelig_inngrep",
            "E",
            "derived",
            "system_a",
            "A",
            bands=bands,
            sides=("L", "R", "S"),
            dtype="str",
            aggregation="block_cause",
            note="is the missing vegetation human-caused (Ja/Nei, N/A when OK)",
        )
        for bands in (P10, P20, ())
    ],
    *[
        Attribute(
            "Kantvegetasjon_A",
            "E",
            "derived",
            "system_a",
            "A",
            sides=("S",),
            bands=bands,
            dtype="str",
            aggregation="kantveg_ratio",
        )
        for bands in (P10, P20, ())
    ],
    # ── L3 derived: width ladder ───────────────────────────────────────────
    Attribute(
        "Bredde", "K", "derived", "width", "B", dtype="str", aggregation="mode"
    ),
    Attribute("w_h", "K", "derived", "width", "B"),
    Attribute("w_eff", "K", "derived", "width", "B", bands=P10),
    Attribute("orient_fakt", "K", "derived", "width", "B", sides=("S",)),
    # ── L3 derived: filtrering ─────────────────────────────────────────────
    *[
        Attribute(n, "K", "derived", "filtrering", "B", bands=P10)
        for n in (
            "v_faktor_SS",
            "v_faktor_NP",
            "N_sub",
            "H_SS",
            "H_PN",
            "Y_SS",
            "Y_X",
            "Y_N",
            "Y_P",
            "Y_N_total",
            "Y_SS_final",
            "Y_P_final",
            "Y_N_final",
            "Y_N_total_final",
        )
    ],
    Attribute("Grofte_bypass", "K", "derived", "filtrering", "B", bands=P10),
    # ── L3 derived: usikkerhet ─────────────────────────────────────────────
    Attribute(
        "usikkerhet_flagg",
        "E",
        "derived",
        "usikkerhet",
        "B",
        sides=("S",),
        dtype="int",
        note="bitmask of the defaults used (report Tabell C.2): 2 arealbruk, 4 løv",
    ),
    # ── L3 derived: temperatur ─────────────────────────────────────────────
    Attribute("rekk", "K", "derived", "temperatur", "B", bands=P10),
    Attribute("Lauv_faktor_temp", "K", "derived", "temperatur", "B", bands=P10),
    Attribute(
        "Skygge_effekt",
        "K",
        "derived",
        "temperatur",
        "B",
        sides=("S",),
        bands=P10,
    ),
    Attribute(
        "Temp_score", "K", "derived", "temperatur", "B", sides=("S",), bands=P10
    ),
    # ── L3 derived: næring ─────────────────────────────────────────────────
    Attribute("Loev_kvalitet", "K", "derived", "naering", "B", bands=P10),
    Attribute("Naering_score", "K", "derived", "naering", "B", bands=P10),
    Attribute("Naering_effekt", "K", "derived", "naering", "B", bands=P10),
    # ── L3 derived: kontinuitet ────────────────────────────────────────────
    Attribute("Kontin_tre", "K", "derived", "kontinuitet", "B", bands=P10),
    Attribute("Kontin_busk", "K", "derived", "kontinuitet", "B", bands=P10),
    Attribute("Kontin", "K", "derived", "kontinuitet", "B", bands=P10),
    # ── L3 derived: scores ─────────────────────────────────────────────────
    *[
        Attribute(
            n,
            "E",
            "derived",
            "scores",
            "B",
            sides=("S",),
            deliver=n in DELIVERED_SCORES,
        )
        for n in (
            "Filter_score",
            "Filter_score_naering",
            "Filter_score_final",
            "Struktur_score",
            "Kontinuitet",
            "Score_total",
        )
    ],
    # ── L4 network ─────────────────────────────────────────────────────────
    Attribute(
        "Oppstuving",
        "E",
        "network",
        "barriers",
        "A",
        sides=("S",),
        dtype="str",
        aggregation="share_ja",
        deliver=True,
    ),
    Attribute(
        "Fragmentering",
        "E",
        "network",
        "barriers",
        "A",
        sides=("S",),
        dtype="str",
        # Worst of the block, not the most common. Fragmentation is a
        # connectivity property: a barrier close upstream fragments the
        # whole block, not the share of it nearest the barrier, so a block
        # that is mostly "God" with a stretch of "Svært dårlig" is not God.
        # Same logic as Barrier, where any barrier counts.
        aggregation="worst",
        deliver=True,
    ),
    *[
        Attribute(n, "K", "network", "connectivity", "B", sides=("S",))
        for n in ("K_lateral", "K_oppstroems", "K_nettverk")
    ],
    Attribute("K_score", "E", "network", "connectivity", "B", sides=("S",)),
    *[
        Attribute(n, "E", "network", "priority", "B", sides=("S",))
        for n in ("f_konnektivitet", "f_arealbruk", "f_helning")
    ],
    *[
        Attribute(n, "K", "network", "priority", "B", sides=("S",))
        for n in (
            "fremmedart_risiko",
            "w_basis",
            "w_opt",
            "gap_rel",
            "prioritet",
        )
    ],
    Attribute(
        "d",
        "K",
        "network",
        "priority",
        "B",
        sides=("S",),
        note="m along the network to the nearest upstream invasive occurrence; 0 in-segment",
    ),
    Attribute(
        "fremmedart_art",
        "K",
        "network",
        "priority",
        "B",
        sides=("S",),
        dtype="str",
        aggregation="mode",
        note="invasive species recorded in the segment's kantsone (Artsdatabanken)",
    ),
]

# Identity columns carried alongside every attribute. The segments are cut
# from Vann-Nett itself, so the ids ARE the official ones.
IDENTITY_COLUMNS = (
    "segment_id",
    "stretch_id",
    "elvID",
)


def select(
    *,
    layer: str | None = None,
    system: str | None = None,
    producer: str | None = None,
) -> list[Attribute]:
    """Registry entries matching every supplied filter."""
    out = REGISTRY
    if layer is not None:
        out = [a for a in out if a.layer == layer]
    if producer is not None:
        out = [a for a in out if a.producer == producer]
    if system is not None:
        out = [a for a in out if system in a.system]
    return list(out)


def columns(**kwargs) -> list[str]:
    """Flat, de-duplicated column names for the selected attributes."""
    seen: dict[str, None] = {}
    for attr in select(**kwargs):
        for col in attr.columns():
            seen.setdefault(col, None)
    return list(seen)


def delivered(rule: str | None = None) -> list[str]:
    """Columns rolled up in the delivery, optionally filtered by rule.

    This is what drives delivery/aggregate.py, so adding an attribute to the
    delivery means setting `deliver=True` here and nothing else.
    """
    return [
        col
        for attr in REGISTRY
        if attr.deliver and (rule is None or attr.aggregation == rule)
        for col in attr.columns()
    ]
