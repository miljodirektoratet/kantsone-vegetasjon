"""Method constants — the specification, not configuration.

Every number here comes from one of the two product-specification workbooks in
``data/``. They are code rather than YAML on purpose: a typo in a coefficient
would silently corrupt every score, so these belong somewhere importable,
greppable, diffable against the Excel and covered by tests.

Sources:
  A = data/Kantsoner_Morfologi_produktspesifikasjon.xlsx
  B = data/Kantsoner_Morfologi_produktspesifikasjon_LB_v03.xlsx (Attributes tab)
"""

from __future__ import annotations

import math

# ── Project-wide ───────────────────────────────────────────────────────────
CRS = "EPSG:25833"  # ETRS89 / UTM 33N
SEGMENT_LENGTH_M = 10.0

# Buffer bands actually required by the two specs. 0-20 gives the area-
# equivalent width w_h (the 0-3 band was the first rung of the earlier ladder
# and is still measured, but nothing reads it any more); everything else is
# 0-10.
BANDS: tuple[tuple[int, int], ...] = ((0, 3), (0, 10), (0, 20))
PRIMARY_BAND = (0, 10)
AREAL_BAND = (0, 20)
WIDEST_BAND_M = max(b[1] for b in BANDS)

# Tree detection only needs to cover the riparian zones the attributes are
# measured on, plus enough margin that a crown centred just outside the widest
# zone still contributes its overhang. The old pipeline used a 150 m envelope
# around the centerline — that was the arealregnskap clip buffer, reused as
# the detection mask — which meant detecting trees over ~19x the pixels whose
# canopies are ever used.
CROWN_MARGIN_M = 15.0  # generous large-tree crown radius
DETECTION_MARGIN_M = WIDEST_BAND_M + CROWN_MARGIN_M  # from the river edge
MIN_TREE_HEIGHT_M = 0.5
# Watershed segmentation tuning. The library's original smoothing (sigma 3.0,
# median 5) merged neighbouring trees into 287 m2 basins — a ~9.5 m crown
# radius, three times a real spruce. Settled visually on Hurumlandet stretch
# 010-17-1_010000145: these give a 49 m2 median crown (~4 m radius) with one
# detected top per polygon.
TREE_GAUSSIAN_SIGMA = 1.0
TREE_MEDIAN_SIZE = 3
TREE_COMPACTNESS = 0.1
# Vp95 / Vhmed count only pixels above this, per the workbook's
# "pikselbasert (høydeverdier > 1 meter)". Distinct from MIN_TREE_HEIGHT_M,
# which is the crown-detection floor and sits lower.
CANOPY_STAT_MIN_HEIGHT_M = 1.0
# The delivered single_trees layer keeps crowns at least this tall.
SINGLE_TREE_MIN_HEIGHT_M = 1.0
SINGLE_TREE_FIELDS = (
    "tree_id",
    "tree_height",
    "height_class",
    "top_x",
    "top_y",
    "top_z",
    "bottom_z",
    "stretch_id",
    "dist_2D",
    "dist_3D",
    "geometry",
)
# A ditch this far off the river's own bearing crosses the bank rather than
# running along it, so its runoff bypasses the buffer.
GROFT_BYPASS_ANGLE_DEG = 45.0

# Canopy height buckets. Vtot sums only trees >= 1 m (V0p5-1 excluded).
HEIGHT_CLASSES = ("0p5-1", "1-5", "5-10", "10pluss")
# Vtot counts trees >= 1 m; the 0.5-1 m shrub class is excluded.
VTOT_CLASSES = ("1-5", "5-10", "10pluss")

# ── elvenett objektType ────────────────────────────────────────────────────
# The source spells these inconsistently, and the misspellings are not rare:
# InnsjoMidtlinje without the o-slash is 96,759 segments (968 km) in an
# earlier national run, and FiktivElv/FiltivElv another 226. Matching on the correct
# spelling alone silently drops them, so every variant seen in the data is
# listed and the sets are matched against, never a literal.
#
# LAKE_OBJEKTTYPER: elvenett's line down the middle of a lake. It has no
# riparian zone of its own -- the kantsone belongs to the shoreline, which
# the lake track measures from the water polygon instead. ElvBekkRegulert
# is here rather than with the rivers because the data says so: 5,950 of
# its 5,959 segments (99.8%) fall inside a Vann-Nett lake polygon, so in
# practice it marks regulated lakes and reservoirs, not regulated rivers.
# ── infrastructure as a cause ──────────────────────────────────────────────
# A bank paved over does not stop being capable of vegetation; it has been
# built on. But PotVeg counts forest, agriculture and soil-covered
# Snaumark, and asphalt is none of those, so a bank built over badly enough
# dropped under the 50% potential test and fell through to "Naturlig
# fravær av kantskog" -- the natural-absence bucket, which is excluded from
# the Kantvegetasjon share entirely. 145,166 banks across the 30 regions,
# 1.17% of all natural-absence classifications, and 9,428 of them state
# Menneskelig_inngrep = "Nei" for ground that is road or building.
#
# The failure was inverted: the more completely a road covered the bank,
# the more certain the code was that the absence was natural. So built
# ground now counts as potential in its own right. The cause reported is
# the existing Bygg / Samferdsel / Konstruert -- the class changes from
# "Naturlig fravær av kantskog" to "Nedbygd kantsone", which is the point;
# the vocabulary did not need a new word.
INFRA_POTENSIELL_MIN_PCT = 50.0

# ── lakes ──────────────────────────────────────────────────────────────────
# A lake has no channel width, and four things key on one: Skygge_effekt
# divides by it, Naering_effekt scales litter per unit length by it, and
# w_basis is derived from it. Rather than invent a width, the width-derived
# terms are dropped for lakes and the ones that remain are restated here.
#
# The recommended kantsone width. The river formula is
# clip(10 + 1.5 x bekk_bredde, 10, 50); a 2 m stream gets 13 m. Pinning
# lakes at the same 10 m floor keeps them comparable with small streams in
# one priority ranking, which is what it is for. 20 m would deliberately
# rank lake shores above small streams; 10 m lets slope, land use and
# connectivity decide instead.
LAKE_W_BASIS_M = 10.0
# Score_total without the Struktur leg, which is temperature and naering
# and therefore wholly width-derived. Stated, not renormalised: 0.5/0.3/0.2
# with the middle removed would give 0.714/0.286, and a weight nobody chose
# is a weight nobody can defend.
LAKE_SCORE_WEIGHTS = {"filter": 0.7, "kontinuitet": 0.3}
# Written onto every segment so the delivery says which track measured it,
# and so the two can be told apart after the fact.
KILDE_LAKE, KILDE_RIVER = "Innsjø", "Elv"
KILDE_COL = "K_S_Kilde"

LAKE_OBJEKTTYPER = (
    "InnsjøMidtlinje",
    "InnsjoMidtlinje",
    "InnsjøMidtlinjeReg",
    "InnsjøRegulert",
    "ElvBekkRegulert",
)
# A constructed line, drawn where no water surface is mapped. The earlier
# reading -- "a connector rather than a watercourse" -- does not survive the
# data: only 0.7% lie inside a lake, 96% carry a Vann-Nett water-body id,
# they run in continuous stretches of 755 m median and up to 6.4 km, and the
# named ones are named bekker. They are small streams grunnkart never drew,
# mostly ditched or culverted through farmland -- which is where the filter
# function of a kantsone matters most. 932 km of them in an earlier national run.
FIKTIV_OBJEKTTYPER = ("ElvelinjeFiktiv", "FiktivElv", "FiltivElv")
# Real watercourses: everything that is neither a lake line nor excluded.
# ElvBekRegulert (one k, 150 segments) is here and not with the lakes --
# none of it falls inside a lake polygon, so unlike ElvBekkRegulert it does
# look like a regulated river.
RIVER_OBJEKTTYPER = (
    "ElvBekk",
    "ElvBekkMidtlinje",
    "ElvBekRegulert",
) + FIKTIV_OBJEKTTYPER

# ── System A ───────────────────────────────────────────────────────────────
# Both K_*_SkogJordbruk and K_*_Vtot are percentages, so the threshold is 50.
KANTVEG_THRESHOLD_PCT = 50.0
# has_veg's threshold at the 0-20 m band, relative to the same (possibly
# per-segment dynamic) threshold used at 0-10 m.
#
# Not a dilution correction. Each band is measured on itself -- the cover
# percentage is vegetation area in that band over the area of that same band
# (primary/vegetation.py) -- so widening the band does not by itself lower
# the number. The bar comes down for a different reason: on the wide rivers
# that get the 20 m band there is often a bare sand or gravel strip along
# the water that grunnkart does not register as such, so the real kantsone
# starts further back than the water's edge. Measured from the edge, that
# strip eats into the band and pushes cover under 50 % even where the
# vegetation behind it is intact.
#
# 0.75 is a blunt margin, wider than the effect it compensates for. In an
# earlier national run, on the stretches that qualify for the 20 m band, the 10-20 m belt is
# about 1 pp *denser* than the 0-10 m belt (41.2 % vs 40.2 %, 991k
# segments), while on narrow stretches it is 1.3 pp sparser -- the sand
# strip is visible, but it is worth ~1 pp, not 25 %. So this will also pass
# segments elsewhere that a 50 % bar would have failed. Worth revisiting
# against a per-segment strip detection rather than a flat discount.
#
# PotVeg's own flat 50 % cutoff is unaffected -- only feeds Kantvegetasjon_A,
# not Potensiell_Kantvegetasjon's PotVeg leg.
KANTVEG_WIDE_THRESHOLD_SCALE = 0.75

# Zones are clipped to Norwegian territory (geometry/landegrense.py). A bank
# left with less than this share of the area it was built with is treated as
# having no data at all, and the both-bank formulas use the other one. Set
# well below any plausible partial overlap: the case this is for measures
# 100 % foreign, and a Norwegian bank beside a generalised border outline
# loses a sliver at most, so nothing lands near the cut. A bank in between
# keeps its Norwegian part and is measured over it, which is why this only
# has to separate "gone" from "there".
GRENSE_KEEP_MIN_FRAC = 0.1

# Manglende kantvegetasjon split in two: which competing land use displaced
# it decides Kantskog_fjernet (agriculture/forestry) vs Nedbygd_kantsone
# (built/infrastructure) -- see derived/kantveg.py. Kantvegetasjon_ikke_mulig
# renamed to make clear it is the *natural* absence case, not a data gap;
# Uklassifisert is the data-gap case, applied where Matched = no.
OK_KANTVEG = "Naturlig kantskog (vegetasjon)"
KANTSKOG_FJERNET = "Kantskog fjernet"
NEDBYGD_KANTSONE = "Nedbygd kantsone"
IKKE_MULIG_KANTVEG = "Naturlig fravær av kantskog"
UKLASSIFISERT_KANTVEG = "Uklassifisert"
# A suspected piped/culverted reach (K_S_Rorlagt_mistanke = Ja). Agreed with
# the client: a bekkelukking is left out of the condition assessment, not
# reported as removed vegetation -- which is what the open-channel rule made
# of nearly all of them (Hurdalsvassdraget: 1,925 of 2,061 "Kantskog fjernet").
IKKE_VURDERT_BEKKELUKKING = "Ikke vurdert (bekkelukking)"

# ── Kantvegetasjon, graded (agreed with the client 2026-10-23) ──────────────
# Each bank of each 10 m segment is graded by how much of its *expected*
# canopy cover is there, instead of yes/no against one threshold. Expected
# cover is the area-weighted mean of these per-class values over the bank's
# 0-10 m zone; ground outside these classes (hei, myr, berg, strand, water)
# expects 0. Values from the reference analysis of kantsoner in verneområder
# and NIBIO naturskog, adjusted with the client (dyrket/bygg/grasmark).
FORVENTET_DEKNING_PCT = {
    "svhoy": 95.0,  # skog, svært høy bonitet
    "hoy": 90.0,  # skog, høy bonitet
    "middels": 85.0,  # skog, middels bonitet
    "lav": 75.0,  # skog, lav bonitet
    "imp": 45.0,  # skog, impediment (not vegetation-capable on its own)
    "dyrket": 95.0,  # dyrket mark
    "bygg": 90.0,  # bebygd, not transport
    "veg": 85.0,  # bebygd, arealbruk TransportTelek
    "gras": 90.0,  # grasmark
}
# Classes that make a bank vegetation-capable; impediment does not.
VEGDYKTIG_KLASSER = (
    "svhoy",
    "hoy",
    "middels",
    "lav",
    "dyrket",
    "bygg",
    "veg",
    "gras",
)
# A bank is assessed when at least this share of its zone is vegetation-
# capable, or when its cover reaches the expected cover anyway.
VEGDYKTIG_MIN_PCT = 50.0
# A reporting unit is classified only when at least this share of its length
# (bekkelukking excluded) was actually measured.
DATADEKNING_MIN = 0.5
# 10 m everywhere: the 0-20 m band on wide stretches is no longer used for
# the official columns (K_S_beregningsbredde is 10 throughout).
ALWAYS_PRIMARY_BAND = True
KANTVEG_GRAD_IKKE_VURDERT = "Ikke vurdert"
KANTVEG_GRAD_IKKE_BEREGNET = "Ikke beregnet"

# Aggregation levels, in roll-up order. Each level is built from the previous.
AGGREGATION_LEVELS = ("100m", "500m", "1000m", "river", "region")
ROLLUP_FACTORS = {"100m": 10, "500m": 5, "1000m": 2}

# ── System B — vegetation width ladder (Attributes rows 11-12) ─────────────
# w_h: the effective riparian width of a bank, the area-equivalent width.
# It is the share of the 0-20 m band that carries crowns, times the 20 m the
# band is wide: a bank with crowns over half of the band has the width of a
# strip of 10 m, wherever in the band the crowns stand. It runs from 0 to
# 20 m without steps (the report asks for a width "uten faste sonegrenser",
# chapter 3.1, and says a strip that starts at 3 m scores lower but not 0).
# The earlier version was a four-rung ladder (0, 3, 10 or 20 m) that gave 0
# whenever the 0-3 m band had under 50 % cover.
WIDTH_BAND_M = 20.0

# Bredde is the delivered class of w_h. The limits are the midpoints between
# the four rungs of the earlier ladder (0, 3, 10 and 20 m), so the labels
# keep their meaning: (upper limit of w_h in m, label), first match wins.
BREDDE_CLASSES = (
    (1.5, "manglende"),
    (6.5, "smale"),
    (15.0, "moderat"),
    (float("inf"), "bred"),
)

# Grunnkart okosystemtypeniva1 groups lakes and rivers separately, and the
# spec's N_SUB and F_AREALBRUK tables are keyed on the Norwegian class names
# from the older okosystemtype_1 (delivered in English then, translated
# below) schema, one level coarser in places. Without this mapping every
# lookup falls through to its default and silently zeroes the nitrogen and
# land-use terms.
# The three semi-natural groups behind Areal_skog / Areal_myr / Areal_Aapen.
# Their sum over both banks, divided by Areal, is k_lateral.
AREAL_CLASS_FROM_GRUNNKART = {
    "skog": ("skog",),
    "myr": ("vatmark",),
    "Aapen": (
        "heiBuskmark",
        "grasmark",
        "liteVegetertMark",
    ),
}

AREALBRUK_FROM_GRUNNKART = {
    "skog": "Skog",
    "bebygdOpparbeidetAreal": "Bebygd og opparbeidet areal",
    "heiBuskmark": "Hei og buskmark",
    "liteVegetertMark": "Lite vegetert mark",
    "dyrketmark": "Dyrket mark",
    "vatmark": "Vatmark",
    "kyststrenderSvabergDyner": "Kyststrender, svaberg og dyner",
    "grasmark": "Grassmark",
    "innsjoerVannmagasiner": "Ferskvann",
    "elverBekkerKanaler": "Ferskvann",
    "hav": "Hav",
}

# Species of tall canopy, from grunnkart's arealdekkeniva2: årsversjon 2025 has
# no treslag field. filtrering reads Lauv and Blandingsskog; Barskog falls to
# its tall_other branch with everything else.
TRESLAG_FROM_AREALDEKKE = {
    f"{ground}{species}": treslag
    for ground in ("skog", "myrTresatt", "jordbrukInnmarksbeiteTresatt")
    for species, treslag in (
        ("Lauv", "Lauv"),
        ("Blanding", "Blandingsskog"),
        ("Barblanding", "Barskog"),
        ("Furu", "Barskog"),
        ("Gran", "Barskog"),
    )
}


# ── System B — filtrering ──────────────────────────────────────────────────
# Vegetation-type modifiers. Keys are the (has_cover, tall_dominant, treslag)
# situations described in the v_faktor rows; see derived/filtrering.py.
V_FAKTOR_SS = {
    "tall_lauv": 1.00,
    "tall_other": 0.95,
    "low": 0.85,
    "Grassmark": 0.75,
    "Hei og buskmark": 0.60,
    "none": 0.0,
}
V_FAKTOR_NP = {
    "tall_lauv": 1.00,
    "tall_other": 0.85,
    "low": 0.70,
    "Grassmark": 0.50,
    "Hei og buskmark": 0.35,
    "none": 0.0,
}

# Nitrogen subsurface factor by arealbruk class.
N_SUB = {
    "Vatmark": 0.80,
    "Kyststrender, svaberg og dyner": 0.80,
    "Skog": 0.65,
    "Hei og buskmark": 0.55,
    "Grassmark": 0.50,
    "Dyrket mark": 0.30,
    "Bebygd og opparbeidet areal": 0.20,
    "Lite vegetert mark": 0.20,
}
N_SUB_DEFAULT = 0.0

# Removal-efficiency curves: Y = ASYMPTOTE * v_faktor * (1 - exp(-RATE * w_eff))
Y_SS_ASYMPTOTE, Y_SS_RATE = 90.9, 0.446  # sediment
Y_X_ASYMPTOTE, Y_X_RATE = 93.2, 0.215  # pesticides (no v_faktor)
Y_N_ASYMPTOTE, Y_N_RATE = 92.0, 0.160  # nitrogen
Y_P_ASYMPTOTE, Y_P_RATE = 89.5, 0.250  # phosphorus

# Y_N_total = SURFACE_WEIGHT * Y_N + SUBSURFACE_WEIGHT * N_sub * 100 (report
# A.2). Filter_score reads Y_N_total_final, the same with the slope modulator
# on the surface term: SURFACE_WEIGHT * Y_N_final + SUBSURFACE_WEIGHT * N_sub *
# 100. Chapter 3.3 puts H_PN on N as on P; the percolation term is N_sub.
Y_N_SURFACE_WEIGHT, Y_N_SUBSURFACE_WEIGHT = 0.4, 0.6

# Slope modulators. H_SS is piecewise around 10 %.
H_SS_LOW_INTERCEPT, H_SS_LOW_SLOPE = 71.02, 2.49
H_SS_HIGH_INTERCEPT, H_SS_HIGH_SLOPE = 134.09, -3.86
H_PN_BREAK_PCT, H_PN_DECAY = 10.0, 0.3

# Grøfte bypass: max(0, 1 - ln(1 + tetthet) / ln(GROFT_LN_BASE))
GROFT_LN_BASE = 11.0
# A grøft only counts when its bearing differs enough from the stream's.
GROFT_ANGLE_MIN_DEG = 45.0

# ── System B — temperatur ──────────────────────────────────────────────────
SHADE_ANGLE_DEG = 37.0
TAN_SHADE_ANGLE = math.tan(math.radians(SHADE_ANGLE_DEG))
# orient_fakt = ORIENT_BASE + ORIENT_AMPLITUDE * sin(azimuth)
ORIENT_BASE, ORIENT_AMPLITUDE = 0.7, 0.3
# Stream width is capped before being used as the shading denominator.
WIDTH_CAP_M = 20.0
# Temp_score scales canopy height against this reference height.
TEMP_REFERENCE_HEIGHT_M = 15.0
LAUV_FAKTOR_TALL, LAUV_FAKTOR_LOW, LAUV_FAKTOR_NONE = 1.0, 0.7, 0.0

# ── System B — næring ──────────────────────────────────────────────────────
LOEV_KVALITET = {
    "tall_lauv": 1.00,
    "tall_blanding": 0.65,
    "tall_other": 0.25,
    "low": 0.40,
    "none": 0.0,
}

# ── System B — kontinuitet ─────────────────────────────────────────────────
KONTIN_TRE_WEIGHT, KONTIN_BUSK_WEIGHT = 0.7, 0.3

# ── System B — scores ──────────────────────────────────────────────────────
FILTER_WEIGHTS = {"SS": 0.25, "P": 0.35, "N": 0.25, "X": 0.15}
FILTER_NAERING_BASE, FILTER_NAERING_GAIN = 0.8, 0.2
STRUKTUR_TEMP_WEIGHT, STRUKTUR_NAERING_WEIGHT = 0.5, 0.5
SCORE_TOTAL_WEIGHTS = {"filter": 0.5, "struktur": 0.3, "kontinuitet": 0.2}

# ── System B — konnektivitet / prioritet ───────────────────────────────────
CONNECTIVITY_DECAY_M = 2500.0  # exp(-d / 2500)
NETWORK_RADIUS_M = 5000.0
F_KONNEKTIVITET_BASE, F_KONNEKTIVITET_GAIN = 1.3, 0.3
# Defaults for what the data does not say, from the report's Tabell C.1: a bank
# with no land use read is åker, the conservative class (highest breddekrav,
# lowest N_sub), and a tall canopy with no species is løv, which the report
# calls representative for the Norwegian lowland. An earlier run did neither: a
# bank with no zone read N_sub 0 and f_arealbruk 1.0, a tall canopy with no
# species filtered as conifer.
AREALBRUK_UNKNOWN = "Dyrket mark"
TRESLAG_UNKNOWN = "Lauv"

# Tabell 11 of the report. Bebygd/vei 1.30 and Myr/naturlig åpne arealer 0.80
# were missing, so Bebygd, Hei og buskmark and Vatmark read the default.
# Grunnkart classes the report has no row for (Lite vegetert mark, Hav,
# Kyststrender, Ferskvann) and a bank with no zone still read the default.
F_AREALBRUK = {
    "Dyrket mark": 1.40,
    "Bebygd og opparbeidet areal": 1.30,
    "Grassmark": 1.20,
    "Skog": 0.90,
    "Hei og buskmark": 0.80,
    "Vatmark": 0.80,
}
F_AREALBRUK_DEFAULT = 1.0
F_HELNING_MAX_PCT = 30.0
W_BASIS_MIN_M, W_BASIS_MAX_M = 10.0, 50.0
W_BASIS_INTERCEPT, W_BASIS_WIDTH_GAIN = 10.0, 1.5
FREMMEDART_IN_SEGMENT = 0.8
FREMMEDART_DECAY_KM = 5.0
# An Artsdatabanken occurrence belongs to a segment when it lies within this
# of that segment's kantsone (the band is 20 m; the points carry their own
# positional uncertainty), and is dropped when its stated uncertainty could
# put it on another reach altogether.
FREMMEDART_MATCH_M = 25.0
FREMMEDART_MAX_UNCERTAINTY_M = 50.0
# Beyond this the upstream decay is under 1 % of the in-segment value.
FREMMEDART_REACH_KM = 25.0
PRIORITET_BASE, PRIORITET_GAIN = 0.2, 0.8

# ── Barriers / continuity (system A) ───────────────────────────────────────
# Oppstuving: walk upstream while width ratio vs. the reference stays high.
OPPSTUVING_WIDTH_RATIO = 1.5
OPPSTUVING_REFERENCE_M = 200.0  # averaged this far downstream of the barrier
OPPSTUVING_MAX_GAP_SEGMENTS = 2

# Fragmentering: distance upstream to the nearest barrier, in metres.
FRAGMENTERING_CLASSES = (
    # (upper bound m, label) — first match wins; None = no barrier upstream
    (2_000, "Svært dårlig"),
    (10_000, "Dårlig"),
    (50_000, "Moderat"),
    (float("inf"), "God"),
)
FRAGMENTERING_NO_BARRIER = "Svært god"
