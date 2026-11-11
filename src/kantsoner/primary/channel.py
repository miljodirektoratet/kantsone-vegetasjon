"""Channel geometry: how wide the river is here, and which way it runs.

`bekk_bredde` sets the baseline recommended buffer width
(`w_basis = clip(10 + 1.5 * bredde, 10, 50)`) and, per stretch, which
analysis band (0-10 or 0-20 m) is the "official" one for that stretch's
segments -- see `beregningsbredde` below. `bekk_azimut` sets the orientation
factor that scales shading — an east-west reach is shaded very differently
from a north-south one at this latitude.
"""

import numpy as np
import pandas as pd
import shapely

from ..io import as_gdf
from ..spec import constants as C

# A stretch whose median channel width reaches the wide band's own width is
# treated as needing that wider ring for its "official" measurement -- a
# narrower ring would clip a real channel this size into land on both sides.
WIDE_STRETCH_THRESHOLD_M = C.AREAL_BAND[1]
# Usikkerhet ("Ja"): the claim/channel used for this segment's calculation
# sits further from the segment's own ELVIS line than plausible drift for its
# width would explain -- most likely the true water is missing from
# grunnkart/FKB and routing locked onto a nearby, unrelated water feature
# instead. A flat floor alone over-flags wide rivers, where a route
# legitimately wanders far from the drawn line while still sitting on the
# right channel; scaling by the channel's own width keeps that from being a
# false positive.
AVSTAND_FLOOR_M = 50.0
AVSTAND_WIDTH_FACTOR = 2.0


def _azimuth(geom):
    """Compass bearing of a segment, degrees clockwise from north.

    `arctan2(d_east, d_north)` per the workbook. Folded onto 0-180: a reach
    running north-east and one running south-west have the same orientation
    to the sun, and only the orientation matters downstream.
    """
    coords = list(geom.coords) if geom.geom_type == "LineString" else None
    if not coords or len(coords) < 2:
        return np.nan
    (x0, y0), (x1, y1) = coords[0], coords[-1]
    deg = np.degrees(np.arctan2(x1 - x0, y1 - y0))
    return float(deg % 180.0)


def extract(segments, claims, water):
    """bekk_bredde, bekk_azimut and the buffer-band/uncertainty columns,
    keyed on segment_id.

    The azimuth comes off the ELVIS segment itself. The width cannot: the
    ELVIS line is coarse and its midpoint is routinely on land. The claim —
    the segment's own piece of the routed channel line — lies inside the
    water by construction, so its midpoint's distance to the nearest bank is
    half the channel. Segments whose stretch never routed have no claim and
    report NaN.
    """
    seg = as_gdf(segments)
    azimut = pd.Series(
        [_azimuth(g) for g in seg.geometry], index=seg.index, dtype="float64"
    )

    bredde = pd.Series(np.nan, index=seg.index, dtype="float64")
    avstand = pd.Series(np.nan, index=seg.index, dtype="float64")
    claims_gdf = as_gdf(claims)
    water_gdf = as_gdf(water)
    if len(claims_gdf) and len(water_gdf):
        banks = shapely.boundary(water_gdf.geometry.values)
        tree = shapely.STRtree(banks)
        mids = shapely.line_interpolate_point(
            claims_gdf.geometry.values, 0.5, normalized=True
        )
        nearest = tree.query_nearest(mids, all_matches=False)
        dist = shapely.distance(mids[nearest[0]], banks[nearest[1]])
        # Twice the distance to the nearest bank is a channel width only
        # when the claim is *in* the channel. Where the route crossed dry
        # ground -- a piped reach, a gap grunnkart never drew -- the
        # nearest bank is some unrelated water far away, and doubling that
        # produced the impossible widths: 41,284 of Hurdalsvassdraget's
        # 131,250 claims have their midpoint on land, median 127 m against
        # 1.9 m for the rest, and they are 61% of everything over 100 m.
        # There is no channel to measure, so the answer is no answer.
        # K_S_Usikkerhet already flags these, and w_basis falls back to its
        # own floor rather than to a number invented from open water.
        wet = shapely.union_all(shapely.make_valid(water_gdf.geometry.values))
        in_channel = shapely.contains(wet, mids[nearest[0]])
        width = np.where(in_channel, dist * 2.0, np.nan)
        by_id = pd.Series(
            width, index=claims_gdf["segment_id"].values[nearest[0]]
        )
        bredde = (
            seg["segment_id"].map(by_id).astype("float64").set_axis(seg.index)
        )

        # How far the claim actually used for this segment's bands sits from
        # the segment's own ELVIS line -- the matching-confidence signal
        # Usikkerhet is built from below.
        claim_by_id = dict(
            zip(
                claims_gdf["segment_id"].values,
                claims_gdf.geometry.values,
                strict=True,
            )
        )
        claim_geom = seg["segment_id"].map(claim_by_id)
        has_claim = claim_geom.notna().to_numpy()
        if has_claim.any():
            seg_mid = shapely.line_interpolate_point(
                seg.geometry.values[has_claim], 0.5, normalized=True
            )
            claim_arr = np.asarray(
                claim_geom.to_numpy()[has_claim], dtype=object
            )
            avstand.iloc[has_claim.nonzero()[0]] = shapely.distance(
                seg_mid, claim_arr
            )

    if "stretch_id" in seg.columns:
        stretch_width = bredde.groupby(seg["stretch_id"]).transform("median")
    else:
        stretch_width = pd.Series(np.nan, index=seg.index, dtype="float64")
    beregningsbredde = np.where(
        (stretch_width >= WIDE_STRETCH_THRESHOLD_M)
        & (not C.ALWAYS_PRIMARY_BAND),
        C.AREAL_BAND[1],
        C.PRIMARY_BAND[1],
    )

    avstand_threshold = np.maximum(
        AVSTAND_FLOOR_M, AVSTAND_WIDTH_FACTOR * bredde.fillna(0.0)
    )
    # Three ways to be uncertain, and the distance test only catches one of
    # them. A fallback segment was measured off an assumed channel rather
    # than mapped water; a segment with no claim at all had nothing measured
    # from. Without these two, "Nei" would read as "we are confident" on
    # exactly the segments there is least reason to be confident about.
    uncertain = avstand > avstand_threshold
    if "fallback" in seg.columns:
        uncertain = uncertain | seg["fallback"].fillna(False).astype(bool)
    uncertain = uncertain | avstand.isna()
    usikkerhet = pd.Series(np.where(uncertain, "Ja", "Nei"), index=seg.index)

    return {
        "K_S_bekk_bredde": bredde.round(2),
        "K_S_bekk_azimut": azimut.round(2),
        "K_S_bekk_bredde_stretch": stretch_width.round(2),
        "K_S_beregningsbredde": pd.Series(
            beregningsbredde, index=seg.index, dtype="int64"
        ),
        "K_S_avstand_beregning": avstand.round(2),
        "K_S_Usikkerhet": usikkerhet,
    }
