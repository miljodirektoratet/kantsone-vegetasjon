"""System A and B formulas, on fixtures rather than GIS data.

These pin the behaviour at the boundaries the specification is explicit
about — the 50 % threshold, the width ladder rungs, the slope hinge — so a
future edit to spec/constants.py cannot quietly move a score.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from kantsoner import derived
from kantsoner.derived import (
    filtrering,
    kontinuitet,
    system_a,
    temperatur,
    width,
)
from kantsoner.spec import constants as C


def make_segments(n: int = 4, **overrides) -> pd.DataFrame:
    """A minimal frame carrying every L2 column the derived layer reads."""
    base = {
        "segment_id": range(1, n + 1),
        "segment_length": [10.0] * n,
        "K_S_bekk_bredde": [5.0] * n,
        "K_S_bekk_azimut": [90.0] * n,
    }
    for side in ("L", "R"):
        base.update(
            {
                f"K_{side}_B0-3_Vtot": [80.0, 20.0, 80.0, 0.0][:n],
                f"K_{side}_B0-10_Vtot": [80.0, 20.0, 60.0, 0.0][:n],
                f"K_{side}_B0-20_Vtot": [80.0, 20.0, 20.0, 0.0][:n],
                f"K_{side}_B0-10_V0p5-1": [5.0] * n,
                f"K_{side}_B0-10_V1-5": [10.0] * n,
                f"K_{side}_B0-10_V5-10": [30.0] * n,
                f"K_{side}_B0-10_V10pluss": [40.0] * n,
                f"K_{side}_B0-10_Vp95": [18.0] * n,
                f"K_{side}_B0-10_Vhmed": [12.0] * n,
                f"K_{side}_B0-10_s": [5.0] * n,
                f"K_{side}_B0-10_SkogJordbruk": [90.0, 90.0, 10.0, 10.0][:n],
                f"K_{side}_B0-10_PotVeg": [90.0, 90.0, 10.0, 10.0][:n],
                # system_a classifies at 0-20 m as well (selected per
                # segment by K_S_beregningsbredde), so these are required.
                f"K_{side}_B0-20_SkogJordbruk": [90.0, 90.0, 10.0, 10.0][:n],
                f"K_{side}_B0-20_PotVeg": [90.0, 90.0, 10.0, 10.0][:n],
                f"K_{side}_B0-10_Arealbruk": ["Skog"] * n,
                f"K_{side}_B0-10_Treslag": ["Lauv"] * n,
                f"K_{side}_B0-20_Areal": [200.0] * n,
                f"K_{side}_B0-20_Areal_skog": [100.0] * n,
                f"K_{side}_B0-20_Areal_myr": [20.0] * n,
                f"K_{side}_B0-20_Areal_Aapen": [10.0] * n,
            }
        )
    base.update(overrides)
    return pd.DataFrame(base)


# ── System A ───────────────────────────────────────────────────────────────
def test_kantvegetasjon_classes():
    out = system_a.compute(make_segments())
    assert list(out["E_S_B0-10_Kantvegetasjon_A"]) == [
        C.OK_KANTVEG,  # Vtot 80 >= 50
        # Vtot 20, but SkogJordbruk 90 -> possible: Manglende-equivalent.
        # None of the 8 arsak candidates are set in this fixture, so
        # agri_total == built_total == 0 and the tie defaults to agri.
        C.KANTSKOG_FJERNET,
        C.OK_KANTVEG,  # Vtot 60 >= 50
        C.IKKE_MULIG_KANTVEG,  # no cover, no potential
    ]


def test_threshold_is_inclusive_at_50():
    df = make_segments(n=1)
    for side in ("L", "R"):
        df[f"K_{side}_B0-10_Vtot"] = [C.KANTVEG_THRESHOLD_PCT]
    out = system_a.compute(df)
    assert out["E_S_B0-10_Kantvegetasjon_A"].iloc[0] == C.OK_KANTVEG


def test_missing_zone_is_not_possible_rather_than_null():
    """A side with no riparian zone must resolve to a class, not NaN."""
    df = make_segments(n=1)
    for side in ("L", "R"):
        df[f"K_{side}_B0-10_Vtot"] = [np.nan]
        df[f"K_{side}_B0-10_SkogJordbruk"] = [np.nan]
        df[f"K_{side}_B0-10_PotVeg"] = [np.nan]
    out = system_a.compute(df)
    assert out["E_S_B0-10_Kantvegetasjon_A"].iloc[0] == C.IKKE_MULIG_KANTVEG


# ── Area-equivalent width ──────────────────────────────────────────────────
def test_width_is_the_share_of_the_band_with_crowns():
    """w_h = 20 m x Vtot(0-20)/100, with no rungs and no zero rule."""
    out = width.compute(make_segments())
    assert list(out["K_L_w_h"]) == [16.0, 4.0, 4.0, 0.0]
    assert list(out["K_L_Bredde"]) == ["bred", "smale", "smale", "manglende"]


def test_a_strip_set_back_from_the_water_has_a_width():
    """Crowns only 10-20 m out used to read 0 (no cover in the first 3 m)."""
    df = make_segments(n=1)
    for side in ("L", "R"):
        df[f"K_{side}_B0-3_Vtot"] = [0.0]
        df[f"K_{side}_B0-10_Vtot"] = [0.0]
        df[f"K_{side}_B0-20_Vtot"] = [50.0]
    assert width.compute(df)["K_L_w_h"].iloc[0] == pytest.approx(10.0)


def test_width_classes_are_the_midpoints_of_the_old_rungs():
    w = pd.Series([0.0, 1.4, 1.5, 6.4, 6.5, 14.9, 15.0, 20.0])
    assert list(width._bredde(w)) == [
        "manglende", "manglende", "smale", "smale",
        "moderat", "moderat", "bred", "bred",
    ]


def test_missing_zone_gives_no_width_not_nan():
    df = make_segments(n=1)
    df["K_L_B0-20_Vtot"] = [np.nan]
    assert width.compute(df)["K_L_w_h"].iloc[0] == 0.0


def test_w_eff_exceeds_w_h_on_a_slope():
    """w_eff = w_h / cos(atan(s/100)) — the along-ground path is longer."""
    df = make_segments(n=1)
    df["K_L_B0-10_s"] = [30.0]
    out = width.compute(df)
    assert out["K_L_B0-10_w_eff"].iloc[0] > out["K_L_w_h"].iloc[0]


# ── Filtrering ─────────────────────────────────────────────────────────────
def test_h_ss_follows_the_report():
    got = filtrering._h_ss(pd.Series([-40.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0]))
    assert got.round(4).tolist() == [0.0, 0.4612, 0.7102, 0.9592, 0.5689, 0.1829, 0.0]


def test_slope_modulator_hinges_at_ten_percent():
    flat = filtrering.compute(make_segments(n=1))["K_L_B0-10_H_PN"].iloc[0]
    df = make_segments(n=1)
    df["K_L_B0-10_s"] = [40.0]
    steep = filtrering.compute(df)["K_L_B0-10_H_PN"].iloc[0]
    assert flat == pytest.approx(1.0)
    assert 0.0 <= steep < 1.0


def test_removal_efficiency_saturates_with_width():
    """Y_SS is asymptotic: doubling an already-wide strip adds little."""
    narrow, wide = make_segments(n=1), make_segments(n=1)
    narrow["K_L_B0-10_w_eff"] = [3.0]
    wide["K_L_B0-10_w_eff"] = [30.0]
    y_narrow = filtrering.compute(narrow)["K_L_B0-10_Y_SS"].iloc[0]
    y_wide = filtrering.compute(wide)["K_L_B0-10_Y_SS"].iloc[0]
    assert y_narrow < y_wide <= C.Y_SS_ASYMPTOTE


def test_no_vegetation_gives_zero_filtration():
    df = make_segments(n=1)
    for side in ("L", "R"):
        for col in ("Vtot", "V0p5-1", "V1-5", "V5-10", "V10pluss"):
            df[f"K_{side}_B0-10_{col}"] = [0.0]
        df[f"K_{side}_B0-20_Vtot"] = [0.0]
        df[f"K_{side}_B0-10_Arealbruk"] = ["Bebygd og opparbeidet areal"]
    out = filtrering.compute(width.compute(df))
    assert out["K_L_B0-10_Y_SS"].iloc[0] == pytest.approx(0.0)


def test_thin_cover_still_filters():
    """Cover under 50 % used to read 'none' (v_faktor 0) next to a width > 0."""
    df = make_segments(n=1)
    for side in ("L", "R"):
        df[f"K_{side}_B0-10_Vtot"] = [20.0]
        df[f"K_{side}_B0-20_Vtot"] = [20.0]
        df[f"K_{side}_B0-10_Treslag"] = ["Lauv"]
    out = filtrering.compute(width.compute(df))
    assert out["K_L_B0-10_v_faktor_SS"].iloc[0] == pytest.approx(1.0)
    assert out["K_L_B0-10_Y_SS"].iloc[0] > 0.0


def test_a_strip_set_back_reads_the_wide_height_classes():
    """No crowns in 0-10 m: tall or low comes from the 0-20 m classes."""
    df = make_segments(n=2)
    for side in ("L", "R"):
        df[f"K_{side}_B0-10_Vtot"] = [0.0, 0.0]
        for col in ("V0p5-1", "V1-5", "V5-10", "V10pluss"):
            df[f"K_{side}_B0-10_{col}"] = [0.0, 0.0]
        df[f"K_{side}_B0-20_Vtot"] = [40.0, 40.0]
        df[f"K_{side}_B0-20_V1-5"] = [30.0, 0.0]
        df[f"K_{side}_B0-20_V10pluss"] = [0.0, 30.0]
        df[f"K_{side}_B0-10_Treslag"] = ["Lauv", "Lauv"]
    got = filtrering._vegetation_situation(df, "L")
    assert list(got) == ["low", "tall_lauv"]


def test_nitrogen_total_carries_n_sub():
    """Filter_score reads Y_N_total_final = 0.4 Y_N_final + 0.6 N_sub 100."""
    from kantsoner.derived import scores

    df = make_segments(n=1)
    out = filtrering.compute(width.compute(df))
    b = "K_L_B0-10"
    assert out[f"{b}_N_sub"].iloc[0] == pytest.approx(0.65)  # Skog
    assert out[f"{b}_Y_N_total_final"].iloc[0] == pytest.approx(
        0.4 * out[f"{b}_Y_N_final"].iloc[0] + 0.6 * 0.65 * 100, abs=2e-3
    )
    # a bank with N_sub 0 scores lower than one with N_sub 0.65, all else equal
    other = df.copy()
    for side in ("L", "R"):
        other[f"K_{side}_B0-10_Arealbruk"] = ["Hav"]
    low = scores.compute(filtrering.compute(width.compute(other)))
    high = scores.compute(out)
    assert low["E_S_Filter_score"].iloc[0] < high["E_S_Filter_score"].iloc[0]


def test_f_arealbruk_has_the_five_rows_of_tabell_11():
    f = C.F_AREALBRUK
    assert f["Dyrket mark"] == 1.40 and f["Bebygd og opparbeidet areal"] == 1.30
    assert f["Grassmark"] == 1.20 and f["Skog"] == 0.90
    assert f["Hei og buskmark"] == 0.80 and f["Vatmark"] == 0.80


def test_unknown_species_on_a_tall_canopy_is_lauv():
    """Report C.1: løv is the default; a known conifer stays a conifer."""
    df = make_segments(n=3)
    for side in ("L", "R"):
        df[f"K_{side}_B0-10_Treslag"] = ["", "<NA>", "Barskog"]
    got = filtrering._vegetation_situation(df, "L")
    assert list(got) == ["tall_lauv", "tall_lauv", "tall_other"]


def test_a_bank_with_a_zone_but_no_land_use_is_aker_not_nothing():
    """Report C.1: åker (N_sub 0.30, f_arealbruk 1.40). A bank that does not
    exist, or has no zone (not measured), stays empty and cannot win the
    comparison between banks."""
    from kantsoner.columns import arealbruk
    from kantsoner.network import priority

    df = make_segments(n=2)
    df["K_L_B0-10_Arealbruk"] = ["", "Skog"]
    df["K_R_B0-10_Arealbruk"] = ["Skog", "Skog"]
    for side in ("L", "R"):
        df[f"K_{side}_B0-20_Atot_m2"] = [2000.0, 2000.0]
    assert list(arealbruk(df, "L")) == ["Dyrket mark", "Skog"]
    assert filtrering.compute(width.compute(df))["K_L_B0-10_N_sub"].iloc[0] == 0.30
    # the blank bank has a zone: the conservative one wins
    assert list(priority._dominant_arealbruk(df)) == ["Dyrket mark", "Skog"]
    # the same blank bank without a zone was not measured: the forest bank stands
    df["K_L_B0-20_Atot_m2"] = [0.0, 2000.0]
    assert list(arealbruk(df, "L")) == ["", "Skog"]
    assert list(priority._dominant_arealbruk(df)) == ["Skog", "Skog"]
    # no bank measured at all: the default, not an empty class
    df["K_R_B0-20_Atot_m2"] = [0.0, 2000.0]
    assert priority._dominant_arealbruk(df).iloc[0] == "Dyrket mark"
    # a bank marked as not there has no land use
    df["K_L_B0-20_Atot_m2"] = [2000.0, 2000.0]
    df["K_L_Ingen_bredd"] = ["Ja", "Nei"]
    assert list(arealbruk(df, "L")) == ["", "Skog"]


def test_nan_area_is_unknown_not_no_zone():
    from kantsoner import sides

    df = make_segments(n=3)
    df["K_L_B0-20_Atot_m2"] = [np.nan, 0.0, 5.0]
    assert list(sides.no_zone(df, "L")) == [False, True, False]


def test_a_bank_with_no_zone_is_skipped_in_the_two_bank_mean_not_counted_as_zero():
    """The bank has no zone area: not measured, so the other bank stands."""
    from kantsoner.derived import scores

    df = make_segments(n=1)
    for side in ("L", "R"):
        df[f"K_{side}_B0-20_Atot_m2"] = [2000.0]
        df[f"K_{side}_B0-10_Kontin_tre"] = [40.0]
        df[f"K_{side}_B0-10_Kontin_busk"] = [0.0]
    df["K_R_B0-20_Atot_m2"] = [0.0]
    df["K_R_B0-10_Kontin_tre"] = [0.0]
    assert scores._mean_sides(df, "Kontin_tre").iloc[0] == pytest.approx(40.0)
    # a measured bank with nothing on it is a real 0
    df["K_R_B0-20_Atot_m2"] = [2000.0]
    assert scores._mean_sides(df, "Kontin_tre").iloc[0] == pytest.approx(20.0)


def test_a_bank_with_no_zone_has_no_width():
    df = make_segments(n=1)
    for side in ("L", "R"):
        df[f"K_{side}_B0-20_Atot_m2"] = [2000.0]
    df["K_R_B0-20_Atot_m2"] = [0.0]
    out = width.compute(df)
    assert out["K_L_w_h"].iloc[0] == 16.0 and np.isnan(out["K_R_w_h"].iloc[0])
    assert out["K_R_Bredde"].iloc[0] == ""


def test_shading_and_slope_use_only_the_measured_bank():
    """One bank without a zone: the other stands for the pair (2 / n_measured)."""
    from kantsoner.derived import temperatur
    from kantsoner.network import priority

    base = make_segments(n=1)
    for side in ("L", "R"):
        base[f"K_{side}_B0-20_Atot_m2"] = [2000.0]
    base["K_R_B0-10_s"] = [0.0]
    half = base.copy()
    half["K_R_B0-20_Atot_m2"] = [0.0]
    two = temperatur.compute(width.compute(base))["K_S_B0-10_Skygge_effekt"].iloc[0]
    one = temperatur.compute(width.compute(half))["K_S_B0-10_Skygge_effekt"].iloc[0]
    assert one == pytest.approx(two)  # same bank twice = the same reach
    # the slope of the measured bank (5 %), not halved by the unmeasured one
    got = priority.compute(half)["E_S_f_helning"].iloc[0]
    assert got == pytest.approx(1.0 / np.cos(np.arctan(0.05)), abs=1e-3)


def test_usikkerhet_flagg_names_the_defaults_used():
    from kantsoner.derived import usikkerhet

    df = make_segments(n=4)
    for side in ("L", "R"):
        df[f"K_{side}_B0-10_Vtot"] = [80.0] * 4
        df[f"K_{side}_B0-10_Arealbruk"] = ["Skog", "", "Skog", ""]
        df[f"K_{side}_B0-10_Treslag"] = ["Lauv", "Lauv", "", ""]
    got = usikkerhet.compute(df)["E_S_usikkerhet_flagg"]
    assert list(got) == [0, usikkerhet.BIT_AREALBRUK, usikkerhet.BIT_LOEV, 6]
    # a bank with no zone was not measured and sets nothing
    for side in ("L", "R"):
        df[f"K_{side}_B0-20_Atot_m2"] = [0.0] * 4
    assert list(usikkerhet.compute(df)["E_S_usikkerhet_flagg"]) == [0] * 4
    # a bank that does not exist sets nothing
    df["K_L_Ingen_bredd"] = ["Ja"] * 4
    df["K_R_Ingen_bredd"] = ["Ja"] * 4
    assert list(usikkerhet.compute(df)["E_S_usikkerhet_flagg"]) == [0] * 4


# ── Temperatur ─────────────────────────────────────────────────────────────
def test_shading_is_capped_at_one():
    df = make_segments(n=1)
    for side in ("L", "R"):
        df[f"K_{side}_B0-10_Vhmed"] = [40.0]  # very tall canopy, narrow stream
    out = temperatur.compute(width.compute(df))
    assert out["K_S_B0-10_Skygge_effekt"].iloc[0] == pytest.approx(1.0)


def test_wide_channel_shades_less_than_narrow():
    narrow, wide = make_segments(n=1), make_segments(n=1)
    narrow["K_S_bekk_bredde"] = [2.0]
    wide["K_S_bekk_bredde"] = [20.0]
    a = temperatur.compute(width.compute(narrow))[
        "K_S_B0-10_Skygge_effekt"
    ].iloc[0]
    b = temperatur.compute(width.compute(wide))["K_S_B0-10_Skygge_effekt"].iloc[
        0
    ]
    assert a >= b


# ── Kontinuitet + scores ───────────────────────────────────────────────────
def test_kontinuitet_weights_trees_over_shrubs():
    out = kontinuitet.compute(make_segments(n=1))
    tre, busk = (
        out["K_L_B0-10_Kontin_tre"].iloc[0],
        out["K_L_B0-10_Kontin_busk"].iloc[0],
    )
    assert out["K_L_B0-10_Kontin"].iloc[0] == pytest.approx(
        C.KONTIN_TRE_WEIGHT * tre + C.KONTIN_BUSK_WEIGHT * busk
    )


def test_score_total_is_the_weighted_sum():
    out = derived.compute_all(make_segments(n=1))
    w = C.SCORE_TOTAL_WEIGHTS
    assert out["E_S_Score_total"].iloc[0] == pytest.approx(
        out["E_S_Filter_score_final"].iloc[0] * w["filter"]
        + out["E_S_Struktur_score"].iloc[0] * w["struktur"]
        + out["E_S_Kontinuitet"].iloc[0] * w["kontinuitet"],
        abs=1e-3,
    )


def test_filter_score_final_collapses_without_grofter():
    """Absent a grøftekart the bypass factor is 1.0, as the spec notes."""
    out = derived.compute_all(make_segments(n=1))
    assert out["E_S_Filter_score_final"].iloc[0] == pytest.approx(
        out["E_S_Filter_score_naering"].iloc[0], abs=1e-4
    )


def test_pipeline_is_order_dependent_but_complete():
    """compute_all must produce every column the registry declares for L3."""
    from kantsoner.spec import attributes as A

    out = derived.compute_all(make_segments())
    missing = [c for c in A.columns(layer="derived") if c not in out.columns]
    assert not missing, f"derived layer did not produce: {missing}"


def test_derived_does_not_mutate_its_input():
    df = make_segments()
    before = df.copy()
    derived.compute_all(df)
    pd.testing.assert_frame_equal(df, before)


# ── Aggregation ────────────────────────────────────────────────────────────
def test_barrier_aggregates_as_a_count_not_a_share():
    """The spec reports barriers as "Antall". Dividing by n_matched would give
    "barriers per segment", which is not a quantity anyone asked for."""
    import numpy as np

    from kantsoner.delivery import aggregate

    block = pd.DataFrame(
        {
            "length": [100.0],
            "length_matched": [100.0],
            "length_uklassifisert": [0.0],
            "length_ok": [60.0],
            "length_kantskog_fjernet": [10.0],
            "length_nedbygd_kantsone": [10.0],
            "length_ikke": [20.0],
            "sum_E_S_Barrier": [3.0],
            # A count is untouched by length-weighting; a share/mean column
            # is value * length, so 50.0 here is "50 %" over 100 m matched.
            "sum_E_S_Flomvern": [50.0],
        }
    )
    out = aggregate._ratios(block, "100m")
    assert out["100m_Barrier"].iloc[0] == 3.0  # count, untouched
    assert out["100m_Flomvern"].iloc[0] == 0.5  # share of matched
    assert out["100m_Barrier_final"].iloc[0] == 2  # 3 barriers -> grade 2
    assert not np.isnan(out["100m_Kantvegetasjon"].iloc[0])


def test_registry_drives_the_delivered_columns():
    from kantsoner.delivery import aggregate
    from kantsoner.spec import attributes as A

    assert set(aggregate.SHARE_COLUMNS) == set(A.delivered("share_ja"))
    assert set(aggregate.COUNT_COLUMNS) == set(A.delivered("count_ja"))
    assert "E_S_Barrier" in aggregate.COUNT_COLUMNS


def test_grade_scales_match_the_workbook():
    """Every boundary from the Aggregation tab, including the class edges —
    these are left-closed, and getting an edge wrong shifts a whole class."""
    from kantsoner.spec import scales

    # _Kantvegetasjon_final: <5 -> 5; 5-20 -> 4; 20-50 -> 3; 50-80 -> 2; 80+ -> 1
    assert list(
        scales.pct_to_grade([0, 4.9, 5, 19.9, 20, 49.9, 50, 79.9, 80, 100])
    ) == [5, 5, 4, 4, 3, 3, 2, 2, 1, 1]
    # _Barrier_final: none -> 5; 1-2 -> 3; 3 -> 2; >=4 -> 1
    assert list(scales.barrier_count_to_grade([0, 1, 2, 3, 4, 10])) == [
        5,
        3,
        3,
        2,
        1,
        1,
    ]
    # Sállir tabell 8: <20 -> 1; 20-40 -> 2; 40-60 -> 3; 60-80 -> 4; 80+ -> 5
    assert list(scales.score_b_to_grade([0, 20, 40, 60, 80, 100])) == [
        1,
        2,
        3,
        4,
        5,
        5,
    ]
    assert list(
        scales.fragmentering_to_grade(["Svært god", "Moderat", "Svært dårlig"])
    ) == [5, 3, 1]


def test_missing_required_input_raises_rather_than_zeroing():
    """`num` defaults a missing column to 0.0 — correct for optional sources,
    but a required one must not silently produce scores computed from zeros."""
    df = make_segments(n=2).drop(columns=["K_L_B0-10_Vhmed", "K_R_B0-10_Vhmed"])
    with pytest.raises(ValueError, match="Vhmed"):
        derived.compute_all(df)
    # explicit opt-out still works, for debugging
    out = derived.compute_all(df, allow_missing=True)
    assert "E_S_Score_total" in out.columns


def test_band_token_matches_the_data():
    """Band tokens follow the workbook: B0-3, not the zero-padded B0-03."""
    from kantsoner.spec import attributes as A

    assert A.band_token((0, 3)) == "B0-3"
    assert A.band_token((0, 10)) == "B0-10"
    assert A.band_token((3, 10)) == "B3-10"


def test_rollup_merges_multipart_block_geometries():
    """Each level rolls up the previous, so a 500 m block is built from 100 m
    geometries that are already MultiLineStrings where the river had a gap."""
    from shapely.geometry import LineString, MultiLineString

    from kantsoner.delivery.aggregate import _merge_lines

    contiguous = _merge_lines(
        [LineString([(0, 0), (1, 0)]), LineString([(1, 0), (2, 0)])]
    )
    assert contiguous.geom_type == "LineString"

    with_gap = _merge_lines(
        [
            MultiLineString([[(0, 0), (1, 0)], [(5, 0), (6, 0)]]),
            LineString([(1, 0), (2, 0)]),
        ]
    )
    assert with_gap.geom_type == "MultiLineString"
    assert _merge_lines([None]) is None


def test_grade_polarity_is_consistent_across_parameters():
    """The 1-5 scale grades deviation from the reference state, and every
    aggregated parameter is a pressure with 0 % undisturbed -- Kantvegetasjon
    is the share of vegetatable length *missing* vegetation. Full vegetation
    and no engineering grade 5 together; none and everywhere, 1."""
    from kantsoner.delivery import aggregate

    pristine = pd.DataFrame(
        {
            "length": [100.0],
            "length_matched": [100.0],
            "length_uklassifisert": [0.0],
            "length_ok": [100.0],
            "length_kantskog_fjernet": [0.0],
            "length_nedbygd_kantsone": [0.0],
            "length_ikke": [0.0],
            "sum_E_S_Forbygning": [0.0],
            "sum_E_S_Flomvern": [0.0],
        }
    )
    out = aggregate._ratios(pristine, "river")
    # Full vegetation, no engineering: every parameter at its best grade.
    assert out["river_Kantvegetasjon_final"].iloc[0] == 5
    assert out["river_Forbygning_final"].iloc[0] == 5
    assert out["river_Flomvern_final"].iloc[0] == 5

    degraded = pristine.assign(
        length_ok=[0.0],
        length_nedbygd_kantsone=[100.0],
        # value * length: 100 % over the full 100 m matched length.
        sum_E_S_Forbygning=[100.0],
        sum_E_S_Flomvern=[100.0],
    )
    out = aggregate._ratios(degraded, "river")
    assert out["river_Kantvegetasjon_final"].iloc[0] == 1
    assert out["river_Forbygning_final"].iloc[0] == 1


def test_vtot_sums_classes_present_on_only_some_segments():
    """Index-aligned Series addition yields NaN where a class is absent, which
    silently zeroed Vtot for every segment lacking one height class."""
    import pandas as pd

    by_class = {
        "1-5": pd.Series({1: 10.0, 2: 5.0}),
        "5-10": pd.Series({2: 7.0}),
        "10pluss": pd.Series({3: 20.0}),
    }
    parts = [
        by_class[hc] for hc in ("1-5", "5-10", "10pluss") if hc in by_class
    ]
    vtot = pd.concat(parts).groupby(level=0).sum()
    assert vtot.to_dict() == {1: 10.0, 2: 12.0, 3: 20.0}
