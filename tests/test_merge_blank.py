"""merge_delivery.blank_unclassified: the score and priority of rows that read
Uklassifisert are emptied, nothing else is touched. Synthetic table, no data."""

import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import merge_delivery as M  # noqa: E402


def _table():
    con = sqlite3.connect(":memory:")
    con.execute(
        'CREATE TABLE "official_centerline_10m" (fid INTEGER PRIMARY KEY, segment_id TEXT, '
        '"E_S_Kantvegetasjon_A" TEXT, "E_S_Score_total" REAL, "K_S_prioritet" REAL, '
        '"E_S_Struktur_score" REAL)'
    )
    rows = [
        (1, "a", "Uklassifisert", 0.0, 1.4, 0.1),  # both banks without a zone
        (2, "b", "Uklassifisert", 0.32, 0.9, 0.2),  # one bank without a zone
        (3, "c", "Naturlig kantskog (vegetasjon)", 0.5, 0.2, 0.3),
        (4, "d", "Ikke vurdert (bekkelukking)", 0.4, 0.3, 0.4),
        (5, "e", None, 0.6, 0.1, 0.5),
    ]
    con.executemany('INSERT INTO "official_centerline_10m" VALUES (?,?,?,?,?,?)', rows)
    return con


def test_blanks_only_uklassifisert_rows_and_only_the_two_columns():
    con = _table()
    n = M.blank_unclassified(con, "official_centerline_10m")
    assert n == 2  # without the aggregation's list only Uklassifisert goes
    got = {
        r[0]: r[1:]
        for r in con.execute(
            'SELECT fid, "E_S_Score_total", "K_S_prioritet", "E_S_Struktur_score" '
            'FROM "official_centerline_10m"'
        )
    }
    assert got[1] == (None, None, 0.1)
    assert got[2] == (None, None, 0.2)  # the other score columns stay
    assert got[3] == (0.5, 0.2, 0.3)
    assert got[4] == (0.4, 0.3, 0.4)  # a culvert row is not Uklassifisert
    assert got[5] == (0.6, 0.1, 0.5)


def test_a_culvert_the_run_did_not_compute_is_blanked_too():
    con = _table()
    n = M.blank_unclassified(con, "official_centerline_10m", not_computed=["d"])
    assert n == 3
    got = {
        r[0]: r[1:]
        for r in con.execute(
            'SELECT fid, "E_S_Score_total", "K_S_prioritet" FROM "official_centerline_10m"'
        )
    }
    assert got[4] == (None, None)  # Ikke vurdert (bekkelukking), Beregnet = no
    assert got[3] == (0.5, 0.2)  # a computed row is not touched


def test_not_computed_reads_beregnet_no_from_the_aggregation(tmp_path):
    db = tmp_path / "agg.gpkg"
    con = sqlite3.connect(db)
    con.execute('CREATE TABLE "official_centerline_10m" (fid INTEGER PRIMARY KEY, segment_id TEXT, Beregnet TEXT)')
    con.executemany(
        'INSERT INTO "official_centerline_10m" (segment_id, Beregnet) VALUES (?,?)',
        [("a", "yes"), ("b", "no"), ("c", "yes"), ("d", "no")],
    )
    con.commit()
    con.close()
    assert sorted(M.not_computed(db, "official_centerline_10m")) == ["b", "d"]


def test_a_table_without_the_columns_raises():
    con = sqlite3.connect(":memory:")
    con.execute('CREATE TABLE "t" (fid INTEGER PRIMARY KEY, "E_S_Kantvegetasjon_A" TEXT)')
    with pytest.raises(sqlite3.OperationalError):
        M.blank_unclassified(con, "t")


def test_the_constant_is_the_runs():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from kantsoner.spec import constants as C

    assert M.UNCLASSIFIED == C.UKLASSIFISERT_KANTVEG
