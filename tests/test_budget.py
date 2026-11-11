"""The stage memory budget the parallel run admits stages through."""

import threading
import time

import pytest

from kantsoner import cli


@pytest.fixture(autouse=True)
def plenty_free(monkeypatch):
    monkeypatch.setattr(cli, "_free_gb", lambda: 1000.0)


def test_stage_estimate_is_capped_at_the_budget():
    assert cli._stage_gb("attributes", 0, 100) == cli.BUDGET_BASE_GB
    assert cli._stage_gb("derived", 1_000_000, 100) == 100


def test_a_stage_waits_for_room_and_gets_it_on_release():
    b = cli._Budget(10, floor=0)
    b.acquire("a", 8)
    got = threading.Event()
    t = threading.Thread(target=lambda: (b.acquire("b", 5), got.set()))
    t.start()
    assert not got.wait(0.3)
    b.release("a")
    assert got.wait(2)
    assert b.used() == 5
    t.join()


def test_an_oversized_stage_runs_alone():
    b = cli._Budget(10, floor=0)
    b.acquire("giant", 50)
    assert b.used() == 50


def test_no_admission_below_the_free_memory_floor(monkeypatch):
    b = cli._Budget(100, floor=25)
    b.acquire("a", 1)
    monkeypatch.setattr(cli, "_free_gb", lambda: 10.0)
    assert not b._fits("b", 1)


def test_a_starved_stage_reserves_its_room(monkeypatch):
    b = cli._Budget(100, floor=0)
    b.acquire("a", 40)
    b.waiting["giant"] = (time.monotonic() - cli.STARVE_S - 1, 70)
    assert not b._fits("giant", 70)  # starved, but still has to fit
    assert not b._fits("small", 10)  # would push the giant further out
    b.release("a")
    b.acquire("b", 20)
    assert b._fits("small", 10)  # 20 + 10 + 70 leaves room for both
    assert not b._fits("medium", 11)
    assert b._fits("giant", 70)


def test_try_acquire_gives_up_and_forgets_the_request():
    b = cli._Budget(10, floor=0)
    b.acquire("a", 8)
    assert not b.try_acquire("b", 5, 0.2)
    assert "b" not in b.waiting
    assert b.used() == 8


def test_surcharge_goes_ahead_untracked_when_memory_is_free(monkeypatch):
    from kantsoner import membudget

    b = cli._Budget(10, floor=5)
    b.acquire("a", 10)
    monkeypatch.setattr(b, "try_acquire", lambda *a: False)
    membudget.bind(b, "r", 10_000, 5)
    try:
        with membudget.surcharge(0.3, "zoning"):  # 3 GB, 1000 free
            assert b.used() == 10
    finally:
        membudget._bound = None


def test_release_of_an_unknown_key_is_harmless():
    b = cli._Budget(10, floor=0)
    b.release("never-acquired")
    assert b.used() == 0
