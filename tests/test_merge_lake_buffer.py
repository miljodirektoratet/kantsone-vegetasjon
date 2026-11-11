"""lake_buffer10 carries the same attributes as buffer10 (plus lake_id)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import merge_delivery as M  # noqa: E402


def test_lake_buffer10_has_every_column_of_buffer10_and_lake_id():
    assert set(M.BUFFER_FIELDS) <= set(M.LAKE_BUFFER_FIELDS)
    assert set(M.LAKE_BUFFER_FIELDS) - set(M.BUFFER_FIELDS) == {"lake_id"}


def test_the_shared_columns_come_in_the_same_order():
    shared = [k for k in M.BUFFER_FIELDS if k in M.LAKE_BUFFER_FIELDS]
    assert shared == [k for k in M.LAKE_BUFFER_FIELDS if k in M.BUFFER_FIELDS]


def test_both_buffers_are_documented_with_the_segment_condition():
    assert list(M.SEGMENT_VALUES) == ["E_S_Kantvegetasjon_fjernet", "E_S_Kantvegetasjon_grad"]
    assert M.SEGMENT_SOURCE == {
        M.BUFFER: "official_centerline_10m",
        M.LAKE_BUFFER: "lake_shore",
    }
