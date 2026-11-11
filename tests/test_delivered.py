"""Which centerline segments reach the delivery (stages._delivered).

The delivery is coupled to Vann-Nett on vannforekomstid, so a segment without
one goes. A lake midline goes only inside a Vann-Nett lake, which is delivered
as a lake under its own -L id; any other lake belongs to a river water body
and its midline goes out with that river. A region with no Vann-Nett lake, or
a run without the lake layer, has no lake_id column at all -- that must not
read as "every midline is inside a lake".
"""

import pandas as pd

from kantsoner import stages
from kantsoner.spec import constants as C

MID, RIVER = "InnsjøMidtlinje", "ElvBekk"


def test_midline_with_a_river_id_stays_without_lake_info():
    gdf = pd.DataFrame(
        {
            "objektType": [MID, MID, RIVER, RIVER],
            "vannforekomstid": ["123-386-R", None, "123-386-R", None],
        }
    )
    assert stages._delivered(gdf).tolist() == [True, False, True, False]


def test_midline_inside_a_vann_nett_lake_goes():
    gdf = pd.DataFrame(
        {
            "objektType": [MID, MID, RIVER],
            "vannforekomstid": ["027-186-R", "027-186-R", "027-186-R"],
            "lake_id": ["027-1524-L", None, None],
        }
    )
    assert stages._delivered(gdf).tolist() == [False, True, True]


def test_lake_shore_segments_go_to_the_lake_files():
    gdf = pd.DataFrame(
        {
            "objektType": [RIVER, None],
            "vannforekomstid": ["027-186-R", None],
            "lake_id": [None, "027-1524-L"],
            C.KILDE_COL: [C.KILDE_RIVER, C.KILDE_LAKE],
        }
    )
    assert stages._delivered(gdf).tolist() == [True, False]
