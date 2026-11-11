"""Band restriction for the per-segment kantsoner.

The kantsoner themselves are built in `geometry.elvis`; what remains here is
the one operation every extractor shares — cutting a batch of zone cells to
the ring of the band a column is named for.
"""

import geopandas as gpd

from ..io import as_gdf


def cut(cells, ring):
    """Restrict per-segment zone cells to one band's ring.

    Cells are built at the widest band, so every extractor must cut them to
    the band its columns are named for. `stretch_id == id` keeps a segment
    from claiming the ring of a neighbouring stretch.
    """
    out = gpd.overlay(
        as_gdf(cells), as_gdf(ring)[["id", "geometry"]], how="intersection"
    )
    return out[out["stretch_id"] == out["id"]].drop(columns="id")
