"""Kantsoner — riparian zone morphology from remote sensing.

Layering; dependencies run one way only::

    spec/     imports nothing
    io        imports only spec
    geometry -> primary -> derived -> network -> delivery

The stage names differ from the package names: `attributes` runs primary/,
`derived` runs derived/ + network/, `aggregation` runs delivery/.

`derived/` must never import geopandas or rasterio; that rule keeps the
expensive layer separate from the cheap one and is enforced by tests/.
"""

import os
from pathlib import Path

# Some machines export PROJ_LIB=/usr/share/proj/, whose proj.db predates what
# the wheels expect — GDAL then silently loses EPSG resolution and any write
# with a declared CRS fails. rasterio bundles a current copy; prefer it.
_bundled = Path(__import__("rasterio").__file__).parent / "proj_data"
if _bundled.exists():
    os.environ["PROJ_LIB"] = os.environ["PROJ_DATA"] = str(_bundled)

__version__ = "1.0.0"
