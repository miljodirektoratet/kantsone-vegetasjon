"""The layering rule, enforced.

`derived/` must never import geopandas or rasterio. That single constraint is
what keeps the expensive layer separate from the cheap one and makes the whole
system-B method testable with a ten-row fixture and no GIS data.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

DERIVED = Path(__file__).resolve().parents[1] / "src" / "kantsoner" / "derived"
BANNED = {"geopandas", "rasterio", "fiona", "pyogrio", "shapely", "osgeo"}


def _imported_modules(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(a.name.split(".")[0] for a in node.names)
        elif (
            isinstance(node, ast.ImportFrom) and node.level == 0 and node.module
        ):
            found.add(node.module.split(".")[0])
    return found


@pytest.mark.parametrize(
    "path", sorted(DERIVED.glob("*.py")), ids=lambda p: p.name
)
def test_derived_imports_no_gis(path: Path) -> None:
    offending = _imported_modules(path) & BANNED
    assert not offending, (
        f"{path.name} imports {sorted(offending)}. L3 must stay pure tabular — "
        "move any geometry work into primary/."
    )


def test_spec_imports_nothing_from_the_package() -> None:
    """spec/ sits at the bottom: it may not import sibling layers."""
    spec = DERIVED.parent / "spec"
    for path in spec.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level > 1:
                pytest.fail(
                    f"{path.name} reaches outside spec/: {ast.dump(node)[:80]}"
                )
