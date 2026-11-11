"""Every function the stages call must actually exist.

`merge_creeks_into_rivers` has now gone missing twice — once because the
original script was never committed, once because a later edit truncated the
module. Both times it stayed undetected until a run reached that line, ten
minutes in. This checks the whole call graph at import time instead.
"""

from __future__ import annotations

import ast
import inspect
import re
from pathlib import Path

import pytest

from kantsoner import stages

SRC = Path(stages.__file__)


def _stage_calls() -> list[tuple[str, str, int]]:
    """(module_alias, attribute, lineno) for every `mod.func(...)` in stages."""
    tree = ast.parse(SRC.read_text(encoding="utf-8"))
    out = []
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
        ):
            out.append((node.func.value.id, node.func.attr, node.lineno))
    return out


@pytest.mark.parametrize(
    "alias,attr,lineno", _stage_calls(), ids=lambda v: str(v)
)
def test_stage_call_target_exists(alias: str, attr: str, lineno: int) -> None:
    target = getattr(stages, alias, None)
    if not inspect.ismodule(target):
        pytest.skip(f"{alias} is not a module-level import")
    assert hasattr(target, attr), (
        f"stages.py:{lineno} calls {alias}.{attr}(), which does not exist in "
        f"{target.__name__}"
    )


def test_every_stage_command_is_implemented() -> None:
    from kantsoner.cli import STAGES

    for name in STAGES:
        assert callable(getattr(stages, name, None)), f"stages.{name} missing"


def test_no_module_imports_a_package_that_does_not_exist() -> None:
    """Function-local imports resolve only when the line runs, so a stale one
    survives both `ruff` and an import sweep. `kantsoner.sources` was merged
    into `kantsoner.io` and a local import in graph.py kept pointing at it —
    which only surfaced four stages into a run."""
    import ast
    import importlib.util
    from pathlib import Path

    root = Path(stages.__file__).parent
    bad = []
    for path in root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.ImportFrom) or not node.level:
                continue
            parts = path.relative_to(root).parts[:-1]
            base = ["kantsoner", *parts][: len(parts) + 1 - (node.level - 1)]
            target = ".".join([*base, node.module or ""]).rstrip(".")
            if importlib.util.find_spec(target) is None:
                bad.append(f"{path.name}:{node.lineno} -> {target}")
    assert not bad, (
        "imports pointing at modules that do not exist: " + "; ".join(bad)
    )


def test_no_frame_receiving_param_is_coerced_to_a_path() -> None:
    """L1/L2 hand intermediates to each other in memory, so a parameter that
    used to be a path can now be a GeoDataFrame. `Path(x)`, `str(x)` as a
    filename, or `list_layers(str(x))` on such a parameter raises only when
    that line runs — four stages into a 20-minute pipeline.

    This has now bitten seven times across four rounds of widening; the fix is
    always `as_gdf` / `input_layers`, so the raw forms are banned outright in
    the modules the stages call with frames.
    """
    from pathlib import Path

    root = Path(stages.__file__).parent
    banned = re.compile(
        r"list_layers\(str\(|Path\((input_path|river_path|segments_path|"
        r"rivers_path|centerlines_path|official_path|buf_0_10_path|"
        r"buf_0_20_path|elvevannforekomster_path|river_polygons_path|"
        r"buffer_left_path|buffer_right_path|vegetation_path|riparian_path)\)"
    )
    offenders = []
    for path in sorted(root.rglob("*.py")):
        for i, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), 1
        ):
            allowed = (
                "OUTPUT_PATH" in line
                or "isinstance" in line  # explicitly guarded
                or "eval_path" in line  # always paths.primary
                or "list_layers(str(p))"
                in line  # on-disk branch of a dict check
            )
            if banned.search(line) and not allowed:
                offenders.append(f"{path.name}:{i}: {line.strip()[:70]}")
    assert not offenders, (
        "these coerce a possibly-in-memory input to a path; use as_gdf() or "
        "input_layers():\n  " + "\n  ".join(offenders)
    )


def test_every_assigned_stage_call_returns_a_value() -> None:
    """`x = mod.fn(...)` in stages.py requires mod.fn to actually return.

    Converting the L1 functions from write-to-disk to return-a-frame missed
    `create_creek_lines`, so the next call received None and failed with
    "None: No such file or directory" — an error that points nowhere near the
    cause."""
    import ast
    from pathlib import Path

    root = Path(stages.__file__).parent
    tree = ast.parse(SRC.read_text(encoding="utf-8"))
    wanted = {
        (n.value.func.value.id, n.value.func.attr)
        for n in ast.walk(tree)
        if isinstance(n, ast.Assign)
        and isinstance(n.value, ast.Call)
        and isinstance(n.value.func, ast.Attribute)
        and isinstance(n.value.func.value, ast.Name)
    }
    missing = []
    for mod, fn in sorted(wanted):
        for path in root.rglob(f"{mod}.py"):
            for node in ast.parse(path.read_text(encoding="utf-8")).body:
                if isinstance(node, ast.FunctionDef) and node.name == fn:
                    returns = any(
                        isinstance(x, ast.Return) and x.value is not None
                        for x in ast.walk(node)
                    )
                    if not returns:
                        missing.append(f"{mod}.{fn}")
    assert not missing, (
        "stages.py assigns from these, but they never return a value: "
        + ", ".join(missing)
    )


def test_band_extractors_cut_zones_to_their_own_band() -> None:
    """Zone cells are built at the widest band, so an extractor whose columns
    are named B0-10 must cut them to the 0-10 ring itself — otherwise it
    measures 0-20 m of ground and labels the result 0-10."""
    import inspect

    from kantsoner.primary import landcover, vegetation

    for mod in (vegetation, landcover):
        assert "segment_zones.cut" in inspect.getsource(mod), (
            f"{mod.__name__} does not cut its zones to a band"
        )
