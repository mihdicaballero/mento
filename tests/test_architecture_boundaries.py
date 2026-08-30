"""Mechanical enforcement of the codes/ layer boundary.

ADR-0002 says the checker does not compute; ADR-0005 adds that the equations do
not convert. Both were rules that depended on review discipline, which is
exactly the kind of boundary that erodes. These tests read the source with `ast`
so a violation fails CI instead of needing to be noticed.
"""

import ast
from pathlib import Path

import pytest

EQUATIONS_ROOT = Path(__file__).resolve().parent.parent / "mento" / "codes"

EQUATION_MODULES = sorted(EQUATIONS_ROOT.glob("*/equations/*.py"))

# Modules that carry units are the boundary; equations must stay below it.
UNIT_MODULES = {"pint", "mento.units"}


def _parse(path: Path) -> ast.Module:
    # utf-8-sig, not utf-8: a file saved with a BOM on Windows would otherwise
    # raise SyntaxError here and fail these tests for the wrong reason.
    return ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))


def test_equation_modules_exist():
    """Guards the other tests: an empty glob would make them vacuously pass."""
    assert EQUATION_MODULES, f"no equation modules found under {EQUATIONS_ROOT}"


@pytest.mark.parametrize("path", EQUATION_MODULES, ids=lambda p: p.name)
def test_equations_do_not_import_units(path: Path) -> None:
    """ADR-0005: equations take and return floats, so they never see a Quantity."""
    offenders = []
    for node in ast.walk(_parse(path)):
        if isinstance(node, ast.Import):
            offenders += [a.name for a in node.names if a.name.split(".")[0] in {"pint"}]
        elif isinstance(node, ast.ImportFrom) and node.module:
            root = node.module.split(".")[0]
            if node.module in UNIT_MODULES or root == "pint":
                offenders.append(node.module)
    assert not offenders, f"{path.name} imports units: {offenders}. Equations work in plain floats (ADR-0005)."


@pytest.mark.parametrize("path", EQUATION_MODULES, ids=lambda p: p.name)
def test_equations_are_free_functions(path: Path) -> None:
    """No `self`, so no element state can leak into an equation."""
    tree = _parse(path)
    classes = [n.name for n in ast.walk(tree) if isinstance(n, ast.ClassDef)]
    assert not classes, f"{path.name} defines classes {classes}; equations are free functions."

    selfish = [
        n.name
        for n in ast.walk(tree)
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and any(a.arg == "self" for a in n.args.args)
    ]
    assert not selfish, f"{path.name} has functions taking self: {selfish}."


@pytest.mark.parametrize("path", EQUATION_MODULES, ids=lambda p: p.name)
def test_every_public_equation_cites_its_clause(path: Path) -> None:
    """A formula whose clause is not written down cannot be reviewed against the code."""
    missing = []
    for node in ast.walk(_parse(path)):
        if isinstance(node, ast.FunctionDef) and not node.name.startswith("_"):
            doc = ast.get_docstring(node) or ""
            if "ACI" not in doc and "EN 199" not in doc and "CIRSOC" not in doc:
                missing.append(node.name)
    assert not missing, f"{path.name}: no code clause cited in the docstring of {missing}."
