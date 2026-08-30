"""Mechanical enforcement of the codes/ layer boundary.

ADR-0002 says the checker does not compute; ADR-0005 adds that the equations do
not convert. Both were rules that depended on review discipline, which is
exactly the kind of boundary that erodes. These tests read the source with `ast`
so a violation fails CI instead of needing to be noticed.
"""

import ast
import subprocess
import sys
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


@pytest.mark.parametrize(
    "statement",
    [
        "from mento import RectangularBeam, Forces, Node",
        "from mento.codes.aci_318_19.equations import shear",
        "from mento.codes.en_1992_2004.equations import flexure",
        "import mento.rebar",
        "import mento.codes",
    ],
    ids=lambda s: s.split()[-1],
)
def test_each_entry_point_imports_in_a_fresh_interpreter(statement: str) -> None:
    """Import cycles hide from the suite, because pytest has already imported
    half the package by the time any test runs.

    A cycle between ``mento.rebar`` and the code modules once shipped green
    here for exactly that reason and only failed for a user typing
    ``from mento import RectangularBeam`` first. Each of these runs cold.
    """
    result = subprocess.run(
        [sys.executable, "-c", statement],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"`{statement}` fails from a cold interpreter:\n{result.stderr}"


def test_importing_the_code_layer_does_not_drag_in_presentation() -> None:
    """``mento.codes`` is calculation; matplotlib and docx belong to the report
    layer and must not be pulled in by importing it."""
    probe = "import sys, mento.codes; print('matplotlib' in sys.modules, 'docx' in sys.modules)"
    result = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "False False", f"mento.codes pulled in presentation: {result.stdout!r}"


MENTO_ROOT = EQUATIONS_ROOT.parent

# Every element: the classes a user builds and calls a check on.
ELEMENT_MODULES = [
    MENTO_ROOT / name
    for name in (
        "beam.py",
        "slab.py",
        "shear_wall.py",
        "punching.py",
        "beam_summary.py",
        "shear_wall_summary.py",
    )
]

PRESENTATION_LIBRARIES = {"matplotlib", "docx", "IPython"}


@pytest.mark.parametrize("path", ELEMENT_MODULES, ids=lambda p: p.name)
def test_elements_do_not_import_presentation_libraries(path: Path) -> None:
    """Phase 3: an element is geometry and orchestration, not a renderer.

    Drawing lives in ``mento.plots``, tables and reports in ``mento.reports``,
    and the element only ever delegates. Importing matplotlib, python-docx or
    IPython here is the signal that presentation has leaked back in.

    A ``TYPE_CHECKING`` import is fine — a return annotation costs nothing at
    runtime — so only real imports count.
    """
    tree = _parse(path)
    type_checking_only = {
        id(node)
        for parent in ast.walk(tree)
        if isinstance(parent, ast.If)
        and (
            (isinstance(parent.test, ast.Name) and parent.test.id == "TYPE_CHECKING")
            or (isinstance(parent.test, ast.Attribute) and parent.test.attr == "TYPE_CHECKING")
        )
        for node in ast.walk(parent)
    }

    offenders = []
    for node in ast.walk(tree):
        if id(node) in type_checking_only:
            continue
        if isinstance(node, ast.Import):
            offenders += [a.name for a in node.names if a.name.split(".")[0] in PRESENTATION_LIBRARIES]
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module.split(".")[0] in PRESENTATION_LIBRARIES:
                offenders.append(node.module)

    assert not offenders, (
        f"{path.name} imports presentation libraries {offenders}. "
        "Move the rendering into mento.plots or mento.reports and delegate to it."
    )


BEAM_CODE_MODULES = [
    EQUATIONS_ROOT / "ACI_318_19_beam.py",
    EQUATIONS_ROOT / "EN_1992_2004_beam.py",
]


@pytest.mark.parametrize("path", BEAM_CODE_MODULES, ids=lambda p: p.name)
def test_code_modules_do_not_build_report_tables(path: Path) -> None:
    """Phase 3: labelled, rounded, unit-annotated tables are presentation.

    They used to make up roughly a third of each code module. They live in
    ``mento/report_tables.py`` now, and the beam decides when to build them —
    the design code does not, which is what lets a check skip them entirely.

    ``ACI_318_19_wall.py`` is deliberately not in this list: the wall report
    has not been moved yet.
    """
    offenders = [
        node.name
        for node in ast.walk(_parse(path))
        if isinstance(node, ast.FunctionDef)
        and (node.name.startswith("_initialize_dicts") or node.name.startswith("_compile_results"))
    ]
    assert not offenders, f"{path.name} still builds report tables: {offenders}"


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
