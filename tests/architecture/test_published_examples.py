"""A ``published_example`` mark is a public claim, so every marked test names its source.

The release workflow counts the tests marked ``published_example`` and mento-web
shows the number as "tests that reproduce a case validated outside mento". A mark
without a traceable source makes that number wrong, and 21 of the 52 tests once
marked asserted mento's own output. So every marked test
carries, in its docstring, a paragraph that starts with ``Source:`` and says which
document the number comes from and where in it it sits::

    Source: CRSI, Design Guide on the ACI 318 Building Code Requirements for
    Structural Concrete, §6.9.2, Example 6.2, p. 6-66, Table 6.24.
    Source: CSI Software Verification, ETABS, "ACI 318-19 Example 001", p. 2.
    Source: ETABS run recorded in BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet
    Flexion, row 29, column BC.

A Calcpad sheet is not a source. It reproduces a calculation with units, which
is how a contributor checks one, but it publishes nothing a reader can hold mento
to: the number has to come from the book, code, design guide or program run the
sheet follows. A ``Source:`` that names a Calcpad sheet fails, and the tests
that once rested on one stay in ``tests/`` as regression tests, unmarked.

The test modules are read with ``ast``, as ``test_architecture_boundaries`` does,
so a mark added without its source fails CI instead of needing to be noticed in
review. A last test checks that reading against what pytest itself collects with
``-m published_example``, which is the command the release workflow runs.
"""

import ast
import dataclasses
import re
import subprocess
import sys
from pathlib import Path

import pytest

TESTS_ROOT = Path(__file__).resolve().parents[1]
TEST_MODULES = sorted(TESTS_ROOT.rglob("test_*.py"))
#: Where the marked tests live, and only they: a reader goes straight to it.
VALIDATION = TESTS_ROOT / "validation"
MARK = "published_example"

# Where in the document the number is. A workbook of program runs needs the
# sheet and the row, column or cell; a publication needs the page, table,
# figure, equation, example or section.
_LOCATOR = re.compile(
    r"\bsheet\b|\brow\b|\bcolumn\b|\bcell\b|\bpage\b|\bp\.\s*\d|\bpp\.\s*\d|§"
    r"|\btable\b|\bfig\.|\bfigure\b|\beq\.|\bequation\b|\bexample\b|\bsection\b",
    re.IGNORECASE,
)
# A Calcpad sheet, by name or by file: a reproduction, not a source.
_CALCPAD = re.compile(r"\bcalcpad\b|\.cpd\b", re.IGNORECASE)


@dataclasses.dataclass(frozen=True)
class MarkedTest:
    module: Path
    name: str
    lineno: int
    docstring: str | None

    @property
    def id(self) -> str:
        return f"{self.module.name}::{self.name}"


def _parse(path: Path) -> ast.Module:
    # utf-8-sig, not utf-8: a file saved with a BOM on Windows would otherwise
    # raise SyntaxError here and fail these tests for the wrong reason.
    return ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))


def _is_the_mark(decorator: ast.expr) -> bool:
    """``pytest.mark.published_example`` or ``pytest.mark.published_example(...)``."""
    node = decorator.func if isinstance(decorator, ast.Call) else decorator
    return (
        isinstance(node, ast.Attribute)
        and node.attr == MARK
        and isinstance(node.value, ast.Attribute)
        and node.value.attr == "mark"
        and isinstance(node.value.value, ast.Name)
        and node.value.value.id == "pytest"
    )


def _marked_tests() -> list[MarkedTest]:
    found: list[MarkedTest] = []
    for path in TEST_MODULES:
        for node in ast.walk(_parse(path)):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if any(_is_the_mark(decorator) for decorator in node.decorator_list):
                found.append(MarkedTest(path, node.name, node.lineno, ast.get_docstring(node)))
    return found


MARKED = _marked_tests()


def _source_paragraph(docstring: str) -> str | None:
    """The paragraph that starts with ``Source:``, joined into one line."""
    lines = docstring.splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith("Source:"):
            paragraph: list[str] = []
            for following in lines[i:]:
                if not following.strip():
                    break
                paragraph.append(following.strip())
            return " ".join(paragraph)
    return None


def test_marked_tests_exist() -> None:
    """Guards the others: with no mark found they would pass vacuously, and the
    release would publish 0."""
    assert MARKED, f"no test marked {MARK} found under {TESTS_ROOT}"


@pytest.mark.parametrize("marked", MARKED, ids=lambda m: m.id)
def test_each_marked_test_cites_its_source(marked: MarkedTest) -> None:
    where = f"{marked.module.name}:{marked.lineno} {marked.name}"
    assert marked.docstring, f"{where} is marked {MARK} but has no docstring: say where its numbers come from"
    source = _source_paragraph(marked.docstring)
    assert source, f"{where} is marked {MARK} but its docstring has no 'Source:' paragraph"
    assert not _CALCPAD.search(source), (
        f"{where}: the 'Source:' paragraph names a Calcpad sheet, which reproduces a calculation "
        f"but does not publish one; cite the book, code, guide or program run it follows: {source!r}"
    )
    assert _LOCATOR.search(source), (
        f"{where}: the 'Source:' paragraph names no place in the document "
        f"(sheet + row/column/cell; page, table, figure, equation, example or section): {source!r}"
    )


def test_the_marked_tests_live_in_validation_and_nothing_else_does() -> None:
    """``tests/validation/`` holds the tests that reproduce a case validated outside mento.

    So a reader goes to one folder for the book, guide and program-run cases, and to
    the rest of ``tests/`` for the ones written against mento itself. A marked
    test elsewhere, or an unmarked one in there, fails.
    """
    misplaced = [m.id for m in MARKED if m.module.parent != VALIDATION]
    assert misplaced == [], f"marked {MARK} outside tests/validation/: {misplaced}"
    marked = {(m.module, m.name) for m in MARKED}
    unmarked = [
        f"{path.name}::{node.name}"
        for path in sorted(VALIDATION.glob("test_*.py"))
        for node in ast.walk(_parse(path))
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name.startswith("test_")
        and (path, node.name) not in marked
    ]
    assert unmarked == [], f"tests/validation/ holds tests not marked {MARK}: {unmarked}"


def test_ast_reading_matches_what_pytest_collects() -> None:
    """The release workflow counts ``pytest --collect-only -m published_example``.

    The marks read with ``ast`` have to be exactly the ones pytest applies, or a
    mark applied some other way (``pytestmark``, a class) would be published
    without ever being asked for its source.
    """
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "--collect-only",
            "-q",
            "--override-ini=addopts=",
            "-p",
            "no:cacheprovider",
            "-m",
            MARK,
            str(TESTS_ROOT),
        ],
        capture_output=True,
        text=True,
        cwd=TESTS_ROOT.parent,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    collected = set()
    for line in result.stdout.splitlines():
        if "::" not in line:
            continue
        module, test_id = line.strip().split("::", 1)
        collected.add(f"{Path(module).name}::{test_id.split('[', 1)[0]}")
    assert collected == {marked.id for marked in MARKED}
