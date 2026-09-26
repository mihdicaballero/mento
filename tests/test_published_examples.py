"""A ``published_example`` mark is a public claim, so every marked test names its source.

The release workflow counts the tests marked ``published_example`` and mento-web
shows the number as "tests that reproduce a case validated outside mento". A mark
without a traceable source makes that number wrong, and 21 of the 52 tests once
marked asserted mento's own output. So every marked test
carries, in its docstring, a paragraph that starts with ``Source:`` and says which
document the number comes from and where in it it sits::

    Source: BEAM-01-Flexure-Rectangle ACI 318-19-v6.xlsm, sheet Flexion, row 29, column BC.
    Source: Calcpad "ACI 318-19 Beam Shear 01 - Imperial.cpd".
    Source: The Concrete Centre, "How to design concrete structures using Eurocode 2",
    3. Slabs, p. 3, Table 5.

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

TESTS_ROOT = Path(__file__).resolve().parent
TEST_MODULES = sorted(TESTS_ROOT.glob("test_*.py"))
MARK = "published_example"

# Where in the document the number is. A Calcpad sheet is one computation, so
# its file name is enough; a workbook needs the sheet and the row, column or
# cell; a publication needs the page, table, figure, equation, example or section.
_LOCATOR = re.compile(
    r"\.cpd\b|\bsheet\b|\brow\b|\bcolumn\b|\bcell\b|\bpage\b|\bp\.\s*\d|\bpp\.\s*\d|§"
    r"|\btable\b|\bfig\.|\bfigure\b|\beq\.|\bequation\b|\bexample\b|\bsection\b",
    re.IGNORECASE,
)


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
    assert _LOCATOR.search(source), (
        f"{where}: the 'Source:' paragraph names no place in the document "
        f"(a .cpd sheet; sheet + row/column/cell; page, table, figure, equation or section): {source!r}"
    )


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
