"""The documentation names tests that exist.

The theory pages end in a validation table that names, row by row, the test that
pins each rule. A test renamed without its citation leaves the page vouching for
a rule with nothing behind it: a2d7857 renamed
``test_footing_bars_are_capped_at_300_mm`` when the cap stopped being a flat
300 mm, fixed the citation in ``one_way_slab.rst`` and left the one in
``footing.rst``. Reading the pages here catches the next one in CI.

A citation is a double-backquoted literal that is a test name, ``test_…``; one
that ends in ``*`` names a family, and needs some test that starts with it. Paths
such as ``tests/test_beam.py`` are not test names and are left alone. The test
modules are read with ``ast``, as ``test_published_examples`` does, so a method
of a test class counts as much as a module-level function.
"""

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs" / "source"
_CITATION = re.compile(r"``(test_\w+)(\*?)``")


def _defined_tests() -> set[str]:
    names: set[str] = set()
    for module in (ROOT / "tests").glob("test_*.py"):
        tree = ast.parse(module.read_text(encoding="utf-8"))
        names.update(
            node.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_")
        )
    return names


def test_every_test_the_documentation_cites_exists() -> None:
    defined = _defined_tests()
    missing = []
    for page in sorted(DOCS.rglob("*.rst")):
        for name, family in _CITATION.findall(page.read_text(encoding="utf-8")):
            found = any(test.startswith(name) for test in defined) if family else name in defined
            if not found:
                missing.append(f"{page.relative_to(ROOT).as_posix()}: {name}{family}")
    assert missing == []


def test_a_family_citation_needs_a_test_that_starts_with_it() -> None:
    """The reading itself: a name, a family, and a path that is neither."""
    text = "``test_one`` and ``test_fam_*`` but not ``tests/test_beam.py``"
    assert _CITATION.findall(text) == [("test_one", ""), ("test_fam_", "*")]
