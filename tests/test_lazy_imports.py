"""The presentation libraries are imported when they are used, not before.

``import mento`` and a whole design used to import matplotlib, seaborn, IPython
and python-docx, because the report and plot modules named them at the top. A
calculation needs none of them, and they were most of the start-up time -- and,
under Pyodide, most of the download. These tests run in a subprocess: the rest
of the suite has long since imported all four, and ``sys.modules`` of the test
process says nothing about a fresh interpreter.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

HEAVY_LIBRARIES = ("matplotlib", "seaborn", "IPython", "docx")

BUILD_AND_DESIGN = """
    import sys

    sys.path.insert(0, {repo!r})

    import mento
    from mento import Concrete_ACI_318_19, Concrete_EN_1992_2004, SteelBar, RectangularBeam, Node, Forces
    from mento import MPa, cm, mm, kN, kNm

    concrete = {concrete}(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="B500S", f_y=500 * MPa)
    beam = RectangularBeam(
        label="101", concrete=concrete, steel_bar=steel, width=20 * cm, height=60 * cm, c_c=25 * mm
    )
    node = Node(section=beam, forces=[Forces(label="C1", V_z=80 * kN, M_y=100 * kNm)])
    node.design()
    node.check()
    node.check_flexure()
    node.check_shear()
    beam.flexure_design, beam.shear_design, beam.reinforcement
"""

REPORT_LOADED = """
    loaded = [name for name in {heavy!r} if name in sys.modules]
    print(",".join(loaded))
"""


def _run(body: str, **fields: object) -> str:
    script = textwrap.dedent(body).format(repo=str(REPO_ROOT), heavy=HEAVY_LIBRARIES, **fields)
    completed = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        encoding="utf-8",
        cwd=REPO_ROOT,
        timeout=300,
    )
    assert completed.returncode == 0, completed.stderr
    return completed.stdout.strip().splitlines()[-1] if completed.stdout.strip() else ""


@pytest.mark.parametrize("concrete", ["Concrete_ACI_318_19", "Concrete_EN_1992_2004"])
def test_a_design_and_a_check_import_no_presentation_library(concrete: str) -> None:
    loaded = _run(BUILD_AND_DESIGN + REPORT_LOADED, concrete=concrete)
    assert loaded == "", f"a design imported {loaded}; import it inside the function that uses it"


def test_importing_every_public_name_imports_no_presentation_library() -> None:
    """The lazy ``__getattr__`` of the package reaches every element and summary."""
    body = """
    import sys

    sys.path.insert(0, {repo!r})

    import mento

    for name in mento.__all__:
        getattr(mento, name)
    """
    loaded = _run(body + REPORT_LOADED)
    assert loaded == "", f"importing the public API imported {loaded}"


def test_plot_is_what_imports_matplotlib() -> None:
    """The other half: the drawing still works, and is the thing that pays for it."""
    body = (
        BUILD_AND_DESIGN
        + """
    assert "matplotlib" not in sys.modules
    import os

    os.environ["MPLBACKEND"] = "Agg"
    beam.plot()
    print("matplotlib" in sys.modules)
    """
    )
    assert _run(body, concrete="Concrete_ACI_318_19") == "True"


def test_plot_without_matplotlib_says_what_to_install() -> None:
    body = (
        BUILD_AND_DESIGN
        + """
    import builtins

    real_import = builtins.__import__

    def no_matplotlib(name, *args, **kwargs):
        if name.split(".")[0] == "matplotlib":
            raise ModuleNotFoundError("No module named 'matplotlib'", name="matplotlib")
        return real_import(name, *args, **kwargs)

    builtins.__import__ = no_matplotlib
    try:
        beam.plot()
    except ImportError as error:
        print(error)
    """
    )
    message = _run(body, concrete="Concrete_ACI_318_19")
    assert "matplotlib" in message and "pip install" in message


def test_views_fall_back_to_print_without_ipython() -> None:
    body = (
        BUILD_AND_DESIGN
        + """
    import builtins

    real_import = builtins.__import__

    def no_ipython(name, *args, **kwargs):
        if name.split(".")[0] == "IPython":
            raise ModuleNotFoundError("No module named 'IPython'", name="IPython")
        return real_import(name, *args, **kwargs)

    builtins.__import__ = no_ipython
    beam.data
    """
    )
    printed = _run(body, concrete="Concrete_ACI_318_19")
    assert "Beam 101" in printed
