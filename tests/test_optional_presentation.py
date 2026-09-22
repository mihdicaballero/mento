"""What happens when a presentation library is missing, run in-process.

``tests/test_lazy_imports.py`` proves in a subprocess that nothing imports these
libraries early. This module covers the other side -- the fallbacks -- by making
the import fail, which a subprocess cannot report coverage for.
"""

from __future__ import annotations

import builtins
import sys
from typing import Any, Callable

import pytest

import mento.codes
from mento import Column, Concrete_ACI_318_19, Forces, PunchingNode, PunchingSlab, RectangularBeam, ShearWall, SteelBar
from mento import MPa, cm, kN, m, mm
from mento.plots import plotting_import_error
from mento.reports import _notebook, summaries
from mento import results
from mento.results import DocumentBuilder, configure_plot_settings


@pytest.fixture()
def without(monkeypatch: pytest.MonkeyPatch) -> Callable[[str], None]:
    """Make ``import <library>`` fail, as it does where the library is not installed."""

    def block(library: str) -> None:
        real_import = builtins.__import__

        def fake_import(name: str, *args: Any, **kwargs: Any) -> Any:
            if name.split(".")[0] == library:
                raise ModuleNotFoundError(f"No module named {library!r}", name=library)
            return real_import(name, *args, **kwargs)

        # The plot modules are cached once anything has drawn; drop them so the
        # import inside plot() runs again and meets the missing library.
        for module in [name for name in sys.modules if name.startswith("mento.plots.")]:
            monkeypatch.delitem(sys.modules, module)
        monkeypatch.setattr(builtins, "__import__", fake_import)

    return block


@pytest.fixture()
def concrete() -> Concrete_ACI_318_19:
    return Concrete_ACI_318_19(name="C25", f_c=25 * MPa)


@pytest.fixture()
def steel_bar() -> SteelBar:
    return SteelBar(name="ADN420", f_y=420 * MPa)


# ============================================================================
# matplotlib
# ============================================================================


def test_a_missing_plotting_library_is_reworded() -> None:
    error = plotting_import_error(ModuleNotFoundError("No module named 'matplotlib'", name="matplotlib.pyplot"))
    assert "pip install matplotlib" in str(error)
    assert error.name == "matplotlib.pyplot"


def test_any_other_import_error_is_left_alone() -> None:
    original = ModuleNotFoundError("No module named 'mento.plots.nope'", name="mento.plots.nope")
    assert plotting_import_error(original) is original


def test_beam_plot_without_matplotlib(
    without: Callable[[str], None], concrete: Concrete_ACI_318_19, steel_bar: SteelBar
) -> None:
    beam = RectangularBeam(
        label="B", concrete=concrete, steel_bar=steel_bar, width=20 * cm, height=50 * cm, c_c=25 * mm
    )
    without("matplotlib")
    with pytest.raises(ImportError, match="pip install matplotlib"):
        beam.plot()


def test_wall_plot_without_matplotlib(
    without: Callable[[str], None], concrete: Concrete_ACI_318_19, steel_bar: SteelBar
) -> None:
    wall = ShearWall(
        label="W", concrete=concrete, steel_bar=steel_bar, thickness=20 * cm, length=3 * m, height=3 * m, c_c=20 * mm
    )
    without("matplotlib")
    with pytest.raises(ImportError, match="pip install matplotlib"):
        wall.plot()


def test_punching_plot_without_matplotlib(
    without: Callable[[str], None], concrete: Concrete_ACI_318_19, steel_bar: SteelBar
) -> None:
    slab = PunchingSlab(concrete=concrete, steel_bar=steel_bar, h=25 * cm, c_c=25 * mm)
    column = Column(shape="rectangular", b=40 * cm, h=40 * cm, position="interior")
    node = PunchingNode(slab=slab, column=column, forces=Forces(label="ELU", V_z=500 * kN))
    without("matplotlib")
    with pytest.raises(ImportError, match="pip install matplotlib"):
        node.plot()


def test_plot_settings_without_seaborn(without: Callable[[str], None]) -> None:
    without("seaborn")
    with pytest.raises(ImportError, match="Plotting needs seaborn"):
        configure_plot_settings()


def test_add_figure_without_matplotlib(without: Callable[[str], None]) -> None:
    builder = DocumentBuilder(title="Figure")
    without("matplotlib")
    # The import is the first thing add_figure does, so the figure is never touched.
    with pytest.raises(ImportError, match="pip install matplotlib"):
        builder.add_figure(object())  # type: ignore[arg-type]


# ============================================================================
# IPython
# ============================================================================


def test_markdown_without_ipython_is_plain_text(without: Callable[[str], None]) -> None:
    without("IPython")
    rendered = _notebook.Markdown("**DCR** = 0.9")
    assert rendered.data == "**DCR** = 0.9"
    assert str(rendered) == "**DCR** = 0.9"


def test_display_without_ipython_prints(without: Callable[[str], None], capsys: pytest.CaptureFixture[str]) -> None:
    without("IPython")
    _notebook.display(_notebook.Markdown("Beam 101"))
    assert capsys.readouterr().out == "Beam 101\n"


def test_markdown_with_ipython_is_ipythons() -> None:
    from IPython.display import Markdown

    assert isinstance(_notebook.Markdown("x"), Markdown)


# ============================================================================
# Names that are built on access
# ============================================================================


def test_the_width_constants_are_still_importable() -> None:
    from docx.shared import Cm

    from mento.results import DETAIL_TABLE_WIDTHS, LIMIT_TABLE_WIDTHS

    assert DETAIL_TABLE_WIDTHS == [Cm(6.1), Cm(2.2), Cm(3.0), Cm(1.4)]
    assert len(LIMIT_TABLE_WIDTHS) == 6
    assert summaries.BEAM_DATA_WIDTHS == [Cm(2), Cm(1), Cm(1), Cm(1)] + [Cm(0.9)] * 11


@pytest.mark.parametrize("module", [results, summaries])
def test_an_unknown_module_attribute_is_still_an_attribute_error(module: Any) -> None:
    with pytest.raises(AttributeError, match="NOT_A_WIDTH"):
        module.NOT_A_WIDTH


def test_a_code_module_is_reachable_from_the_package(monkeypatch: pytest.MonkeyPatch) -> None:
    # A submodule is an attribute of its package only once it has been imported,
    # which in this process it long has; take it away, as in a fresh interpreter.
    monkeypatch.delattr(mento.codes, "ACI_318_19_beam", raising=False)
    assert mento.ACI_318_19_beam.__name__ == "mento.codes.ACI_318_19_beam"
    assert mento.EN_1992_2004_beam.__name__ == "mento.codes.EN_1992_2004_beam"
