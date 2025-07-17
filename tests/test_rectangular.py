import pytest
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Rectangle, FancyBboxPatch
import numpy as np
from typing import Generator

from mento.material import SteelBar, Concrete
from mento.units import inch, cm, mm, MPa
from mento.rectangular import RectangularSection  # Assuming this is the module name
from mento.settings import GLOBAL_BEAM_SETTINGS
from mento.results import CUSTOM_COLORS

@pytest.fixture()
def concrete_c25() -> Concrete:
    return Concrete("C25")

@pytest.fixture()
def steel_b500s() -> SteelBar:
    return SteelBar(name="B500S", f_y=500 * MPa)

@pytest.fixture()
def basic_rectangular_section(concrete_c25: Concrete, steel_b500s: SteelBar) -> RectangularSection:
    """Fixture for a basic RectangularSection instance."""
    return RectangularSection(
        label="R101",
        concrete=concrete_c25,
        steel_bar=steel_b500s,
        c_c=25 * mm,
        width=30 * cm, # Example dimensions
        height=50 * cm
    )

@pytest.fixture()

def imperial_rectangular_section(concrete_c25: Concrete, steel_b500s: SteelBar) -> Generator[RectangularSection, None, None]:
    """Fixture for a RectangularSection with imperial units set in concrete."""
    # Temporarily override unit system for concrete for this fixture
    # In a real app, you'd likely configure this via settings or concrete constructor
    concrete_c25.unit_system = "imperial"
    section = RectangularSection(
        label="R102",
        concrete=concrete_c25,
        steel_bar=steel_b500s,
        c_c=1 * inch, # Example imperial cover
        width=12 * inch, # Example imperial width
        height=24 * inch
    )
    yield section # Use yield to ensure cleanup happens after the test
    concrete_c25.unit_system = "metric" # Reset for other tests if fixture is reused

# --- Tests for RectangularSection ---

def test_rectangular_section_initialization(basic_rectangular_section: RectangularSection) -> None:
    """Test RectangularSection specific initialization."""
    assert basic_rectangular_section.width == 30 * cm
    assert basic_rectangular_section.height == 50 * cm
    assert basic_rectangular_section.label == "R101"
    assert basic_rectangular_section.c_c == 25 * mm

def test_rectangular_section_geometric_properties(basic_rectangular_section: RectangularSection) -> None:
    """Test calculation of geometric properties (A_x, I_y, I_z)."""
    expected_A_x = (30 * cm) * (50 * cm)
    assert basic_rectangular_section.A_x == expected_A_x.to("cm**2")

    expected_I_y = (30 * cm) * (50 * cm)**3 / 12
    assert basic_rectangular_section.I_y == expected_I_y.to("cm**4")

    expected_I_z = (50 * cm) * (30 * cm)**3 / 12
    assert basic_rectangular_section.I_z == expected_I_z.to("cm**4")


def test_rectangular_section_settings_property(basic_rectangular_section: RectangularSection) -> None:
    """Test that the settings property returns GLOBAL_BEAM_SETTINGS."""
    assert basic_rectangular_section.settings is GLOBAL_BEAM_SETTINGS

# --- Test for the plot method ---

def test_rectangular_section_plot_smoke_test(basic_rectangular_section: RectangularSection) -> None:
    """
    Smoke test for the plot method.
    Simply calls the method and ensures no exceptions are raised.
    Closes the plot to prevent it from showing during tests.
    """
    try:
        basic_rectangular_section.plot()
        # Close the plot to prevent it from popping up
        plt.close()
    except Exception as e:
        pytest.fail(f"Plot method raised an exception: {e}")

def test_rectangular_section_plot_components(basic_rectangular_section: RectangularSection) -> None:
    """
    Test that the plot method adds the expected matplotlib components
    (main rectangle, stirrup rectangles, annotations).
    """
    basic_rectangular_section.plot()
    ax = basic_rectangular_section._ax # Access the private _ax created by plot

    assert ax is not None, "Plot method should create and assign an Axes object."

    # Check for patches (rectangles, FancyBboxPatch)
    # Expecting: 1 main Rectangle, 2 FancyBboxPatches for stirrup
    # + potentially other internal patches depending on matplotlib's rendering
    # It's safer to check for specific types and properties if possible.

    # Filter for Rectangle and FancyBboxPatch
    rectangles = [p for p in ax.patches if isinstance(p, Rectangle)]
    fancy_bboxes = [p for p in ax.patches if isinstance(p, FancyBboxPatch)]

    assert len(rectangles) >= 1, "Expected at least one Rectangle patch (the main section)."
    assert len(fancy_bboxes) == 2, "Expected exactly two FancyBboxPatch objects for the stirrup."

    # Check for annotations (dimensions)
    # Annotations are stored in ax.texts for text, ax.lines or ax.artists for arrows
    # matplotlib.axes.Axes.annotate returns an Annotation object.
    # Text labels are found in ax.texts
    assert len(ax.texts) == 4, "Expected 4 text annotations for dimensions (width, height)."

    # You could add more specific checks, e.g.:
    # - Check the coordinates of the main rectangle:
    main_rect = next((p for p in rectangles if p.get_x() == 0 and p.get_y() == 0), None)
    assert main_rect is not None
    assert np.isclose(main_rect.get_width(), basic_rectangular_section.width.magnitude)
    assert np.isclose(main_rect.get_height(), basic_rectangular_section.height.magnitude)
    assert main_rect.get_edgecolor() == to_rgba(CUSTOM_COLORS["dark_gray"])

    # - Check stirrup patch properties (example, will need exact dimensions based on calculation)
    #   This gets more complex quickly, focusing on existence and basic properties is often enough.

    plt.close() # Close the figure created by the plot method

def test_rectangular_section_plot_imperial_units_text(imperial_rectangular_section: RectangularSection) -> None:
    """
    Test that the plot method displays dimensions in imperial units
    when the concrete's unit system is set to 'imperial'.
    """
    imperial_rectangular_section.plot()
    ax = imperial_rectangular_section._ax
    assert ax is not None, "The plot method did not assign an Axes object to _ax."

    # Extract text labels
    text_labels = [t.get_text() for t in ax.texts]

    # Check if the labels contain "inch" (represented by "in" in pint's compact format)
    # Adjust based on how pint formats "inch" in your setup (e.g., "in", "inch", etc.)
    # "{:~P}" often uses common abbreviations.
    assert any("in" in label for label in text_labels), "Expected imperial unit 'in' in dimension text."
    assert "{:.0f~P}".format(imperial_rectangular_section.width.to("inch")) in text_labels
    assert "{:.0f~P}".format(imperial_rectangular_section.height.to("inch")) in text_labels

    plt.close()

def test_rectangular_section_plot_metric_units_text(basic_rectangular_section: RectangularSection) -> None:
    """
    Test that the plot method displays dimensions in metric units (cm)
    when the concrete's unit system is not 'imperial'.
    """
    basic_rectangular_section.plot()  # Ensure the plot method is called
    ax = basic_rectangular_section._ax
    assert ax is not None, "The plot method did not assign an Axes object to _ax."

    # Extract text labels
    text_labels = [t.get_text() for t in ax.texts]

    # Check if the labels contain "cm"
    assert any("cm" in label for label in text_labels), "Expected metric unit 'cm' in dimension text."
    assert "{:.0f~P}".format(basic_rectangular_section.width.to("cm")) in text_labels
    assert "{:.0f~P}".format(basic_rectangular_section.height.to("cm")) in text_labels

    plt.close()
