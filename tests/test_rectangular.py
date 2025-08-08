import pytest
from typing import Generator

from mento.material import SteelBar, Concrete
from mento.units import inch, cm, mm, MPa
from mento.rectangular import RectangularSection  # Assuming this is the module name


@pytest.fixture()
def concrete_c25() -> Concrete:
    return Concrete("C25")


@pytest.fixture()
def steel_b500s() -> SteelBar:
    return SteelBar(name="B500S", f_y=500 * MPa)


@pytest.fixture()
def basic_rectangular_section(
    concrete_c25: Concrete, steel_b500s: SteelBar
) -> RectangularSection:
    """Fixture for a basic RectangularSection instance."""
    return RectangularSection(
        label="R101",
        concrete=concrete_c25,
        steel_bar=steel_b500s,
        c_c=25 * mm,
        width=30 * cm,  # Example dimensions
        height=50 * cm,
    )


@pytest.fixture()
def imperial_rectangular_section(
    concrete_c25: Concrete, steel_b500s: SteelBar
) -> Generator[RectangularSection, None, None]:
    """Fixture for a RectangularSection with imperial units set in concrete."""
    # Temporarily override unit system for concrete for this fixture
    # In a real app, you'd likely configure this via settings or concrete constructor
    concrete_c25.unit_system = "imperial"
    section = RectangularSection(
        label="R102",
        concrete=concrete_c25,
        steel_bar=steel_b500s,
        c_c=1 * inch,  # Example imperial cover
        width=12 * inch,  # Example imperial width
        height=24 * inch,
    )
    yield section  # Use yield to ensure cleanup happens after the test
    concrete_c25.unit_system = "metric"  # Reset for other tests if fixture is reused


# --- Tests for RectangularSection ---


def test_rectangular_section_initialization(
    basic_rectangular_section: RectangularSection,
) -> None:
    """Test RectangularSection specific initialization."""
    assert basic_rectangular_section.width == 30 * cm
    assert basic_rectangular_section.height == 50 * cm
    assert basic_rectangular_section.label == "R101"
    assert basic_rectangular_section.c_c == 25 * mm


def test_rectangular_section_geometric_properties(
    basic_rectangular_section: RectangularSection,
) -> None:
    """Test calculation of geometric properties (A_x, I_y, I_z)."""
    expected_A_x = (30 * cm) * (50 * cm)
    assert basic_rectangular_section.A_x == expected_A_x.to("cm**2")

    expected_I_y = (30 * cm) * (50 * cm) ** 3 / 12
    assert basic_rectangular_section.I_y == expected_I_y.to("cm**4")

    expected_I_z = (50 * cm) * (30 * cm) ** 3 / 12
    assert basic_rectangular_section.I_z == expected_I_z.to("cm**4")
