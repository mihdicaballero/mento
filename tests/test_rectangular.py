import pytest
from typing import Generator

from mento.material import SteelBar, Concrete
from mento.units import inch, cm, mm, MPa
from mento.rectangular import RectangularSection  # Assuming this is the module name


@pytest.fixture()
def basic_rectangular_section(concrete_c25: Concrete, steel_b500s: SteelBar) -> RectangularSection:
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


@pytest.mark.parametrize(
    ("width", "height", "expected_exception", "message"),
    [
        pytest.param(0 * cm, 50 * cm, ValueError, "width", id="zero-width"),
        pytest.param(-30 * cm, 50 * cm, ValueError, "width", id="negative-width"),
        pytest.param(30 * MPa, 50 * cm, TypeError, "width", id="width-wrong-units"),
        pytest.param(30, 50 * cm, TypeError, "width", id="width-not-quantity"),
        pytest.param(
            float("nan") * cm,
            50 * cm,
            ValueError,
            "width",
            id="nan-width",
        ),
        pytest.param(30 * cm, 0 * cm, ValueError, "height", id="zero-height"),
        pytest.param(30 * cm, -50 * cm, ValueError, "height", id="negative-height"),
        pytest.param(30 * cm, 50 * MPa, TypeError, "height", id="height-wrong-units"),
        pytest.param(30 * cm, 50, TypeError, "height", id="height-not-quantity"),
        pytest.param(
            30 * cm,
            float("inf") * cm,
            ValueError,
            "height",
            id="infinite-height",
        ),
    ],
)
def test_rectangular_section_rejects_invalid_dimensions(
    concrete_c25: Concrete,
    steel_b500s: SteelBar,
    width,
    height,
    expected_exception,
    message,
) -> None:
    with pytest.raises(expected_exception, match=message):
        RectangularSection(
            concrete=concrete_c25,
            steel_bar=steel_b500s,
            c_c=25 * mm,
            width=width,
            height=height,
        )


@pytest.mark.parametrize(
    "c_c",
    [
        pytest.param(15 * cm, id="half-width"),
        pytest.param(16 * cm, id="greater-than-half-width"),
    ],
)
def test_rectangular_section_rejects_excessive_cover(
    concrete_c25: Concrete,
    steel_b500s: SteelBar,
    c_c,
) -> None:
    with pytest.raises(ValueError, match="c_c"):
        RectangularSection(
            concrete=concrete_c25,
            steel_bar=steel_b500s,
            c_c=c_c,
            width=30 * cm,
            height=50 * cm,
        )


def test_rectangular_section_accepts_cover_below_limit(
    concrete_c25: Concrete,
    steel_b500s: SteelBar,
) -> None:
    section = RectangularSection(
        concrete=concrete_c25,
        steel_bar=steel_b500s,
        c_c=14.9 * cm,
        width=30 * cm,
        height=50 * cm,
    )

    assert section.c_c == 14.9 * cm


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
