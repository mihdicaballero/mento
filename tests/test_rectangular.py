import pytest
from mento.material import SteelBar, Concrete
from mento.units import ksi, inch
from mento.rectangular import RectangularSection  # Assuming this is the module name


@pytest.fixture()
def setup_section() -> RectangularSection:
    concrete = Concrete("C25")
    steel_bar = SteelBar(name="ADN 420", f_y=60 * ksi)
    width = 10 * inch
    height = 16 * inch
    section = RectangularSection(
        label="V101",
        concrete=concrete,
        steel_bar=steel_bar,
        width=width,
        height=height,
        c_c=1.5 * inch,
    )
    return section


def test_width(setup_section: RectangularSection) -> None:
    section = setup_section
    expected_width = 10 * inch
    assert section.width.magnitude == expected_width.magnitude
    assert str(section.width.units) == str(expected_width.units)


def test_height(setup_section: RectangularSection) -> None:
    section = setup_section
    expected_height = 16 * inch
    assert section.height.magnitude == expected_height.magnitude
    assert str(section.height.units) == str(expected_height.units)


def test_area(setup_section: RectangularSection) -> None:
    section = setup_section
    assert section.A_x.to("inch**2").magnitude == pytest.approx(10 * 16, rel=1e-3)


def test_moment_of_inertia_y(setup_section: RectangularSection) -> None:
    section = setup_section
    assert section.I_y.to("inch**4").magnitude == pytest.approx(
        10 * (16**3) / 12, rel=1e-3
    )


def test_moment_of_inertia_z(setup_section: RectangularSection) -> None:
    section = setup_section
    assert section.I_z.to("inch**4").magnitude == pytest.approx(
        16 * (10**3) / 12, rel=1e-3
    )
