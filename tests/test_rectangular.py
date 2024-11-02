import pytest
from mento.material import SteelBar, Concrete
from mento.units import ksi, inch, cm
from mento.rectangular import RectangularSection  # Assuming this is the module name


@pytest.fixture()
def setup_section() -> RectangularSection:
    concrete = Concrete('C25')
    steel_bar = SteelBar(name="ADN 420", f_y=60 * ksi)
    width = 10 * inch
    height = 16 * inch
    section = RectangularSection(concrete=concrete, steel_bar=steel_bar, width=width, height=height)
    return section

def test_width(setup_section: RectangularSection) -> None:
    section = setup_section
    expected_width = 10 * inch
    assert section.width.magnitude == expected_width.to(cm).magnitude
    assert str(section.width.units) == str(expected_width.to(cm).units)


def test_height(setup_section: RectangularSection) -> None:
    section = setup_section
    expected_height = 16 * inch
    assert section.height.magnitude == expected_height.to(cm).magnitude
    assert str(section.height.units) == str(expected_height.to(cm).units)


def test_area(setup_section: RectangularSection) -> None:
    section = setup_section
    expected_area = (10 * inch) * (16 * inch)
    assert section.A_x.magnitude == expected_area.to(cm**2).magnitude
    assert str(section.A_x.units) == str(expected_area.to(cm**2).units)


def test_moment_of_inertia_y(setup_section: RectangularSection) -> None:
    section = setup_section
    expected_I_y = ((10 * inch) * (16 * inch) ** 3) / 12
    assert section.I_y.magnitude == expected_I_y.to(cm**4).magnitude
    assert str(section.I_y.units) == str(expected_I_y.to(cm**4).units)


def test_moment_of_inertia_z(setup_section: RectangularSection) -> None:
    section = setup_section
    expected_I_z = ((16 * inch) * (10 * inch) ** 3) / 12
    assert section.I_z.magnitude == expected_I_z.to(cm**4).magnitude
    assert str(section.I_z.units) == str(expected_I_z.to(cm**4).units)

def test_default_d_value(setup_section: RectangularSection) -> None:
    section = setup_section
    # Here we check the initial value calculation for `d`, assuming `c_c`, `_stirrup_d_b`, and `_long_d_b` are set.
    expected_d = section._height - (section.c_c + section._stirrup_d_b + section._long_d_b / 2)
    assert section._d == expected_d