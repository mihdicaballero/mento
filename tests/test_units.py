import pytest
from pint import UnitRegistry

# Initialize UnitRegistry fixture
@pytest.fixture
def unit_registry() -> UnitRegistry:
    return UnitRegistry(system='mks')


def test_metric_units(unit_registry: UnitRegistry) -> None:
    ureg = unit_registry
    assert (2 * ureg.centimeter).to(ureg.meter).magnitude == 0.02
    assert (3 * ureg.megapascal).to(ureg.pascal).magnitude == 3_000_000
    assert (4 * ureg.kilogram).to(ureg.gram).magnitude == 4000
    assert (1 * ureg.millimeter).to(ureg.meter).magnitude == 0.001


def test_imperial_units(unit_registry: UnitRegistry) -> None:
    ureg = unit_registry
    assert (1 * ureg.psi).to(ureg.pascal).magnitude == pytest.approx(6894.76, 0.001)
    assert (1 * ureg.inch).to(ureg.foot).magnitude == 1 / 12
    assert (1 * ureg.kip).to(ureg.pound_force).magnitude == pytest.approx(1000, 0.001)


def test_derived_units(unit_registry: UnitRegistry) -> None:
    ureg = unit_registry
    N = 3 * ureg.kilonewton
    A = 3 * ureg.meter ** 2
    stress = (N / A).to(ureg.pascal)
    assert stress.magnitude == 1000


def test_conversion_and_formatting(unit_registry: UnitRegistry) -> None:
    ureg = unit_registry
    a = 3 * ureg.meter
    N = 3 * ureg.kilonewton
    assert N.to('kip').magnitude == pytest.approx(0.674, 0.001)
    assert '{:~P}'.format(a) == '3 m'


def test_frequency_calculation(unit_registry: UnitRegistry) -> None:
    ureg = unit_registry
    wavelength = 1550 * ureg.nanometer
    frequency = (ureg.speed_of_light / wavelength).to('Hz')
    assert frequency.magnitude == pytest.approx(1.934e14, 1e12)

# Run the tests using pytest framework
