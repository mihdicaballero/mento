import pytest
from pint import Quantity
from pint.facets.plain import PlainQuantity
from typing import Generator

from mento.forces import Forces
from mento.units import kN, kNm, kip

# Helper to check if a quantity matches another, allowing for slight float differences
# Option A: Explicitly allow PlainQuantity
def assert_quantity_equal(q1: Quantity | PlainQuantity, q2: Quantity | PlainQuantity, rtol: float = 1e-9, atol: float = 1e-12) -> None:
    """Asserts two Pint Quantities are numerically and unit-wise equal."""
    assert q1.units == q2.units, f"Units mismatch: {q1.units} (actual) vs {q2.units} (expected)"
    assert pytest.approx(q1.magnitude, rel=rtol, abs=atol) == q2.magnitude, \
        f"Magnitudes mismatch: {q1.magnitude} (actual) vs {q2.magnitude} (expected) with units {q1.units}"

# --- Fixtures ---

@pytest.fixture(autouse=True)
def reset_forces_id_counter() -> Generator[None, None, None]:
    """Fixture to reset the Forces._last_id counter before each test."""
    original_last_id = Forces._last_id
    Forces._last_id = 0
    yield # Allow the test to run
    Forces._last_id = original_last_id # Restore if needed, though usually not critical after tests

@pytest.fixture
def default_metric_forces() -> Forces:
    """A Forces instance with default metric values."""
    return Forces()

@pytest.fixture
def custom_metric_forces() -> Forces:
    """A Forces instance with custom metric values."""
    return Forces(label="Custom Metric", N_x=100 * kN, V_z=50 * kN, M_y=200 * kNm)

@pytest.fixture
def custom_imperial_forces() -> Forces:
    """A Forces instance with custom imperial values."""
    # Note: When creating, you can use any compatible units. The class will convert for display.
    return Forces(label="Custom Imperial", N_x=20 * kip, V_z=10 * kip, M_y=150 * kNm, unit_system="imperial")


# --- Tests for Initialization (`__init__`) ---

def test_forces_default_initialization(default_metric_forces: Forces) -> None:
    """Test Forces initialization with default values."""
    assert default_metric_forces.label is None
    assert_quantity_equal(default_metric_forces._N_x, 0 * kN)
    assert_quantity_equal(default_metric_forces._V_z, 0 * kN)
    assert_quantity_equal(default_metric_forces._M_y, 0 * kNm)
    assert default_metric_forces.unit_system == "metric"
    assert default_metric_forces.id == 1 # First instance after reset

def test_forces_custom_initialization_metric(custom_metric_forces: Forces) -> None:
    """Test Forces initialization with custom metric values."""
    assert custom_metric_forces.label == "Custom Metric"
    assert_quantity_equal(custom_metric_forces._N_x, 100 * kN)
    assert_quantity_equal(custom_metric_forces._V_z, 50 * kN)
    assert_quantity_equal(custom_metric_forces._M_y, 200 * kNm)
    assert custom_metric_forces.unit_system == "metric"
    assert custom_metric_forces.id == 1 # First instance for this test (due to reset fixture)

def test_forces_custom_initialization_imperial(custom_imperial_forces: Forces) -> None:
    """Test Forces initialization with custom imperial values."""
    assert custom_imperial_forces.label == "Custom Imperial"
    assert_quantity_equal(custom_imperial_forces._N_x, 20 * kip)
    assert_quantity_equal(custom_imperial_forces._V_z, 10 * kip)
    # 150 kNm to ft*kip: 150 kN*m * (0.224809 kip/kN) * (3.28084 ft/m) = 110.63 ft*kip approx
    assert_quantity_equal(custom_imperial_forces._M_y, 150 * kNm) # Stored internally as kNm
    assert custom_imperial_forces.unit_system == "imperial"
    assert custom_imperial_forces.id == 1

def test_forces_id_incrementing() -> None:
    """Test that the _id increments correctly for multiple instances."""
    f1 = Forces()
    f2 = Forces()
    f3 = Forces()
    assert f1.id == 1
    assert f2.id == 2
    assert f3.id == 3
    assert Forces._last_id == 3

def test_forces_id_is_read_only(default_metric_forces: Forces) -> None:
    """Test that the 'id' property is read-only."""
    with pytest.raises(AttributeError, match="property 'id' of 'Forces' object has no setter"):
        default_metric_forces.id = 99 # type: ignore [misc]

# --- Tests for Properties (`N_x`, `V_z`, `M_y`) ---

def test_N_x_property_metric(custom_metric_forces: Forces) -> None:
    """Test N_x property returns correct value in metric."""
    assert_quantity_equal(custom_metric_forces.N_x, 100 * kN)

def test_V_z_property_metric(custom_metric_forces: Forces) -> None:
    """Test V_z property returns correct value in metric."""
    assert_quantity_equal(custom_metric_forces.V_z, 50 * kN)

def test_M_y_property_metric(custom_metric_forces: Forces) -> None:
    """Test M_y property returns correct value in metric."""
    assert_quantity_equal(custom_metric_forces.M_y, 200 * kNm)

def test_N_x_property_imperial(custom_imperial_forces: Forces) -> None:
    """Test N_x property returns correct value in imperial."""
    # 20 kip (already kip)
    assert_quantity_equal(custom_imperial_forces.N_x, 20 * kip)

def test_V_z_property_imperial(custom_imperial_forces: Forces) -> None:
    """Test V_z property returns correct value in imperial."""
    # 10 kip (already kip)
    assert_quantity_equal(custom_imperial_forces.V_z, 10 * kip)

def test_M_y_property_imperial(custom_imperial_forces: Forces) -> None:
    """Test M_y property returns correct value in imperial."""
    # 150 kNm converted to ft*kip (110.63 ft*kip approx)
    assert_quantity_equal(custom_imperial_forces.M_y, (150 * kNm).to("ft*kip"))


# --- Tests for `get_forces()` method ---

def test_get_forces_metric(custom_metric_forces: Forces) -> None:
    """Test get_forces returns correct dictionary in metric."""
    forces_dict = custom_metric_forces.get_forces()
    assert_quantity_equal(forces_dict["N_x"], 100 * kN)
    assert_quantity_equal(forces_dict["V_z"], 50 * kN)
    assert_quantity_equal(forces_dict["M_y"], 200 * kNm)

def test_get_forces_imperial(custom_imperial_forces: Forces) -> None:
    """Test get_forces returns correct dictionary in imperial."""
    forces_dict = custom_imperial_forces.get_forces()
    assert_quantity_equal(forces_dict["N_x"], 20 * kip)
    assert_quantity_equal(forces_dict["V_z"], 10 * kip)
    assert_quantity_equal(forces_dict["M_y"], (150 * kNm).to("ft*kip"))


# --- Tests for `set_forces()` method ---

def test_set_forces_updates_values(default_metric_forces: Forces) -> None:
    """Test set_forces updates the internal force values."""
    default_metric_forces.set_forces(N_x=20 * kN, V_z=5 * kN, M_y=30 * kNm)
    assert_quantity_equal(default_metric_forces._N_x, 20 * kN)
    assert_quantity_equal(default_metric_forces._V_z, 5 * kN)
    assert_quantity_equal(default_metric_forces._M_y, 30 * kNm)

def test_set_forces_with_default_values(default_metric_forces: Forces) -> None:
    """Test set_forces with some default values to ensure they remain unchanged."""
    default_metric_forces.set_forces(N_x=20 * kN) # Only set N_x
    assert_quantity_equal(default_metric_forces._N_x, 20 * kN)
    assert_quantity_equal(default_metric_forces._V_z, 0 * kN) # Should remain default
    assert_quantity_equal(default_metric_forces._M_y, 0 * kNm) # Should remain default


# --- Tests for `compare_to()` method ---

def test_compare_to_V_z_true(custom_metric_forces: Forces) -> None:
    """Test compare_to for V_z when self > other."""
    other_forces = Forces(V_z=25 * kN)
    assert custom_metric_forces.compare_to(other_forces, by="V_z") is True

def test_compare_to_V_z_false(custom_metric_forces: Forces) -> None:
    """Test compare_to for V_z when self <= other."""
    other_forces_equal = Forces(V_z=50 * kN)
    other_forces_greater = Forces(V_z=75 * kN)
    assert custom_metric_forces.compare_to(other_forces_equal, by="V_z") is False
    assert custom_metric_forces.compare_to(other_forces_greater, by="V_z") is False

def test_compare_to_N_x(custom_metric_forces: Forces) -> None:
    """Test compare_to for N_x."""
    other_forces = Forces(N_x=50 * kN)
    assert custom_metric_forces.compare_to(other_forces, by="N_x") is True
    other_forces_larger = Forces(N_x=150 * kN)
    assert custom_metric_forces.compare_to(other_forces_larger, by="N_x") is False

def test_compare_to_M_y(custom_metric_forces: Forces) -> None:
    """Test compare_to for M_y."""
    other_forces = Forces(M_y=150 * kNm)
    assert custom_metric_forces.compare_to(other_forces, by="M_y") is True
    other_forces_larger = Forces(M_y=250 * kNm)
    assert custom_metric_forces.compare_to(other_forces_larger, by="M_y") is False

def test_compare_to_different_units(custom_metric_forces: Forces) -> None:
    """Test compare_to with Forces objects having different internal units."""
    # custom_metric_forces.V_z is 50 kN
    other_forces_imperial = Forces(V_z=10 * kip) # 10 kip = 44.48 kN
    assert custom_metric_forces.compare_to(other_forces_imperial, by="V_z") is True # 50 kN > 44.48 kN

def test_compare_to_invalid_attribute(custom_metric_forces: Forces) -> None:
    """Test compare_to raises ValueError for invalid attribute."""
    with pytest.raises(ValueError, match="Comparison attribute must be one of 'N_x', 'V_z', or 'M_y'"):
        custom_metric_forces.compare_to(Forces(), by="invalid_attr")


# --- Tests for `__str__()` method ---

def test_str_representation_metric(custom_metric_forces: Forces) -> None:
    """Test __str__ representation for metric forces."""
    s = str(custom_metric_forces)
    # The ID will depend on the order of tests due to the fixture setup, so we only check format
    assert f"Force ID: {custom_metric_forces.id}, Label: Custom Metric, N_x: 100.00 kN, V_z: 50.00 kN, M_y: 200.00 kN·m" in s

def test_str_representation_imperial(custom_imperial_forces: Forces) -> None:
    """Test __str__ representation for imperial forces."""
    s = str(custom_imperial_forces)
    assert f"Force ID: {custom_imperial_forces.id}, Label: Custom Imperial, N_x: 20.00 kip, V_z: 10.00 kip, M_y: {custom_imperial_forces.M_y:.2f~P}" in s
    # M_y needs careful checking due to potential floating point and pint's default format.
    # The {custom_imperial_forces.M_y:.1f~P} format specifier ensures it's formatted as pint would do for `str()`
    # and matches the expected output including units with unicode mu for 'micro'.