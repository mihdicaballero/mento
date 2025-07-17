import pytest
from pint import Quantity
from typing import List

from mento.material import SteelBar, Concrete
from mento.units import MPa, mm, kN, m
from mento.section import Section  # Assuming this is the module name
from mento.forces import Forces # Assuming mento.forces.Forces is the correct path


@pytest.fixture()
def setup_basic_section() -> Section:
    concrete = Concrete("C25")
    steel_bar = SteelBar(name="B500S", f_y=500 * MPa)
    section = Section(label="V101", concrete=concrete, steel_bar=steel_bar, c_c = 25 * mm)
    return section

@pytest.fixture()
def concrete_c25() -> Concrete:
    """Fixture for a C25 concrete."""
    return Concrete("C25")

@pytest.fixture()
def steel_b500s() -> SteelBar:
    """Fixture for B500S steel."""
    return SteelBar(name="B500S", f_y=500 * MPa)

@pytest.fixture()
def basic_section(concrete_c25: Concrete, steel_b500s: SteelBar) -> Section:
    """Fixture for a basic Section instance."""
    return Section(label="V101", concrete=concrete_c25, steel_bar=steel_b500s, c_c=25 * mm)

@pytest.fixture()
def sample_forces() -> List[Forces]:
    """Fixture for a list of sample Forces objects."""
    # Assuming Forces has a constructor like Forces(N, Vx, Vy, Mx, My, Mz) or similar
    # Adjust this based on the actual Forces class definition
    return [
        Forces(N_x=100 * kN, V_z=5 * kN, M_y=20 * kN * m),
        Forces(N_x=50 * kN, V_z=5 * kN, M_y=10 * kN * m)
    ]

# --- Tests for Section Initialization and Attributes ---

def test_section_initialization(basic_section: Section, concrete_c25: Concrete, steel_b500s: SteelBar) -> None:
    """Test that the Section is initialized correctly with basic attributes."""
    assert basic_section.label == "V101"
    assert basic_section.concrete is concrete_c25 # Check if it's the exact same object
    assert basic_section.steel_bar is steel_b500s # Check if it's the exact same object
    assert basic_section.c_c == 25 * mm
    assert isinstance(basic_section.c_c, Quantity)

def test_section_default_node_and_label() -> None:
    """Test that optional attributes have correct default values."""
    concrete = Concrete("C30")
    steel = SteelBar(name="A400", f_y=400 * MPa)
    section = Section(concrete=concrete, steel_bar=steel, c_c=30 * mm)
    assert section.label is None
    assert section.node is None

def test_section_id_increment(basic_section: Section) -> None:
    """Test that the section ID increments correctly for new instances."""
    section1 = basic_section # This is the first instance created by the fixture
    
    # Create new instances to ensure _last_id is correctly managed
    concrete_temp = Concrete("C25")
    steel_temp = SteelBar(name="B500S", f_y=500 * MPa)

    section2 = Section(
        label="V102", concrete=concrete_temp, steel_bar=steel_temp, c_c=25 * mm
    )
    section3 = Section(
        label="V103", concrete=concrete_temp, steel_bar=steel_temp, c_c=25 * mm
    )

    # The fixture's instance will have an ID. Subsequent instances will increment from there.
    # The exact value depends on how many times Section has been instantiated in tests.
    # A more robust test would be to test relative increment.
    assert section2.id == section1.id + 1
    assert section3.id == section2.id + 1

def test_section_id_read_only(basic_section: Section) -> None:
    """Test that the 'id' property is read-only."""
    with pytest.raises(AttributeError):
        basic_section.id = 999 # type: ignore [misc] # Intentional type ignore for testing attribute error

# --- Tests for Placeholder Methods (and future implementations) ---

def test_check_shear_placeholder(basic_section: Section, sample_forces: List[Forces]) -> None:
    """Test the placeholder for check_shear method."""
    # This test simply ensures the method can be called without error for now.
    # When implementation starts, this test will need to be updated.
    try:
        basic_section.check_shear(sample_forces)
    except NotImplementedError:
        pytest.fail("check_shear should currently be a pass or raise NotImplementedError if intended.")
    except Exception as e:
        pytest.fail(f"check_shear raised an unexpected error: {e}")

def test_design_shear_placeholder(basic_section: Section, sample_forces: List[Forces]) -> None:
    """Test the placeholder for design_shear method."""
    try:
        basic_section.design_shear(sample_forces)
    except NotImplementedError:
        pytest.fail("design_shear should currently be a pass or raise NotImplementedError if intended.")
    except Exception as e:
        pytest.fail(f"design_shear raised an unexpected error: {e}")

def test_check_flexure_placeholder(basic_section: Section, sample_forces: List[Forces]) -> None:
    """Test the placeholder for check_flexure method."""
    try:
        basic_section.check_flexure(sample_forces)
    except NotImplementedError:
        pytest.fail("check_flexure should currently be a pass or raise NotImplementedError if intended.")
    except Exception as e:
        pytest.fail(f"check_flexure raised an unexpected error: {e}")

def test_design_flexure_placeholder(basic_section: Section, sample_forces: List[Forces]) -> None:
    """Test the placeholder for design_flexure method."""
    try:
        basic_section.design_flexure(sample_forces)
    except NotImplementedError:
        pytest.fail("design_flexure should currently be a pass or raise NotImplementedError if intended.")
    except Exception as e:
        pytest.fail(f"design_flexure raised an unexpected error: {e}")

# Add similar placeholder tests for shear_results_detailed, flexure_results_detailed, etc.
def test_shear_results_detailed_placeholder(basic_section: Section, sample_forces: List[Forces]) -> None:
    """Test placeholder for shear_results_detailed."""
    try:
        basic_section.shear_results_detailed(sample_forces[0]) # Pass one force
        basic_section.shear_results_detailed() # Test without force
    except NotImplementedError:
        pytest.fail("shear_results_detailed should currently be a pass or raise NotImplementedError.")
    except Exception as e:
        pytest.fail(f"shear_results_detailed raised an unexpected error: {e}")