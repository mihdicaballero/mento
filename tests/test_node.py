import pytest
import pandas as pd
from unittest.mock import (
    MagicMock,
    call,
)  # Import MagicMock and call for robust mocking
from typing import List, Generator

from mento.section import Section  # Assuming Section is correctly imported
from mento.forces import Forces  # Assuming Forces is correctly imported
from mento.node import Node
from mento.units import kN, m  # For creating sample Forces


# --- Fixtures ---


@pytest.fixture
def mock_section() -> MagicMock:
    """
    A mock Section object to pass to Node.
    Its methods will be mocked to track calls.
    We also mock `node` attribute, `label`, and `results` property.
    """
    mock_sec = MagicMock(spec=Section)  # spec=Section ensures mock has Section's attributes/methods
    mock_sec.node = None  # Node will set this
    mock_sec.label = "TestSection"  # For __repr__ test
    mock_sec.results = None  # Mock the property to return None (matching current implementation)
    mock_sec.check.return_value = None  # check() returns None
    mock_sec.design.return_value = None  # design() returns None
    mock_sec.check_flexure.return_value = pd.DataFrame({"Result": ["Flexure Check OK"]})
    mock_sec.design_flexure.return_value = pd.DataFrame({"Result": ["Flexure Design Done"]})
    mock_sec.check_shear.return_value = pd.DataFrame({"Result": ["Shear Check OK"]})
    mock_sec.design_shear.return_value = pd.DataFrame({"Result": ["Shear Design Done"]})
    mock_sec.shear_results_detailed.return_value = None  # Assuming no specific return
    mock_sec.shear_results_detailed_doc.return_value = None
    mock_sec.flexure_results_detailed.return_value = None
    mock_sec.flexure_results_detailed_doc.return_value = None
    return mock_sec


@pytest.fixture
def sample_force_single() -> Forces:
    """A single sample Forces object."""
    return Forces(N_x=100 * kN, V_z=10 * kN, M_y=20 * kN * m)


@pytest.fixture
def sample_forces_list() -> List[Forces]:
    """A list of sample Forces objects."""
    return [
        Forces(N_x=50 * kN, V_z=5 * kN, M_y=10 * kN * m),
        Forces(N_x=75 * kN, V_z=8 * kN, M_y=15 * kN * m),
    ]


@pytest.fixture(autouse=True)
def reset_node_id_counter() -> Generator[None, None, None]:
    """Fixture to reset the Node._last_id counter before each test."""
    original_last_id = Node._last_id
    Node._last_id = 0
    yield  # Allow the test to run
    Node._last_id = original_last_id  # Restore original state if needed, though not critical here


# --- Tests for Initialization ---


def test_node_initialization_with_single_force(mock_section: MagicMock, sample_force_single: Forces) -> None:
    """Test Node initialization with a single Forces object."""
    node = Node(mock_section, sample_force_single)
    assert node.section is mock_section
    assert node.forces == [sample_force_single]
    assert node.section.node is node  # Check if section.node is set to this node


def test_node_initialization_with_list_of_forces(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test Node initialization with a list of Forces objects."""
    node = Node(mock_section, sample_forces_list)
    assert node.section is mock_section
    assert node.forces == sample_forces_list
    assert node.section.node is node


def test_node_initialization_invalid_section_type(sample_force_single: Forces) -> None:
    """Test Node initialization with an invalid section type."""
    with pytest.raises(TypeError, match="section must be an instance of Section"):
        Node("not a section", sample_force_single)  # type: ignore


def test_node_initialization_invalid_forces_type(mock_section: MagicMock) -> None:
    """Test Node initialization with an invalid forces type."""
    with pytest.raises(TypeError, match="forces must be an instance of Forces or a list of Forces"):
        Node(mock_section, "not a force")  # type: ignore

    with pytest.raises(TypeError, match="forces must be an instance of Forces or a list of Forces"):
        Node(mock_section, [Forces(N_x=1 * kN), "not a force"])  # type: ignore


# --- Tests for ID Management ---


def test_node_id_increment(mock_section: MagicMock, sample_force_single: Forces) -> None:
    """Test that Node IDs increment correctly."""
    node1 = Node(mock_section, sample_force_single)
    node2 = Node(mock_section, sample_force_single)  # Create another node

    assert node1.id == 1
    assert node2.id == 2
    assert Node._last_id == 2  # Check the class counter


def test_node_id_read_only(mock_section: MagicMock, sample_force_single: Forces) -> None:
    """Test that the 'id' property is read-only."""
    node = Node(mock_section, sample_force_single)
    with pytest.raises(AttributeError):
        node.id = 999  # type: ignore [misc]


# --- Tests for Force Management ---


def test_add_single_force(mock_section: MagicMock, sample_force_single: Forces) -> None:
    """Test adding a single force to the node."""
    node = Node(mock_section, [])  # Initialize with empty forces
    node.add_forces(sample_force_single)
    assert node.forces == [sample_force_single]


def test_add_list_of_forces(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test adding a list of forces to the node."""
    node = Node(mock_section, Forces(N_x=10 * kN))  # Start with one force
    node.add_forces(sample_forces_list)
    assert len(node.forces) == 1 + len(sample_forces_list)
    assert node.forces[1:] == sample_forces_list


def test_add_forces_invalid_type(mock_section: MagicMock) -> None:
    """Test adding forces with invalid type."""
    node = Node(mock_section, [])
    with pytest.raises(TypeError, match="forces must be an instance of Forces or a list of Forces"):
        node.add_forces("invalid_force")  # type: ignore

    with pytest.raises(TypeError, match="forces must be an instance of Forces or a list of Forces"):
        node.add_forces([Forces(N_x=1 * kN), "invalid"])  # type: ignore


def test_get_forces_list(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test getting the list of forces."""
    node = Node(mock_section, sample_forces_list)
    retrieved_forces = node.get_forces_list()
    assert retrieved_forces == sample_forces_list
    assert retrieved_forces is node.forces  # Should return the actual list, not a copy


def test_clear_forces(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test clearing all forces from the node."""
    node = Node(mock_section, sample_forces_list)
    assert len(node.forces) > 0
    node.clear_forces()
    assert node.forces == []


# --- Tests for check() and design() methods ---


def test_check_delegation(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test that check() calls section.check() with correct forces."""
    node = Node(mock_section, sample_forces_list)
    result = node.check()
    mock_section.check.assert_called_once_with(sample_forces_list)
    assert result is None  # section.check() currently returns None


def test_design_delegation(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test that design() calls section.design() with correct forces."""
    node = Node(mock_section, sample_forces_list)
    result = node.design()
    mock_section.design.assert_called_once_with(sample_forces_list)
    assert result is None  # section.design() currently returns None


def test_check_after_adding_forces(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test check() after dynamically adding forces to the node."""
    node = Node(mock_section, sample_forces_list[0])
    node.add_forces(sample_forces_list[1])

    result = node.check()
    assert result is None
    # Verify check was called with the updated forces list
    mock_section.check.assert_called_once()


def test_design_after_clearing_forces(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test design() after clearing and re-adding forces."""
    node = Node(mock_section, sample_forces_list)

    node.clear_forces()
    assert len(node.forces) == 0

    result = node.design()
    assert result is None
    mock_section.design.assert_called_once_with([])


def test_check_and_design_sequence(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test calling check() followed by design() in sequence."""
    node = Node(mock_section, sample_forces_list)

    check_result = node.check()
    assert check_result is None

    design_result = node.design()
    assert design_result is None

    # Verify both were called with the same forces
    mock_section.check.assert_called_once_with(sample_forces_list)
    mock_section.design.assert_called_once_with(sample_forces_list)


# --- Tests for Delegation to Section Methods ---


def test_check_flexure_delegation(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test that check_flexure calls section.check_flexure with correct forces."""
    node = Node(mock_section, sample_forces_list)
    result_df = node.check_flexure()
    mock_section.check_flexure.assert_called_once_with(sample_forces_list)
    assert isinstance(result_df, pd.DataFrame)
    assert result_df.equals(pd.DataFrame({"Result": ["Flexure Check OK"]}))


def test_design_flexure_delegation(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test that design_flexure calls section.design_flexure with correct forces."""
    node = Node(mock_section, sample_forces_list)
    result_df = node.design_flexure()
    mock_section.design_flexure.assert_called_once_with(sample_forces_list)
    assert isinstance(result_df, pd.DataFrame)
    assert result_df.equals(pd.DataFrame({"Result": ["Flexure Design Done"]}))


def test_check_shear_delegation(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test that check_shear calls section.check_shear with correct forces."""
    node = Node(mock_section, sample_forces_list)
    result_df = node.check_shear()
    mock_section.check_shear.assert_called_once_with(sample_forces_list)
    assert isinstance(result_df, pd.DataFrame)
    assert result_df.equals(pd.DataFrame({"Result": ["Shear Check OK"]}))


def test_design_shear_delegation(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test that design_shear calls section.design_shear with correct forces."""
    node = Node(mock_section, sample_forces_list)
    result_df = node.design_shear()
    mock_section.design_shear.assert_called_once_with(sample_forces_list)
    assert isinstance(result_df, pd.DataFrame)
    assert result_df.equals(pd.DataFrame({"Result": ["Shear Design Done"]}))


def test_shear_results_detailed_delegation(mock_section: MagicMock, sample_force_single: Forces) -> None:
    """Test shear_results_detailed calls section.shear_results_detailed with correct force."""
    node = Node(mock_section, [sample_force_single])
    node.shear_results_detailed(sample_force_single)
    mock_section.shear_results_detailed.assert_called_once_with(sample_force_single)
    # Also test with None
    node.shear_results_detailed(None)
    mock_section.shear_results_detailed.assert_has_calls([call(sample_force_single), call(None)])


def test_shear_results_detailed_doc_delegation(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test shear_results_detailed_doc calls section.shear_results_detailed_doc."""
    node = Node(mock_section, sample_forces_list)
    node.shear_results_detailed_doc()
    mock_section.shear_results_detailed_doc.assert_called_once_with(None)  # Default to None


def test_flexure_results_detailed_delegation(mock_section: MagicMock, sample_force_single: Forces) -> None:
    """Test flexure_results_detailed calls section.flexure_results_detailed."""
    node = Node(mock_section, [sample_force_single])
    node.flexure_results_detailed(sample_force_single)
    mock_section.flexure_results_detailed.assert_called_once_with(sample_force_single)


def test_flexure_results_detailed_doc_delegation(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test flexure_results_detailed_doc calls section.flexure_results_detailed_doc."""
    node = Node(mock_section, sample_forces_list)
    node.flexure_results_detailed_doc()
    mock_section.flexure_results_detailed_doc.assert_called_once_with(None)


# --- Tests for results property ---


def test_results_property(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test that results property delegates to section.results."""
    node = Node(mock_section, sample_forces_list)
    result = node.results
    # Since section.results currently returns None, node.results should also return None
    assert result is None


def test_results_property_read_only(mock_section: MagicMock, sample_forces_list: List[Forces]) -> None:
    """Test that results property is read-only (cannot be set)."""
    node = Node(mock_section, sample_forces_list)
    with pytest.raises(AttributeError):
        node.results = pd.DataFrame()  # type: ignore [misc]


# --- Tests for __repr__ ---


def test_node_repr_with_forces(mock_section: MagicMock, sample_force_single: Forces) -> None:
    """Test the __repr__ method when forces are present."""
    node = Node(mock_section, sample_force_single)
    expected_repr_part1 = f"Node ID: {node.id} - Section label: {mock_section.label}\nForces Applied:\n"
    # Assuming sample_force_single.__str__ or __repr__ returns a sensible string
    expected_repr_part2 = f" - {str(sample_force_single)}"
    assert expected_repr_part1 in repr(node)
    assert expected_repr_part2 in repr(node)


def test_node_repr_without_forces(mock_section: MagicMock) -> None:
    """Test the __repr__ method when no forces are present."""
    node = Node(mock_section, [])
    expected_repr = f"Node ID: {node.id} - Section label: {mock_section.label}\nForces Applied:No forces applied"
    assert repr(node) == expected_repr
