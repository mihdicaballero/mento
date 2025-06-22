import pytest
from mento.node import Node
from mento.forces import Forces
from mento.material import Concrete_ACI_318_19, SteelBar
from mento.beam import RectangularBeam
from mento.units import psi, kN, kNm, inch


@pytest.fixture
def basic_section() -> RectangularBeam:
    concrete = Concrete_ACI_318_19(name="fc 4000", f_c=4000 * psi)
    steel = SteelBar(name="fy 60000", f_y=60000 * psi)
    return RectangularBeam(
        label="Test Beam",
        concrete=concrete,
        steel_bar=steel,
        width=12 * inch,
        height=24 * inch,
    )


@pytest.fixture
def basic_force() -> Forces:
    return Forces(label="Test Force", V_z=40 * kN, M_y=400 * kNm)


def test_node_initialization(
    basic_section: RectangularBeam, basic_force: Forces
) -> None:
    """Test that a Node can be initialized with a section and forces"""
    node = Node(section=basic_section, forces=basic_force)

    assert node.id == 1  # First node should have ID 1
    assert node.section == basic_section
    assert len(node.get_forces_list()) == 1
    assert node.get_forces_list()[0] == basic_force


def test_node_with_multiple_forces(
    basic_section: RectangularBeam, basic_force: Forces
) -> None:
    """Test that a Node can be initialized with multiple forces"""
    force2 = Forces(label="Test Force 2", V_z=50 * kN, M_y=500 * kNm)
    node = Node(section=basic_section, forces=[basic_force, force2])

    assert len(node.get_forces_list()) == 2
    assert node.get_forces_list()[0] == basic_force
    assert node.get_forces_list()[1] == force2


def test_add_forces(basic_section: RectangularBeam, basic_force: Forces) -> None:
    """Test that forces can be added to a Node after initialization"""
    node = Node(section=basic_section, forces=basic_force)
    force2 = Forces(label="Test Force 2", V_z=50 * kN, M_y=500 * kNm)

    node.add_forces(force2)
    assert len(node.get_forces_list()) == 2
    assert node.get_forces_list()[1] == force2


def test_add_multiple_forces(
    basic_section: RectangularBeam, basic_force: Forces
) -> None:
    """Test that multiple forces can be added at once"""
    node = Node(section=basic_section, forces=basic_force)
    force2 = Forces(label="Test Force 2", V_z=50 * kN, M_y=500 * kNm)
    force3 = Forces(label="Test Force 3", V_z=60 * kN, M_y=600 * kNm)

    node.add_forces([force2, force3])
    assert len(node.get_forces_list()) == 3
    assert node.get_forces_list()[1] == force2
    assert node.get_forces_list()[2] == force3


def test_clear_forces(basic_section: RectangularBeam, basic_force: Forces) -> None:
    """Test that forces can be cleared from a Node"""
    node = Node(section=basic_section, forces=basic_force)
    assert len(node.get_forces_list()) == 1

    node.clear_forces()
    assert len(node.get_forces_list()) == 0


def test_node_id_increment(basic_section: RectangularBeam, basic_force: Forces) -> None:
    """Test that Node IDs increment correctly"""
    node1 = Node(section=basic_section, forces=basic_force)
    node2 = Node(section=basic_section, forces=basic_force)
    node3 = Node(section=basic_section, forces=basic_force)

    assert node1.id == 6
    assert node2.id == 7
    assert node3.id == 8


def test_repr_method(basic_section: RectangularBeam, basic_force: Forces) -> None:
    """Test the string representation of a Node"""
    node = Node(section=basic_section, forces=basic_force)
    repr_str = repr(node)

    assert f"Node ID: {node.id}" in repr_str
    assert "Test Beam" in repr_str  # From section label
    assert "Test Force" in repr_str  # From force label


def test_section_node_reference(
    basic_section: RectangularBeam, basic_force: Forces
) -> None:
    """Test that the Node is properly referenced in its Section"""
    node = Node(section=basic_section, forces=basic_force)
    assert basic_section.node == node
