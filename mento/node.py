from typing import Union, Optional, List
import pandas as pd

from mento.section import Section
from mento.forces import Forces


class Node:
    """
    Node Class
    The `Node` class represents a structural node that is associated with a section
    (e.g., a beam) and can have forces applied to it. It provides methods for managing
    forces, performing structural checks, and designing for flexure and shear.
    Attributes:
        section (Section): The section (e.g., beam) associated with the node.
        forces (List[Forces]): A list of forces applied to the node.
    Methods:
        add_forces(forces: Union[Forces, List[Forces]]) -> None:
            Adds one or more forces to the node.
        get_forces_list() -> List[Forces]:
            Returns the list of forces applied to this node.
        clear_forces() -> None:
            Removes all forces from the node.
        check_flexure() -> pd.DataFrame:
            Performs a flexural check on the section using the applied forces.
        design_flexure() -> pd.DataFrame:
            Designs the section for flexure based on the applied forces.
        check_shear() -> pd.DataFrame:
            Performs a shear check on the section using the applied forces.
        design_shear() -> pd.DataFrame:
            Designs the section for shear based on the applied forces.
        shear_results_detailed(force: Optional[Forces] = None) -> None:
            Provides detailed shear results for a specific force or all forces.
        shear_results_detailed_doc(force: Optional[Forces] = None) -> None:
            Provides detailed shear results documentation for a specific force or all forces.
        flexure_results_detailed(force: Optional[Forces] = None) -> None:
            Provides detailed flexure results for a specific force or all forces.
        flexure_results_detailed_doc(force: Optional[Forces] = None) -> None:
            Provides detailed flexure results documentation for a specific force or all forces.
        results -> None:
            A property that provides beam results for Jupyter Notebook.
        id -> int:
            A read-only property to access the private `_id` of the node.
    """

    _id: int
    _last_id: int = 0  # Class variable to keep track of last assigned ID

    def __init__(self, section: Section, forces: Union[Forces, List[Forces]]) -> None:
        # Initialize the private ID automatically
        Node._last_id += 1  # Increment the class variable for the next ID
        self._id = Node._last_id  # Assign the next available ID

        if not isinstance(section, Section):
            raise TypeError("section must be an instance of Section")

        # Ensure forces is either a single Forces object or a list of Forces
        if isinstance(forces, Forces):
            self.forces = [forces]  # Wrap single Forces object in a list
        elif isinstance(forces, list) and all(
            isinstance(force, Forces) for force in forces
        ):
            self.forces = forces  # Assign the list directly if valid
        else:
            raise TypeError("forces must be an instance of Forces or a list of Forces")

        self.section = section
        # Set the node reference in the section (beam)
        self.section.node = self  # Associate the Node with the RectangularBeam

    def add_forces(self, forces: Union[Forces, List[Forces]]) -> None:
        if isinstance(forces, Forces):
            self.forces.append(forces)  # Append single Forces object
        elif isinstance(forces, list) and all(
            isinstance(force, Forces) for force in forces
        ):
            self.forces.extend(forces)  # Extend list with multiple Forces objects
        else:
            raise TypeError("forces must be an instance of Forces or a list of Forces")

    def get_forces_list(self) -> List[Forces]:
        """Returns the list of forces applied to this node."""
        return self.forces

    def clear_forces(self) -> None:
        """Remove all forces from the node."""
        self.forces.clear()

    def check_flexure(self) -> pd.DataFrame:
        return self.section.check_flexure(self.forces)

    def design_flexure(self) -> pd.DataFrame:
        return self.section.design_flexure(self.forces)

    def check_shear(self) -> pd.DataFrame:
        return self.section.check_shear(self.forces)

    def design_shear(self) -> pd.DataFrame:
        return self.section.design_shear(self.forces)

    def shear_results_detailed(self, force: Optional[Forces] = None) -> None:
        return self.section.shear_results_detailed(force)

    def shear_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        return self.section.shear_results_detailed_doc(force)

    def flexure_results_detailed(self, force: Optional[Forces] = None) -> None:
        return self.section.flexure_results_detailed(force)

    def flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        return self.section.flexure_results_detailed_doc(force)

    def __repr__(self) -> str:
        # Get the section label (assuming Section has a `_id` attribute)
        section_label = f"Section label: {self.section.label}"

        # Get the list of forces applied
        forces_list = [
            str(force) for force in self.forces
        ]  # Assuming Forces has a `__repr__` or `__str__` method
        forces_str = (
            "\n  - " + "\n  - ".join(forces_list)
            if forces_list
            else "No forces applied"
        )

        # Combine section label and forces into a single string
        return f"Node ID: {self.id} - {section_label}\nForces Applied:{forces_str}"

    # Beam results for Jupyter Notebook
    @property
    def results(self) -> None:
        return self.section.results

    @property
    def id(self) -> int:
        """Read-only property to access the private _id."""
        return self._id
