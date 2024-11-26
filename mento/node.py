from typing import List
from mento.section import Section
from mento.forces import Forces
from typing import Union


class Node:
    def __init__(self, section: Section, forces: Union[Forces, List[Forces]]) -> None:
        if not isinstance(section, Section):
            raise TypeError("section must be an instance of Section")
        
        # Ensure forces is either a single Forces object or a list of Forces
        if isinstance(forces, Forces):
            self.forces = [forces]  # Wrap single Forces object in a list
        elif isinstance(forces, list) and all(isinstance(force, Forces) for force in forces):
            self.forces = forces  # Assign the list directly if valid
        else:
            raise TypeError("forces must be an instance of Forces or a list of Forces")
        
        self.section = section
        # Set the node reference in the section (beam)
        self.section.node = self  # Associate the Node with the RectangularBeam

    def add_forces(self, forces: Union[Forces, List[Forces]]) -> None:
        if isinstance(forces, Forces):
            self.forces.append(forces)  # Append single Forces object
        elif isinstance(forces, list) and all(isinstance(force, Forces) for force in forces):
            self.forces.extend(forces)  # Extend list with multiple Forces objects
        else:
            raise TypeError("forces must be an instance of Forces or a list of Forces")

    def get_forces_list(self) -> List[Forces]:
        """Returns the list of forces applied to this node."""
        return self.forces

    def reset_forces(self) -> None:
        """Reset all forces in the node to zero."""   
        if isinstance(self.forces, Forces):
            # Replace the single Forces object with a zero-initialized Forces object
            self.forces = [Forces()]
        elif isinstance(self.forces, list):
            # Replace each force in the list with a zero-initialized Forces object
            self.forces = [Forces() for _ in self.forces]

