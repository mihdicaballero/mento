from typing import List
from mento.section import Section
from mento.forces import Forces


class Node:
    def __init__(self, section: Section, forces_list: List[Forces]) -> None:
        if not isinstance(section, Section):
            raise TypeError("section must be an instance of Section")
        if not all(isinstance(forces, Forces) for forces in forces_list):
            raise TypeError("All items in forces_list must be instances of class Forces,")
        self.section = section
        self.forces_list = forces_list
        # Set the node reference in the section (beam)
        self.section.node = self  # Associate the Node with the RectangularBeam

    def add_forces(self, forces: Forces) -> None:
        if not isinstance(forces, Forces):
            raise TypeError("forces must be an instance of class Forces,")
        self.forces_list.append(forces)

    def get_forces(self) -> List[Forces]:
        """Returns the list of forces applied to this node."""
        return self.forces_list
