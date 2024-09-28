from section import Section
from forces import Forces
from devtools import debug
from typing import List

class Node:
    def __init__(self, section: Section, forces_list: List[Forces]) -> None:
        if not isinstance(section, Section):
            raise TypeError("section must be an instance of Section")
        if not all(isinstance(forces, Forces) for forces in forces_list):
            raise TypeError("All items in forces_list must be instances of Forces")
        self.section = section
        self.forces_list = forces_list

    def add_forces(self, forces: Forces) -> None:
        if not isinstance(forces, Forces):
            raise TypeError("forces must be an instance of Forces")
        self.forces_list.append(forces)

    def get_total_moment(self) -> float:
        return sum(forces.get_moment() for forces in self.forces_list)


def main() -> None:
    # Ejemplo de uso
    section = Section(name="V10")
    forces1 = Forces(My=100)
    forces2 = Forces(My=200)

    node = Node(section=section, forces_list=[forces1, forces2])

    debug(f"Total Moment: {node.get_total_moment()}")

    # Añadir otra instancia de Forces
    forces3 = Forces(My=50)
    node.add_forces(forces3)
    debug(f"Updated Total Moment: {node.get_total_moment()}")

if __name__ == "__main__":
    main()
