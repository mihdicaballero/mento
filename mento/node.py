from section import Section
from forces import Forces
from units import kN, m
from devtools import debug
from typing import List

class Node:
    def __init__(self, section: Section, forces_list: List[Forces]) -> None:
        if not isinstance(section, Section):
            raise TypeError("section must be an instance of Section")
        if not all(isinstance(forces, Forces) for forces in forces_list):
            raise TypeError("All items in forces_list must be instances of class Forces,")
        self.section = section
        self.forces_list = forces_list

    def add_forces(self, forces: Forces) -> None:
        if not isinstance(forces, Forces):
            raise TypeError("forces must be an instance of class Forces,")
        self.forces_list.append(forces)


def main() -> None:
    # Ejemplo de uso
    section = Section(name="V10")
    forces1 = Forces(My=100*kN*m)
    forces2 = Forces(My=200*kN*m)

    node = Node(section=section, forces_list=[forces1, forces2])

    debug(node.forces_list)

    # Añadir otra instancia de Forces
    forces3 = Forces(My=50*kN*m)
    node.add_forces(forces3)
    debug(node.forces_list)

if __name__ == "__main__":
    main()
