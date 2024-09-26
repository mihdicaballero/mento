from section import Section
from forces import Forces

class Node:
    def __init__(self, section, forces_list):
        if not isinstance(section, Section):
            raise TypeError("section must be an instance of Section")
        if not all(isinstance(forces, Forces) for forces in forces_list):
            raise TypeError("All items in forces_list must be instances of Forces")
        self.section = section
        self.forces_list = forces_list

    def add_forces(self, forces):
        if not isinstance(forces, Forces):
            raise TypeError("forces must be an instance of Forces")
        self.forces_list.append(forces)

    def get_total_moment(self):
        return sum(forces.get_moment() for forces in self.forces_list)


def main():
    # Ejemplo de uso
    section = Section(material="Concrete", SteelBar_long="Longitudinal Steel", SteelBar_trans="Transverse Steel", cc=25)
    forces1 = Forces(My=100)
    forces2 = Forces(My=200)

    node = Node(section=section, forces_list=[forces1, forces2])

    print(f"Node Section Material: {node.section.material}")
    print(f"Total Moment: {node.get_total_moment()}")

    # Añadir otra instancia de Forces
    forces3 = Forces(My=50)
    node.add_forces(forces3)
    print(f"Updated Total Moment: {node.get_total_moment()}")

if __name__ == "__main__":
    main()