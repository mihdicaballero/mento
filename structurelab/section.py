from dataclasses import dataclass
import forallpeople
forallpeople.environment('structural', top_level=True)
import material
from settings import Settings



# Definir algunas unidades adicionales útiles
cm = 1e-2 * m  # type: ignore

@dataclass
class Section:
    def __init__(self, name:str):
        self.name=name
        self.settings=Settings()


@dataclass
class ConcreteSection(Section):
    def __init__(self, name: str, concrete, steelBar: str):  # type: ignore
        super().__init__(name)
        self.concrete = concrete
        self.steelBar = steelBar
        self.cc = self.settings.cc



@dataclass
class RectangularConcreteSection(ConcreteSection):
    def __init__(self, name: str, concrete, steelBar: str, width: float, depth: float):
        super().__init__(name, concrete, steelBar)
        self.width = width
        self.depth = depth

    def get_area(self):
        return self.width * self.depth



@dataclass
class Beam(RectangularConcreteSection):
    def __init__(self, name: str, concrete: material.Concrete, steelBar: str, width: float, depth: float):  # type: ignore
        super().__init__(name, concrete, steelBar, width, depth)

    def design_flexure_ACI_318_19():
        # Determination of Balanced ratio
        print("HOLA")
        pass

    '''
    def desgin_flexure():
        if isinstance(self.concrete, concreteACI):
            self.design_flexure_ACI()
        else:
            raise ValueError("Tipo de hormigón no soportado")
    '''
        



def main():
    # Ejemplo de uso
    concrete=material.create_concrete(name="H30",f_c=30*MPa, design_code="ACI 318-19") # type: ignore
    section = Beam(
        name="Sección Rectangular",
        concrete=concrete,
        steelBar="Barras Longitudinales",
        width=400 * mm,  # type: ignore
        depth=500 * mm,  # type: ignore
    )

    print(f"Nombre de la sección: {section.name}")
    print(f"Área de la sección: {section.get_area()}")


if __name__ == "__main__":
    main()