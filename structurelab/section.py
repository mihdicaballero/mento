from dataclasses import dataclass
import forallpeople
forallpeople.environment('structural', top_level=True)

# Definir algunas unidades adicionales útiles
cm = 1e-2 * m  # type: ignore

@dataclass
class Section:
    name: str

@dataclass
class ConcreteSection(Section):
    concrete: str
    steelBar_long: str
    steelBar_trans: str
    cc: float  # Sin valor predeterminado aquí

    def __init__(self, name: str, concrete: str, steelBar_long: str, steelBar_trans: str, cc: float):  # type: ignore
        super().__init__(name)
        self.concrete = concrete
        self.steelBar_long = steelBar_long
        self.steelBar_trans = steelBar_trans
        self.cc = cc

@dataclass
class RectangularConcreteSection(ConcreteSection):
    width: float
    depth: float

    def __init__(self, name: str, concrete: str, steelBar_long: str, steelBar_trans: str, cc: float, width: float, depth: float):  # type: ignore
        super().__init__(name, concrete, steelBar_long, steelBar_trans, cc)
        self.width = width
        self.depth = depth

    def get_area(self):
        return self.width * self.depth

# Ejemplo de uso
section = RectangularConcreteSection(
    name="Sección Rectangular",
    concrete="Concreto",
    steelBar_long="Barras Longitudinales",
    steelBar_trans="Barras Transversales",
    cc=25 * mm,  # type: ignore
    width=300 * mm,  # type: ignore
    depth=500 * mm,  # type: ignore
)

print(f"Nombre de la sección: {section.name}")
print(f"Área de la sección: {section.get_area()}")
