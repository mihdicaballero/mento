from dataclasses import dataclass
from mento.concrete.section import ConcreteSection
from mento.material import Concrete, SteelBar, create_concrete
from mento. settings import Settings
from typing import Optional
from mento.units import psi, ksi, inch
from devtools import debug

@dataclass
class RectangularConcreteSection(ConcreteSection):
    def __init__(self, name: str, concrete: Concrete, steel_bar: SteelBar, width: float, 
                 height: float, settings: Optional[Settings] = None):
        super().__init__(name, concrete, steel_bar, settings)
        self._width = width
        self._height = height
        # TODO(Mihdi): This implementation is not handling edge cases properly
        self.d = 0.9*self._height # Initial value 
        self._A_x = self._width * self._height
        self._I_y = self._width*self._height**3/12
        self._I_z = self._height*self._width**3/12

    @property
    def width(self) -> float:
        "Beam width."
        return self._width

    @property
    def height(self) -> float:
        "Beam height."
        return self._height
    
    @property
    def A_x(self) -> float:
        "Cross section area."
        return self._A_x
    
    @property
    def I_y(self) -> float:
        "Moment of inertia about the Y axis."
        return self._I_y
    
    @property
    def I_z(self) -> float:
        "Moment of inertia about the Z axis."
        return self._I_z
    
def main() -> None:
    concrete = create_concrete(name="C4",f_c=4000*psi, design_code="ACI 318-19") 
    steelBar = SteelBar(name="ADN 420", f_y=60*ksi) 
    section = RectangularConcreteSection(name="V-10x16",concrete=concrete,steel_bar=steelBar,
                                         width=10*inch, height=16*inch)
    debug(section.width, section.height)
    debug(section.A_x, section.I_y, section.I_z)

if __name__ == "__main__":
    main()
