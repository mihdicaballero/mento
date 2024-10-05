from dataclasses import dataclass
from mento.concrete.section import ConcreteSection
from mento.material import SteelBar, Concrete
from mento.settings import Settings
from mento.units import ksi, inch
from devtools import debug
from pint import Quantity
from pint.facets.plain import PlainQuantity
from typing import Optional
    
@dataclass
class RectangularConcreteSection(ConcreteSection):
    def __init__(self, concrete: Concrete, steel_bar: SteelBar, width: PlainQuantity, height: PlainQuantity, 
                 settings: Optional[Settings] = None) -> None:
        super().__init__(concrete, steel_bar, settings)
        self._width = width
        self._height = height
        # TODO: This implementation is not handling edge cases properly
        self.d = 0.9*self._height # Initial value 
        self._A_x = self._width * self._height
        self._I_y = self._width*self._height**3/12
        self._I_z = self._height*self._width**3/12
        if settings:
            self.settings.update(settings.settings)  # Update with any provided settings

    @property
    def width(self) -> PlainQuantity:
        "Beam width."
        return self._width.to('cm')

    @property
    def height(self) -> PlainQuantity:
        "Beam height."
        return self._height.to('cm')
    
    @property
    def A_x(self) -> Quantity:
        "Cross section area."
        return self._A_x.to('cm**2')
    
    @property
    def I_y(self) -> Quantity:
        "Moment of inertia about the Y axis."
        return self._I_y.to('cm**4')
    
    @property
    def I_z(self) -> Quantity:
        "Moment of inertia about the Z axis."
        return self._I_z.to('cm**4')
    
def main() -> None:
    concrete = Concrete('C25') 
    steel_bar = SteelBar(name="ADN 420", f_y=60*ksi) 
    section = RectangularConcreteSection(concrete=concrete, steel_bar=steel_bar, width=10*inch, height=16*inch)
    debug(section.width, section.height)
    debug(section.A_x, section.I_y, section.I_z)
    debug(section.get_settings)

if __name__ == "__main__":
    main()
