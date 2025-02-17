from dataclasses import dataclass
from mento.section import Section
from mento.material import SteelBar, Concrete
from pint.facets.plain import PlainQuantity
from typing import Optional, Dict, Any
    
@dataclass
class RectangularSection(Section):
    def __init__(self, concrete: Concrete, steel_bar: SteelBar, width: PlainQuantity, height: PlainQuantity, 
                 settings: Optional[Dict[str, Any]] = None) -> None:
        super().__init__(concrete, steel_bar, settings)

        if settings:
            self.settings.update(settings)  # Update with any provided settings

        self._width = width
        self._height = height
        self._A_x = self._width * self._height
        self._I_y = self._width*self._height**3/12
        self._I_z = self._height*self._width**3/12


    @property
    def width(self) -> PlainQuantity:
        "Beam width."
        return self._width.to('cm')

    @property
    def height(self) -> PlainQuantity:
        "Beam height."
        return self._height.to('cm')
    
    @property
    def A_x(self) -> PlainQuantity:
        "Cross section area."
        return self._A_x.to('cm**2')
    
    @property
    def I_y(self) -> PlainQuantity:
        "Moment of inertia about the Y axis."
        return self._I_y.to('cm**4')
    
    @property
    def I_z(self) -> PlainQuantity:
        "Moment of inertia about the Z axis."
        return self._I_z.to('cm**4')
