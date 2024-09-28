from dataclasses import dataclass
from mento.concrete.section import ConcreteSection
from mento.material import Concrete, SteelBar
from mento. settings import Settings
from typing import Optional

@dataclass
class RectangularConcreteSection(ConcreteSection):
    def __init__(self, name: str, concrete: Concrete, steel_bar: SteelBar, width: float, 
                 depth: float, settings: Optional[Settings] = None):
        super().__init__(name, concrete, steel_bar, settings)
        self._width = width
        self._depth = depth
        self.d = 0.9*self._depth # Initial value

    def get_area(self) -> float:
        return self._width * self._depth