from dataclasses import dataclass
from structurelab.concrete.section import ConcreteSection

@dataclass
class RectangularConcreteSection(ConcreteSection):
    def __init__(self, name: str, concrete, steelBar: str, width: float, depth: float):
        super().__init__(name, concrete, steelBar)
        self._width = width
        self._depth = depth

    def get_area(self):
        return self._width * self._depth