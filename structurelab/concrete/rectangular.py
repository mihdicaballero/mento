import sys
import os

# Add the project root to PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from dataclasses import dataclass
from structurelab.concrete.section import ConcreteSection

@dataclass
class RectangularConcreteSection(ConcreteSection):
    def __init__(self, name: str, concrete, steelBar: str, width: float, depth: float, settings=None):
        super().__init__(name, concrete, steelBar, settings)
        self._width = width
        self._depth = depth

    def get_area(self):
        return self._width * self._depth