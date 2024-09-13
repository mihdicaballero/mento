import sys
import os

# Add the project root to PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from dataclasses import dataclass
from structurelab.section import Section


@dataclass
class ConcreteSection(Section):
    def __init__(self, name: str, concrete, steelBar: str):  # type: ignore
        super().__init__(name)
        self.concrete = concrete
        self.steelBar = steelBar
        self.cc = self._settings.get_setting('clear_cover')