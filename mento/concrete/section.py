import sys
import os

# Add the project root to PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from dataclasses import dataclass
from mento.section import Section


@dataclass
class ConcreteSection(Section):
    def __init__(self, name: str, concrete, steelBar: str, settings=None):  # type: ignore
        super().__init__(name, settings)
        self.concrete = concrete
        self.steelBar = steelBar
        