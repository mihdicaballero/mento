from dataclasses import dataclass
from mento.section import Section
from mento.material import Concrete, SteelBar
from mento.settings import Settings
from typing import Optional


@dataclass
class ConcreteSection(Section):
    def __init__(self, name: str, concrete: Concrete, steel_bar: SteelBar, settings: Optional[Settings] = None):
        super().__init__(name, settings)
        self.concrete = concrete
        self.steel_bar = steel_bar
        self.cc = self._settings.get_setting('clear_cover')
        self.stirrup_d_b = self._settings.get_setting('stirrup_diameter')
        