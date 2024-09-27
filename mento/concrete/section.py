from dataclasses import dataclass
from mento.section import Section


@dataclass
class ConcreteSection(Section):
    def __init__(self, name: str, concrete, steelBar: str, settings=None):
        super().__init__(name, settings)
        self.concrete = concrete
        self.steel_bar = steelBar
        self.cc = self._settings.get_setting('clear_cover')
        self.stirrup_d_b = self._settings.get_setting('stirrup_diameter')
        