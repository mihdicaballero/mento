from dataclasses import dataclass
from mento.section import Section
from mento.material import Concrete, SteelBar
from mento.settings import Settings
from devtools import debug
from mento.units import MPa
from typing import Optional


@dataclass
class ConcreteSection(Section):
    def __init__(self, concrete: Concrete, steel_bar: SteelBar, settings: Optional[Settings] = None) -> None:
        super().__init__(settings)
        self.concrete = concrete
        self.steel_bar = steel_bar
        self.cc = self.settings.get_setting('clear_cover')
        self.stirrup_d_b = self.settings.get_setting('stirrup_diameter')
        if settings:
            self.settings.update(settings.settings)  # Update with any provided settings
    

def main() -> None:
    concrete = Concrete('C25')
    steel_bar = SteelBar(name='B500S', f_y=500*MPa)
    concrete_section = ConcreteSection(concrete= concrete, steel_bar=steel_bar) 
    debug(concrete_section)
    debug(concrete_section.get_settings)
    custom_settings = {'clear_cover': 20}
    concrete_section.update_settings(custom_settings)
    debug(concrete_section.get_settings)

if __name__ == "__main__":
    main()