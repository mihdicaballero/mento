from dataclasses import dataclass
from devtools import debug
from typing import Optional, TYPE_CHECKING, Dict, Any

from mento.units import MPa
from mento.section import Section
from mento.material import Concrete, SteelBar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

@dataclass
class ConcreteSection(Section):
    def __init__(self, concrete: Concrete, steel_bar: SteelBar, settings: Optional[Dict[str, Any]] = None) -> None:
        super().__init__(settings)
        self.concrete = concrete
        self.steel_bar = steel_bar
        # Update default settings, if entered by user
        if settings:
            self.settings.update(settings)  # Update with any provided settings

        self.c_c: PlainQuantity = self.settings.get_setting('clear_cover')
        self._stirrup_d_b: PlainQuantity = self.settings.get_setting('stirrup_diameter_ini')
        self._long_d_b: PlainQuantity = self.settings.get_setting('longitudinal_diameter_ini')
    

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