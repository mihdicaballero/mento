from mento.settings import Settings
from dataclasses import dataclass, field
from typing import Optional, TYPE_CHECKING, Dict, Any
from devtools import debug

from mento.units import MPa
from mento.material import Concrete, SteelBar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

@dataclass
class Section:
    _id: int
    settings: Settings = field(default_factory=Settings)
    _last_id: int = field(default=0, init=False, repr=False)  # Class variable to keep track of last assigned ID

    def __init__(self, concrete: Concrete, steel_bar: SteelBar, settings: Optional[Dict[str, Any]] = None) -> None:
        # Initialize the private ID automatically
        Section._last_id += 1  # Increment the class variable for the next ID
        self._id = Section._last_id  # Assign the next available ID
        self.concrete = concrete
        self.steel_bar = steel_bar
        
        # Initialize with default settings
        self.settings = Settings() 
        if settings:
            self.settings.update(settings)  # Update with any provided settings
        self.c_c: PlainQuantity = self.settings.get_setting('clear_cover')
        self._stirrup_d_b: PlainQuantity = self.settings.get_setting('stirrup_diameter_ini')
        self._long_d_b: PlainQuantity = self.settings.get_setting('longitudinal_diameter_ini')
        

    @property
    def id(self) -> int:
        """Read-only property to access the private _id."""
        return self._id

    @property
    def get_settings(self) -> Dict[str, Any]:
        """Returns the complete current settings as a dictionary."""
        return self.settings.settings  # Access the settings dictionary from the Settings instance
       
    def update_settings(self, new_settings: Dict[str, Any]) -> None:
        """Updates settings with new values."""
        self.settings.update(new_settings)


def main() -> None:
    concrete = Concrete('C25')
    steel_bar = SteelBar(name='B500S', f_y=500*MPa)
    section = Section(concrete= concrete, steel_bar=steel_bar) 
    debug(section.get_settings, section.id)
    custom_settings = {'clear_cover': 20}
    section.update_settings(custom_settings)
    debug(section.settings.default_settings)
    debug(section.get_settings)
    debug(section)
    debug(section.get_settings)
    custom_settings = {'clear_cover': 20}
    section.update_settings(custom_settings)
    debug(section.get_settings)

if __name__ == "__main__":
    main()
        

