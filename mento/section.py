from mento.settings import Settings
from dataclasses import dataclass, field
from typing import Optional, Dict, Any
from devtools import debug

@dataclass
class Section:
    _id: int
    settings: Settings = field(default_factory=Settings)
    _last_id: int = field(default=0, init=False, repr=False)  # Class variable to keep track of last assigned ID

    def __init__(self, settings: Optional[Dict[str, Any]] = None) -> None:
        # Initialize the private ID automatically
        Section._last_id += 1  # Increment the class variable for the next ID
        self._id = Section._last_id  # Assign the next available ID
        
        # Initialize with default settings
        self.settings = Settings() 
        if settings:
            self.settings.update(settings)  # Update with any provided settings

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
    section = Section()
    debug(section.get_settings, section.id)
    custom_settings = {'clear_cover': 20}
    section.update_settings(custom_settings)
    debug(section.settings.default_settings)
    debug(section.get_settings)

if __name__ == "__main__":
    main()
        

