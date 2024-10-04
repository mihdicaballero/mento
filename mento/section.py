from mento.settings import Settings
from dataclasses import dataclass, field
from typing import Optional, Dict, Any
from devtools import debug

@dataclass
class Section: 
    
    _last_id: int = field(default=0, init=False, repr=False)  # Class variable to keep track of last assigned ID

    def __init__(self, label: str, settings: Optional[Settings] = None):
        # Initialize the private ID automatically
        Section._last_id += 1  # Increment the class variable for the next ID
        self._id = Section._last_id  # Assign the next available ID

        self._label: str = label
        self._settings: Settings = settings if settings is not None else Settings()

    @property
    def id(self) -> int:
        """Read-only property to access the private _id."""
        return self._id
    
    @property
    def label(self) -> str:
        """Returns the label of the section."""
        return self._label
    
    @property
    def settings(self) -> Dict[str, Any]:
        """Returns the current settings as a dictionary."""
        return self._settings.settings

    def update_settings(self, new_settings: Dict[str, Any]) -> None:
        """Updates settings with new values."""
        self._settings.update(new_settings)

def main() -> None:
    section = Section(label="V-10x16")
    debug(section.label, section.settings, section.id)
    custom_settings = {'clear_cover': 20}
    section.update_settings(custom_settings)
    debug(section.settings)

if __name__ == "__main__":
    main()
        

