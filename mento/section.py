from mento.settings import Settings
from dataclasses import dataclass
from typing import Optional, Dict, Any

@dataclass
class Section: 
    def __init__(self, name: str, settings: Optional[Settings] = None):
        self._name: str = name
        self._settings: Settings = settings if settings is not None else Settings()

    def get_name(self) -> str:
        """Returns the name of the section."""
        return self._name

    def get_settings(self) -> Dict[str, Any]:
        """Returns the current settings as a dictionary."""
        return self._settings.settings

    def update_settings(self, new_settings: Dict[str, Any]) -> None:
        """Updates settings with new values."""
        self._settings.update(new_settings)
        

