from mento.settings import Settings
from dataclasses import dataclass




@dataclass
class Section:
    def __init__(self, name:str, settings=None):
        self._name=name
        self._settings = settings if settings is not None else Settings()
    
    def get_name(self):
        return self._name
    
    def get_settings(self):
        """Returns the current settings."""
        return self._settings.settings
    
    def update_settings(self, new_settings: dict):
        """Updates settings with new values."""
        self._settings.update(new_settings)
        

