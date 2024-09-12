from dataclasses import dataclass
import forallpeople
forallpeople.environment('structural', top_level=True)

from structurelab.settings import Settings



@dataclass
class Section:
    def __init__(self, name:str, settings=None):
        self._name=name
        self._settings = settings if settings is not None else Settings()
    
    def get_name(self):
        return self._name