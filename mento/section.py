from dataclasses import dataclass, field
from typing import Optional, TYPE_CHECKING, Dict, Any, List

from mento.settings import Settings
from mento.material import Concrete, SteelBar
from mento.forces import Forces

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity
    from mento.node import Node
    
@dataclass
class Section:
    _id: int
    settings: Settings = field(default_factory=Settings)
    _last_id: int = field(default=0, init=False, repr=False)  # Class variable to keep track of last assigned ID
    node: Optional["Node"] = None  # Use forward declaration with Optional["Node"]

    def __init__(self, concrete: Concrete, steel_bar: SteelBar, settings: Optional[Dict[str, Any]] = None) -> None:
        # Initialize the private ID automatically
        Section._last_id += 1  # Increment the class variable for the next ID
        self._id = Section._last_id  # Assign the next available ID
        self.concrete = concrete
        self.steel_bar = steel_bar
        
        # Initialize with default settings
        self.settings = Settings(concrete) 
        if settings:
            self.settings.update(settings)  # Update with any provided settings
        self.c_c: PlainQuantity = self.settings.get_setting('clear_cover')
        self._stirrup_d_b: PlainQuantity = self.settings.get_setting('stirrup_diameter_ini')
        

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

    def _check_shear(self, forces: List[Forces]) -> None:
        pass

    def _design_shear(self, forces: List[Forces]) -> None:
        pass

    def _check_flexure(self, forces: List[Forces]) -> None:
        pass

    def _design_flexure(self, forces: List[Forces]) -> None:
        pass

    def _shear_results_detailed(self, force: Optional[Forces] = None) -> None:
        pass

    def _shear_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        pass

    def _flexure_results_detailed(self, force: Optional[Forces] = None) -> None:
        pass

    def _flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        pass
    
    # Beam results for Jupyter Notebook
    @property
    def results(self) -> None:
        return self.results

