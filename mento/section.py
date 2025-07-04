from dataclasses import dataclass, field
from typing import TYPE_CHECKING, List, Optional

from mento.material import Concrete, SteelBar
from mento.forces import Forces
from mento.units import mm

from pint import Quantity

if TYPE_CHECKING:
    from mento.node import Node


@dataclass
class Section:
    _id: int
    c_c: Quantity = field(default = 25*mm)
    _last_id: int = field(
        default=0, init=False, repr=False
    )  # Class variable to keep track of last assigned ID
    node: Optional["Node"] = None  # Use forward declaration with Optional["Node"]
    label: Optional[str] = None

    def __init__(
        self,
        label: Optional[str],
        concrete: Concrete,
        steel_bar: SteelBar,
    ) -> None:
        # Initialize the private ID automatically
        Section._last_id += 1  # Increment the class variable for the next ID
        self._id = Section._last_id  # Assign the next available ID
        self.concrete = concrete
        self.steel_bar = steel_bar
        self.label = label
         
        # self._stirrup_d_b: Quantity = self.settings.stirrup_diameter_ini
        # self._layers_spacing: Quantity = self.settings.layers_spacing

    @property
    def id(self) -> int:
        """Read-only property to access the private _id."""
        return self._id

    def check_shear(self, forces: List[Forces]) -> None:
        pass

    def design_shear(self, forces: List[Forces]) -> None:
        pass

    def check_flexure(self, forces: List[Forces]) -> None:
        pass

    def design_flexure(self, forces: List[Forces]) -> None:
        pass

    def shear_results_detailed(self, force: Optional[Forces] = None) -> None:
        pass

    def shear_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        pass

    def flexure_results_detailed(self, force: Optional[Forces] = None) -> None:
        pass

    def flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        pass

    # Beam results for Jupyter Notebook
    @property
    def results(self) -> None:
        return self.results
