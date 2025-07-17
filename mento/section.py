from dataclasses import dataclass, field
from typing import TYPE_CHECKING, List, Optional

from mento.material import Concrete, SteelBar
from mento.forces import Forces

from pint import Quantity

if TYPE_CHECKING:
    from mento.node import Node


@dataclass
class Section:
    """
    Represents a section composed of concrete and steel reinforcement.

    Attributes:
        concrete (Concrete): The concrete material used in the section.
        steel_bar (SteelBar): The steel reinforcement used in the section.
        c_c (Quantity): The concrete cover or another relevant quantity.
        node (Optional[Node]): The node associated with this section, if any.
        label (Optional[str]): An optional label for the section.

    Methods:
        id: Read-only property to access the section's unique ID.
        check_shear(forces: List[Forces]): Checks the section for shear under given forces.
        design_shear(forces: List[Forces]): Designs the section for shear under given forces.
        check_flexure(forces: List[Forces]): Checks the section for flexure under given forces.
        design_flexure(forces: List[Forces]): Designs the section for flexure under given forces.
        shear_results_detailed(force: Optional[Forces] = None): Provides detailed shear results.
        shear_results_detailed_doc(force: Optional[Forces] = None): Provides detailed shear results for documentation.
        flexure_results_detailed(force: Optional[Forces] = None): Provides detailed flexure results.
        flexure_results_detailed_doc(force: Optional[Forces] = None): Provides detailed flexure results for documentation.
        results: Property to access beam results for Jupyter Notebook.
    """

    concrete: Concrete = field(kw_only=True)
    steel_bar: SteelBar = field(kw_only=True)
    c_c: Quantity = field(kw_only=True)
    _id: int = field(init=False)
    _last_id: int = field(default=0, init=False, repr=False)
    node: Optional["Node"] = field(default=None, init=False)
    label: Optional[str] = field(default=None, kw_only=True)

    def __post_init__(self) -> None:
        Section._last_id += 1
        self._id = Section._last_id

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
