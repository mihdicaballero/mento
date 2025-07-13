from dataclasses import dataclass, field
from typing import TYPE_CHECKING, List, Optional

from mento.material import Concrete, SteelBar
from mento.forces import Forces

from pint import Quantity

if TYPE_CHECKING:
    from mento.node import Node


@dataclass
class Section:
    concrete: Concrete
    steel_bar: SteelBar
    c_c: Quantity
    _id: int = field(init=False)
    _last_id: int = field(
        default=0, init=False, repr=False
    )
    node: Optional["Node"] = field(default=None, init=False)
    label: str = field(default=None)

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
