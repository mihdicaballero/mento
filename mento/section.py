import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, ClassVar, List, Optional

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

    #: How the element reaches the ground, in the terms the minimum
    #: reinforcement clauses distinguish. ``"free"`` is anything spanning
    #: between supports; ``"soil"`` is an element bearing directly on the
    #: ground, which ACI 318-19 §9.6.1.1(b) exempts from the flexural minimum
    #: written for a beam and EN 1992-1-1 reinforces to the halved geometric
    #: minimum of a foundation instead.
    #:
    #: A ClassVar and not a field on purpose: it says what kind of element this
    #: is, so it belongs to the class -- ``Footing`` sets it -- rather than
    #: being a constructor argument any beam could be handed.
    support: ClassVar[str] = "free"

    concrete: Concrete = field(kw_only=True)
    steel_bar: SteelBar = field(kw_only=True)
    c_c: Quantity = field(kw_only=True)
    _id: int = field(init=False)
    _last_id: int = field(default=0, init=False, repr=False)
    node: Optional["Node"] = field(default=None, init=False)
    label: Optional[str] = field(default=None, kw_only=True)

    def __post_init__(self) -> None:
        # Concrete cover must be a physical length.
        if not isinstance(self.c_c, Quantity) or not self.c_c.check("[length]"):
            raise TypeError("c_c must be a length Quantity.")

        cover_mm = self.c_c.to("mm").magnitude
        if not math.isfinite(cover_mm) or cover_mm < 0:
            raise ValueError("c_c must be greater than or equal to zero.")

        # Assign the section ID only after validating the cover.
        Section._last_id += 1
        self._id = Section._last_id

    @property
    def id(self) -> int:
        """Read-only property to access the private _id."""
        return self._id

    def check(self, forces: List[Forces]) -> None:
        pass

    def design(self, forces: List[Forces]) -> None:
        pass

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

    def _get_units_row_shear(self) -> None:
        pass

    # Beam results for Jupyter Notebook
    @property
    def results(self) -> None:
        pass
