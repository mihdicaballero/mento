from typing import Optional
from dataclasses import dataclass, field
from pint import Quantity
import matplotlib.pyplot as plt

from mento.section import Section


@dataclass
class RectangularSection(Section):
    """
    Represents a rectangular cross-section for structural analysis and design.
    This class extends the generic `Section` class to provide properties and methods specific to rectangular sections,
    including geometric properties (area, moments of inertia), access to design settings, and plotting capabilities.

    Attributes:
        width (Quantity): The width of the rectangular section.
        height (Quantity): The height of the rectangular section.

    Properties:
        settings (BeamSettings): Access to global design rules and settings.
        A_x (Quantity): Cross-sectional area, returned in cm².
        I_y (Quantity): Moment of inertia about the Y axis, returned in cm⁴.
        I_z (Quantity): Moment of inertia about the Z axis, returned in cm⁴.
    Methods:
        plot(): Plots the rectangular section with dimensions and stirrup representation, including rounded corners and thickness.
    """

    # Fields unique to RectangularSection
    width: Quantity = field(kw_only=True)
    height: Quantity = field(kw_only=True)
    _ax: Optional[plt.Axes] = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        super().__post_init__()
        self._A_x = self.width * self.height
        self._I_y = self.width * self.height**3 / 12
        self._I_z = self.height * self.width**3 / 12

    @property
    def A_x(self) -> Quantity:
        "Cross section area."
        return self._A_x.to("cm**2")

    @property
    def I_y(self) -> Quantity:
        "Moment of inertia about the Y axis."
        return self._I_y.to("cm**4")

    @property
    def I_z(self) -> Quantity:
        "Moment of inertia about the Z axis."
        return self._I_z.to("cm**4")
