import math
from typing import Optional
from dataclasses import dataclass, field
from pint import Quantity
from matplotlib.axes import Axes

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
    _ax: Optional[Axes] = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        # Width must be a physical length.
        if not isinstance(self.width, Quantity) or not self.width.check("[length]"):
            raise TypeError("width must be a length Quantity.")

        width_mm = self.width.to("mm").magnitude
        if not math.isfinite(width_mm) or width_mm <= 0:
            raise ValueError("width must be greater than zero.")

        # Height must be a physical length.
        if not isinstance(self.height, Quantity) or not self.height.check("[length]"):
            raise TypeError("height must be a length Quantity.")

        height_mm = self.height.to("mm").magnitude
        if not math.isfinite(height_mm) or height_mm <= 0:
            raise ValueError("height must be greater than zero.")

        # Validate the general section properties.
        super().__post_init__()

        # Cover must fit inside the smallest section dimension.
        cover_mm = self.c_c.to("mm").magnitude
        max_cover_mm = min(width_mm, height_mm) / 2

        if cover_mm >= max_cover_mm:
            raise ValueError(
                "c_c must be less than half of the smallest section dimension."
            )

        # Calculate the rectangular section properties.
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
