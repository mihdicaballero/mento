from typing import Optional
from dataclasses import dataclass, field
from pint import Quantity
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch


from mento.section import Section
from mento.results import CUSTOM_COLORS
from mento.settings import BeamSettings, GLOBAL_BEAM_SETTINGS


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
        self._stirrup_d_b = self.settings.stirrup_diameter_ini

    @property
    def settings(self) -> BeamSettings:
        """Access global design rules."""
        return GLOBAL_BEAM_SETTINGS

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

    def plot(self) -> None:
        """
        Plots the rectangular section with a dark gray border, light gray hatch, and dimensions.
        Also plots the stirrup with rounded corners and thickness.
        """

        # Create figure and axis
        fig, self._ax = plt.subplots()

        # Create a rectangle patch for the section
        rect = Rectangle(
            (0, 0),
            self.width.magnitude,
            self.height.magnitude,
            linewidth=1.3,
            edgecolor=CUSTOM_COLORS["dark_gray"],
            facecolor=CUSTOM_COLORS["light_gray"],
        )
        self._ax.add_patch(rect)

        # Calculate stirrup dimensions
        c_c = self.c_c.to("cm").magnitude
        stirrup_width = self.width.to("cm").magnitude - 2 * c_c
        stirrup_height = self.height.to("cm").magnitude - 2 * c_c
        stirrup_thickness = self._stirrup_d_b.to("cm").magnitude

        # Create rounded corners for the stirrup
        inner_radius = stirrup_thickness * 2
        outer_radius = stirrup_thickness * 3

        # Create the outer rounded rectangle for the stirrup
        outer_rounded_rect = FancyBboxPatch(
            (c_c, c_c),  # Bottom-left corner
            stirrup_width,  # Width
            stirrup_height,  # Height
            boxstyle=f"Round, pad=0, rounding_size={outer_radius}",  # Rounded corners
            edgecolor=CUSTOM_COLORS["dark_blue"],
            facecolor="white",
            linewidth=1,
        )
        self._ax.add_patch(outer_rounded_rect)

        # Create the inner rounded rectangle for the stirrup (offset by thickness)
        inner_rounded_rect = FancyBboxPatch(
            (c_c + stirrup_thickness, c_c + stirrup_thickness),  # Bottom-left corner
            stirrup_width - 2 * stirrup_thickness,  # Width
            stirrup_height - 2 * stirrup_thickness,  # Height
            boxstyle=f"Round, pad=0, rounding_size={inner_radius}",  # Rounded corners
            edgecolor=CUSTOM_COLORS["dark_blue"],
            facecolor=CUSTOM_COLORS["light_gray"],
            linewidth=1,
        )
        self._ax.add_patch(inner_rounded_rect)

        # Set plot limits with some padding
        padding = max(self.width.magnitude, self.height.magnitude) * 0.2
        self._ax.set_xlim(-padding, self.width.magnitude + padding)
        self._ax.set_ylim(-padding, self.height.magnitude + padding)

        # Text and dimension offsets
        dim_offset = 2.5
        text_offset = dim_offset + 2
        # Add width dimension
        self._ax.annotate(
            "",  # No text here, text is added separately
            xy=(0, -dim_offset),  # Start of arrow (left side)
            xytext=(self.width.magnitude, -dim_offset),  # End of arrow (right side)
            arrowprops={
                "arrowstyle": "<->",
                "lw": 1,
                "color": CUSTOM_COLORS["dark_blue"],
            },
        )
        if self.concrete.unit_system == "imperial":
            # Example: format to 2 decimal places, then use pint's compact (~P) format
            width = "{:.0f~P}".format(self.width.to("inch"))
            height = "{:.0f~P}".format(self.height.to("inch"))
        else:
            width = "{:.0f~P}".format(self.width.to("cm"))
            height = "{:.0f~P}".format(self.height.to("cm"))
        # Add width dimension text below the arrow
        self._ax.text(
            self.width.magnitude / 2,  # Center of the arrow
            -text_offset,  # Slightly below the arrow
            width,
            ha="center",
            va="top",
            color=CUSTOM_COLORS["dark_gray"],
        )

        # Add height dimension
        self._ax.annotate(
            "",  # No text here, text is added separately
            xy=(-dim_offset, 0),  # Start of arrow (bottom)
            xytext=(-dim_offset, self.height.magnitude),  # End of arrow (top)
            arrowprops={
                "arrowstyle": "<->",
                "lw": 1,
                "color": CUSTOM_COLORS["dark_blue"],
            },
        )
        # Add height dimension text to the left of the arrow
        self._ax.text(
            -text_offset,  # Slightly to the left of the arrow
            self.height.magnitude / 2,  # Center of the arrow
            height,
            ha="right",
            va="center",
            color=CUSTOM_COLORS["dark_gray"],
            rotation=90,  # Rotate text vertically
        )

        # Set aspect of the plot to be equal
        self._ax.set_aspect("equal")
        # Remove axes for better visualization
        self._ax.axis("off")
