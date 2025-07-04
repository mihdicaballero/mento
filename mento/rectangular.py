from typing import Optional
from dataclasses import dataclass
from pint.facets.plain import PlainQuantity
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch


from mento.section import Section
from mento.material import SteelBar, Concrete
from mento.results import CUSTOM_COLORS


@dataclass
class RectangularSection(Section):
    label: Optional[str] = None

    def __init__(
        self,
        label: Optional[str],
        concrete: Concrete,
        steel_bar: SteelBar,
        width: PlainQuantity,
        height: PlainQuantity,
    ) -> None:
        super().__init__(label, concrete, steel_bar)

        self._width = width
        self._height = height
        self._A_x = self._width * self._height
        self._I_y = self._width * self._height**3 / 12
        self._I_z = self._height * self._width**3 / 12

        self._ax: Optional[plt.Axes] = None

    @property
    def width(self) -> PlainQuantity:
        "Beam width."
        return self._width.to("cm")

    @property
    def height(self) -> PlainQuantity:
        "Beam height."
        return self._height.to("cm")

    @property
    def A_x(self) -> PlainQuantity:
        "Cross section area."
        return self._A_x.to("cm**2")

    @property
    def I_y(self) -> PlainQuantity:
        "Moment of inertia about the Y axis."
        return self._I_y.to("cm**4")

    @property
    def I_z(self) -> PlainQuantity:
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
            width = "{:~P}".format(self.width.to("inch"))
            height = "{:~P}".format(self.height.to("inch"))
        else:
            width = "{:~P}".format(self.width.to("cm"))
            height = "{:~P}".format(self.height.to("cm"))
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
