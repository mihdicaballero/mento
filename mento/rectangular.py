from dataclasses import dataclass
from pint.facets.plain import PlainQuantity
from typing import Optional, Dict, Any
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from mento.section import Section
from mento.material import SteelBar, Concrete
from mento.results import CUSTOM_COLORS
    
@dataclass
class RectangularSection(Section):
    def __init__(self, concrete: Concrete, steel_bar: SteelBar, width: PlainQuantity, height: PlainQuantity, 
                 settings: Optional[Dict[str, Any]] = None) -> None:
        super().__init__(concrete, steel_bar, settings)

        if settings:
            self.settings.update(settings)  # Update with any provided settings

        self._width = width
        self._height = height
        self._A_x = self._width * self._height
        self._I_y = self._width*self._height**3/12
        self._I_z = self._height*self._width**3/12


    @property
    def width(self) -> PlainQuantity:
        "Beam width."
        return self._width.to('cm')

    @property
    def height(self) -> PlainQuantity:
        "Beam height."
        return self._height.to('cm')
    
    @property
    def A_x(self) -> PlainQuantity:
        "Cross section area."
        return self._A_x.to('cm**2')
    
    @property
    def I_y(self) -> PlainQuantity:
        "Moment of inertia about the Y axis."
        return self._I_y.to('cm**4')
    
    @property
    def I_z(self) -> PlainQuantity:
        "Moment of inertia about the Z axis."
        return self._I_z.to('cm**4')

    def plot(self) -> None:
        """
        Plots the rectangular section with a dark gray border, light gray hatch, and dimensions.
        """
    
        # Create figure and axis
        fig, ax = plt.subplots()

        # Create a rectangle patch
        rect = Rectangle((0, 0), self.width.magnitude, self.height.magnitude, 
                         linewidth=1, edgecolor=CUSTOM_COLORS['dark_gray'], facecolor=CUSTOM_COLORS['light_gray'])

        # Add the rectangle to the axis
        ax.add_patch(rect)

        # Set plot limits with some padding
        padding = max(self.width.magnitude, self.height.magnitude) * 0.2
        ax.set_xlim(-padding, self.width.magnitude + padding)
        ax.set_ylim(-padding, self.height.magnitude + padding)

        # Text and dimension offsets
        dim_offset = 2.5
        text_offset = dim_offset + 2
        # Add width dimension
        ax.annotate(
            '',  # No text here, text is added separately
            xy=(0, -dim_offset),  # Start of arrow (left side)
            xytext=(self.width.magnitude, -dim_offset),  # End of arrow (right side)
            arrowprops={"arrowstyle": '<->', "lw": 1, "color": CUSTOM_COLORS['dark_blue']}
        )
        if self.concrete.unit_system == "imperial":
            width = f'{self.width.to('inch')}'
            height = f'{self.height.to('inch')}'
        else:
            width = f'{self.width.to('cm')}'
            height = f'{self.height.to('cm')}'
        # Add width dimension text below the arrow
        ax.text(
            self.width.magnitude / 2,  # Center of the arrow
            -text_offset,  # Slightly below the arrow
            width, 
            ha='center', 
            va='top',
            color =  CUSTOM_COLORS['dark_gray'],
        )

        # Add height dimension
        ax.annotate(
            '',  # No text here, text is added separately
            xy=(-dim_offset, 0),  # Start of arrow (bottom)
            xytext=(-dim_offset, self.height.magnitude),  # End of arrow (top)
            arrowprops={"arrowstyle": '<->', "lw": 1, "color": CUSTOM_COLORS['dark_blue']}
        )
        # Add height dimension text to the left of the arrow
        ax.text(
            -text_offset,  # Slightly to the left of the arrow
            self.height.magnitude / 2,  # Center of the arrow
            height, 
            ha='right', 
            va='center',
            color =  CUSTOM_COLORS['dark_gray'],
            rotation=90  # Rotate text vertically
        )

        # Set aspect of the plot to be equal
        ax.set_aspect('equal')
        # Remove axes for better visualization
        ax.axis('off')

        # Show plot
        plt.show()

