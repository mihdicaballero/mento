from __future__ import annotations
from typing import Any, TYPE_CHECKING
from dataclasses import dataclass
import math
import numpy as np

from mento.beam import RectangularBeam
from mento.units import m, mm, cm, inch

if TYPE_CHECKING:
    from pint import Quantity


@dataclass
class OneWaySlab(RectangularBeam):
    """
    One-way slab section. Inherits all flexure and shear checks from RectangularBeam,
    but uses bar diameters + spacing instead of number of bars for longitudinal reinforcement.
    """

    def __post_init__(self) -> None:
        super().__post_init__()

        self.mode = "slab"

        self._stirrup_d_b = 0 * mm
        self._stirrup_s_l = 0 * mm
        self._A_v = 0 * mm**2 / m

        # Allow slabs to place as many bars as minimum spacing permits
        self.settings.max_bars_per_layer = max(self.settings.max_bars_per_layer, 200)
        # Force same diameter for all bars within a layer
        self.settings.max_diameter_diff = 0 * self.settings.max_diameter_diff

        self._initialize_longitudinal_rebar_attributes()

    def _initialize_longitudinal_rebar_attributes(self) -> None:
        """Initialize all rebar-related attributes with default values."""

        # Unit system default starter rebar for initial effective height
        # In this case is 0 mm since the slab might not have bottom or top rebar
        if self.concrete.unit_system == "metric":
            self._n1_b, self._d_b1_b = 0, 0 * mm
            self._n1_t, self._d_b1_t = 0, 0 * mm
        else:
            self._n1_b, self._d_b1_b = 0, 0 * inch
            self._n1_t, self._d_b1_t = 0, 0 * inch

        # Bottom rebar defaults (diameter and spacing)
        self._s_b1_b = 0 * mm
        self._d_b2_b, self._s_b2_b = 0 * mm, 0 * mm
        self._d_b3_b, self._s_b3_b = 0 * mm, 0 * mm
        self._d_b4_b, self._s_b4_b = 0 * mm, 0 * mm
        self._n2_b = 0  # No position 2 in slabs
        self._n4_b = 0  # No position 4 in slabs

        # Top rebar defaults (diameter and spacing)
        self._s_b1_t = 0 * mm
        self._d_b2_t, self._s_b2_t = 0 * mm, 0 * mm
        self._d_b3_t, self._s_b3_t = 0 * mm, 0 * mm
        self._d_b4_t, self._s_b4_t = 0 * mm, 0 * mm
        self._n2_t = 0  # No position 2 in slabs
        self._n4_t = 0  # No position 4 in slabs

        # Update dependent attributes
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def set_slab_transverse_rebar(
        self,
        d_b: Quantity = 0 * mm,
        s_long: Quantity = 0 * cm,
        s_trans: Quantity = 0 * cm,
    ) -> None:
        """Sets the transverse rebar in the slab section.

        A zero diameter or a zero spacing in either direction means no stirrups,
        which clears the transverse reinforcement.
        """
        self._stirrup_s_l = s_long
        if d_b == 0 * mm or s_long == 0 * cm or s_trans == 0 * cm:
            # No stirrups: same state as a section without transverse rebar
            self._A_v = 0 * mm**2 / m
            self._stirrup_n = 0
            self._stirrup_d_b = 0 * mm
            self._stirrup_s_w = 0 * mm
        else:
            n_legs_per_unit_width = self.width / s_trans
            A_db = (d_b**2) * math.pi / 4  # Area of one stirrup leg per unit width
            self._A_v = A_db * n_legs_per_unit_width / s_long  # Legs area per unit length
            self._stirrup_d_b = d_b
            self._stirrup_s_w = s_trans
            # Both shear checks read _stirrup_n to decide whether the section
            # carries shear reinforcement at all, and the drawings read it as a
            # number of stirrups. Leaving it at zero was why a slab given
            # transverse rebar was still checked as if it had none. A strip has
            # no closed stirrup to count -- the legs at s_trans need not divide
            # into it a whole number of times -- so this is the equivalent
            # number of rows; _A_v above is what the capacity is computed from.
            self._stirrup_n = max(1, round(n_legs_per_unit_width.to("dimensionless").magnitude / 2))

        # Update effective heights
        self._update_effective_heights()

    def _leg_spacing_across_width(self) -> Quantity:
        """The transverse spacing the strip was detailed with.

        A beam derives this from the geometry of its stirrup cage. A slab does
        not have one: the legs sit on a grid, and how far apart they are across
        the strip is chosen, not implied. Deriving it the beam's way put the
        legs of a designed slab at a spacing it had never been given.
        """
        return self._stirrup_s_w

    def _apply_transverse_design(self, design: Any) -> None:
        """A slab is detailed by a spacing in each direction, not by a stirrup count."""
        self.set_slab_transverse_rebar(
            d_b=design["d_b"],
            s_long=design["s_l"],
            s_trans=design["s_w"],
        )

    def set_slab_longitudinal_rebar_bot(
        self,
        d_b1: Quantity = 0 * mm,
        s_b1: Quantity = 0 * mm,
        d_b3: Quantity = 0 * mm,
        s_b3: Quantity = 0 * mm,
    ) -> None:
        """Update the bottom rebar configuration and recalculate attributes.

        Args:
            d_b1: Diameter of position 1 bottom rebar (default: 0 mm).
            s_b1: Spacing of position 1 bottom rebar (default: 0 mm).
            d_b3: Diameter of position 3 bottom rebar (default: 0 mm).
            s_b3: Spacing of position 3 bottom rebar (default: 0 mm).
        """
        self._d_b1_b = d_b1 if d_b1 != 0 * mm else self._d_b1_b
        self._s_b1_b = s_b1 if s_b1 != 0 * mm else self._s_b1_b
        self._d_b3_b = d_b3 if d_b3 != 0 * mm else self._d_b3_b
        self._s_b3_b = s_b3 if s_b3 != 0 * mm else self._s_b3_b
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def set_slab_longitudinal_rebar_top(
        self,
        d_b1: Quantity = 0 * mm,
        s_b1: Quantity = 0 * mm,
        d_b3: Quantity = 0 * mm,
        s_b3: Quantity = 0 * mm,
    ) -> None:
        """Update the top rebar configuration and recalculate attributes.

        Args:
            d_b1: Diameter of position 1 top rebar (default: 0 mm).
            s_b1: Spacing of position 1 top rebar (default: 0 mm).
            d_b3: Diameter of position 3 top rebar (default: 0 mm).
            s_b3: Spacing of position 3 top rebar (default: 0 mm).

        """
        self._d_b1_t = d_b1 if d_b1 != 0 * mm else self._d_b1_t
        self._s_b1_t = s_b1 if s_b1 != 0 * mm else self._s_b1_t
        self._d_b3_t = d_b3 if d_b3 != 0 * mm else self._d_b3_t
        self._s_b3_t = s_b3 if s_b3 != 0 * mm else self._s_b3_t
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def _calculate_longitudinal_rebars(self) -> None:
        """Calculate the total rebar area for a slab, given spacing and slab width."""

        # Helper function: Calculate number of bars given spacing and slab width
        def calculate_bars(spacing: Quantity, width: Quantity) -> int:
            if spacing == 0 * mm:
                return 0  # Default to 0 bars if spacing is zero
            # Round up to ensure full bars (e.g., 3.2 bars → 4 bars)
            return int(np.ceil((width.to("cm").magnitude / spacing.to("cm").magnitude)))  # Returns dimensionless count

        # --- BOTTOM REBAR ---
        self._n1_b = calculate_bars(self._s_b1_b, self.width)
        self._n3_b = calculate_bars(self._s_b3_b, self.width)

        # --- TOP REBAR ---
        self._n1_t = calculate_bars(self._s_b1_t, self.width)
        self._n3_t = calculate_bars(self._s_b3_t, self.width)
