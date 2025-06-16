from __future__ import annotations
from typing import Optional, Dict, Any, TYPE_CHECKING
from dataclasses import dataclass
import math
import numpy as np
from devtools import debug

from mento.beam import RectangularBeam
from mento.forces import Forces
from mento.node import Node
from mento.material import (
    Concrete,
    SteelBar,
    Concrete_ACI_318_19,
    # Concrete_EN_1992_2004,
)
from mento.units import m, mm, cm, inch, MPa, kNm

if TYPE_CHECKING:
    from pint import Quantity


@dataclass
class OneWaySlab(RectangularBeam):
    def __init__(
        self,
        label: Optional[str],
        concrete: Concrete,
        steel_bar: SteelBar,
        thickness: Quantity = 0 * m,
        width: Quantity = 0 * m,
        settings: Optional[Dict[str, Any]] = None,
    ):
        # Initialize as a rectangular beam
        super().__init__(
            label=label,
            concrete=concrete,
            steel_bar=steel_bar,
            width=width,
            height=thickness,
            settings=settings,
        )
        self._stirrup_d_b = 0 * mm  # Converts stirrup to 0 mm for slabs

    def _initialize_longitudinal_rebar_attributes(self) -> None:
        """Initialize all rebar-related attributes with default values."""
        # Bottom rebar defaults (diameter and spacing)
        self._d_b1_b, self._s_b1_b = 0 * mm, 0 * mm
        self._d_b2_b, self._s_b2_b = 0 * mm, 0 * mm
        self._d_b3_b, self._s_b3_b = 0 * mm, 0 * mm
        self._d_b4_b, self._s_b4_b = 0 * mm, 0 * mm

        # Top rebar defaults (diameter and spacing)
        self._d_b1_t, self._s_b1_t = 0 * mm, 0 * mm
        self._d_b2_t, self._s_b2_t = 0 * mm, 0 * mm
        self._d_b3_t, self._s_b3_t = 0 * mm, 0 * mm
        self._d_b4_t, self._s_b4_t = 0 * mm, 0 * mm

        # Default starter rebar (for effective depth calculation only)
        if self.concrete.unit_system == "metric":
            self._d_b1_b, self._d_b1_t = (
                8 * mm,
                8 * mm,
            )  # Diameter only, spacing left as 0
        else:
            self._d_b1_b, self._d_b1_t = (3 * inch) / 8, (3 * inch) / 8

        # Update dependent attributes
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def set_slab_transverse_rebar(
        self,
        d_b: Quantity = 0 * mm,
        s_long: Quantity = 0 * cm,
        s_trans: Quantity = 0 * cm,
    ) -> None:
        """Sets the transverse rebar in the slab section."""
        self._stirrup_s_l = s_long
        n_legs_per_unit_width = self.width / s_trans
        A_db = (d_b**2) * math.pi / 4  # Area of one stirrup leg per unit width
        self._A_v = A_db * n_legs_per_unit_width / s_long  # Legs area per unit length

        # Update effective heights
        self._update_effective_heights()

    def set_slab_longitudinal_rebar_bot(
        self,
        d_b1: Quantity = 0 * mm,
        s_b1: Quantity = 0 * mm,
        d_b2: Quantity = 0 * mm,
        s_b2: Quantity = 0 * mm,
        d_b3: Quantity = 0 * mm,
        s_b3: Quantity = 0 * mm,
        d_b4: Quantity = 0 * mm,
        s_b4: Quantity = 0 * mm,
    ) -> None:
        """Update the bottom rebar configuration and recalculate attributes.

        Args:
            d_b1: Diameter of position 1 bottom rebar (default: 0 mm).
            s_b1: Spacing of position 1 bottom rebar (default: 0 mm).
            ... (similar for positions 2-4).
        """
        self._d_b1_b = d_b1 if d_b1 != 0 * mm else self._d_b1_b
        self._s_b1_b = s_b1 if s_b1 != 0 * mm else self._s_b1_b
        self._d_b2_b = d_b2 if d_b2 != 0 * mm else self._d_b2_b
        self._s_b2_b = s_b2 if s_b2 != 0 * mm else self._s_b2_b
        self._d_b3_b = d_b3 if d_b3 != 0 * mm else self._d_b3_b
        self._s_b3_b = s_b3 if s_b3 != 0 * mm else self._s_b3_b
        self._d_b4_b = d_b4 if d_b4 != 0 * mm else self._d_b4_b
        self._s_b4_b = s_b4 if s_b4 != 0 * mm else self._s_b4_b
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def set_slab_longitudinal_rebar_top(
        self,
        d_b1: Quantity = 0 * mm,
        s_b1: Quantity = 0 * mm,
        d_b2: Quantity = 0 * mm,
        s_b2: Quantity = 0 * mm,
        d_b3: Quantity = 0 * mm,
        s_b3: Quantity = 0 * mm,
        d_b4: Quantity = 0 * mm,
        s_b4: Quantity = 0 * mm,
    ) -> None:
        """Update the top rebar configuration and recalculate attributes."""
        self._d_b1_t = d_b1 if d_b1 != 0 * mm else self._d_b1_t
        self._s_b1_t = s_b1 if s_b1 != 0 * mm else self._s_b1_t
        self._d_b2_t = d_b2 if d_b2 != 0 * mm else self._d_b2_t
        self._s_b2_t = s_b2 if s_b2 != 0 * mm else self._s_b2_t
        self._d_b3_t = d_b3 if d_b3 != 0 * mm else self._d_b3_t
        self._s_b3_t = s_b3 if s_b3 != 0 * mm else self._s_b3_t
        self._d_b4_t = d_b4 if d_b4 != 0 * mm else self._d_b4_t
        self._s_b4_t = s_b4 if s_b4 != 0 * mm else self._s_b4_t
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def _calculate_longitudinal_rebars(self) -> None:
        """Calculate the total rebar area for a slab, given spacing and slab width."""

        # Helper function: Calculate number of bars given spacing and slab width
        def calculate_bars(spacing: Quantity, width: Quantity) -> int:
            if spacing == 0 * mm:
                return 0
            # Round up to ensure full bars (e.g., 3.2 bars → 4 bars)
            return int(
                np.ceil((width / spacing).magnitude)
            )  # Returns dimensionless count

        # --- BOTTOM REBAR ---
        self._n1_b = calculate_bars(self._s_b1_b, self.width)
        self._n2_b = calculate_bars(self._s_b2_b, self.width)
        self._n3_b = calculate_bars(self._s_b3_b, self.width)
        self._n4_b = calculate_bars(self._s_b4_b, self.width)

        # --- TOP REBAR ---
        self._n1_t = calculate_bars(self._s_b1_t, self.width)
        self._n2_t = calculate_bars(self._s_b2_t, self.width)
        self._n3_t = calculate_bars(self._s_b3_t, self.width)
        self._n4_t = calculate_bars(self._s_b4_t, self.width)


if __name__ == "__main__":
    # Example usage
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    slab = OneWaySlab(
        label="Slab 01",
        concrete=concrete,
        steel_bar=steelBar,
        thickness=200 * mm,
        width=1 * m,
    )
    # Set only position 1 bottom rebar (diameter=10mm, spacing=150mm)
    slab.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=20 * cm)
    # Set top rebar positions 1 and 2 (leave others unchanged)
    slab.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=15 * cm)
    debug(
        slab._A_s_bot, slab._A_s_top, slab._available_s_bot, slab._available_s_top
    )  # Debugging output for areas
    f1 = Forces(label="C1", M_y=20 * kNm)
    f2 = Forces(label="C1", M_y=-20 * kNm)
    node_1 = Node(slab, [f1, f2])
    print(node_1.check_flexure())
    node_1.flexure_results_detailed()
    # slab.plot()
