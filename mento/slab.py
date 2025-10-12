from __future__ import annotations
from typing import TYPE_CHECKING, Any
from dataclasses import dataclass
import math
import numpy as np
import pandas as pd
from devtools import debug

from mento.beam import RectangularBeam
from mento.forces import Forces
from mento.node import Node
from mento.material import (
    SteelBar,
    Concrete_ACI_318_19,
    # Concrete_EN_1992_2004,
)
from mento.rebar import Rebar
from mento.units import m, mm, cm, inch, MPa, kNm, kN, kip, ksi

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
        # Initialize with 12 mm or #4 bars in top & bottom.
        if self.concrete.unit_system == "metric":
            self._d_b1_b, self._d_b1_t = (
                10 * mm,
                10 * mm,
            )  # Diameter only, spacing left as 0
        else:
            self._d_b1_b, self._d_b1_t = (4 * inch) / 8, (4 * inch) / 8

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
                np.ceil((width.to("cm").magnitude / spacing.to("cm").magnitude))
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

    def _update_effective_heights(self) -> None:
        """Update effective heights and depths for moment and shear calculations."""
        # Consider average effective height so it works for analysis in X and Y directions for
        # a two-way slab.
        self._c_mec_bot = self.c_c + self._bot_rebar_centroid + self._d_b1_b / 2
        self._c_mec_top = self.c_c + self._top_rebar_centroid + self._d_b1_t / 2
        self._d_bot = self.height - self._c_mec_bot
        self._d_top = self.height - self._c_mec_top
        # Use bottom or top effective height
        self._d_shear = min(self._d_bot, self._d_top)


def slab_metric() -> None:
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    slab = OneWaySlab(
        label="Slab 01",
        concrete=concrete,
        steel_bar=steelBar,
        width=1 * m,
        height=20 * cm,
        c_c=25 * mm,
    )
    # Set only position 1 bottom rebar (diameter=10mm, spacing=150mm)
    # slab.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=20 * cm)
    # Set top rebar positions 1 and 2 (leave others unchanged)
    # slab.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=15 * cm)
    # debug(
    #     slab._A_s_bot, slab._A_s_top, slab._available_s_bot, slab._available_s_top
    # )  # Debugging output for areas
    f1 = Forces(label="C1", M_y=80 * kNm)
    # f2 = Forces(label="C1", V_z=30 * kN, M_y=-20 * kNm)
    node_1 = Node(slab, [f1])
    print(node_1.design_flexure())
    print(slab._n1_b, slab._d_b1_b, slab._n2_b, slab._d_b2_b)
    # print(slab.flexure_design_results_bot)
    # print(node_1.check_flexure())
    # node_1.flexure_results_detailed()
    # print(node_1.check_shear())
    # node_1.shear_results_detailed()
    # slab.plot()


def slab_imperial() -> None:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4 * ksi)
    steelBar = SteelBar(name="ADN 420", f_y=60 * ksi)
    slab = OneWaySlab(
        label="Slab 01",
        concrete=concrete,
        steel_bar=steelBar,
        width=12 * inch,
        height=7 * inch,
        c_c=0.75 * inch,
    )
    # Set only position 1 bottom rebar
    slab.set_slab_longitudinal_rebar_bot(d_b1=0.5 * inch, s_b1=12 * inch)
    # Set top rebar positions 1 and 2
    slab.set_slab_longitudinal_rebar_top(d_b1=0.5 * inch, s_b1=12 * inch)
    debug(
        slab._A_s_bot,
        slab._A_s_top,
        slab._n1_b,
        slab._available_s_bot,
        slab._available_s_top,
    )  # Debugging output for areas
    f1 = Forces(label="C1", M_y=20 * kNm)
    f2 = Forces(label="C1", V_z=1.52 * kip, M_y=-20 * kNm)
    node_1 = Node(slab, [f1, f2])
    # print(node_1.check_flexure())
    # node_1.flexure_results_detailed()
    node_1.check_shear()
    node_1.shear_results_detailed()


if __name__ == "__main__":
    slab_metric()
