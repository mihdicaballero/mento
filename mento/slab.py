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

        self._stirrup_d_b = 0 * mm
        self._stirrup_s_l = 0 * mm
        self._A_v = 0 * mm**2/m

        # Allow slabs to place as many bars as minimum spacing permits
        self.settings.max_bars_per_layer = max(self.settings.max_bars_per_layer, 200)
        # Force same diameter for all bars within a layer
        self.settings.max_diameter_diff = 0 * self.settings.max_diameter_diff

        self._initialize_longitudinal_rebar_attributes()

    def _create_rebar_designer(self) -> Rebar:
        """Use the slab-tailored rebar optimizer."""
        return SlabRebar(self)

    def _apply_longitudinal_design_bot(self, design: dict) -> None:
        merged = self._merge_layer_design(design)
        super()._apply_longitudinal_design_bot(merged)

    def _apply_longitudinal_design_top(self, design: dict) -> None:
        merged = self._merge_layer_design(design)
        super()._apply_longitudinal_design_top(merged)

    def _merge_layer_design(self, design: dict) -> dict:
        """Collapse multi-group beam layouts into slab-friendly layers."""
        merged = dict(design)

        n1 = int(design.get("n_1", 0))
        n2 = int(design.get("n_2", 0))
        total_layer1 = n1 + n2
        d1 = design.get("d_b1")
        if total_layer1 == 0:
            merged["n_1"] = 0
            merged["d_b1"] = None
        else:
            merged["n_1"] = total_layer1
            merged["d_b1"] = d1
        merged["n_2"] = 0
        merged["d_b2"] = None

        n3 = int(design.get("n_3", 0))
        n4 = int(design.get("n_4", 0))
        total_layer2 = n3 + n4
        if total_layer2 == 0:
            merged["n_3"] = 0
            merged["d_b3"] = None
        else:
            merged["n_3"] = total_layer1 if total_layer1 > 0 else total_layer2
            d3 = design.get("d_b3")
            merged["d_b3"] = d3 if d3 is not None else d1
        merged["n_4"] = 0
        merged["d_b4"] = None

        merged["total_bars"] = merged["n_1"] + merged["n_3"]
        return merged

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
                12 * mm,
                12 * mm,
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


class SlabRebar(Rebar):
    """Rebar optimizer that enforces slab-specific layout rules."""

    @staticmethod
    def _same_diameter(a: Any, b: Any) -> bool:
        if b is None:
            return a is None
        if a is None:
            return False
        try:
            return (a - b).magnitude == 0
        except AttributeError:
            return a == b

    def _normalize_slab_rows(self, df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df

        normalized_rows = []
        for _, row in df.iterrows():
            n1 = int(row.get("n_1", 0))
            n2 = int(row.get("n_2", 0))
            total_layer1 = n1 + n2
            if total_layer1 == 0:
                continue

            d1 = row.get("d_b1")
            d2 = row.get("d_b2")
            if d2 is not None and not self._same_diameter(d2, d1):
                continue

            n3 = int(row.get("n_3", 0))
            n4 = int(row.get("n_4", 0))
            total_layer2 = n3 + n4
            if total_layer2 not in (0, total_layer1):
                continue

            d3 = row.get("d_b3")
            d4 = row.get("d_b4")
            if d3 is not None and not self._same_diameter(d3, d1):
                continue
            if d4 is not None and not self._same_diameter(d4, d1):
                continue

            new_row = row.copy()
            new_row["n_1"] = total_layer1
            new_row["d_b1"] = d1
            new_row["n_2"] = 0
            new_row["d_b2"] = None
            if total_layer2 == 0:
                new_row["n_3"] = 0
                new_row["d_b3"] = None
            else:
                new_row["n_3"] = total_layer1
                new_row["d_b3"] = d1
            new_row["n_4"] = 0
            new_row["d_b4"] = None
            new_row["total_bars"] = new_row["n_1"] + new_row["n_3"]
            normalized_rows.append(new_row)

        if not normalized_rows:
            return df

        normalized_df = pd.DataFrame(normalized_rows)
        normalized_df = normalized_df.drop_duplicates(
            subset=["n_1", "d_b1", "n_3", "d_b3", "total_as"]
        )
        if "functional" in normalized_df.columns:
            normalized_df.sort_values(by=["functional"], inplace=True)
        normalized_df.reset_index(drop=True, inplace=True)
        return normalized_df

    def longitudinal_rebar_ACI_318_19(
        self, A_s_req: "Quantity", A_s_max: "Quantity" | None = None
    ) -> pd.DataFrame:
        super().longitudinal_rebar_ACI_318_19(A_s_req, A_s_max)
        normalized = self._normalize_slab_rows(self._long_combos_df)
        self._long_combos_df = normalized
        return normalized.head(10)


def slab_metric() -> None:
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    slab = OneWaySlab(
        label="Slab 01",
        concrete=concrete,
        steel_bar=steelBar,
        width = 1 * m,
        height = 20 * cm,
        c_c = 25 * mm
    )
    # Set only position 1 bottom rebar (diameter=10mm, spacing=150mm)
    slab.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=20 * cm)
    # Set top rebar positions 1 and 2 (leave others unchanged)
    slab.set_slab_longitudinal_rebar_top(d_b1=12 * mm, s_b1=15 * cm)
    debug(
        slab._A_s_bot, slab._A_s_top, slab._available_s_bot, slab._available_s_top
    )  # Debugging output for areas
    f1 = Forces(label="C1", M_y=20 * kNm)
    f2 = Forces(label="C1", V_z=30 * kN, M_y=-20 * kNm)
    node_1 = Node(slab, [f1, f2])
    print(node_1.check_flexure())
    node_1.flexure_results_detailed()
    # node_1.check_shear()
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
        height = 7 * inch,
        c_c = 0.75 * inch
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
