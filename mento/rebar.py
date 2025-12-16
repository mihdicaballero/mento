from __future__ import annotations
from typing import Any, Dict, TYPE_CHECKING, Tuple
import math
import pandas as pd
import numpy as np

from mento.units import psi, mm, cm, inch, MPa
from mento.material import Concrete_ACI_318_19

if TYPE_CHECKING:
    from mento.beam import RectangularBeam
    from pint import Quantity
    from pandas import DataFrame


class Rebar:
    def __init__(self, beam: RectangularBeam):
        """
        Initializes the Rebar object with the associated beam and settings.
        """

        self.mode = getattr(beam, "mode", "beam")  # "beam" or "slab"

        self.beam = beam
        self._long_combos_df: DataFrame = pd.DataFrame()
        self._trans_combos_df: DataFrame = pd.DataFrame()
        # Precompute spacing limits as floats in mm
        self._clear_limit_mm = self.beam.settings.clear_spacing.to("mm").magnitude
        self._vibrator_mm = self.beam.settings.vibrator_size.to("mm").magnitude
        self._clear_spacing = self.beam.settings.clear_spacing.to("mm")
        # Unit system default rebar
        if self.beam.concrete.unit_system == "metric":
            self.rebar_diameters = [
                6 * mm,
                8 * mm,
                10 * mm,
                12 * mm,
                16 * mm,
                20 * mm,
                25 * mm,
                32 * mm,
            ]
            self.rebar_areas = {d: (math.pi * d**2) / 4 for d in self.rebar_diameters}
        else:
            self.rebar_diameters = [
                3 * inch / 8,
                4 * inch / 8,
                5 * inch / 8,
                6 * inch / 8,
                7 * inch / 8,
                8 * inch / 8,
                1.128 * inch,
                1.27 * inch,
                1.41 * inch,
                1.693 * inch,
            ]
            rebar_areas_list = [(d**2 * np.pi / 4) for d in self.rebar_diameters]
            self.rebar_areas = dict(zip(self.rebar_diameters, rebar_areas_list))

    @property
    def longitudinal_rebar_design(self) -> DataFrame:
        return self._long_combos_df.iloc[0]

    @property
    def transverse_rebar_design(self) -> DataFrame:
        return self._trans_combos_df.iloc[0]

    ##########################################################
    # TRANSVERSE REBAR DESIGN
    ##########################################################

    def calculate_max_spacing_ACI_318_19(self, V_s_req: Quantity, A_cv: Quantity) -> Tuple[Quantity, Quantity]:
        """
        Calculate the maximum allowable spacing across the length and width of the beam
        based on design requirements.

        Parameters:
        -----------
        V_s_req : float
            Required shear force for the rebar.
        A_cv : float
            Effective shear area of the concrete section.

        Returns:
        --------
        tuple
            (s_max_l, s_max_w): The maximum spacing across the length and width of the beam.
        """

        f_c = self.beam.concrete.f_c
        if isinstance(self.beam.concrete, Concrete_ACI_318_19):
            lambda_factor = self.beam.concrete.lambda_factor

        # Determine maximum spacing based on V_s_req condition
        # Maximum spacing for lower shear demand
        if self.beam.concrete.unit_system == "metric":
            if V_s_req <= 0.083 * lambda_factor * math.sqrt(f_c / MPa) * MPa * A_cv:
                s_max_l = min(self.beam._d_shear / 2, 60 * cm)
                s_max_w = min(self.beam._d_shear, 60 * cm)
            else:
                # Maximum spacing for higher shear demand
                s_max_l = min(self.beam._d_shear / 4, 30 * cm)
                s_max_w = min(self.beam._d_shear / 2, 30 * cm)
        else:
            if V_s_req <= 4 * lambda_factor * math.sqrt(f_c / psi) * psi * A_cv:
                s_max_l = min(self.beam._d_shear / 2, 24 * inch)
                s_max_w = min(self.beam._d_shear, 24 * inch)
            else:
                s_max_l = min(self.beam._d_shear / 4, 12 * inch)
                s_max_w = min(self.beam._d_shear / 2, 12 * inch)

        return s_max_l, s_max_w

    def calculate_max_spacing_EN_1992_2004(self, alpha: float) -> Tuple[Quantity, Quantity]:
        """
        Calculate the maximum allowable spacing across the length and width of the beam
        based on design requirements for EN 1992-2004.

        Parameters:
        -----------
        alpha: stirrups angle

        Returns:
        --------
        tuple
            (s_max_l, s_max_w): The maximum spacing along the length and width of the beam.
        """
        # Maximum of 40 cm for stirrup spacing
        s_max_l = min(0.75 * self.beam._d_shear * (1 + 1 / math.tan(alpha)), 40 * cm)
        s_max_w = min(0.75 * self.beam._d_shear, 60 * cm)

        return s_max_l, s_max_w

    def transverse_rebar_ACI_318_19(self, V_s_req: Quantity) -> Any:
        if self.beam.concrete.unit_system == "metric":
            valid_diameters = self.rebar_diameters[2:5]  # Minimum 10 mm
        else:
            valid_diameters = self.rebar_diameters[0:3]

        A_cv = self.beam.width * self.beam._d_shear
        s_max_l, s_max_w = self.calculate_max_spacing_ACI_318_19(V_s_req, A_cv)

        return valid_diameters, s_max_l, s_max_w

    def transverse_rebar_CIRSOC_201_25(self, V_s_req: Quantity) -> Any:
        valid_diameters = self.rebar_diameters[0:4]  # Minimum 6 mm, maximum 12 mm

        A_cv = self.beam.width * self.beam._d_shear
        s_max_l, s_max_w = self.calculate_max_spacing_ACI_318_19(V_s_req, A_cv)

        return valid_diameters, s_max_l, s_max_w

    def transverse_rebar_EN_1992_2004(self, alpha: float) -> Any:
        valid_diameters = self.rebar_diameters[0:4]  # Minimum 6 mm, maximum 12 mm

        s_max_l, s_max_w = self.calculate_max_spacing_EN_1992_2004(alpha)

        return valid_diameters, s_max_l, s_max_w

    def transverse_rebar(self, A_v_req: Quantity, V_s_req: Quantity, alpha: float) -> DataFrame:
        """
        Computes the required transverse reinforcement based on ACI 318-19.

        Args:
            A_v_req: Required area for transverse reinforcement.

        Returns:
            A dictionary containing the transverse rebar design.
        """

        # Prepare the list for valid combinations
        valid_combinations = []

        # Get code specific limitations
        if self.beam.concrete.design_code == "ACI 318-19":
            valid_diameters, s_max_l, s_max_w = self.transverse_rebar_ACI_318_19(V_s_req)
        elif self.beam.concrete.design_code == "CIRSOC 201-25":
            valid_diameters, s_max_l, s_max_w = self.transverse_rebar_CIRSOC_201_25(V_s_req)
        elif self.beam.concrete.design_code == "EN 1992-2004":
            valid_diameters, s_max_l, s_max_w = self.transverse_rebar_EN_1992_2004(alpha)

        # Iterate through available diameters
        for d_b in valid_diameters:
            # Start with 2 legs = 1 stirrup
            n_legs = 2

            # Start with maximum allowed spacing s_max_l
            if self.beam.concrete.unit_system == "metric":
                s_l: Quantity = math.floor(s_max_l.to("cm").magnitude) * cm
            else:
                s_l = math.floor(s_max_l.to("inch").magnitude) * inch

            while True:
                # Calculate spacing based on current legs
                n_stirrups = math.ceil(n_legs / 2)  # Number of stirrups based on number of legs
                n_legs_actual = n_stirrups * 2  # Ensure legs are even
                # Consider 1 leg less for spacing laong width
                s_w = (self.beam.width - 2 * self.beam.c_c - self.beam._stirrup_d_b) / (n_legs_actual - 1)  # noqa: F841

                A_db = self.rebar_areas[d_b]  # Area of a stirrup bar
                A_vs = n_legs_actual * A_db  # Area of vertical stirrups
                A_v: Quantity = A_vs / s_l  # Area of vertical stirrups per unit length

                # Store the valid combination if spacing is also valid
                if self.beam.concrete.unit_system == "metric":
                    # Check if the calculated A_v meets or exceeds the required A_v
                    if A_v >= A_v_req:
                        valid_combinations.append(
                            {
                                "n_stir": int(n_stirrups),
                                "d_b": d_b,
                                "s_l": s_l.to("cm"),  # spacing along length
                                "s_w": s_w.to("cm"),  # spacing along width
                                "A_v": A_v.to("cm**2/m"),
                                "s_max_l": s_max_l.to("cm"),
                                "s_max_w": s_max_w.to("cm"),
                            }
                        )
                        # Stop checking larger diameters
                        break

                    # If A_v is insufficient, reduce s_l by 1 cm
                    s_l -= 1 * cm
                    # If s_l becomes less than 5 cm, increase the number of legs by 2 and reset s_l to s_max_l
                    if s_l < 5 * cm:  # If spacing is less than 5 cm, increase 1 stirrup
                        n_legs += 2
                        s_l = math.floor(s_max_l.to("cm").magnitude) * cm  # Reset s_l to the max allowed spacing
                else:
                    # Check if the calculated A_v meets or exceeds the required A_v
                    if A_v >= A_v_req:
                        valid_combinations.append(
                            {
                                "n_stir": int(n_stirrups),
                                "d_b": d_b,
                                "s_l": s_l.to("inch"),  # spacing along length
                                "s_w": s_w.to("inch"),  # spacing along width
                                "A_v": A_v.to("inch**2/ft"),
                                "s_max_l": s_max_l.to("inch"),
                                "s_max_w": s_max_w.to("inch"),
                            }
                        )
                        # Stop checking larger diameters
                        break

                    # If A_v is insufficient, reduce s_l by 1 inch
                    s_l -= 1 * inch
                    # If s_l becomes less than 2 inch, increase the number of legs by 2 and reset s_l to s_max_l
                    if s_l < 2 * inch:  # If spacing is less than 2 inch, increase 1 stirrup
                        n_legs += 2
                        s_l = math.floor(s_max_l.to("inch").magnitude) * inch  # Reset s_l to the max allowed spacing

        # Create a DataFrame with all valid combinations
        df_combinations = pd.DataFrame(valid_combinations)

        # Sort combinations by the total rebar area required (ascending)
        # Sort by 'A_v' first, then by 'n_stir' to prioritize fewer bars
        df_combinations.sort_values(by=["n_stir", "A_v"], inplace=True)
        df_combinations.reset_index(drop=True, inplace=True)
        self._trans_combos_df = df_combinations
        return df_combinations

    ##########################################################
    # LONGITUDINAL REBAR DESIGN
    ##########################################################

    def longitudinal_rebar_ACI_318_19(
        self,
        A_s_req: Quantity,
        A_s_max: Quantity | None = None,
        mech_cover: Quantity | None = None,
    ) -> DataFrame:
        """
        Computes the required longitudinal reinforcement based on ACI 318-19.

        Args:
            A_s_req: Required longitudinal rebar area.
            A_s_max: Optional maximum allowable longitudinal rebar area. If not
                provided, the limit defaults to ``10 * A_s_req``.

        Returns:
            A DataFrame containing the best combinations of rebar details. If no
            combination satisfies ``A_s_req``, returns the combination with the
            maximum possible area.
        """
        self.A_s_req = A_s_req
        self.original_mech_cover = mech_cover
        effective_width = self.beam.width - 2 * (self.beam.c_c + self.beam._stirrup_d_b)

        # --- Early exit: no steel required -------------------------------------
        if A_s_req.to("cm**2").magnitude == 0:
            df = pd.DataFrame(
                [
                    {
                        "n_1": 0,
                        "d_b1": 0 * mm,
                        "n_2": 0,
                        "d_b2": None,
                        "n_3": 0,
                        "d_b3": None,
                        "n_4": 0,
                        "d_b4": None,
                        "total_as": 0 * cm**2,
                        "total_bars": 0,
                        "clear_spacing": self.beam.settings.clear_spacing.to("mm"),
                    }
                ]
            )
            self._long_combos_df = df
            return df

        # Variables to track the combinations
        valid_combinations = []
        best_fallback_combination = None  # To store the best fallback design
        max_fallback_area = 0 * cm**2  # To track the maximum area in fallback cases
        # Create a list of rebar diameters that are equal to or greater than the minimum diameter
        if self.beam.concrete.unit_system == "metric":
            self.min_long_rebar = 10 * mm
        else:
            self.min_long_rebar = 3 * inch / 8

        # Filter valid rebar diameters based on the minimum longitudinal diameter per design code and beam settings
        valid_rebar_diameters = [
            d
            for d in self.rebar_diameters
            if d >= self.min_long_rebar
            and d >= self.beam.settings.minimum_longitudinal_diameter
            and d <= self.beam.settings.max_longitudinal_diameter
        ]

        for d_b1 in valid_rebar_diameters:  # Without taking Ø6 as a possible solution
            for d_b2 in [d for d in valid_rebar_diameters if d <= d_b1]:
                # Quick upper-bound check for this diameter pair
                area1 = self.rebar_areas[d_b1]  # cm²
                area2 = self.rebar_areas[d_b2]

                max_bars = self.beam.settings.max_bars_per_layer

                # Max for layer 1 (n1 fixed, n2 up to max_bars)
                A_layer1_max = 2 * area1 + max_bars * area2  # n1=2 fixed

                if self.mode == "slab":
                    # Slab: mirror layer 2 => max total As = 2 * layer1
                    A_total_max = 2 * A_layer1_max
                else:
                    # Beam: rough upper bound (could be improved)
                    A_total_max = 2 * A_layer1_max
                # Max limit for As (same logic as later)
                max_limit = 10 * A_s_req
                if A_s_max is not None:
                    max_limit = min(max_limit, A_s_max)

                # If even the maximum possible As with these diameters
                # is less than required, skip all n2/n3/n4 loops.
                if A_total_max < min(A_s_req, max_limit):
                    continue

                for d_b3 in [d for d in valid_rebar_diameters if d <= d_b2]:
                    for d_b4 in [d for d in valid_rebar_diameters if d <= d_b3]:
                        # Condition 5: |d_b1 - d_b2| and |d_b3 - d_b4| must not exceed max_diameter_diff
                        # Apply diameter difference condition across all combinations
                        # Ensure that no two bars exceed `max_diameter_diff`
                        diameters = [d_b1, d_b2, d_b3, d_b4]
                        diameters = [d for d in diameters if d is not None]  # Filter out None values for unused bars

                        # Apply the diameter difference check across all bars
                        if not self._check_diameter_differences(diameters):
                            continue  # Skip this combination if any diameter pair exceeds max_diameter_diff`

                        n1 = 0 if self.A_s_req == 0 * cm**2 else 2  # This is a fixed value for every beam

                        # Iterate over possible numbers of bars in each group
                        for n2 in range(0, self.beam.settings.max_bars_per_layer + 1):  # n2 can be 0 or more
                            if n1 + n2 > self.beam.settings.max_bars_per_layer:
                                continue  # Skip if the total bars in layer 1 exceed the limit

                            spacing_ok = True

                            # ALWAYS compute spacing for this layer-1 combo
                            if not self._check_spacing(n1, n2, d_b1, d_b2, effective_width):
                                continue

                            # Calculate area for layer 1
                            A_s_layer_1 = n1 * self.rebar_areas[d_b1] + (
                                n2 * self.rebar_areas[d_b2] if n2 > 0 else 0 * cm**2
                            )

                            # Condition 6 and 7: Check clear spacing in layer 1
                            if n2 > 0:
                                spacing_ok = self._check_spacing(n1, n2, d_b1, d_b2, effective_width)
                                if not spacing_ok:
                                    continue

                            max_limit = 10 * A_s_req
                            if A_s_max is not None:
                                max_limit = min(max_limit, A_s_max)
                            max_limit = max(max_limit, n1 * self.rebar_areas[self.min_long_rebar])

                            if A_s_layer_1 > max_limit:
                                break  # further n2 will only increase area

                            # Check if total area from layer 1 is enough for required A_s
                            # And also less than the maximum limit
                            if A_s_layer_1 >= A_s_req and A_s_layer_1 <= max_limit:
                                total_as = A_s_layer_1  # Only consider layer 1
                                total_bars = n1 + n2  # Total bars only in layer 1
                                valid_combinations.append(
                                    {
                                        "n_1": n1,
                                        "d_b1": d_b1,
                                        "n_2": n2,
                                        "d_b2": d_b2 if n2 > 0 else None,  # Display as None if n2 is 0
                                        "n_3": 0,  # No bars in layer 2
                                        "d_b3": None,
                                        "n_4": 0,  # No bars in layer 2
                                        "d_b4": None,
                                        "total_as": total_as.to("cm**2"),
                                        "total_bars": total_bars,
                                        "clear_spacing": self._clear_spacing.to("mm"),
                                    }
                                )
                            else:
                                # Track the combination with the maximum possible area (fallback)
                                if A_s_layer_1 > max_fallback_area and A_s_layer_1 <= max_limit:
                                    max_fallback_area = A_s_layer_1
                                    best_fallback_combination = {
                                        "n_1": n1,
                                        "d_b1": d_b1,
                                        "n_2": n2,
                                        "d_b2": d_b2 if n2 > 0 else None,
                                        "n_3": 0,
                                        "d_b3": None,
                                        "n_4": 0,
                                        "d_b4": None,
                                        "total_as": A_s_layer_1.to("cm**2"),
                                        "total_bars": n1 + n2,
                                        "clear_spacing": self._clear_spacing.to("mm"),
                                    }

                            # =============================================================
                            # --- Layer 2 combinations (beam vs slab logic) ---
                            # =============================================================

                            if self.mode == "slab":
                                # ---------------------------------------------------------
                                # In slab mode:
                                #  - Only two cases are considered:
                                #       (1) one single layer
                                #       (2) second layer identical to the first (mirror)
                                #  - n3 = n1, n4 = n2 when a second layer exists
                                # ---------------------------------------------------------
                                for has_second_layer in [False, True]:
                                    if not has_second_layer:
                                        n3, n4 = 0, 0
                                    else:
                                        n3, n4 = n1, n2

                                    # --- Compute total reinforcement in layer 2 ------------
                                    A_s_layer_2 = n3 * self.rebar_areas[d_b3] + (
                                        n4 * self.rebar_areas[d_b4] if n4 > 0 else 0 * cm**2
                                    )

                                    # Condition 4: area of layer 1 ≥ area of layer 2
                                    if A_s_layer_1 < A_s_layer_2:
                                        continue

                                    # Condition 6 and 7: check clear spacing in layer 2
                                    if n4 > 0 and not self._check_spacing(n3, n4, d_b3, d_b4, effective_width):
                                        continue

                                    # --- Compute total reinforcement and evaluate -----------
                                    total_as = A_s_layer_1 + A_s_layer_2
                                    if total_as >= A_s_req and total_as <= max_limit:
                                        total_bars = n1 + n2 + n3 + n4
                                        valid_combinations.append(
                                            {
                                                "n_1": n1,
                                                "d_b1": d_b1,
                                                "n_2": n2,
                                                "d_b2": d_b2 if n2 > 0 else None,
                                                "n_3": n3,
                                                "d_b3": d_b3 if n3 > 0 else None,
                                                "n_4": n4,
                                                "d_b4": d_b4 if n4 > 0 else None,
                                                "total_as": total_as.to("cm**2"),
                                                "total_bars": total_bars,
                                                "clear_spacing": self._clear_spacing.to("mm"),
                                            }
                                        )
                                    else:
                                        # Track fallback combination with maximum As
                                        if total_as > max_fallback_area and total_as <= max_limit:
                                            max_fallback_area = total_as
                                            best_fallback_combination = {
                                                "n_1": n1,
                                                "d_b1": d_b1,
                                                "n_2": n2,
                                                "d_b2": d_b2 if n2 > 0 else None,
                                                "n_3": n3,
                                                "d_b3": d_b3 if n3 > 0 else None,
                                                "n_4": n4,
                                                "d_b4": d_b4 if n4 > 0 else None,
                                                "total_as": total_as.to("cm**2"),
                                                "total_bars": n1 + n2 + n3 + n4,
                                                "clear_spacing": self._clear_spacing.to("mm"),
                                            }

                            else:
                                # ---------------------------------------------------------
                                # Normal beam logic (original)
                                # ---------------------------------------------------------
                                for n3 in [0, 2]:
                                    for n4 in range(0, self.beam.settings.max_bars_per_layer + 1):
                                        # Ensure layer 2 bars are not more than layer 1 bars
                                        if n3 + n4 > n1 + n2:
                                            continue
                                        if n3 == 0 and n4 > 0:
                                            continue
                                        if n3 + n4 > self.beam.settings.max_bars_per_layer:
                                            continue

                                        A_s_layer_2 = n3 * self.rebar_areas[d_b3] + (
                                            n4 * self.rebar_areas[d_b4] if n4 > 0 else 0 * cm**2
                                        )

                                        if A_s_layer_1 < A_s_layer_2:
                                            continue
                                        if n4 > 0 and not self._check_spacing(n3, n4, d_b3, d_b4, effective_width):
                                            continue

                                        total_as = A_s_layer_1 + A_s_layer_2
                                        if total_as >= A_s_req and total_as <= max_limit:
                                            total_bars = n1 + n2 + n3 + n4
                                            valid_combinations.append(
                                                {
                                                    "n_1": n1,
                                                    "d_b1": d_b1,
                                                    "n_2": n2,
                                                    "d_b2": d_b2 if n2 > 0 else None,
                                                    "n_3": n3,
                                                    "d_b3": d_b3 if n3 > 0 else None,
                                                    "n_4": n4,
                                                    "d_b4": d_b4 if n4 > 0 else None,
                                                    "total_as": total_as.to("cm**2"),
                                                    "total_bars": total_bars,
                                                    "clear_spacing": self._clear_spacing.to("mm"),
                                                }
                                            )
                                        else:
                                            if total_as > max_fallback_area and total_as <= max_limit:
                                                max_fallback_area = total_as
                                                best_fallback_combination = {
                                                    "n_1": n1,
                                                    "d_b1": d_b1,
                                                    "n_2": n2,
                                                    "d_b2": d_b2 if n2 > 0 else None,
                                                    "n_3": n3,
                                                    "d_b3": d_b3 if n3 > 0 else None,
                                                    "n_4": n4,
                                                    "d_b4": d_b4 if n4 > 0 else None,
                                                    "total_as": total_as.to("cm**2"),
                                                    "total_bars": n1 + n2 + n3 + n4,
                                                    "clear_spacing": self._clear_spacing.to("mm"),
                                                }

        # Convert valid combinations to DataFrame
        df = pd.DataFrame(valid_combinations)
        # Drop duplicate rows based on the specified columns
        df = df.drop_duplicates(subset=["n_1", "d_b1", "n_2", "d_b2", "n_3", "d_b3", "n_4", "d_b4"])

        # If no valid combinations satisfy A_s_req, use the best fallback combination
        if df.empty and best_fallback_combination is not None:
            df = pd.DataFrame([best_fallback_combination])

        # Only calculate penalties if we have valid combinations
        if not df.empty:
            modified_df = self._calculate_penalties_long_rebar(df)
            # Sort by 'Functional' to sort by the best options
            modified_df.sort_values(by=["functional"], inplace=True)
            modified_df.reset_index(drop=True, inplace=True)
            self._long_combos_df = modified_df
            return modified_df.head(10)
        else:
            # Return empty DataFrame with expected structure if no combinations found
            self._long_combos_df = df
            return df

    def longitudinal_rebar_EN_1992_2004(self, A_s_req: Quantity) -> None:
        self.longitudinal_rebar_ACI_318_19(A_s_req)  # TODO WE HAVE TO CHANGE THIS

    def _check_spacing(
        self,
        n1: int,
        n2: int,
        d_b1: Quantity,
        d_b2: Quantity,
        effective_width: Quantity,
    ) -> bool:
        """
        Checks the clear spacing between rebars in a layer.

        Parameters:
            n1 (int): Number of bars in the first group of the layer.
            n2 (int): Number of bars in the second group of the layer.
            d_b1 (Quantity): Diameter of bars in the first group of the layer.
            d_b2 (Quantity): Diameter of bars in the second group of the layer.
            effective_width (Quantity): The effective width available for bar placement.

        Returns:
            bool: True if the clear spacing satisfies the design limits, False otherwise.
        """

        # Convert everything once to mm (floats)
        eff_mm = effective_width.to("mm").magnitude
        d1_mm = d_b1.to("mm").magnitude
        d2_mm = d_b2.to("mm").magnitude

        total_bars = n1 + n2

        total_bar_width_mm = n1 * d1_mm + n2 * d2_mm
        clear_mm = (eff_mm - total_bar_width_mm) / (total_bars - 1)

        # Store as Quantity for reporting
        self._clear_spacing = clear_mm * mm

        # Determine the maximum clear spacing limit
        # Max required clear spacing in mm (precomputed)
        max_clear_spacing_mm = max(
            self._clear_limit_mm,
            self._vibrator_mm,
            max(d1_mm, d2_mm),
        )

        # Check if the clear spacing is within limits
        return clear_mm >= max_clear_spacing_mm

    def _check_diameter_differences(self, diameters: list) -> bool:
        """
        Checks that all diameter differences between bars in the list
        do not exceed max_diameter_diff.
        """
        for i, d1 in enumerate(diameters):
            for d2 in diameters[i + 1 :]:
                if abs(d1 - d2) > self.beam.settings.max_diameter_diff:
                    return False
        return True

    def _calculate_penalties_long_rebar(
        self,
        df: pd.DataFrame,
        alpha: float = 3.5,
        beta: float = 0.30,
        gamma: float = 0.25,
        delta: float = 1,
        epsilon: float = 0,  # No penalty for beams, just slabs
    ) -> pd.DataFrame:
        """
        Calculate penalties for rebar configurations and add them as columns to the DataFrame.

        Args:
            df (pd.DataFrame): The input DataFrame containing rebar configurations.
            alpha (float): Weight for the number of bars penalty.
            beta (float): Weight for the diameter difference penalty.
            gamma (float): Weight for the layer penalty.

        Returns:
            pd.DataFrame: The modified DataFrame with penalty columns and the final functional.
        """

        # Adjust penalty weights depending on element type
        if getattr(self, "mode", "beam") == "slab":
            # Slabs prefer many small bars → penalize large diameters and spacing more
            alpha *= 0.8  # reduce area weight (slabs have smaller demand differences)
            beta *= 0.1  # stronger sensitivity to number of bars
            gamma *= 0.5  # less penalty for multi-diameter use
            delta *= 0.5  # allow two layers if needed
            epsilon = 0.75  # strong penalty on large bar sizes

        # Calculate minimum bars and minimum area of steel
        min_bars = df["total_bars"].min()
        min_as = df["total_as"].min()

        # Diameter difference penalty function
        def diameter_difference_penalty(row: pd.Series) -> float:
            """
            Calculate the penalty for variation in bar diameters.

            Args:
                row (pd.Series): A row of the DataFrame.

            Returns:
                float: The standard deviation of the diameters (penalty).
            """
            diameters = []
            if row["n_1"] > 0:
                diameters.extend([row["d_b1"].magnitude] * row["n_1"])
            if row["n_2"] > 0:
                diameters.extend([row["d_b2"].magnitude] * row["n_2"])
            if row["n_3"] > 0:
                diameters.extend([row["d_b3"].magnitude] * row["n_3"])
            if row["n_4"] > 0:
                diameters.extend([row["d_b4"].magnitude] * row["n_4"])
            return np.std(diameters) if diameters else 0.0

        # Layer penalty function
        def layer_penalty(row: pd.Series) -> int:
            """
            Calculate the penalty for using additional layers of reinforcement.

            Args:
                row (pd.Series): A row of the DataFrame.

            Returns:
                int: 1 if additional layers are used, 0 otherwise.
            """
            # Check if any bars are in n_3, n_4, or future layers
            if row["n_3"] > 0 or row["n_4"] > 0:
                return 1  # Penalize if additional layers are used
            return 0  # No penalty if only the first layer is used

        # Calculate penalties and add them as columns
        df["area_penalty"] = df.apply(lambda row: (alpha * row["total_as"] / min_as), axis=1)
        # Prefer moderate bar counts, where very high or very low will be penalized
        df["bars_penalty"] = beta * ((df["total_bars"] - min_bars) / min_bars) ** 2
        df["diameter_penalty"] = gamma * df.apply(diameter_difference_penalty, axis=1)
        df["layer_penalty"] = delta * df.apply(layer_penalty, axis=1)

        # Diameter size penalty, for very large or very small
        min_d = df[["d_b1", "d_b2", "d_b3", "d_b4"]].stack().dropna().apply(lambda x: x.magnitude).min()
        df["diameter_size_penalty"] = epsilon * (
            df.apply(
                lambda r: (
                    max(
                        [
                            d.magnitude
                            for d, n in zip(
                                [r["d_b1"], r["d_b2"], r["d_b3"], r["d_b4"]],
                                [r["n_1"], r["n_2"], r["n_3"], r["n_4"]],
                            )
                            if n > 0 and d is not None
                        ]
                    )
                    / min_d
                    - 1
                ),
                axis=1,
            )
        )
        # Slab penalty for large spacing
        if getattr(self, "mode", "beam") == "slab":
            max_spacing_allowed = 300  # mm
            df["spacing_penalty"] = df["clear_spacing"].apply(
                lambda s: 0
                if s.magnitude <= max_spacing_allowed
                else (s.magnitude - max_spacing_allowed) / max_spacing_allowed
            )
        else:
            df["spacing_penalty"] = 0

        # Calculate the final functional
        df["functional"] = df.apply(
            lambda row: (
                row["area_penalty"]
                + row["bars_penalty"]
                + row["diameter_penalty"]
                + row["layer_penalty"]
                + row["diameter_size_penalty"]
                + row["spacing_penalty"]
            ),
            axis=1,
        )

        return df

    def longitudinal_rebar(self, A_s_req: Quantity, A_s_max: Quantity | None = None) -> Dict[str, Any]:
        """
        Selects the appropriate longitudinal rebar method based on the design
        code.

        Args:
            A_s_req: Required longitudinal rebar area.
            A_s_max: Optional maximum allowable longitudinal rebar area.
        """
        if self.beam.concrete.design_code == "ACI 318-19" or self.beam.concrete.design_code == "CIRSOC 201-25":
            return self.longitudinal_rebar_ACI_318_19(A_s_req, A_s_max)
        elif self.beam.concrete.design_code == "EN 1992-2004":
            return self.longitudinal_rebar_EN_1992_2004(A_s_req)
