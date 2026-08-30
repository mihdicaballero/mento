from __future__ import annotations
from typing import Any, Dict, TYPE_CHECKING, Tuple
import math
import pandas as pd
import numpy as np

from mento.codes.aci_318_19.equations import shear as aci_shear_eq
from mento.codes.en_1992_2004.equations import shear as en_shear_eq
from mento.precompute import CANONICAL, DISPLAY, section_floats
from mento.units import mm, cm, inch

if TYPE_CHECKING:
    from mento.beam import RectangularBeam
    from pint import Quantity
    from pandas import DataFrame

# `cm**2` raises a pint Unit to a power, which is far from free. The rebar search
# builds one area Quantity per candidate it keeps, so the unit is built once here
# instead of on every call.
_CM2 = cm**2


def max_stirrup_spacing_ACI_318_19(beam: RectangularBeam, V_s_req: float, A_cv: float) -> Tuple[float, float]:
    """Stirrup spacing limits of ACI 318-19 Table 9.7.6.2.2, for a beam.

    Floats in the beam's own unit system, in and out: the shear check runs
    entirely in floats, so a pint signature here would put the boundary back in
    the middle of it. :meth:`Rebar.calculate_max_spacing_ACI_318_19` wraps this
    for the design path, which still speaks pint.

    A module-level function rather than a ``Rebar`` method because it needs
    nothing from the bar catalogue: building a whole ``Rebar`` for it cost more
    than the check it serves.
    """
    sec = section_floats(beam)
    return aci_shear_eq.max_stirrup_spacing(
        V_s_req,
        sec.f_c,
        beam.concrete.lambda_factor,
        A_cv,
        sec.d_shear,
        is_imperial=sec.is_imperial,
    )


def max_stirrup_spacing_EN_1992_2004(beam: RectangularBeam, alpha: float) -> Tuple[float, float]:
    """Stirrup spacing limits of EN 1992-1-1 §9.2.2(6) and (8), for a beam.

    Floats in mm, in and out, for the same reason as the ACI counterpart.
    :meth:`Rebar.calculate_max_spacing_EN_1992_2004` wraps it for the design
    path, which still speaks pint.
    """
    return en_shear_eq.max_stirrup_spacing(section_floats(beam).d_shear, alpha)


class RebarDesignInfeasibleError(Exception):
    """Raised when no valid rebar combination fits the section geometry and
    the code-imposed limits (A_s_req, A_s_max, bar diameter, spacing, layers).

    Typical trigger: very narrow sections combined with high A_s_req (small
    b + high fy or high Mu). Callers should catch this and either fall back
    to "best-effort" behavior or surface the infeasibility to the user.
    """


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
        if self._long_combos_df.empty:
            raise RebarDesignInfeasibleError(
                "No valid longitudinal rebar combination found — the required "
                "steel area cannot be fit in the section given the geometry "
                "(width, max diameter, layers, spacing limits)."
            )
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

        Parameters
        ----------
        V_s_req : float
            Required shear force for the rebar.
        A_cv : float
            Effective shear area of the concrete section.

        Returns
        -------
        tuple
            (s_max_l, s_max_w): The maximum spacing across the length and width of the beam.
        """

        sec = section_floats(self.beam)
        canonical = CANONICAL[sec.is_imperial]
        display = DISPLAY[sec.is_imperial]
        s_max_l, s_max_w = max_stirrup_spacing_ACI_318_19(
            self.beam,
            V_s_req.to(canonical["force"]).magnitude,
            A_cv.to(canonical["area"]).magnitude,
        )
        return (
            (s_max_l * canonical["length"]).to(display["length"]),
            (s_max_w * canonical["length"]).to(display["length"]),
        )

    def calculate_max_spacing_EN_1992_2004(self, alpha: float) -> Tuple[Quantity, Quantity]:
        """
        Calculate the maximum allowable spacing across the length and width of the beam
        based on design requirements for EN 1992-2004.

        Parameters
        ----------
        alpha: stirrups angle

        Returns
        -------
        tuple
            (s_max_l, s_max_w): The maximum spacing along the length and width of the beam.
        """
        s_max_l, s_max_w = max_stirrup_spacing_EN_1992_2004(self.beam, alpha)
        return (s_max_l * mm).to(cm), (s_max_w * mm).to(cm)

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

    def min_legs_along_width(self, d_b: Quantity, s_max_w: Quantity) -> int:
        """
        Smallest even number of stirrup legs whose transverse spacing fits within ``s_max_w``.

        The legs are spread evenly across the section, so the ``n_legs - 1`` gaps between
        them have to cover the distance separating the centres of the outermost pair,
        ``width - 2 * c_c - d_b``. Both ACI 318-19 Table 9.7.6.2.2 and EN 1992-1-1 9.2.2(8)
        cap that gap, which is what forces a wide section to carry more than the two legs of
        a single stirrup regardless of how much area the shear demand asks for.

        Parameters
        ----------
        d_b : Quantity
            Diameter of the stirrup bar being tried.
        s_max_w : Quantity
            Maximum spacing allowed across the width of the beam.

        Returns
        -------
        int
            Number of legs, always even and never below 2.
        """
        outer_span = self.beam.width - 2 * self.beam.c_c - d_b
        if outer_span <= 0 * cm or s_max_w <= 0 * cm:
            return 2
        # The tolerance keeps a span that divides exactly from rounding up a whole gap.
        n_gaps = math.ceil((outer_span / s_max_w).to("dimensionless").magnitude - 1e-9)
        n_legs = max(2, n_gaps + 1)
        # Legs come in pairs, one closed stirrup each.
        return n_legs + (n_legs % 2)

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
            # Start from the fewest legs that keep the transverse spacing within s_max_w,
            # rather than from a single stirrup: on a wide section two legs never comply.
            n_legs = self.min_legs_along_width(d_b, s_max_w)

            # Start with maximum allowed spacing s_max_l
            if self.beam.concrete.unit_system == "metric":
                s_l: Quantity = math.floor(s_max_l.to("cm").magnitude) * cm
            else:
                s_l = math.floor(s_max_l.to("inch").magnitude) * inch

            while True:
                # Calculate spacing based on current legs
                n_stirrups = math.ceil(n_legs / 2)  # Number of stirrups based on number of legs
                n_legs_actual = n_stirrups * 2  # Ensure legs are even
                # n_legs - 1 gaps span the distance between the outermost leg centres.
                # This uses the candidate diameter: self.beam._stirrup_d_b still holds
                # whatever the previous pass left behind, which is not what is being tried.
                s_w = (self.beam.width - 2 * self.beam.c_c - d_b) / (n_legs_actual - 1)

                A_db = self.rebar_areas[d_b]  # Area of a stirrup bar
                A_vs = n_legs_actual * A_db  # Area of vertical stirrups
                A_v: Quantity = A_vs / s_l  # Area of vertical stirrups per unit length

                # Store the valid combination if spacing is also valid
                if self.beam.concrete.unit_system == "metric":
                    # Check if the calculated A_v meets or exceeds the required A_v, and
                    # that the legs are close enough together across the width.
                    if A_v >= A_v_req and s_w <= s_max_w:
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
                    # Check if the calculated A_v meets or exceeds the required A_v, and
                    # that the legs are close enough together across the width.
                    if A_v >= A_v_req and s_w <= s_max_w:
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
        max_fallback_cm2 = 0.0  # To track the maximum area in fallback cases
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

        # --- Float view of the search space ------------------------------------
        # The loops below evaluate tens of thousands of candidate layouts. Running
        # that arithmetic on pint Quantities was the dominant cost of the whole
        # flexure design, so units are stripped once here and re-applied only to
        # the handful of combinations actually kept (ADR-0005). Areas are in cm²
        # and diameters in mm throughout this block.
        areas_cm2 = [self.rebar_areas[d].to("cm**2").magnitude for d in valid_rebar_diameters]
        diams_mm = [d.to("mm").magnitude for d in valid_rebar_diameters]
        A_s_req_cm2 = A_s_req.to("cm**2").magnitude
        eff_width_mm = effective_width.to("mm").magnitude
        max_diam_diff_mm = self.beam.settings.max_diameter_diff.to("mm").magnitude
        max_bars = self.beam.settings.max_bars_per_layer

        # n1 is fixed at 2, and A_s_req > 0 is guaranteed by the early exit above,
        # so both area limits are loop-invariant and are computed once.
        n1 = 2
        skip_limit_cm2 = 10 * A_s_req_cm2
        if A_s_max is not None:
            skip_limit_cm2 = min(skip_limit_cm2, A_s_max.to("cm**2").magnitude)
        max_limit_cm2 = max(skip_limit_cm2, n1 * self.rebar_areas[self.min_long_rebar].to("cm**2").magnitude)

        # valid_rebar_diameters is ascending, so "every diameter <= d_bN" is a
        # prefix of it and the nested loops can walk indices instead of
        # re-filtering the list with a pint comparison on each pass.
        for i1 in range(len(valid_rebar_diameters)):
            area1, d1_mm = areas_cm2[i1], diams_mm[i1]
            for i2 in range(i1 + 1):
                area2, d2_mm = areas_cm2[i2], diams_mm[i2]

                # Quick upper-bound check for this diameter pair.
                # Max for layer 1 (n1 fixed, n2 up to max_bars); layer 2 can at
                # best mirror it, for both beam and slab mode.
                A_layer1_max = 2 * area1 + max_bars * area2  # n1=2 fixed
                A_total_max = 2 * A_layer1_max

                # If even the maximum possible As with these diameters
                # is less than required, skip all n2/n3/n4 loops.
                if A_total_max < min(A_s_req_cm2, skip_limit_cm2):
                    continue

                for i3 in range(i2 + 1):
                    area3 = areas_cm2[i3]
                    for i4 in range(i3 + 1):
                        area4, d4_mm = areas_cm2[i4], diams_mm[i4]

                        # Condition 5: no two bar diameters may differ by more
                        # than max_diameter_diff. The four are in descending
                        # order here, so the widest pair is (d_b1, d_b4).
                        if d1_mm - d4_mm > max_diam_diff_mm:
                            continue

                        d_b1 = valid_rebar_diameters[i1]
                        d_b2 = valid_rebar_diameters[i2]
                        d_b3 = valid_rebar_diameters[i3]
                        d_b4 = valid_rebar_diameters[i4]

                        # Iterate over possible numbers of bars in each group
                        for n2 in range(0, max_bars + 1):  # n2 can be 0 or more
                            if n1 + n2 > max_bars:
                                continue  # Skip if the total bars in layer 1 exceed the limit

                            clear_mm = self._layer_clear_spacing_mm(n1, n2, d1_mm, d2_mm, eff_width_mm)
                            if clear_mm is None:
                                continue
                            self._clear_spacing = clear_mm * mm

                            # Calculate area for layer 1
                            A_s_layer_1 = n1 * area1 + (n2 * area2 if n2 > 0 else 0.0)

                            if A_s_layer_1 > max_limit_cm2:
                                break  # further n2 will only increase area

                            # Check if total area from layer 1 is enough for required A_s
                            # And also less than the maximum limit
                            if A_s_layer_1 >= A_s_req_cm2 and A_s_layer_1 <= max_limit_cm2:
                                # Only consider layer 1 — no bars in layer 2
                                valid_combinations.append(
                                    self._long_combo(n1, d_b1, n2, d_b2, 0, None, 0, None, A_s_layer_1, clear_mm)
                                )
                            else:
                                # Track the combination with the maximum possible area (fallback)
                                if A_s_layer_1 > max_fallback_cm2 and A_s_layer_1 <= max_limit_cm2:
                                    max_fallback_cm2 = A_s_layer_1
                                    best_fallback_combination = self._long_combo(
                                        n1, d_b1, n2, d_b2, 0, None, 0, None, A_s_layer_1, clear_mm
                                    )

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
                                    A_s_layer_2 = n3 * area3 + (n4 * area4 if n4 > 0 else 0.0)

                                    # --- Compute total reinforcement and evaluate -----------
                                    total_as = A_s_layer_1 + A_s_layer_2
                                    if total_as >= A_s_req_cm2 and total_as <= max_limit_cm2:
                                        valid_combinations.append(
                                            self._long_combo(n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4, total_as, clear_mm)
                                        )
                                    else:
                                        # Track fallback combination with maximum As
                                        if total_as > max_fallback_cm2 and total_as <= max_limit_cm2:
                                            max_fallback_cm2 = total_as
                                            best_fallback_combination = self._long_combo(
                                                n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4, total_as, clear_mm
                                            )

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
                                        A_s_layer_2 = n3 * area3 + (n4 * area4 if n4 > 0 else 0.0)

                                        total_as = A_s_layer_1 + A_s_layer_2
                                        if total_as >= A_s_req_cm2 and total_as <= max_limit_cm2:
                                            valid_combinations.append(
                                                self._long_combo(
                                                    n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4, total_as, clear_mm
                                                )
                                            )
                                        else:
                                            if total_as > max_fallback_cm2 and total_as <= max_limit_cm2:
                                                max_fallback_cm2 = total_as
                                                best_fallback_combination = self._long_combo(
                                                    n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4, total_as, clear_mm
                                                )

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

    def longitudinal_rebar_EN_1992_2004(
        self,
        A_s_req: Quantity,
        A_s_max: Quantity | None = None,
        mech_cover: Quantity | None = None,
    ) -> None:
        # The bar-selection strategy (fit the area with the fewest, most uniform
        # bars that still respect spacing and layer limits) is geometry, not code
        # provisions, so EN reuses the ACI selector. Only the areas fed into it
        # come from EN 1992-2004.
        self.longitudinal_rebar_ACI_318_19(A_s_req, A_s_max, mech_cover)

    def _layer_clear_spacing_mm(
        self,
        n1: int,
        n2: int,
        d1_mm: float,
        d2_mm: float,
        eff_width_mm: float,
    ) -> float | None:
        """
        Clear spacing between the bars of one layer, in mm.

        Everything is a plain float: this runs inside the innermost loop of the
        longitudinal rebar search, where pint arithmetic dominated the cost.

        Parameters:
            n1 (int): Number of bars in the first group of the layer.
            n2 (int): Number of bars in the second group of the layer.
            d1_mm (float): Diameter of the first group of bars, in mm.
            d2_mm (float): Diameter of the second group of bars, in mm.
            eff_width_mm (float): Width available for bar placement, in mm.

        Returns:
            float | None: The clear spacing in mm, or None when it falls below
            the design limits (clear spacing, vibrator size, largest diameter).
        """
        clear_mm = (eff_width_mm - (n1 * d1_mm + n2 * d2_mm)) / (n1 + n2 - 1)

        # Determine the maximum clear spacing limit; the first two terms are
        # precomputed in __init__.
        max_clear_spacing_mm = max(self._clear_limit_mm, self._vibrator_mm, d1_mm, d2_mm)

        return clear_mm if clear_mm >= max_clear_spacing_mm else None

    def _long_combo(
        self,
        n1: int,
        d_b1: Quantity,
        n2: int,
        d_b2: Quantity | None,
        n3: int,
        d_b3: Quantity | None,
        n4: int,
        d_b4: Quantity | None,
        total_as_cm2: float,
        clear_mm: float,
    ) -> Dict[str, Any]:
        """
        Builds one row of the longitudinal rebar combination table.

        This is where the search loop's floats become Quantities again: it runs
        once per candidate kept, not once per candidate evaluated. A group's
        diameter is reported as None when the group holds no bars.
        """
        return {
            "n_1": n1,
            "d_b1": d_b1,
            "n_2": n2,
            "d_b2": d_b2 if n2 > 0 else None,
            "n_3": n3,
            "d_b3": d_b3 if n3 > 0 else None,
            "n_4": n4,
            "d_b4": d_b4 if n4 > 0 else None,
            "total_as": total_as_cm2 * _CM2,
            "total_bars": n1 + n2 + n3 + n4,
            "clear_spacing": clear_mm * mm,
        }

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

        # These penalties were computed with df.apply(axis=1), which builds a
        # Series for every row; on the few hundred candidates a design produces
        # that machinery, not the arithmetic, was the cost. The columns are
        # pulled out as plain lists once and scored in a single pass.
        groups = list(
            zip(
                df["n_1"].tolist(),
                df["d_b1"].tolist(),
                df["n_2"].tolist(),
                df["d_b2"].tolist(),
                df["n_3"].tolist(),
                df["d_b3"].tolist(),
                df["n_4"].tolist(),
                df["d_b4"].tolist(),
            )
        )

        diameter_penalties = []
        layer_penalties = []
        max_d_per_row = []
        min_d = None
        for n_1, d_1, n_2, d_2, n_3, d_3, n_4, d_4 in groups:
            # One entry per bar, so the spread below is weighted by bar count.
            diameters = []
            if n_1 > 0:
                diameters.extend([d_1.magnitude] * n_1)
            if n_2 > 0:
                diameters.extend([d_2.magnitude] * n_2)
            if n_3 > 0:
                diameters.extend([d_3.magnitude] * n_3)
            if n_4 > 0:
                diameters.extend([d_4.magnitude] * n_4)

            # Penalty for variation in bar diameters
            diameter_penalties.append(np.std(diameters) if diameters else 0.0)
            # Penalty for using a second layer of reinforcement
            layer_penalties.append(1 if (n_3 > 0 or n_4 > 0) else 0)

            # Every row reaching this point has n_1 == 2, so `diameters` is
            # never empty (the A_s_req == 0 case returns before scoring).
            max_d_per_row.append(max(diameters))
            row_min = min(diameters)
            min_d = row_min if min_d is None else min(min_d, row_min)

        # Calculate penalties and add them as columns
        min_as_mag = min_as.magnitude
        df["area_penalty"] = [alpha * q.magnitude / min_as_mag for q in df["total_as"]]
        # Prefer moderate bar counts, where very high or very low will be penalized
        df["bars_penalty"] = beta * ((df["total_bars"] - min_bars) / min_bars) ** 2
        df["diameter_penalty"] = [gamma * p for p in diameter_penalties]
        df["layer_penalty"] = [delta * p for p in layer_penalties]

        # Diameter size penalty, for very large or very small
        df["diameter_size_penalty"] = [epsilon * (d / min_d - 1) for d in max_d_per_row]

        # Slab penalty for large spacing
        if getattr(self, "mode", "beam") == "slab":
            max_spacing_allowed = 300  # mm
            df["spacing_penalty"] = [
                (
                    0.0
                    if s.magnitude <= max_spacing_allowed
                    else (s.magnitude - max_spacing_allowed) / max_spacing_allowed
                )
                for s in df["clear_spacing"]
            ]
        else:
            df["spacing_penalty"] = 0

        # Calculate the final functional
        df["functional"] = (
            df["area_penalty"]
            + df["bars_penalty"]
            + df["diameter_penalty"]
            + df["layer_penalty"]
            + df["diameter_size_penalty"]
            + df["spacing_penalty"]
        )

        return df

    def longitudinal_rebar(
        self,
        A_s_req: Quantity,
        A_s_max: Quantity | None = None,
        mech_cover: Quantity | None = None,
    ) -> Dict[str, Any]:
        """
        Selects the appropriate longitudinal rebar method based on the design
        code.

        Args:
            A_s_req: Required longitudinal rebar area.
            A_s_max: Optional maximum allowable longitudinal rebar area.
            mech_cover: Optional mechanical cover to the bar centroid, used as
                the starting geometry for the layer layout.
        """
        if self.beam.concrete.design_code == "ACI 318-19" or self.beam.concrete.design_code == "CIRSOC 201-25":
            return self.longitudinal_rebar_ACI_318_19(A_s_req, A_s_max, mech_cover)
        elif self.beam.concrete.design_code == "EN 1992-2004":
            return self.longitudinal_rebar_EN_1992_2004(A_s_req, A_s_max, mech_cover)
