from __future__ import annotations
from typing import Any, Dict, TYPE_CHECKING, Tuple
import math
import pandas as pd

from mento.units import psi, mm, cm

if TYPE_CHECKING:
    from mento.beam import RectangularBeam
    from pint.facets.plain import PlainQuantity
    from pandas import DataFrame

class Rebar:
    def __init__(self, beam: RectangularBeam):
        """
        Initializes the Rebar object with the associated beam and settings.
        """
        self.beam = beam
        self.min_clear_spacing = self.beam.settings.get_setting('clear_spacing')
        self.vibrator_size = self.beam.settings.get_setting('vibrator_size')
        self.max_diameter_diff = self.beam.settings.get_setting('max_diameter_diff')
        self.max_bars_per_layer = self.beam.settings.get_setting('max_bars_per_layer')
        self.min_long_rebar = self.beam.settings.get_setting('minimum_longitudinal_diameter')
        self.rebar_diameters = [6*mm, 8*mm, 10*mm, 12*mm, 16*mm, 20*mm, 25*mm, 32*mm]
        self.rebar_areas = {d: (math.pi * d ** 2) / 4 for d in self.rebar_diameters}
        self._long_combos_df: DataFrame = None
        self._trans_combos_df: DataFrame = None
        self._clear_spacing: PlainQuantity = 0*mm

    @property
    def longitudinal_rebar_design(self) -> DataFrame:
        # If combinations_df has not been calculated yet, calculate it
        if self._long_combos_df is None:
            raise ValueError("First run a longitudinal rebar design method, like "
                              "\"beam_longitudinal_rebar_ACI_318_19\".")
        elif self._long_combos_df.empty:
            raise ValueError("No valid combinations found.")
        return self._long_combos_df.iloc[0]
    
    @property
    def transverse_rebar_design(self) -> DataFrame:
        # If combinations_df has not been calculated yet, calculate it
        if self._trans_combos_df is None:
            raise ValueError("First run a transverse rebar design method, like "
                              "\"beam_transverse_rebar_ACI_318_19\".")
        elif self._trans_combos_df.empty:
            raise ValueError("No valid combinations found.")
        return self._trans_combos_df.iloc[0]

    
    def longitudinal_rebar_ACI_318_19(self, A_s_req: PlainQuantity) -> DataFrame:
        """
        Computes the required longitudinal reinforcement based on ACI 318-19.
        
        Args:
            A_s_req: Required longitudinal rebar area.

        Returns:
            A DataFrame containing the best combinations of rebar details.
        """
        self.A_s_req = A_s_req
        effective_width = self.beam.width - 2 * (self.beam.c_c + self.beam._stirrup_d_b)

        # Variables to track the combinations
        valid_combinations = []
        # Create a list of rebar diameters that are equal to or greater than the minimum diameter
        valid_rebar_diameters = [d for d in self.rebar_diameters if d >= self.min_long_rebar]

        for d_b1 in valid_rebar_diameters: # Without taking Ø6 as a possible solution
                for d_b2 in [d for d in valid_rebar_diameters if d <= d_b1]:
                    for d_b3 in [d for d in valid_rebar_diameters if d <= d_b2]:
                        for d_b4 in [d for d in valid_rebar_diameters if d <= d_b3]:

                            # Condition 5: |d_b1 - d_b2| and |d_b3 - d_b4| must not exceed max_diameter_diff
                            # Ensure all diameter combinations satisfy the ordering constraint
                            if not (d_b1 >= d_b2 >= d_b3 >= d_b4):
                                continue  # Skip combinations that do not meet the ordering condition
                            # Apply diameter difference condition across all combinations
                            # Ensure that no two bars exceed `max_diameter_diff`
                            diameters = [d_b1, d_b2, d_b3, d_b4]
                            diameters = [d for d in diameters if d is not None]  #Filter out None values for unused bars
                            
                            # Apply the diameter difference check across all bars
                            if not self._check_diameter_differences(diameters):
                                continue  # Skip this combination if any diameter pair exceeds max_diameter_diff`

                            n1 = 2 # This is a fixed value for every beam 
                            
                            # Iterate over possible numbers of bars in each group
                            for n2 in range(0, self.max_bars_per_layer + 1):  # n2 can be 0 or more
                                # Check spacing for the first set of bars 
                                self._check_spacing(n1, n2, d_b1, d_b2, effective_width)

                                if n1 + n2 > self.max_bars_per_layer:
                                    continue  # Skip if the total bars in layer 1 exceed the limit

                                # Calculate area for layer 1
                                A_s_layer_1 = (
                                    n1 * self.rebar_areas[d_b1] +
                                    (n2 * self.rebar_areas[d_b2] if n2 > 0 else 0*cm**2)
                                )

                                # Condition 6 and 7: Check clear spacing in layer 1
                                if n2 > 0 and not self._check_spacing(n1, n2, d_b1, d_b2, effective_width):
                                    continue

                                A_s_max = max(1.25*A_s_req, n1 * self.rebar_areas[self.min_long_rebar])
                                # Check if total area from layer 1 is enough for required A_s
                                # And also less than 25% greater than A_s_req
                                if A_s_layer_1 >= A_s_req and A_s_layer_1 <= A_s_max:
                                    total_as = A_s_layer_1  # Only consider layer 1
                                    total_bars = n1 + n2    # Total bars only in layer 1
                                    valid_combinations.append({
                                        'n_1': n1,
                                        'd_b1': d_b1,
                                        'n_2': n2,
                                        'd_b2': d_b2 if n2 > 0 else None,  # Display as None if n2 is 0
                                        'n_3': 0,  # No bars in layer 2
                                        'd_b3': None,
                                        'n_4': 0,  # No bars in layer 2
                                        'd_b4': None,
                                        'total_as': total_as.to('cm**2'),
                                        'total_bars': total_bars,
                                        'clear_spacing': self._clear_spacing.to('mm')
                                    })

                                # Now check combinations where bars are added in layer 2 (n3 and n4)
                                for n3 in [0, 2]:  # n3 can be 0 or fixed at 2 if present
                                    for n4 in range(0, self.max_bars_per_layer + 1):
                                        # Ensure layer 2 bars are not more than layer 1 bars
                                        if n3 + n4 > n1 + n2:
                                            continue  # Skip if layer 2 bars exceed layer 1 bars
                                        if n3 == 0 and n4 > 0:
                                            continue  # If n3 is 0, n4 must also be 0
                                        if n3 + n4 > self.max_bars_per_layer:
                                            continue  # Skip if the total bars in layer 2 exceed the limit

                                        # Layer 2 area calculation handling zero values
                                        A_s_layer_2 = (
                                            n3 * self.rebar_areas[d_b3] +
                                            (n4 * self.rebar_areas[d_b4] if n4 > 0 else 0 * cm ** 2)
                                        )
                                        # Condition 4: Area of layer 1 must be >= area of layer 2
                                        if A_s_layer_1 < A_s_layer_2:
                                            continue

                                        # Condition 6 and 7: Check clear spacing in layer 2
                                        if n4 > 0 and not self._check_spacing(n3, n4, d_b3, d_b4, effective_width):
                                            continue

                                        # Check if total area is enough for required A_s
                                        total_as = A_s_layer_1 + A_s_layer_2
                                        if total_as >= A_s_req and total_as <= A_s_max:
                                            total_bars = n1 + n2 + n3 + n4  # Count the total number of bars
                                            valid_combinations.append({
                                                'n_1': n1,
                                                'd_b1': d_b1,
                                                'n_2': n2,
                                                'd_b2': d_b2 if n2 > 0 else None,
                                                'n_3': n3,
                                                'd_b3': d_b3 if n3 > 0 else None,
                                                'n_4': n4,
                                                'd_b4': d_b4 if n4 > 0 else None,
                                                'total_as': total_as.to('cm**2'),
                                                'total_bars': total_bars,
                                                'clear_spacing': self._clear_spacing.to('mm')
                                            })




        # If no valid combination is found, raise an error
        if not valid_combinations:
            raise ValueError("Cannot fit the required reinforcement within "
                             "the beam width considering clear cover and spacing.")

        # Convert valid combinations to DataFrame
        df = pd.DataFrame(valid_combinations)
        # Drop duplicate rows based on the specified columns
        df = df.drop_duplicates(subset=['n_1', 'd_b1', 'n_2', 'd_b2', 'n_3', 'd_b3', 'n_4', 'd_b4'])

        # Sort by 'total_as' first, then by 'total_bars' to prioritize fewer bars
        df.sort_values(by=['total_bars', 'total_as'], inplace=True)
        df.reset_index(drop=True, inplace=True)
        self._long_combos_df = df

        return df.head(7)
        
        # return df
    def _check_spacing(self, n1:int, n2:int, d_b1: PlainQuantity, d_b2:PlainQuantity,
                       effective_width: PlainQuantity) -> bool:
        """
        Checks the clear spacing between rebars in a layer.

        Parameters:
            n1 (int): Number of bars in the first group of the layer.
            n2 (int): Number of bars in the second group of the layer.
            d_b1 (PlainQuantity): Diameter of bars in the first group of the layer.
            d_b2 (PlainQuantity): Diameter of bars in the second group of the layer.
            effective_width (PlainQuantity): The effective width available for bar placement.

        Returns:
            bool: True if the clear spacing satisfies the design limits, False otherwise.
        """
        def layer_clear_spacing(n_a: int, d_a: PlainQuantity, n_b: int, d_b: PlainQuantity) -> PlainQuantity:
            """
            Helper function to calculate clear spacing for a given layer.

            Parameters:
                n_a (int): Number of bars in the first group of the layer.
                d_a (PlainQuantity): Diameter of bars in the first group of the layer.
                n_b (int): Number of bars in the second group of the layer.
                d_b (PlainQuantity): Diameter of bars in the second group of the layer.

            Returns:
                PlainQuantity: Clear spacing for the given layer.
            """
            total_bars = n_a + n_b
            if total_bars <= 1:
                return effective_width - max(d_a, d_b)  # Clear space for one bar
            total_bar_width = n_a * d_a + n_b * d_b
            return (effective_width - total_bar_width) / (total_bars - 1)

        # Calculate the clear spacing for the layer
        self._clear_spacing = layer_clear_spacing(n1, d_b1, n2, d_b2)

        # Determine the maximum clear spacing limit
        max_clear_spacing = max(self.min_clear_spacing, self.vibrator_size, max(d_b1, d_b2))

        # Check if the clear spacing is within limits
        return self._clear_spacing >= max_clear_spacing

    def _check_diameter_differences(self, diameters: list) -> bool:
        """
        Checks that all diameter differences between bars in the list
        do not exceed max_diameter_diff.
        """
        for i, d1 in enumerate(diameters):
            for d2 in diameters[i + 1:]:
                if abs(d1 - d2) > self.max_diameter_diff:
                    return False
        return True


    def longitudinal_rebar_EN_1992(self, A_s_req: PlainQuantity) -> None:
        pass


    def transverse_rebar_ACI_318_19(self, A_v_req: PlainQuantity, 
                                         V_s_req: PlainQuantity) -> DataFrame:
        """
        Computes the required transverse reinforcement based on ACI 318-19.

        Args:
            A_v_req: Required area for transverse reinforcement.
            V_s_req: Shear demand.

        Returns:
            A dictionary containing the transverse rebar design.
        """

        A_cv = self.beam.width*self.beam.d

        s_max_l, s_max_w = self.calculate_max_spacing_ACI(V_s_req, A_cv)

        # Prepare the list for valid combinations
        valid_combinations = []     

        # Iterate through available diameters from 10 mm to 16 mm
        for d_b in self.rebar_diameters[2:5]:
            # Start with 2 legs = 1 stirrup
            n_legs = 2

            # Start with maximum allowed spacing s_max_l
            s_l: PlainQuantity = math.floor(s_max_l.to('cm').magnitude)*cm

            while True:
                # Calculate spacing based on current legs
                n_stirrups = math.ceil(n_legs / 2)  # Number of stirrups based on number of legs
                n_legs_actual = n_stirrups * 2      # Ensure legs are even
                # Consider 1 leg less for spacing laong width
                s_w = (self.beam.width - 2 * self.beam.c_c - self.beam._stirrup_d_b) / (n_legs_actual - 1)  # noqa: F841

                A_db = self.rebar_areas[d_b]  # Area of a stirrup bar
                A_vs = n_legs_actual * A_db  # Area of vertical stirrups
                A_v: PlainQuantity = A_vs / s_l  # Area of vertical stirrups per unit length

                # Check if the calculated A_v meets or exceeds the required A_v
                if A_v >= A_v_req:
                    # Store the valid combination if spacing is also valid
                    valid_combinations.append({
                        'n_stir': int(n_stirrups),
                        'd_b': d_b,
                        's_l': s_l.to('cm'),  # spacing along length
                        's_w': s_w.to('cm'), # spacing along width
                        'A_v': A_v.to('cm**2/m'),
                        's_max_l': s_max_l.to('cm'),
                        's_max_w': s_max_w.to('cm')
                    })
                    # Stop checking larger diameters
                    break  
                    
                # If A_v is insufficient, reduce s_l by 1 cm
                s_l -= 1 * cm
                # If s_l becomes less than 5 cm, increase the number of legs by 2 and reset s_l to s_max_l
                if s_l < 5 * cm:  # If spacing is less than 5 cm, increase 1 stirrup
                    n_legs += 2
                    s_l = math.floor(s_max_l.to('cm').magnitude)*cm # Reset s_l to the max allowed spacing
                # Break the loop if legs exceed the limit of legs or max stirrup diameter
                if n_legs > 6:
                    break

        # Create a DataFrame with all valid combinations
        df_combinations = pd.DataFrame(valid_combinations)

        # Sort combinations by the total rebar area required (ascending)
        df_combinations.sort_values(by='A_v', inplace=True)
        df_combinations.reset_index(drop=True, inplace=True)
        self._trans_combos_df = df_combinations
        return df_combinations
    
    def transverse_rebar_EN_1992_2004(self, A_v_req: PlainQuantity, V_s_req: PlainQuantity) -> None:
        pass
    
    # Factory method to select the transverse rebar method
    def longitudinal_rebar(self, A_s_req: PlainQuantity) -> Dict[str, Any]:
        """
        Selects the appropriate longitudinal rebar method based on the design code.
        """
        if self.beam.concrete.design_code=="ACI 318-19":
            return self.longitudinal_rebar_ACI_318_19(A_s_req)
        # elif self.beam.concrete.design_code=="EN 1992":
        #     return self.beam_longitudinal_rebar_EN_1992(A_s_req)
        # elif self.beam.concrete.design_code=="EHE-08":
        #     return self.beam_longitudinal_rebar_EHE_08(A_s_req)
        else:
            raise ValueError(f"Longitudinal design method not implemented \
                             for concrete type: {type(self.beam.concrete).__name__}")
    
    # Factory method to select the longitudinal rebar method
    def transverse_rebar(self, A_v_req: PlainQuantity, V_s_req: PlainQuantity) -> DataFrame:
        """
        Selects the appropriate transverse rebar method based on the design code.
        """
        if self.beam.concrete.design_code=="ACI 318-19":
            return self.transverse_rebar_ACI_318_19(A_v_req, V_s_req)
        # elif self.beam.concrete.design_code=="EN 1992":
        #     return self.beam_transverse_rebar_EN_1992(A_v_req, V_s_req)
        # elif self.beam.concrete.design_code=="EHE-08":
        #     return self.beam_transverse_rebar_EHE_08(A_v_req, V_s_req)
        else:
            raise ValueError(f"Shear design method not implemented \
                             for concrete type: {type(self.beam.concrete).__name__}")
        
    def calculate_max_spacing_ACI(self, V_s_req: PlainQuantity,
                                   A_cv: PlainQuantity) -> Tuple[PlainQuantity, PlainQuantity]:
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

        f_c=self.beam.concrete.f_c
        lambda_factor = self.beam.settings.get_setting('lambda')
        
        # Determine maximum spacing based on V_s_req condition
        if V_s_req <= 4 * lambda_factor * math.sqrt(f_c / psi) * psi * A_cv:
            # Maximum spacing for lower shear demand
            s_max_l = min(self.beam.d / 2, 60 * cm)
            s_max_w = min(self.beam.d, 60 * cm)
        else:
            # Maximum spacing for higher shear demand
            s_max_l = min(self.beam.d / 4, 30 * cm)
            s_max_w = min(self.beam.d / 2, 30 * cm)
        
        return s_max_l, s_max_w  
    
    def calculate_max_spacing_EN_1992_2004(self, alpha: float) -> Tuple[PlainQuantity, PlainQuantity]:
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
        s_max_l = 0.75*self.beam.d*(1+1 / math.tan(alpha))
        s_max_w = min(0.75*self.beam.d, 60 * cm)
           
        return s_max_l, s_max_w
