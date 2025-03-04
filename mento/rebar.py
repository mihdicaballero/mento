from __future__ import annotations
from typing import Any, Dict, TYPE_CHECKING, Tuple
import math
import pandas as pd
import numpy as np

from mento.units import psi, mm, cm, inch, MPa

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
        self._long_combos_df: DataFrame = None
        self._trans_combos_df: DataFrame = None
        self._clear_spacing: PlainQuantity = 0*mm
        # Unit system default rebar
        if self.beam.concrete.unit_system == "metric":
            self.rebar_diameters = [6*mm, 8*mm, 10*mm, 12*mm, 16*mm, 20*mm, 25*mm, 32*mm]
            self.rebar_areas = {d: (math.pi * d ** 2) / 4 for d in self.rebar_diameters}
        else:
            self.rebar_diameters = [3*inch/8, 4*inch/8, 5*inch/8, 6*inch/8, 7*inch/8, 8*inch/8, 1.128*inch,
                                    1.27*inch, 1.41*inch, 1.693*inch]
            rebar_areas_list = [0.1104*inch**2, 0.1963*inch**2, 0.3068*inch**2, 0.4418*inch**2, 0.6013*inch**2,
                                0.7854*inch**2, 0.9940*inch**2, 1.27*inch**2, 1.56*inch**2, 2.25*inch**2]
            self.rebar_areas = dict(zip(self.rebar_diameters, rebar_areas_list))

    @property
    def longitudinal_rebar_design(self) -> DataFrame:
        return self._long_combos_df.iloc[0]
    
    @property
    def transverse_rebar_design(self) -> DataFrame:
        return self._trans_combos_df.iloc[0]

    
    def longitudinal_rebar_ACI_318_19(self, A_s_req: PlainQuantity) -> DataFrame:
        """
        Computes the required longitudinal reinforcement based on ACI 318-19.
        
        Args:
            A_s_req: Required longitudinal rebar area.

        Returns:
            A DataFrame containing the best combinations of rebar details.
            If no combination satisfies A_s_req, returns the combination with the maximum possible area.
        """
        self.A_s_req = A_s_req
        effective_width = self.beam.width - 2 * (self.beam.c_c + self.beam._stirrup_d_b)

        # Variables to track the combinations
        valid_combinations = []
        best_fallback_combination = None  # To store the best fallback design
        max_fallback_area = 0 * cm**2  # To track the maximum area in fallback cases
        # Create a list of rebar diameters that are equal to or greater than the minimum diameter
        if self.beam.concrete.unit_system == "metric":
            self.min_long_rebar = 10*mm
        else: 
            self.min_long_rebar = 3*inch/8
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

                                A_s_max = max(10*A_s_req, n1 * self.rebar_areas[self.min_long_rebar])
                                # Check if total area from layer 1 is enough for required A_s
                                # And also less than A_s_max
                                if A_s_layer_1 >= A_s_req and A_s_layer_1 <= A_s_max:
                                # if A_s_layer_1 >= A_s_req:
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
                                else:
                                    # Track the combination with the maximum possible area (fallback)
                                    if A_s_layer_1 > max_fallback_area:
                                        max_fallback_area = A_s_layer_1
                                        best_fallback_combination = {
                                            'n_1': n1,
                                            'd_b1': d_b1,
                                            'n_2': n2,
                                            'd_b2': d_b2 if n2 > 0 else None,
                                            'n_3': 0,
                                            'd_b3': None,
                                            'n_4': 0,
                                            'd_b4': None,
                                            'total_as': A_s_layer_1.to('cm**2'),
                                            'total_bars': n1 + n2,
                                            'clear_spacing': self._clear_spacing.to('mm')
                                        }

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
                                        # if total_as >= A_s_req: 
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
                                        else:
                                            # Track the combination with the maximum possible area (fallback)
                                            if total_as > max_fallback_area:
                                                max_fallback_area = total_as
                                                best_fallback_combination = {
                                                    'n_1': n1,
                                                    'd_b1': d_b1,
                                                    'n_2': n2,
                                                    'd_b2': d_b2 if n2 > 0 else None,
                                                    'n_3': n3,
                                                    'd_b3': d_b3 if n3 > 0 else None,
                                                    'n_4': n4,
                                                    'd_b4': d_b4 if n4 > 0 else None,
                                                    'total_as': total_as.to('cm**2'),
                                                    'total_bars': n1 + n2 + n3 + n4,
                                                    'clear_spacing': self._clear_spacing.to('mm')
                                                }

        # Convert valid combinations to DataFrame
        df = pd.DataFrame(valid_combinations)
        # Drop duplicate rows based on the specified columns
        df = df.drop_duplicates(subset=['n_1', 'd_b1', 'n_2', 'd_b2', 'n_3', 'd_b3', 'n_4', 'd_b4'])
        # If no valid combinations satisfy A_s_req, use the best fallback combination
        if df.empty and best_fallback_combination is not None:
            df = pd.DataFrame([best_fallback_combination])

        modified_df: pd.DataFrame = self._calculate_penalties_long_rebar(df)
        # Sort by 'Functional' to sort by the best options
        modified_df.sort_values(by=['functional'], inplace=True)
        # modified_df.drop(columns=['area_penalty','bars_penalty', 'diameter_penalty',
        #                           'layer_penalty', 'combined_penalty', 'functional'], inplace=True)
        modified_df.reset_index(drop=True, inplace=True)
        self._long_combos_df = modified_df

        return modified_df.head(10)


    def longitudinal_rebar_EN_1992(self, A_s_req: PlainQuantity) -> None:
        pass
  
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
    
    def _calculate_penalties_long_rebar(self, df: pd.DataFrame,
                                        alpha: float =  3.5,
                                        beta: float =   0.3,
                                        gamma: float =  0.25,
                                        delta: float =  1) -> pd.DataFrame:
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

        # Calculate minimum bars and minimum area of steel
        min_bars = df['total_bars'].min()
        min_as = df['total_as'].min()

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
            if row['n_1'] > 0:
                diameters.extend([row['d_b1'].magnitude] * row['n_1'])
            if row['n_2'] > 0:
                diameters.extend([row['d_b2'].magnitude] * row['n_2'])
            if row['n_3'] > 0:
                diameters.extend([row['d_b3'].magnitude] * row['n_3'])
            if row['n_4'] > 0:
                diameters.extend([row['d_b4'].magnitude] * row['n_4'])
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
            if row['n_3'] > 0 or row['n_4'] > 0:
                return 1  # Penalize if additional layers are used
            return 0  # No penalty if only the first layer is used

        
        # Calculate penalties and add them as columns
        df['area_penalty'] = df.apply(lambda row: (alpha*row['total_as']/min_as), axis=1)
        df['bars_penalty'] = beta * (df['total_bars'] - min_bars)/min_bars
        df['diameter_penalty'] = gamma * df.apply(diameter_difference_penalty, axis=1)
        df['layer_penalty'] = delta * df.apply(layer_penalty, axis=1)

        # Calculate the final functional
        df['functional'] = df.apply(
            lambda row: (
                row['area_penalty'] +
                row['bars_penalty'] +
                row['diameter_penalty'] +
                row['layer_penalty']
            ),
            axis=1
        )

        return df

    def transverse_rebar_ACI_318_19(self, V_s_req: PlainQuantity) -> Any:

        if self.beam.concrete.unit_system == "metric":
            valid_diameters = self.rebar_diameters[2:5] #Minimum 10 mm
        else: 
            valid_diameters = self.rebar_diameters[0:3]
        
        A_cv = self.beam.width*self.beam._d_shear
        s_max_l, s_max_w = self.calculate_max_spacing_ACI_318_19(V_s_req, A_cv)

        return valid_diameters, s_max_l, s_max_w
    
    def transverse_rebar_CIRSOC_201_25(self, V_s_req: PlainQuantity) -> Any:

        valid_diameters = self.rebar_diameters[0:5] #Minimum 6 mm
        
        A_cv = self.beam.width*self.beam._d_shear
        s_max_l, s_max_w = self.calculate_max_spacing_ACI_318_19(V_s_req, A_cv)

        return valid_diameters, s_max_l, s_max_w
    
    
    def transverse_rebar(self, A_v_req: PlainQuantity, V_s_req: PlainQuantity) -> DataFrame:
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
        if self.beam.concrete.design_code=="ACI 318-19":
            valid_diameters, s_max_l, s_max_w = self.transverse_rebar_ACI_318_19(V_s_req)
        elif self.beam.concrete.design_code=="CIRSOC 201-25":
            valid_diameters, s_max_l, s_max_w = self.transverse_rebar_CIRSOC_201_25(V_s_req)
        else:
            raise ValueError(f"Shear design method not implemented \
                             for concrete type: {type(self.beam.concrete).__name__}")

        # Iterate through available diameters
        for d_b in valid_diameters:
            # Start with 2 legs = 1 stirrup
            n_legs = 2

            # Start with maximum allowed spacing s_max_l
            if self.beam.concrete.unit_system == "metric":
                s_l: PlainQuantity = math.floor(s_max_l.to('cm').magnitude)*cm
            else:
                s_l = math.floor(s_max_l.to('inch').magnitude)*inch

            while True:
                # Calculate spacing based on current legs
                n_stirrups = math.ceil(n_legs / 2)  # Number of stirrups based on number of legs
                n_legs_actual = n_stirrups * 2      # Ensure legs are even
                # Consider 1 leg less for spacing laong width
                s_w = (self.beam.width - 2 * self.beam.c_c - self.beam._stirrup_d_b) / (n_legs_actual - 1)  # noqa: F841

                A_db = self.rebar_areas[d_b]  # Area of a stirrup bar
                A_vs = n_legs_actual * A_db  # Area of vertical stirrups
                A_v: PlainQuantity = A_vs / s_l  # Area of vertical stirrups per unit length

                # Store the valid combination if spacing is also valid
                if self.beam.concrete.unit_system == "metric":
                    # Check if the calculated A_v meets or exceeds the required A_v
                    if A_v >= A_v_req:
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
                else:
                    # Check if the calculated A_v meets or exceeds the required A_v
                    if A_v >= A_v_req:
                        valid_combinations.append({
                            'n_stir': int(n_stirrups),
                            'd_b': d_b,
                            's_l': s_l.to('inch'),  # spacing along length
                            's_w': s_w.to('inch'), # spacing along width
                            'A_v': A_v.to('inch**2/ft'),
                            's_max_l': s_max_l.to('inch'),
                            's_max_w': s_max_w.to('inch')
                        })
                        # Stop checking larger diameters
                        break  
                        
                    # If A_v is insufficient, reduce s_l by 1 inch
                    s_l -= 1 * inch
                    # If s_l becomes less than 2 inch, increase the number of legs by 2 and reset s_l to s_max_l
                    if s_l < 2 * inch:  # If spacing is less than 2 inch, increase 1 stirrup
                        n_legs += 2
                        s_l = math.floor(s_max_l.to('inch').magnitude)*inch # Reset s_l to the max allowed spacing
                    # Break the loop if legs exceed the limit of legs or max stirrup diameter
                    if n_legs > 6:
                        break

        # Create a DataFrame with all valid combinations
        df_combinations = pd.DataFrame(valid_combinations)

        # Sort combinations by the total rebar area required (ascending)
        # Sort by 'A_v' first, then by 'n_stir' to prioritize fewer bars
        df_combinations.sort_values(by=['n_stir','A_v'], inplace=True)
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
        if self.beam.concrete.design_code=="ACI 318-19" or self.beam.concrete.design_code=="CIRSOC 201-25":
            return self.longitudinal_rebar_ACI_318_19(A_s_req)
        # elif self.beam.concrete.design_code=="EN 1992":
        #     return self.beam_longitudinal_rebar_EN_1992(A_s_req)
        # elif self.beam.concrete.design_code=="EHE-08":
        #     return self.beam_longitudinal_rebar_EHE_08(A_s_req)
        else:
            raise ValueError(f"Longitudinal design method not implemented \
                             for concrete type: {type(self.beam.concrete).__name__}")
    
    # Factory method to select the longitudinal rebar method
    # def transverse_rebar(self, A_v_req: PlainQuantity, V_s_req: PlainQuantity) -> DataFrame:
    #     """
    #     Selects the appropriate transverse rebar method based on the design code.
    #     """
    #     if self.beam.concrete.design_code=="ACI 318-19":
    #         return self.transverse_rebar_ACI_318_19(A_v_req, V_s_req)
    #     # elif self.beam.concrete.design_code=="EN 1992":
    #     #     return self.beam_transverse_rebar_EN_1992(A_v_req, V_s_req)
    #     # elif self.beam.concrete.design_code=="EHE-08":
    #     #     return self.beam_transverse_rebar_EHE_08(A_v_req, V_s_req)
    #     else:
    #         raise ValueError(f"Shear design method not implemented \
    #                          for concrete type: {type(self.beam.concrete).__name__}")
        
    def calculate_max_spacing_ACI_318_19(self, V_s_req: PlainQuantity,
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
                s_max_l = min(self.beam._d_shear / 2, 24*inch)
                s_max_w = min(self.beam._d_shear, 24*inch)
            else:
                s_max_l = min(self.beam._d_shear / 4, 12*inch)
                s_max_w = min(self.beam._d_shear / 2, 12*inch)
        
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
        s_max_l = 0.75*self.beam._d_shear*(1+1 / math.tan(alpha))
        s_max_w = min(0.75*self.beam._d_shear, 60 * cm)
           
        return s_max_l, s_max_w


