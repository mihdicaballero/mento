from __future__ import annotations
import math
from mento.units import psi, mm, inch, m, cm
from typing import Any, Dict, TYPE_CHECKING
import pandas as pd

if TYPE_CHECKING:
    from mento.concrete.beam import RectangularConcreteBeam
    from pint.facets.plain import PlainQuantity
    from pandas import DataFrame

class Rebar:
    def __init__(self, beam: RectangularConcreteBeam):
        """
        Initializes the Rebar object with the associated beam and settings.
        """
        self.beam = beam
        self.clear_spacing = self.beam.settings.get_setting('clear_spacing')
        self.rebar_diameters = [6*mm, 8*mm, 10*mm, 12*mm, 16*mm, 20*mm, 25*mm, 32*mm]
        self.rebar_areas = {d: (math.pi * d ** 2) / 4 for d in self.rebar_diameters}
        self._beam_long_combos_df: DataFrame = None
        self._beam_trans_combos_df: DataFrame = None

    @property
    def beam_longitudinal_rebar_design(self) -> DataFrame:
        # If combinations_df has not been calculated yet, calculate it
        if self._beam_long_combos_df is None:
            raise ValueError("First run a longitudinal rebar design method, like "
                              "\"beam_longitudinal_rebar_ACI_318_19\".")
        elif self._beam_long_combos_df.empty:
            raise ValueError("No valid combinations found.")
        return self._beam_long_combos_df.iloc[0]
    
    def beam_longitudinal_rebar_ACI_318_19(self, A_s_req: PlainQuantity) -> DataFrame:
        """
        Computes the required longitudinal reinforcement based on ACI 318-19.
        
        Args:
            A_s_req: Required longitudinal rebar area.

        Returns:
            A dictionary containing the best combination of rebar details.
        """
        self.A_s_req = A_s_req
        effectivewidth = self.beam.width - 2 * (self.beam.cc + self.beam.stirrup_d_b)

        # Variables to track the combinations
        valid_combinations = []
        layer = 1 # Assuming a single layer for simplicity
        
        # Try all combinations of rebar sizes and numbers
        for diameter in self.rebar_diameters:
            rebar_area = self.rebar_areas[diameter]  # Convert diameter to rebar area
            num_bars = 2  # Start with 2 bars minimum
            total_as = 0 * cm**2
            
            while num_bars * diameter + (num_bars - 1) * self.clear_spacing <= effectivewidth:
                total_as = (num_bars * rebar_area)

                if total_as >= A_s_req:
                    available_spacing: PlainQuantity = (effectivewidth - (num_bars * diameter)) / (num_bars - 1)
                    valid_combinations.append({
                        'layer_1': layer,
                        'num_bars_1': num_bars,
                        'diameter_1': diameter.to('mm'),
                        'total_as': total_as.to('cm**2'),
                        'available_spacing_1': available_spacing.to('mm')
                    })
                    break  # Exit the loop when the required area is satisfied
                        
                num_bars += 1

        # If no valid combination is found, raise a ValueError
        if not valid_combinations:
            raise ValueError("Cannot fit the required reinforcement within "
                            "the beam width considering clear cover and spacing.")

        # Convert the list of valid combinations to a pandas DataFrame
        df = pd.DataFrame(valid_combinations)

        # Sort the DataFrame by 'total_as' (total area of steel)
        df.sort_values(by='total_as', inplace=True)
        df.reset_index(drop=True, inplace=True)
        self._beam_long_combos_df = df
        
        return df
    
    def beam_longitudinal_rebar_EN_1992(self, A_s_req: PlainQuantity) -> None:
        pass

    def beam_longitudinal_rebar_EHE_08(self, A_s_req: PlainQuantity) -> None:
        pass

    def beam_transverse_rebar_ACI_318_19(self, A_v_req: PlainQuantity, 
                                         V_s_req: PlainQuantity) -> Dict[str, PlainQuantity]:
        """
        Computes the required transverse reinforcement based on ACI 318-19.

        Args:
            A_v_req: Required area for transverse reinforcement.
            V_s_req: Shear demand.

        Returns:
            A dictionary containing the transverse rebar design.
        """
        f_c=self.beam.concrete.f_c
        lambda_factor = self.beam.settings.get_setting('lambda')
        A_cv = self.beam.width*self.beam.d

        # Check if V_s_req <= 4 * lambda * sqrt(f_c) * A_cv
        if V_s_req <= 4 * lambda_factor * math.sqrt(f_c / psi) * psi * A_cv:
            # Maximum spacing across the length of the beam
            s_max_l = min(self.beam.d / 2, 24 * inch)
        else:
            # Maximum spacing across the length of the beam
            s_max_l = min(self.beam.d / 4, 12 * inch)

        # # Check if V_s_req <= 4 * lambda * sqrt(f_c) * A_cv
        # if V_s_req <= 4 * lambda_factor * math.sqrt(f_c / psi) * psi * A_cv:
        #     # Maximum spacing across the length of the beam
        #     s_max_l = min(self.beam.d / 2, 60*cm) 
        #     # Maximum spacing across the width of the beam
        #     s_max_w = min(self.beam.d, 60*cm) 
        #     # Spacing along length
        #     s = math.floor(self.beam.d / (2*cm)) * cm
        # else:
        #     # Maximum spacing across the length of the beam
        #     s_max_l = min(self.beam.d / 4, 30*cm) 
        #     # Maximum spacing across the width of the beam
        #     s_max_w = min(self.beam.d / 2, 30*cm) 
        #     # Spacing along length
        #     s = math.floor(self.beam.d/(4*inch))*inch

        # Prepare the list for valid combinations
        valid_combinations = []     

        # Iterate through available diameters
        for d_b in self.rebar_diameters:
            # Start with 2 legs
            n_legs = 2
            while n_legs <= 4:
                # Calculate spacing based on current legs
                n_stirrups = math.ceil(n_legs / 2)  # Number of stirrups based on number of legs
                n_legs_actual = n_stirrups * 2      # Ensure legs are even
                s_w = (self.beam.width - 2 * self.beam.cc - self.beam.stirrup_d_b) / (n_legs_actual - 1)

                # Ensure spacing along the width is within maximum allowed spacing
                if s_w <= s_max_l:
                    # Calculate the required area of transverse reinforcement
                    A_vs_req = A_v_req * (self.beam.d / n_legs_actual)  # A_v based on spacing along length
                    A_db = self.rebar_areas[d_b]  # Area of a stirrup bar
                    A_vs = n_legs_actual * A_db  # Area of vertical stirrups
                    A_v = A_vs / (self.beam.d / n_legs_actual)  # Area of vertical stirrups per unit length

                    # Store valid combination
                    valid_combinations.append({
                        'n_stirrups': n_stirrups,
                        'n_legs': n_legs_actual,
                        'd_b': d_b,
                        'spacing_length': s_max_l,  # spacing along length
                        'spacing_width': s_w,
                        'A_v': A_v,
                    })

                    # Check if the spacing is less than 4 cm and increment legs if necessary
                    if s_w < 4 * cm:
                        n_legs += 1  # Increase legs for next iteration
                    else:
                        break  # If spacing is acceptable, break out of the loop for this diameter
                else:
                    break  # If spacing exceeds max spacing, break out of the loop for this diameter

        # Create a DataFrame with all valid combinations
        df_combinations = pd.DataFrame(valid_combinations)

        # Sort combinations by the total rebar area required (ascending)
        df_combinations = df_combinations.sort_values(by='A_v')
        self._beam_trans_combos_df = df_combinations
        return df_combinations

        # # Ensure that the calculated spacing is within the maximum allowed spacing
        # s = min(s, s_max_l)
        # # Required transverse rebar per stirrup distance
        # A_vs_req = A_v_req*s
        # # Number of legs of the stirrups across the beam width
        # n_legs_req = math.ceil((self.beam.width - 2 * self.beam.cc - self.beam.stirrup_d_b) / s_max_w) + 1
        # n_stirrups = math.ceil(n_legs_req/2)
        # n_legs = n_stirrups*2
        # # Spacing along the width of the beam
        # s_w = (self.beam.width - 2 * self.beam.cc - self.beam.stirrup_d_b) / (n_legs - 1)
        # # Ensure that spacing along the width is within the maximum allowed spacing
        # s_w = min(s_w, s_max_w)
        # # Minimum bar diameter (in inches)
        # d_b_min = math.sqrt((4 * A_vs_req/(m**2)) / (math.pi * n_legs))*m

        # # Find the smallest available bar diameter greater than or equal to d_bmin
        # d_b = max(3/8*inch, min(filter(lambda db: db >= d_b_min, self.rebar_diameters)))
        # # Area of a stirrup bar 
        # A_db = self.rebar_areas[d_b]  # Convert diameter to rebar area
        # # Vertical stirrups angle
        # # alpha = 90  # degrees
        # # Area of vertical stirrups
        # A_vs = n_legs * A_db
        # # Area of vertical stirrups per unit length (in^2/ft)
        # A_v = A_vs / s
        # return {
        #         'n_stirrups': n_stirrups,
        #         'd_b': d_b,
        #         'spacing': s,
        #         'A_v': A_v,
        #     }
    
    def beam_transverse_rebar_EN_1992(self, A_v_req: PlainQuantity, V_s_req: PlainQuantity) -> None:
        pass

    def beam_transverse_rebar_EHE_08(self, A_v_req: PlainQuantity, V_s_req: PlainQuantity) -> None:
        pass
    
    # Factory method to select the transverse rebar method
    def beam_longitudinal_rebar(self, A_s_req: PlainQuantity) -> Dict[str, Any]:
        """
        Selects the appropriate longitudinal rebar method based on the design code.
        """
        if self.beam.concrete.design_code=="ACI 318-19":
            return self.beam_longitudinal_rebar_ACI_318_19(A_s_req)
        # elif self.beam.concrete.design_code=="EN 1992":
        #     return self.beam_longitudinal_rebar_EN_1992(A_s_req)
        # elif self.beam.concrete.design_code=="EHE-08":
        #     return self.beam_longitudinal_rebar_EHE_08(A_s_req)
        else:
            raise ValueError(f"Longitudinal design method not implemented \
                             for concrete type: {type(self.beam.concrete).__name__}")
    
    # Factory method to select the longitudinal rebar method
    def beam_transverse_rebar(self, A_v_req: PlainQuantity, V_s_req: PlainQuantity) -> Dict[str, PlainQuantity]:
        """
        Selects the appropriate transverse rebar method based on the design code.
        """
        if self.beam.concrete.design_code=="ACI 318-19":
            return self.beam_transverse_rebar_ACI_318_19(A_v_req, V_s_req)
        # elif self.beam.concrete.design_code=="EN 1992":
        #     return self.beam_transverse_rebar_EN_1992(A_v_req, V_s_req)
        # elif self.beam.concrete.design_code=="EHE-08":
        #     return self.beam_transverse_rebar_EHE_08(A_v_req, V_s_req)
        else:
            raise ValueError(f"Shear design method not implemented \
                             for concrete type: {type(self.beam.concrete).__name__}")
        
        