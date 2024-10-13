from __future__ import annotations
from typing import Any, Dict, TYPE_CHECKING
import math
import pandas as pd

from mento.units import psi, mm, cm

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
        self._long_combos_df: DataFrame = None
        self._trans_combos_df: DataFrame = None

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
            A dictionary containing the best combination of rebar details.
        """
        self.A_s_req = A_s_req
        effectivewidth = self.beam.width - 2 * (self.beam.c_c + self.beam._stirrup_d_b)

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
        self._long_combos_df = df
        
        return df
    
    def longitudinal_rebar_EN_1992(self, A_s_req: PlainQuantity) -> None:
        pass

    def longitudinal_rebar_EHE_08(self, A_s_req: PlainQuantity) -> None:
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
        f_c=self.beam.concrete.f_c
        lambda_factor = self.beam.settings.get_setting('lambda')
        A_cv = self.beam.width*self.beam.d

        # Check if V_s_req <= 4 * lambda * sqrt(f_c) * A_cv
        if V_s_req <= 4 * lambda_factor * math.sqrt(f_c / psi) * psi * A_cv:
            # Maximum spacing across the length of the beam
            s_max_l = min(self.beam.d / 2, 60*cm)
            # Maximum spacing across the width of the beam
            s_max_w = min(self.beam.d, 60*cm) 
        else:
            # Maximum spacing across the length of the beam
            s_max_l = min(self.beam.d / 4, 30*cm)
            # Maximum spacing across the width of the beam
            s_max_w = min(self.beam.d / 2, 30*cm)  # noqa: F841

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
                        'n_stir': n_stirrups,
                        'd_b': d_b,
                        's_l': s_l.to('cm'),  # spacing along length
                        'A_v': A_v.to('cm**2/m'),
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
    
    def transverse_rebar_EN_1992(self, A_v_req: PlainQuantity, V_s_req: PlainQuantity) -> None:
        pass

    def transverse_rebar_EHE_08(self, A_v_req: PlainQuantity, V_s_req: PlainQuantity) -> None:
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
        
        