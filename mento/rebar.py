from __future__ import annotations
import math
from mento.units import psi, mm, inch, m, cm
from typing import Any, Dict, TYPE_CHECKING

if TYPE_CHECKING:
    from beam import Beam  # Import only for type checking to avoid circular import issues

class Rebar:
    def __init__(self, beam: Beam):
        """
        Initializes the Rebar object with the associated beam and settings.
        """
        self.beam = beam
        self.settings = beam.get_settings()
        self.clear_spacing = self.settings.get('clear_spacing')
        self.rebar_diameters = [6*mm, 8*mm, 10*mm, 12*mm, 16*mm, 20*mm, 25*mm, 32*mm]
        self.rebar_areas = {d: (math.pi * d ** 2) / 4 for d in self.rebar_diameters}

    def beam_longitudinal_rebar_ACI_318_19(self, A_s_req: float) -> Dict[str, Any]:
        """
        Computes the required longitudinal reinforcement based on ACI 318-19.
        
        Args:
            A_s_req: Required steel area in square inches.

        Returns:
            A dictionary containing the best combination of rebar details.
        """
        self.A_s_req = A_s_req
        effective_width = self.beam._width - 2 * (self.beam.cc + self.beam.stirrup_d_b)

        # Variables to track the best combination
        best_combination = None
        min_total_area = float('inf')
        layer = 1 # Assuming a single layer for simplicity
        
        # Try all combinations of rebar sizes and numbers
        for diameter in self.rebar_diameters:
            rebar_area = self.rebar_areas[diameter]  # Convert diameter to rebar area
            num_bars = 2  # Start with 2 bars minimum
            total_as = 0 * cm**2
            
            while num_bars * diameter + (num_bars - 1) * self.clear_spacing <= effective_width:
                total_as = num_bars * rebar_area
                
                if total_as >= A_s_req:
                    if total_as < min_total_area:
                        min_total_area = total_as
                        available_spacing = (effective_width - (num_bars * diameter)) / (num_bars - 1)
                        best_combination = {
                            'layer_1': layer,
                            'num_bars_1': num_bars,
                            'diameter_1': diameter,
                            'total_as': total_as,
                            'available_spacing_1': available_spacing
                        }
                    break  # Exit the loop when the required area is satisfied
                
                num_bars += 1

        # If no valid combination is found, raise a ValueError
        if best_combination is None:
            raise ValueError("Cannot fit the required reinforcement within "
                             "the beam width considering clear cover and spacing.")
        
                
        return best_combination
    
    def beam_longitudinal_rebar_EN_1992(self, A_s_req: float) -> None:
        pass

    def beam_longitudinal_rebar_EHE_08(self, A_s_req: float) -> None:
        pass

    def beam_transverse_rebar_ACI_318_19(self, A_v_req: float, V_s_req: float) -> Dict[str, Any]:
        """
        Computes the required transverse reinforcement based on ACI 318-19.

        Args:
            A_v_req: Required area for transverse reinforcement.
            V_s_req: Shear demand.

        Returns:
            A dictionary containing the transverse rebar design.
        """
        concrete_properties=self.beam.concrete.get_properties()
        f_c=concrete_properties["f_c"]
        lambda_factor = self.beam._settings.get_setting('lambda')
        A_cv = self.beam._width*self.beam.d

        # Check if V_s_req <= 4 * lambda * sqrt(f_c) * A_cv
        if V_s_req <= 4 * lambda_factor * math.sqrt(f_c / psi) * psi * A_cv:
            # Maximum spacing across the length of the beam
            s_max_l = min(self.beam.d / 2, 24*inch) 
            # Maximum spacing across the width of the beam
            s_max_w = min(self.beam.d, 24*inch) 
            # Spacing along length
            s = math.floor(self.beam.d / 2)*inch
        else:
            # Maximum spacing across the length of the beam
            s_max_l = min(self.beam.d / 4, 12*inch) 
            # Maximum spacing across the width of the beam
            s_max_w = min(self.beam.d / 2, 12*inch) 
            # Spacing along length
            s = math.floor(self.beam.d / 4) *inch

        # Ensure that the calculated spacing is within the maximum allowed spacing
        s = min(s, s_max_l)
        # Required transverse rebar per stirrup distance
        A_vs_req = A_v_req*s
        # Number of legs of the stirrups across the beam width
        n_legs_req = math.ceil((self.beam._width - 2 * self.beam.cc - self.beam.stirrup_d_b) / s_max_w) + 1
        n_stirrups = math.ceil(n_legs_req/2)
        n_legs = n_stirrups*2
        # Spacing along the width of the beam
        s_w = (self.beam._width - 2 * self.beam.cc - self.beam.stirrup_d_b) / (n_legs - 1)
        # Ensure that spacing along the width is within the maximum allowed spacing
        s_w = min(s_w, s_max_w)
        # Minimum bar diameter (in inches)
        d_b_min = math.sqrt((4 * A_vs_req.value) / (math.pi * n_legs))*m

        # Find the smallest available bar diameter greater than or equal to d_bmin
        d_b = max(3/8*inch, min(filter(lambda db: db >= d_b_min, self.rebar_diameters)))
        # Area of a stirrup bar 
        A_db = self.rebar_areas[d_b]  # Convert diameter to rebar area
        # Vertical stirrups angle
        # alpha = 90  # degrees
        # Area of vertical stirrups
        A_vs = n_legs * A_db
        # Area of vertical stirrups per unit length (in^2/ft)
        A_v = A_vs / s
        return {
                'n_stirrups': n_stirrups,
                'd_b': d_b,
                'spacing': s,
                'A_v': A_v,
            }
    
    def beam_transverse_rebar_EN_1992(self, A_v_req: float, V_s_req: float) -> None:
        pass

    def beam_transverse_rebar_EHE_08(self, A_v_req: float, V_s_req: float) -> None:
        pass
    
    # Factory method to select the transverse rebar method
    def beam_longitudinal_rebar(self, A_s_req: float) -> Dict[str, Any]:
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
    def beam_transverse_rebar(self, A_v_req: float, V_s_req: float) -> Dict[str, Any]:
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