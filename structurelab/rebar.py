import forallpeople
forallpeople.environment('structural',top_level=True)
cm = 1e-2*m # type: ignore
import math
from structurelab import material

class Rebar:
    rebar_diameters = [6*mm, 8*mm, 10*mm, 12*mm, 16*mm, 20*mm, 25*mm, 32*mm] #type: ignore
    rebar_areas = {d: (math.pi * d ** 2) / 4 for d in rebar_diameters}
# Add beam class of rebar
    def __init__(self, beam):
            self.beam = beam
            self.settings = beam.get_settings()
            self.cc = self.settings.get('clear_cover')
            self.clear_spacing = self.settings.get('clear_spacing')
            self.def_stirrup_db = self.settings.get('stirrup_diameter')

    def beam_longitudinal_rebar_ACI_318_19(self, A_s_req:float):
        self.A_s_req = A_s_req
        effective_width = self.beam._width - 2 * (self.cc + self.def_stirrup_db)
        # The method finds the best combination of rebar for the required As
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
                    if best_combination is None or total_as < min_total_area:
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
            raise ValueError("Cannot fit the required reinforcement within the beam width considering clear cover and spacing.")
        
                
        return best_combination
    
    def beam_longitudinal_rebar_EN_1992(self):
        pass
    def beam_longitudinal_rebar_EHE_08(self):
        pass

    def beam_transverse_rebar_ACI_318_19(self, A_v_req:float, V_s_req:float, lambda_factor:float, f_c:float, d:float):
        self.def_stirrup_db = self.settings.get('stirrup_diameter')
        best_combination = None
        min_A_v = float('inf')
        # Check if V_s_req <= 4 * lambda * sqrt(f_c) * A_cv
        A_cv = self.beam._width*d
        if V_s_req <= 4 * lambda_factor * math.sqrt(f_c / psi) * psi * A_cv: #type: ignore
            # Maximum spacing across the length of the beam
            s_max_l = min(d / 2, 24*inch)  #type: ignore
            # Maximum spacing across the width of the beam
            s_max_w = min(d, 24*inch)  #type: ignore
            # Spacing along length
            s = math.floor(d / 2)*inch #type: ignore
        else:
            # Maximum spacing across the length of the beam
            s_max_l = min(d / 4, 12*inch)  #type: ignore
            # Maximum spacing across the width of the beam
            s_max_w = min(d / 2, 12*inch)  #type: ignore
            # Spacing along length
            s = math.floor(d / 4) *inch #type: ignore

        # Ensure that the calculated spacing is within the maximum allowed spacing
        s = min(s, s_max_l)
        # Required transverse rebar per stirrup distance
        A_vs_req = A_v_req*s
        # Number of legs of the stirrups across the beam width
        n_legs_req = math.ceil((self.beam._width - 2 * self.cc - self.def_stirrup_db) / s_max_w) + 1
        n_stirrups = math.ceil(n_legs_req/2)
        n_legs = n_stirrups*2
        # Spacing along the width of the beam
        s_w = (self.beam._width - 2 * self.cc - self.def_stirrup_db) / (n_legs - 1)
        # Ensure that spacing along the width is within the maximum allowed spacing
        s_w = min(s_w, s_max_w)
        # Minimum bar diameter (in inches)
        d_b_min = math.sqrt((4 * A_vs_req.value) / (math.pi * n_legs))*m #type: ignore

        # Find the smallest available bar diameter greater than or equal to d_bmin
        d_b = max(3/8*inch, min(filter(lambda db: db >= d_b_min, self.rebar_diameters))) #type: ignore
        # Area of a stirrup bar 
        A_db = self.rebar_areas[d_b]  # Convert diameter to rebar area
        # Vertical stirrups angle
        alpha = 90  # degrees
        # Area of vertical stirrups
        A_vs = n_legs * A_db
        # Area of vertical stirrups per unit length (in^2/ft)
        A_v = A_vs / s
        best_combination = {
                            'n_stirrups': n_stirrups,
                            'd_b': d_b,
                            's': s,
                            'A_v': A_v,
                        }
        return best_combination
    
    def beam_transverse_rebar_EN_1992(self):
        pass
    
    def beam_transverse_rebar_EHE_08(self):
        pass
    
    # Factory method to select the transverse rebar method
    def beam_transverse_rebar(self,  A_v_req:float, V_s_req:float, lambda_factor:float, f_c:float, d:float):
        if self.beam.concrete.design_code=="ACI 318-19":
            return self.beam_transverse_rebar_ACI_318_19(A_v_req, V_s_req, lambda_factor, f_c, d)
        elif self.beam.concrete.design_code=="EN 1992":
            return self.beam_transverse_rebar_EN_1992()
        elif self.beam.concrete.design_code=="EHE-08":
            return self.beam_transverse_rebar_EHE_08()
        else:
            raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")
    
    # Factory method to select the longitudinal rebar method
    def beam_longitudinal_rebar(self,  A_s_req:float):
        if self.beam.concrete.design_code=="ACI 318-19":
            return self.beam_longitudinal_rebar_ACI_318_19(A_s_req)
        elif self.beam.concrete.design_code=="EN 1992":
            return self.beam_longitudinal_rebar_EN_1992()
        elif self.beam.concrete.design_code=="EHE-08":
            return self.beam_longitudinal_rebar_EHE_08()
        else:
            raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")