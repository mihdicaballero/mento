import sys
import os
from devtools import debug

# Add the project root to PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from dataclasses import dataclass
from structurelab.concrete.rectangular import RectangularConcreteSection
from structurelab import material
from structurelab.rebar import Rebar
from structurelab.units import MPa, ksi, psi, kip, mm, inch, kN, m, cm
from structurelab.results import Formatter
from IPython.display import Markdown, display
import numpy as np
import math

@dataclass
class Beam(RectangularConcreteSection):
    def __init__(self, name: str, concrete: material.Concrete, steelBar: material.SteelBar, 
                 width: float, depth: float,settings=None):  
        super().__init__(name, concrete, steelBar, width, depth, settings)
        self.shear_design_results = None

    def __determine_maximum_flexural_reinforcement_ratio_ACI_318_19(self):
        # Determination of maximum reinforcement ratio
        concrete_properties=self.concrete.get_properties()
        beta_1=concrete_properties["beta_1"]
        f_c=concrete_properties["f_c"]
        epsilon_c=concrete_properties["epsilon_c"]
        rebar_properties=self.steelBar.get_properties()
        f_y=rebar_properties["f_y"]
        epsilon_y=rebar_properties['epsilon_ty']

        epsilon_min_rebar_ACI_318_19=epsilon_y+epsilon_c # ESTO CHEQUEARLO BIEN, CREO QUE ES ASI, PERO REVISAR

        rho_max=0.85*beta_1*f_c/f_y*(epsilon_c/(epsilon_c+epsilon_min_rebar_ACI_318_19))
        return rho_max

    def __calculate_phi_ACI_318_19(self, epsilon_mas_deformado:float):
        # CREO QUE NO LA USO PARA NADA, PERO OJO, NO SE. 
        rebar_properties=self._steelBar.get_properties()
        concrete_properties=self._concrete.get_properties()
        epsilon_c=concrete_properties["epsilon_c"]
        rebar_properties=self._steelBar.get_properties()
        epsilon_y=rebar_properties["epsilon_y"]

        if epsilon_mas_deformado<=epsilon_y:
            return 0.65
        elif epsilon_mas_deformado<=epsilon_y+epsilon_c:
            return (0.9-0.65)*(epsilon_mas_deformado-epsilon_y)/epsilon_c+0.65
        else:
            return 0.9

    def design_flexure_ACI_318_19(self, M_u:float):
        max_rho=self.__determine_maximum_flexural_reinforcement_ratio_ACI_318_19()
        self._settings.load_aci_318_19_settings()
        phi = self._settings.get_setting('phi_t')
        setting_flexural_min_reduction = self._settings.get_setting('flexural_min_reduction')
        concrete_properties=self.concrete.get_properties()
        f_c=concrete_properties["f_c"]
        rebar_properties=self.steelBar.get_properties()
        f_y=rebar_properties["f_y"]
        d=0.9*self._depth # Asumption
        A=-phi*f_y**2/(1.7*f_c*self._width)
        B=phi*f_y*d
        C=-M_u
        A_s_calc=(-B+np.sqrt(B**2 - 4*A*C))/(2*A)
        A_s_min=max((0.25*np.sqrt(f_c / MPa)*MPa/f_y*self._depth*self._width/cm**2) , (1.4*MPa/f_y*self._depth*self._width/cm**2))*cm**2

        if A_s_calc>A_s_min:
            self.__A_s_calculated=A_s_calc
        elif  4*A_s_calc/3 > A_s_min:
            self.__A_s_calculated=A_s_min
        else: 
            if setting_flexural_min_reduction=="True":
                self.__A_s_calculated=4*A_s_calc/3
            else:
                self.__A_s_calculated=A_s_min
            
        rho_max=self.__determine_maximum_flexural_reinforcement_ratio_ACI_318_19()
        A_s_max=rho_max*self._depth*self._width

        if self.__A_s_calculated > A_s_max:
            debug("The section cannot be reinforced; increase the dimensions or use compressive reinforcement.")
            return
        else:
            result={
                'As_min_code':A_s_min,
                'As_required':A_s_calc,
                'As_max':A_s_max,
                'As_adopted':self.__A_s_calculated
            }
            return result

    def design_flexure_EN_1992(M_u):
        pass

    def design_flexure_EHE_08(M_u):
        pass

    def design_flexure(self, M_u:float):
        if self.concrete.design_code=="ACI 318-19":
            return self.design_flexure_ACI_318_19(M_u)
        elif self.concrete.design_code=="EN 1992":
            return self.design_flexure_EN_1992(M_u)
        elif self.concrete.design_code=="EHE-08":
            return self.design_flexure_EHE_08(M_u)
        else:
            raise ValueError("Concrete design code not supported")


    def check_flexure_ACI_318_19(self, M_u:float, A_s, d_b=0.5, n_bars=2):
        pass

    def check_shear_ACI_318_19(self, V_u:float, N_u:float, A_s:float, d_b:float, s:float, n_legs:float):
        concrete_properties=self.concrete.get_properties()
        f_c=concrete_properties["f_c"]
        rebar_properties=self.steelBar.get_properties()
        f_y=rebar_properties["f_y"]
        cc = self._settings.get_setting('clear_cover')
        self._settings.load_aci_318_19_settings()
        phi_v = self._settings.get_setting('phi_v')
        f_yt = min(f_y,60*ksi)
        lambda_factor = self._settings.get_setting('lambda')

        # Minimum shear reinforcement calculation
        # 'Minimum reinforcement should be placed if the factored shear Vu is greater than half the shear capacity of the concrete,
        # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
        # Rebar needed, V_u > φ_v*V_c/2
        A_v_min = max((0.75 * math.sqrt(f_c / psi) * psi/ f_yt) * self._width , (50 * psi/f_yt) * self._width)  
        
        # Shear reinforcement calculations
        A_db = (d_b ** 2) * math.pi / 4  # Area of one stirrup leg
        A_vs = n_legs * A_db  # Total area of stirrups
        A_v = A_vs / s  # Stirrup area per unit length

        d_bs_ini = self._settings.get_setting('longitudinal_diameter') # Longitudinal diameter if none is defined
        d_stirrup = d_b
        d = self._depth-(cc+d_stirrup+d_bs_ini/2)

        V_s = A_v * f_yt * d  # Shear contribution of reinforcement
        phi_V_s = phi_v * V_s  # Reduced shear contribution of reinforcement

        # Effective shear area and longitudinal reinforcement ratio
        A_cv = self._width * d  # Effective shear area
        A_g = self._width * self._depth  # Gross area
        rho_w = A_s / A_cv  # Longitudinal reinforcement ratio
        
        # Size modification factor
        lambda_s = math.sqrt(2 / (1 + d / (10*inch)))

        # Concrete shear strength calculation
        sigma_Nu = min(N_u / (6 * A_g), 0.05 * f_c)  # Axial stress influence

        # Concrete shear capacity depending on whether min rebar is present or not
        if A_v < A_v_min:
            k_c_min_rebar = 8 * lambda_s * lambda_factor * rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + sigma_Nu
        else:           
            k_c_min_rebar = max(2 * lambda_factor * math.sqrt(f_c / psi) * psi + sigma_Nu,
                            8 * lambda_factor * rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + sigma_Nu)

        V_cmin = 0*kip  
        V_cmax = (5 * lambda_factor * math.sqrt(f_c / psi) * psi) * A_cv  # Maximum concrete shear strength
        
        # Calculate actual concrete shear strength
        V_c = min(V_cmax, max(V_cmin, k_c_min_rebar * A_cv))
        phi_V_c = phi_v * V_c  # Reduced concrete shear strength
        
        # Maximum total shear capacity
        V_max = V_c + (8 * lambda_factor * math.sqrt(f_c / psi) * psi) * A_cv
        phi_V_max = phi_v * V_max  # Reduced maximum shear capacity

        if V_u < phi_V_c/2:
            A_v_min = 0*inch
            max_shear_ok = True
        elif phi_V_c/2 < V_u < phi_V_max:
            A_v_min = max((0.75 * math.sqrt(f_c / psi) * psi/ f_yt) * self._width , (50 * psi/f_yt) * self._width)  
            max_shear_ok = True
        else:
            max_shear_ok = False 

        # Total shear strength
        phi_V_n = phi_v * (V_c + V_s)  # Total reduced shear strength (concrete + rebar)

        # Required shear reinforcing nominal strength
        V_s_req = V_u-phi_V_c
        # Required shear reinforcing
        A_v_req = max(V_s_req/(phi_v*f_yt*d), A_v_min)

        # Check results
        result = {
            'A_v_min': A_v_min,  # Minimum shear reinforcement area
            'A_v_req': A_v_req, # Required shear reinforcing area
            'A_v': A_v,  # Provided stirrup reinforcement per unit length
            'phi_V_c': phi_V_c,  # Concrete contribution to shear capacity
            'phi_V_s': phi_V_s,  # Reinforcement contribution to shear capacity
            'phi_V_n': phi_V_n,  # Total shear capacity
            'phi_V_max': phi_V_max,  # Maximum shear capacity
            'shear_ok': V_u <= phi_V_n,  # Check if applied shear is within total capacity
            'max_shear_ok': max_shear_ok,  # Check if applied shear is within max shear capacity
            "FUv" : V_u / phi_V_n 
        }

        return result
    
    def design_shear_ACI_318_19(self, V_u:float, N_u:float, A_s:float):
        concrete_properties=self.concrete.get_properties()
        f_c=concrete_properties["f_c"]
        rebar_properties=self.steelBar.get_properties()
        f_y=rebar_properties["f_y"]
        cc = self._settings.get_setting('clear_cover')
        self._settings.load_aci_318_19_settings()
        phi_v = self._settings.get_setting('phi_v')
        f_yt = min(f_y,60*ksi)
        lambda_factor = self._settings.get_setting('lambda')
        # Inputs into class
        self.V_u = V_u
        self.N_u = N_u
        self.A_s = A_s 
        # Effective height
        d_bs_ini = self._settings.get_setting('longitudinal_diameter') # Longitudinal diameter if none is defined
        d_stirrup_ini = self._settings.get_setting('stirrup_diameter')
        self.d = self._depth-(cc+d_stirrup_ini+d_bs_ini/2)
        # Effective shear area and longitudinal reinforcement ratio
        A_cv = self._width * self.d  # Effective shear area
        A_g = self._width * self._depth  # Gross area
        rho_w = A_s / A_cv  # Longitudinal reinforcement ratio
        
        # Concrete shear strength calculation
        sigma_Nu = min(N_u / (6 * A_g), 0.05 * f_c)  # Axial stress influence

        # Concrete shear capacity assuming that the beam is provided with minimum shear rebar:
        k_c_min_rebar = max(2 * lambda_factor * math.sqrt(f_c / psi) * psi + sigma_Nu,
                            8 * lambda_factor * rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + sigma_Nu)

        V_cmin = 0*kip 
        V_cmax = (5 * lambda_factor * math.sqrt(f_c / psi) * psi) * A_cv  # Maximum concrete shear strength
        
        # Calculate actual concrete shear strength
        V_c = min(V_cmax, max(V_cmin, k_c_min_rebar * A_cv))
        phi_V_c = phi_v * V_c  # Reduced concrete shear strength
        
        # Maximum total shear capacity
        V_max = V_c + (8 * lambda_factor * math.sqrt(f_c / psi) * psi) * A_cv
        phi_V_max = phi_v * V_max  # Reduced maximum shear capacity

        # Minimum shear reinforcement calculation
        # Minimum reinforcement should be placed if the factored shear Vu is greater than half the shear capacity of the concrete, 
        # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
        if V_u < phi_V_c/2:
            A_v_min = 0*inch
            max_shear_ok = True
        elif phi_V_c/2 < V_u < phi_V_max:
            A_v_min = max((0.75 * math.sqrt(f_c / psi) * psi/ f_yt) * self._width , (50 * psi/f_yt) * self._width) 
            max_shear_ok = True
        else:
            max_shear_ok = False 

        # Shear reinforcement calculations
        # Required shear reinforcing nominal strength
        V_s_req = V_u-phi_V_c
        # Required shear reinforcing
        A_v_req = max(V_s_req/(phi_v*f_yt*self.d), A_v_min)

        # Design input for rebar detailing
        self.shear_rebar_input = {
            'A_v_req': A_v_req,
            'V_s_req': V_s_req,
            'd': self.d  # effective depth         
        }

        A_v_design = Rebar(self).beam_transverse_rebar(self.shear_rebar_input)
        # Store the final rebar design results as a property of the Beam instance
        self.transverse_rebar = {
            'n_stirrups': A_v_design['n_stirrups'],
            'd_b': A_v_design['d_b'],
            'spacing': A_v_design['spacing'],
            'A_v': A_v_design['A_v']
        }

        A_v =  A_v_design['A_v']
        V_s = A_v * f_yt * self.d  # Shear contribution of reinforcement
        phi_V_s = phi_v * V_s  # Reduced shear contribution of reinforcement
        # Total shear strength
        phi_V_n = phi_v * (V_c + V_s)  # Total reduced shear strength (concrete + rebar)

        # Design results
        self.shear_design_results = {
            'A_v_min': A_v_min,  # Minimum shear reinforcement area
            'A_v_req': A_v_req, # Required shear reinforcing area
            'A_v': A_v,  # Provided stirrup reinforcement per unit length
            'phi_V_c': phi_V_c,  # Concrete contribution to shear capacity
            'phi_V_s': phi_V_s,  # Reinforcement contribution to shear capacity
            'phi_V_n': phi_V_n,  # Total shear capacity
            'phi_V_max': phi_V_max,  # Maximum shear capacity
            'shear_ok': V_u <= phi_V_n,  # Check if applied shear is within total capacity
            'max_shear_ok': max_shear_ok,  # Check if applied shear is within max shear capacity
            "FUv" : V_u / phi_V_n 
        }
        return self.shear_design_results
    
    def check_shear_EN_1992(self, V_u: float, N_u: float, A_s: float):
        pass
 
    def design_shear_EN_1992(self, V_u: float, N_u: float, A_s: float):
        pass

    def design_shear_EHE_08(self, V_u: float, N_u: float, A_s: float):
        pass

    def check_shear_EHE_08(self, V_u: float, N_u: float, A_s: float):
        pass
  
    # Factory method to select the shear design method
    def design_shear(self, V_u: float, N_u: float, A_s: float):
        if self.concrete.design_code=="ACI 318-19":
            return self.design_shear_ACI_318_19(V_u, N_u, A_s)
        elif self.concrete.design_code=="EN 1992":
            return self.design_shear_EN_1992(V_u, N_u, A_s)
        elif self.concrete.design_code=="EHE-08":
            return self.design_shear_EHE_08(V_u, N_u, A_s)
        else:
            raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")
    
    # Factory method to select the shear check method
    def check_shear(self, V_u: float, N_u: float, A_s: float, d_b:float, s:float, n_legs:int):
        if self.concrete.design_code=="ACI 318-19":
            return self.check_shear_ACI_318_19(V_u, N_u, A_s, d_b, s, n_legs)
        elif self.concrete.design_code=="EN 1992":
            return self.check_shear_EN_1992(V_u, N_u, A_s, d_b, s, n_legs)
        elif self.concrete.design_code=="EHE-08":
            return self.check_shear_EHE_08(V_u, N_u, A_s, d_b, s, n_legs)
        else:
            raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")
        
    # Beam results forJjupyter Notebook
    @property
    def shear_results(self):
        if not self.shear_design_results:
            raise ValueError("Shear design has not been performed yet. Call design_shear first.")

        FU_v = self.shear_design_results['FUv']
        # Create FUFormatter instance and format FU value
        formatter = Formatter()
        formatted_FU = formatter.FU(FU_v)
        rebar_v = f"{self.transverse_rebar['n_stirrups']}eØ{self.transverse_rebar['d_b']}/{self.transverse_rebar['spacing']}"
        # Print results
        markdown_content = f"Armadura transversal {rebar_v}, $A_v$={round(self.transverse_rebar['A_v'].value*10000,2)}" \
                         f"cm²/m, $V_u$={round(self.V_u,0)}, $\\phi V_n$={round(self.shear_design_results['phi_V_n'],0)} → {formatted_FU}"
        # Display the content
        display(Markdown(markdown_content))

        return markdown_content


def main():
    concrete=material.create_concrete(name="H30",f_c=30*MPa, design_code="ACI 318-19") 
    steelBar=material.SteelBar(name="ADN 420", f_y=420*MPa) 
    section = Beam(
        name="V-40x50",
        concrete=concrete,
        steelBar=steelBar,
        width=400 * mm,  
        depth=500 * mm,  
    )

    debug(f"Nombre de la sección: {section.get_name()}")
    resultados=section.design_flexure(500*kN*m)  
    debug(resultados)


def shear():
    # Define custom settings
    custom_settings = {
        'clear_cover': 1.5*inch, 
        'stirrup_diameter': 0.5*inch, 
        'longitudinal_diameter': 1*inch 
        }
    concrete=material.create_concrete(name="C4",f_c=4000*psi, design_code="ACI 318-19") 
    steelBar=material.SteelBar(name="ADN 420", f_y=60*ksi) 
    section = Beam(name="V-10x16",concrete=concrete,steelBar=steelBar,width=10*inch, depth=16*inch)
    section.update_settings(custom_settings)
    V_u=37.727*kip 
    N_u = 0*kip 
    A_s=0.847*inch**2 
    results=section.check_shear(V_u, N_u, A_s, d_b=12*mm, s=6*inch, n_legs=2) 
    debug(results)
    section.design_shear(V_u, N_u, A_s) 
    debug(section.shear_results)
    debug(section.shear_design_results)
    debug(section.transverse_rebar)


if __name__ == "__main__":
    shear()