import sys
import os

# Add the project root to PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from dataclasses import dataclass
from structurelab.concrete.rectangular import RectangularConcreteSection
from structurelab import material
import numpy as np
import math
import forallpeople
forallpeople.environment('structural', top_level=True)
# Definir algunas unidades adicionales útiles
cm = 1e-2 * m  # type: ignore

@dataclass
class Beam(RectangularConcreteSection):
    def __init__(self, name: str, concrete: material.Concrete, steelBar: material.SteelBar, 
                 width: float, depth: float,settings=None):  # type: ignore
        super().__init__(name, concrete, steelBar, width, depth, settings)

    def __determine_maximum_flexural_reinforcement_ratio_ACI_318_19(self):
        # Determination of maximum reinforcement ratio
        concrete_properties=self.concrete.get_properties()
        beta_1=concrete_properties["beta_1"]
        f_c=concrete_properties["f_c"]
        epsilon_c=concrete_properties["epsilon_c"]
        rebar_properties=self.steelBar.get_properties()
        f_y=rebar_properties["f_y"]
        epsilon_ty=rebar_properties['epsilon_ty']

        epsilon_min_rebar_ACI_318_19=epsilon_ty+epsilon_c # ESTO CHEQUEARLO BIEN, CREO QUE ES ASI, PERO REVISAR

        rho_max=0.85*beta_1*f_c/f_y*(epsilon_c/(epsilon_c+epsilon_min_rebar_ACI_318_19))
        return rho_max

    def __calculate_phi_ACI_318_19(self, epsilon_mas_deformado:float):
        # CREO QUE NO LA USO PARA NADA, PERO OJO, NO SE. 
        rebar_properties=self._steelBar.get_properties()
        concrete_properties=self._concrete.get_properties()
        epsilon_c=concrete_properties["epsilon_c"]
        rebar_properties=self._steelBar.get_properties()
        epsilon_ty=rebar_properties["epsilon_ty"]

        if epsilon_mas_deformado<=epsilon_ty:
            return 0.65
        elif epsilon_mas_deformado<=epsilon_ty+epsilon_c:
            return (0.9-0.65)*(epsilon_mas_deformado-epsilon_ty)/epsilon_c+0.65
        else:
            return 0.9

    def design_flexure_ACI_318_19(self, M_u:float):
        self._settings.load_aci_318_19_settings()
        phi = self._settings.get_setting('phi_t')
        setting_flexural_min_reduction = self._settings.get_setting('flexural_min_reduction')
        concrete_properties=self.concrete.get_properties()
        f_c=concrete_properties['f_c']
        beta_1=concrete_properties['beta_1']
        rebar_properties=self.steelBar.get_properties()
        f_y=rebar_properties['f_y']
        epsilon_ty=rebar_properties['epsilon_ty']
        E_s=rebar_properties['E_s']

        d=0.9*self._depth # Asumption, this is the main difference between design and check.
        b=self._width

        
        # Determination of minimum reinforcement
        A_s_min=max((0.25*np.sqrt(f_c / MPa)*MPa/f_y*self._depth*self._width/cm**2) , (1.4*MPa/f_y*self._depth*self._width/cm**2))*cm**2# type: ignore

        # Determination of maximum reinforcement
        rho_max=self.__determine_maximum_flexural_reinforcement_ratio_ACI_318_19()
        A_s_max=rho_max*d*b

        # Determination of required reinforcement
        R_n=M_u/(phi*b*d**2)
        A_s_calc=0.85*f_c*b*d/f_y*(1-np.sqrt(1-2*R_n/(0.85*f_c)))

        if A_s_calc>A_s_min:
            self._A_s_calculated=A_s_calc
        elif  4*A_s_calc/3 > A_s_min:
            self._A_s_calculated=A_s_min
        else: 
            if setting_flexural_min_reduction=='True':
                self._A_s_calculated=4*A_s_calc/3
            else:
                self._A_s_calculated=A_s_min
        if self._A_s_calculated <= A_s_max: 
            self._A_s_comp=0
            result={
                'As_min_code':A_s_min,
                'As_required':A_s_calc,
                'As_max':A_s_max,
                'As_adopted':self._A_s_calculated,
                'As_compression':self._A_s_comp
            }
            return result
        else:
            rho=0.85*beta_1*f_c/f_y*(0.003/(epsilon_ty+0.006))
            M_n_t=rho*f_y*(d-0.59*rho*f_y*d/f_c)*b*d
            M_n_prima=M_u/phi-M_n_t
            c_t=0.003*d/(epsilon_ty+0.006)
            # HAY QUE VER DONDE ESTA DEFINIDO EL d_prima por ahora asumo
            d_prima=5*cm
            f_s_prima=min(0.003*E_s*(1-d_prima/c_t),f_y)
            A_s_prima=M_n_prima/(f_s_prima*(d-d_prima))
            A_s=rho*b*d+A_s_prima
            self._A_s_calculated=A_s
            self._A_s_comp=A_s_prima
            result={
                'As_min_code':A_s_min,
                'As_required':self._A_s_calculated,
                'As_max':A_s_max,
                'As_adopted':self._A_s_calculated,
                'As_compression':self._A_s_comp
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
        f_yt = min(f_y,60*ksi) #type: ignore
        lambda_factor = self._settings.get_setting('lambda')

        # Minimum shear reinforcement calculation
        # 'Minimum reinforcement should be placed if the factored shear Vu is greater than half the shear capacity of the concrete,
        # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
        # Rebar needed, V_u > φ_v*V_c/2
        A_v_min = max((0.75 * math.sqrt(f_c / psi) * psi/ f_yt) * self._width , (50 * psi/f_yt) * self._width)  # type: ignore
        
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
        lambda_s = math.sqrt(2 / (1 + d / (10*inch))) #type: ignore

        # Concrete shear strength calculation
        sigma_Nu = min(N_u / (6 * A_g), 0.05 * f_c)  # Axial stress influence

        # Concrete shear capacity depending on whether min rebar is present or not
        if A_v < A_v_min:
            k_c_min_rebar = 8 * lambda_s * lambda_factor * rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + sigma_Nu
        else:           
            k_c_min_rebar = max(2 * lambda_factor * math.sqrt(f_c / psi) * psi + sigma_Nu, #type: ignore
                            8 * lambda_factor * rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + sigma_Nu) #type: ignore

        V_cmin = 0*kip  # type: ignore
        V_cmax = (5 * lambda_factor * math.sqrt(f_c / psi) * psi) * A_cv  # Maximum concrete shear strength #type: ignore
        
        # Calculate actual concrete shear strength
        V_c = min(V_cmax, max(V_cmin, k_c_min_rebar * A_cv))
        phi_V_c = phi_v * V_c  # Reduced concrete shear strength
        
        # Maximum total shear capacity
        V_max = V_c + (8 * lambda_factor * math.sqrt(f_c / psi) * psi) * A_cv #type: ignore
        phi_V_max = phi_v * V_max  # Reduced maximum shear capacity

        if V_u < phi_V_c/2:
            A_v_min = 0*inch #type: ignore
            max_shear_ok = True
        elif phi_V_c/2 < V_u < phi_V_max:
            A_v_min = max((0.75 * math.sqrt(f_c / psi) * psi/ f_yt) * self._width , (50 * psi/f_yt) * self._width)  # type: ignore
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
        f_yt = min(f_y,60*ksi) #type: ignore
        lambda_factor = self._settings.get_setting('lambda')
        # Effective height
        d_bs_ini = self._settings.get_setting('longitudinal_diameter') # Longitudinal diameter if none is defined
        d_stirrup_ini = self._settings.get_setting('stirrup_diameter')
        d = self._depth-(cc+d_stirrup_ini+d_bs_ini/2)
        # Effective shear area and longitudinal reinforcement ratio
        A_cv = self._width * d  # Effective shear area
        A_g = self._width * self._depth  # Gross area
        rho_w = A_s / A_cv  # Longitudinal reinforcement ratio
        
        # Concrete shear strength calculation
        sigma_Nu = min(N_u / (6 * A_g), 0.05 * f_c)  # Axial stress influence

        # Concrete shear capacity assuming that the beam is provided with minimum shear rebar:
        k_c_min_rebar = max(2 * lambda_factor * math.sqrt(f_c / psi) * psi + sigma_Nu, #type: ignore
                            8 * lambda_factor * rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + sigma_Nu) #type: ignore

        V_cmin = 0*kip  # type: ignore
        V_cmax = (5 * lambda_factor * math.sqrt(f_c / psi) * psi) * A_cv  # Maximum concrete shear strength #type: ignore
        
        # Calculate actual concrete shear strength
        V_c = min(V_cmax, max(V_cmin, k_c_min_rebar * A_cv))
        phi_V_c = phi_v * V_c  # Reduced concrete shear strength
        
        # Maximum total shear capacity
        V_max = V_c + (8 * lambda_factor * math.sqrt(f_c / psi) * psi) * A_cv #type: ignore
        phi_V_max = phi_v * V_max  # Reduced maximum shear capacity

        # Minimum shear reinforcement calculation
        # Minimum reinforcement should be placed if the factored shear Vu is greater than half the shear capacity of the concrete, 
        # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
        if V_u < phi_V_c/2:
            A_v_min = 0*inch #type: ignore
            max_shear_ok = True
        elif phi_V_c/2 < V_u < phi_V_max:
            A_v_min = max((0.75 * math.sqrt(f_c / psi) * psi/ f_yt) * self._width , (50 * psi/f_yt) * self._width)  # type: ignore
            max_shear_ok = True
        else:
            max_shear_ok = False 

        # Shear reinforcement calculations
        # Required shear reinforcing nominal strength
        V_s_req = V_u-phi_V_c
        # Required shear reinforcing
        A_v_req = max(V_s_req/(phi_v*f_yt*d), A_v_min)

        # CALL REBAR CLASS AND DESIGN THE STIRRUPS OF THE SECTION and get A_v provided
        # The result should be A_v =  0.0654498*inch and '1eØ0.5/6'
        A_v =  0.0654498*inch #type: ignore
        V_s = A_v * f_yt * d  # Shear contribution of reinforcement
        phi_V_s = phi_v * V_s  # Reduced shear contribution of reinforcement
        # Total shear strength
        phi_V_n = phi_v * (V_c + V_s)  # Total reduced shear strength (concrete + rebar)

        # Design results
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


def flexure():
    # Example 6.9.2 CRSI Design Guide - p6-64 (355 pdf)
    concrete=material.create_concrete(name="fc4000",f_c=4000*psi, design_code="ACI 318-19") # type: ignore
    steelBar=material.SteelBar(name="G60", f_y=60000*psi) # type: ignore
    section = Beam(
        name="B-28x24",
        concrete=concrete,
        steelBar=steelBar,
        width=28 * inch,  # type: ignore
        depth=24 * inch,  # type: ignore
    )

    print(f"Nombre de la sección: {section.get_name()}")
    resultados=section.design_flexure(156.1*kip*ft)  # type: ignore
    print(resultados)


def shear():
    # Define custom settings
    custom_settings = {
        'clear_cover': 1.5*inch, # type: ignore
        'stirrup_diameter': 0.5*inch, # type: ignore
        'longitudinal_diameter': 1*inch # type: ignore
        }
    concrete=material.create_concrete(name="C4",f_c=4000*psi, design_code="ACI 318-19") # type: ignore
    steelBar=material.SteelBar(name="ADN 420", f_y=60*ksi) # type: ignore
    section = Beam(
        name="V-10x16",
        concrete=concrete,
        steelBar=steelBar,
        width=10*inch,  # type: ignore
        depth=16*inch,  # type: ignore
    )
    section.update_settings(custom_settings)
    V_u=37.727*kip # type: ignore
    N_u = 0*kip # type: ignore
    A_s=0.847*inch**2 # type: ignore
    results=section.check_shear_ACI_318_19(V_u, N_u, A_s, d_b=0.5*inch, s=6*inch, n_legs=2) # type: ignore
    print(results)
    results=section.design_shear_ACI_318_19(V_u, N_u, A_s) # type: ignore
    print(results)

if __name__ == "__main__":
    flexure()
    shear()