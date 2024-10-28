from devtools import debug
from dataclasses import dataclass
from IPython.display import Markdown, display
from typing import Optional, Dict, Any, cast, List
from pint.facets.plain import PlainQuantity
import numpy as np
import pandas as pd
from pandas import DataFrame
import math
import warnings

from mento.concrete.rectangular import RectangularConcreteSection
from mento.material import Concrete, SteelBar, Concrete_ACI_318_19 
from mento.rebar import Rebar
from mento import MPa, ksi, psi, kip, mm, inch, kN, m, cm
from mento.results import Formatter, TablePrinter, DocumentBuilder
from mento.forces import Forces  

@dataclass
class RectangularConcreteBeam(RectangularConcreteSection):
    label: Optional[str] = None

    def __init__(self, label: Optional[str], concrete: Concrete, steel_bar: SteelBar, 
                 width: PlainQuantity, height: PlainQuantity, settings: Optional[Dict[str, Any]] = None):   
        super().__init__(concrete, steel_bar, width, height, settings)
        if settings:
            self.settings.update(settings)  # Update with any provided settings

        self.label = label
        self.shear_design_results: DataFrame = None
        self._stirrup_s_l: PlainQuantity = 0*cm
        self._stirrup_n: int = 0
        self._stirrup_A_v: PlainQuantity = 0*cm**2/m
        self._phi_V_n: PlainQuantity = 0*kN
        self._phi_V_s: PlainQuantity = 0*kN
        self._phi_V_c: PlainQuantity = 0*kN
        self._phi_V_max: PlainQuantity = 0*kN
        self._A_s_req_bot: PlainQuantity = 0*cm**2
        self._V_u: PlainQuantity = 0*kN
        self._N_u: PlainQuantity = 0*kN
        self._A_cv: PlainQuantity = 0*cm**2
        self._k_c_min: PlainQuantity = 0*MPa
        self._sigma_u: PlainQuantity = 0*MPa
        self._A_s: PlainQuantity = 0*cm**2
        self._A_v_req: PlainQuantity = 0*cm**2
    
    @property
    def d(self) -> PlainQuantity:
        "Effective height."
        return self._d
    
    def set_transverse_rebar(self, n_stirrups: int = 0, d_b:PlainQuantity = 0*mm, s_l:PlainQuantity = 0*cm) -> None:
        """Sets the transverse rebar in the object."""
        self._stirrup_n = n_stirrups
        self._stirrup_d_b = d_b
        self._stirrup_s_l = s_l
        # Update effective height d with new values
        self._d = self._height -(self.c_c+self._stirrup_d_b+self._long_d_b/2) # Initial value 

    def set_longitudinal_rebar_bot(self, n1: int, d_b1: PlainQuantity, n2: int, d_b2: PlainQuantity, 
                                n3: int, d_b3: PlainQuantity, n4: int, d_b4: PlainQuantity, 
                                total_as: PlainQuantity, total_bars: int, clear_spacing: PlainQuantity) -> None:
        """Sets the longitudinal rebar in the object."""
        self._longitudinal_n1 = n1
        self._longitudinal_d_b1 = d_b1
        self._longitudinal_n2 = n2
        self._longitudinal_d_b2 = d_b2 if n2 > 0 else None
        self._longitudinal_n3 = n3
        self._longitudinal_d_b3 = d_b3 if n3 > 0 else None
        self._longitudinal_n4 = n4
        self._longitudinal_d_b4 = d_b4 if n4 > 0 else None
        self._total_as_bot = total_as
        self._total_bars_bot = total_bars
        self._clear_spacing_bot = clear_spacing
    
    def __maximum_flexural_reinforcement_ratio(self) -> float:
        if self.concrete.design_code=="ACI 318-19":
            concrete_aci = cast(Concrete_ACI_318_19, self.concrete)  # Cast to the specific subclass
            # Determination of maximum reinforcement ratio
            epsilon_min_rebar_ACI_318_19=self.steel_bar.epsilon_y+concrete_aci.epsilon_c # TODO: REVISAR

            rho_max=0.85*concrete_aci.beta_1/self.steel_bar.f_y*(concrete_aci.epsilon_c/(concrete_aci.epsilon_c+epsilon_min_rebar_ACI_318_19))
            
            return rho_max
        else: 
            return 0

    def __calculate_phi_ACI_318_19(self, epsilon_mas_deformado:float) -> float:
        concrete_properties=self.concrete.get_properties()
        epsilon_c=concrete_properties["epsilon_c"]

        if epsilon_mas_deformado<=self.steel_bar.epsilon_y:
            return 0.65
        elif epsilon_mas_deformado<=self.steel_bar.epsilon_y+epsilon_c:
            return (0.9-0.65)*(epsilon_mas_deformado-self.steel_bar.epsilon_y)/epsilon_c+0.65
        else:
            return 0.9

    def design_flexure_ACI_318_19(self, M_u:float)-> Dict[str, Any]:
        self.settings.load_aci_318_19_settings()
        phi = self.settings.get_setting('phi_t')
        setting_flexural_min_reduction = self.settings.get_setting('flexural_min_reduction')
        concrete_properties=self.concrete.get_properties()
        f_c=concrete_properties['f_c']
        beta_1=concrete_properties['beta_1']
        
        # Determination of minimum reinforcement
        A_s_min=max((3*np.sqrt(f_c / psi)*psi/self.steel_bar.f_y*self.d*self._width),
                     (200*psi/self.steel_bar.f_y*self.d*self._width))

        # Determination of maximum reinforcement
        rho_max=self.__maximum_flexural_reinforcement_ratio()
        A_s_max=rho_max*self.d*self._width

        # Determination of required reinforcement
        R_n=M_u/(phi*self._width*self.d**2)
        A_s_req=0.85*f_c*self._width*self.d/self.steel_bar.f_y*(1-np.sqrt(1-2*R_n/(0.85*f_c)))

        if A_s_req>A_s_min:
            self._A_s_calculated=A_s_req
        elif  4*A_s_req/3 > A_s_min:
            self._A_s_calculated=A_s_min
        else: 
            if setting_flexural_min_reduction=='True':
                self._A_s_calculated=4*A_s_req/3
            else:
                self._A_s_calculated=A_s_min
        if self._A_s_calculated <= A_s_max: 
            self._A_s_comp=0
            result={
                'As_min_code':A_s_min,
                'As_required':A_s_req,
                'As_max':A_s_max,
                'As_adopted':self._A_s_calculated,
                'As_compression':self._A_s_comp
            }
            return result
        else:
            rho=0.85*beta_1*f_c/self.steel_bar.f_y*(0.003/(self.steel_bar.epsilon_y+0.006))
            M_n_t=rho*self.steel_bar.f_y*(self.d-0.59*rho*self.steel_bar.f_y*self.d/f_c)*self.width*self.d
            M_n_prima=M_u/phi-M_n_t
            c_t=0.003*self.d/(self.steel_bar.epsilon_y+0.006)
            # TODO: HAY QUE VER DONDE ESTA DEFINIDO EL d_prima por ahora asumo
            d_prima=5*cm
            f_s_prima=min(0.003*self.steel_bar.E_s*(1-d_prima/c_t),self.steel_bar.f_y)
            A_s_prima=M_n_prima/(f_s_prima*(self.d-d_prima))
            A_s=rho*self.width*self.d+A_s_prima
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


    def design_flexure_EN_1992(self, M_u: float) -> None:
        pass

    def design_flexure_EHE_08(self, M_u: float) -> None:
        pass

    def design_flexure(self, M_u:float) -> Dict[str, Any]:
        if self.concrete.design_code=="ACI 318-19":
            return self.design_flexure_ACI_318_19(M_u)
        # elif self.concrete.design_code=="EN 1992":
        #     return self.design_flexure_EN_1992(M_u)
        # elif self.concrete.design_code=="EHE-08":
        #     return self.design_flexure_EHE_08(M_u)
        else:
            raise ValueError(f"Longitudinal design method not implemented \
                    for concrete type: {type(self.concrete).__name__}")


    # def check_flexure_ACI_318_19(self, M_u:float, A_s, d_b=0.5, n_bars=2) -> None:
    #     pass

    def check_shear_ACI_318_19(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> DataFrame:
        # Set the initial variables
        self._N_u = Force.N_x
        self._V_u = Force.V_z
        d_bs = self._stirrup_d_b
        s_l = self._stirrup_s_l
        n_legs = self._stirrup_n*2
        self._A_s = A_s

        f_c=self.concrete.f_c
        self.settings.load_aci_318_19_settings()
        phi_v = self.settings.get_setting('phi_v')
        self.lambda_factor = self.settings.get_setting('lambda')
        f_yt = self.calculate_f_yt()

        # Minimum shear reinforcement calculation
        self._A_v_min = self.calculate_A_v_min(f_c)
         
        # Shear reinforcement calculations
        A_db = (d_bs ** 2) * math.pi / 4  # Area of one stirrup leg
        A_vs = n_legs * A_db  # Total area of stirrups
        A_v = A_vs / s_l  # Stirrup area per unit length
        self._stirrup_A_v = A_v

        V_s: PlainQuantity = A_v * f_yt * self.d  # Shear contribution of reinforcement
        self._phi_V_s = phi_v * V_s  # Reduced shear contribution of reinforcement
        # Effective shear area and longitudinal reinforcement ratio
        self._A_cv = self.width * self.d  # Effective shear area
        A_g = self.A_x  # Gross area
        self._rho_w = self._A_s.to('cm**2') / self._A_cv.to('cm**2')  # Longitudinal reinforcement ratio
     
        # Size modification factor for Imperial system
        self._lambda_s = math.sqrt(2 / (1 + self.d / (10*inch)))

        # Concrete shear strength calculation
        self._sigma_Nu = min(self._N_u / (6 * A_g), 0.05 * f_c)  # Axial stress influence

        # Concrete shear capacity depending on whether min rebar is present or not
        if A_v < self._A_v_min:
            self._k_c_min = 8 * self._lambda_s * self.lambda_factor * self._rho_w ** (1/3)\
                  * math.sqrt(f_c / psi) * psi + self._sigma_Nu
        else:           
            self._k_c_min = max(2 * self.lambda_factor * math.sqrt(f_c / psi) * psi + self._sigma_Nu,
                            8 * self.lambda_factor * self._rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + self._sigma_Nu)

        V_cmin = 0*kip  
        V_cmax = self.calculate_V_cmax()
        
        # Calculate actual concrete shear strength
        self.V_c = min(V_cmax, max(V_cmin, self._k_c_min * self._A_cv))
        self._phi_V_c = phi_v * self.V_c  # Reduced concrete shear strength
        
        # Maximum total shear capacity
        V_max = self.calculate_V_max()
        self._phi_V_max = phi_v * V_max  # Reduced maximum shear capacity

        if self._V_u < self._phi_V_c/2:
            self._A_v_min = 0*inch
            self._max_shear_ok = True
        elif self._phi_V_c/2 < self._V_u < self._phi_V_max:
            self._A_v_min = self.calculate_A_v_min(f_c)
            self._max_shear_ok = True
        else:
            self._A_v_min = self.calculate_A_v_min(f_c)
            self._max_shear_ok = False 

        # Total shear strength
        self._phi_V_n = phi_v * (self.V_c + V_s)  # Total reduced shear strength (concrete + rebar)

        # Required shear reinforcing nominal strength
        V_s_req: PlainQuantity = self._V_u-self._phi_V_c
        # Required shear reinforcing
        self._A_v_req = max(V_s_req/(phi_v*f_yt*self.d), self._A_v_min)
        self._FUv = (self._V_u.to('kN') / self._phi_V_n.to('kN'))

        # Rebar spacing checks
        section_rebar = Rebar(self)
        n_legs_actual = self._stirrup_n * 2      # Ensure legs are even
        self._s_l = self._stirrup_s_l
        self._s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs_actual - 1) 
        self._s_max_l, self._s_max_w = section_rebar.calculate_max_spacing(V_s_req, self._A_cv)

        # Check results
        results: Dict[str,Any] = {
            'Label': self.label, #Beam label
            'Av,min': self._A_v_min.to('cm ** 2 / m'),  # Minimum shear reinforcement area
            'Av,req': self._A_v_req.to('cm ** 2 / m'), # Required shear reinforcing area
            'Av': A_v.to('cm ** 2 / m'),  # Provided stirrup reinforcement per unit length
            'Vu': self._V_u.to('kN'), # Max Vu for the design
            'ØVc': self._phi_V_c.to('kN'),  # Concrete contribution to shear capacity
            'ØVs': self._phi_V_s.to('kN'),  # Reinforcement contribution to shear capacity
            'ØVn': self._phi_V_n.to('kN'),  # Total shear capacity
            'ØVmax': self._phi_V_max.to('kN'),  # Maximum shear capacity
            'Vu<ØVmax': self._max_shear_ok,  # Check if applied shear is within max shear capacity
            'Vu<ØVn': self._V_u <= self._phi_V_n,  # Check if applied shear is within total capacity
            "DCR" :  self._FUv
        }
        self._initialize_dicts_ACI_318_19()

        return pd.DataFrame([results], index=[0])
    
    def design_shear_ACI_318_19(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> DataFrame:
        self._N_u = Force.N_x
        self._V_u = Force.V_z
        f_c = self.concrete.f_c
        f_yt = self.steel_bar.f_yt
        self.settings.load_aci_318_19_settings()
        phi_v = self.settings.get_setting('phi_v')
        self.lambda_factor = self.settings.get_setting('lambda')
        # Inputs into class
        self._A_s = A_s 
        # Effective shear area and longitudinal reinforcement ratio
        self._A_cv = self.width * self.d  # Effective shear area
        A_g = self.A_x # Gross area
        self._rho_w = self._A_s / self._A_cv  # Longitudinal reinforcement ratio
        
        # Size modification factor for Imperial system
        self._lambda_s = math.sqrt(2 / (1 + self.d / (10*inch)))

        # Concrete shear strength calculation
        self._sigma_Nu = min(self._N_u / (6 * A_g), 0.05 * f_c)  # Axial stress influence

        # Concrete shear capacity assuming that the beam is provided with minimum shear rebar:
        self._k_c_min = max(2 * self.lambda_factor * math.sqrt(f_c / psi) * psi + self._sigma_Nu,
                            8 * self.lambda_factor * self._rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + self._sigma_Nu)

        V_cmin = 0*kip 
        V_cmax = self.calculate_V_cmax()
        
        # Calculate actual concrete shear strength
        self.V_c = min(V_cmax, max(V_cmin, self._k_c_min * self._A_cv))
        self._phi_V_c = phi_v * self.V_c  # Reduced concrete shear strength
        
        # Maximum total shear capacity for Imperial system
        V_max = self.calculate_V_max()
        self._phi_V_max = phi_v * V_max  # Reduced maximum shear capacity

        # Minimum shear reinforcement calculation
        # Minimum reinforcement should be placed if the factored shear Vu 
        # is greater than half the shear capacity of the concrete, 
        # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
        if self._V_u < self._phi_V_c/2:
            self._A_v_min = 0*inch
            self._max_shear_ok = True
        elif self._phi_V_c/2 < self._V_u < self._phi_V_max:
            self._A_v_min = self.calculate_A_v_min(f_c)
            self._max_shear_ok = True
        else:
            self._A_v_min = self.calculate_A_v_min(f_c)
            self._max_shear_ok = False 

        # Shear reinforcement calculations
        # Required shear reinforcing nominal strength
        V_s_req: PlainQuantity = self._V_u-self._phi_V_c
        # Required shear reinforcing
        self._A_v_req = max(V_s_req/(phi_v*self.steel_bar.f_yt*self.d), self._A_v_min)

        section_rebar = Rebar(self)
        self.shear_design_results = section_rebar.transverse_rebar(self._A_v_req, V_s_req)
        # Store the final rebar design results as a property of the Beam instance
        best_design = section_rebar.transverse_rebar_design
        self._s_l = best_design['s_l']
        self._s_w = best_design['s_w']
        self._s_max_l = best_design['s_max_l']
        self._s_max_w = best_design['s_max_w']
        self.set_transverse_rebar(best_design['n_stir'], best_design['d_b'], best_design['s_l'])
        self._stirrup_A_v = best_design['A_v']

        A_v: PlainQuantity =  best_design['A_v']
        V_s = A_v * f_yt * self.d  # Shear contribution of reinforcement
        self._phi_V_s = phi_v * V_s  # Reduced shear contribution of reinforcement
        # Total shear strength
        self._phi_V_n = phi_v * (self.V_c + V_s)  # Total reduced shear strength (concrete + rebar)
        self._FUv = (self._V_u.to('kN') / self._phi_V_n.to('kN'))

        # Design results
        results = {
            'Label': self.label, #Beam label
            'Av,min': self._A_v_min.to('cm ** 2 / m'),  # Minimum shear reinforcement area
            'Av,req': self._A_v_req.to('cm ** 2 / m'), # Required shear reinforcing area
            'Av': A_v.to('cm ** 2 / m'),  # Provided stirrup reinforcement per unit length
            'Vu': self._V_u.to('kN'), # Max Vu for the design
            'ØVc': self._phi_V_c.to('kN'),  # Concrete contribution to shear capacity
            'ØVs': self._phi_V_s.to('kN'),  # Reinforcement contribution to shear capacity
            'ØVn': self._phi_V_n.to('kN'),  # Total shear capacity
            'ØVmax': self._phi_V_max.to('kN'),  # Maximum shear capacity
            'Vu<ØVmax': self._max_shear_ok,  # Check if applied shear is within max shear capacity
            'Vu<ØVn': self._V_u <= self._phi_V_n,  # Check if applied shear is within total capacity
            "DCR" :  self._FUv
        }
        self._initialize_dicts_ACI_318_19()

        return pd.DataFrame([results], index=[0])
    
    def check_shear_EN_1992(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> None:
        return None
 
    def design_shear_EN_1992(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> None:
        return None

    def design_shear_EHE_08(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> None:
        return None

    def check_shear_EHE_08(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> None:
        return None
  
    # Factory method to select the shear design method
    def design_shear(self, Force: Forces, A_s: PlainQuantity = 0*cm**2) -> DataFrame:
        if self.concrete.design_code=="ACI 318-19":
            return self.design_shear_ACI_318_19(Force, A_s)
        # elif self.concrete.design_code=="EN 1992":
        #     return self.design_shear_EN_1992(V_u, N_u, A_s)
        # elif self.concrete.design_code=="EHE-08":
        #     return self.design_shear_EHE_08(V_u, N_u, A_s)
        else:
            raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")
    
    # Factory method to select the shear check method
    def check_shear(self, Force: Forces = Forces(), A_s: PlainQuantity = 0*cm**2) -> DataFrame:
        if self.concrete.design_code=="ACI 318-19":
            return self.check_shear_ACI_318_19(Force, A_s)
        # elif self.concrete.design_code=="EN 1992":
        #     return self.check_shear_EN_1992(V_u, N_u, A_s, d_b, s, n_legs)
        # elif self.concrete.design_code=="EHE-08":
        #     return self.check_shear_EHE_08(V_u, N_u, A_s, d_b, s, n_legs)
        else:
            raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")

    def calculate_A_v_min(self, f_c: PlainQuantity) -> PlainQuantity:
        """Calculate the minimum shear reinforcement based on unit system."""
        # 'Minimum reinforcement should be placed if the factored shear Vu 
        # is greater than half the shear capacity of the concrete,
        # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
        # Rebar needed, V_u > φ_v*V_c/2 for Imperial system
        f_yt = self.calculate_f_yt()
        
        if self.concrete.unit_system == "metric":
            A_v_min = max(
                (0.062 * math.sqrt(f_c / MPa) * MPa / f_yt) * self.width,
                (0.35 * MPa / f_yt) * self.width
            )
        else:
            A_v_min = max(
                (0.75 * math.sqrt(f_c / psi) * psi / f_yt) * self.width,
                (50 * psi / f_yt) * self.width
            )
        
        return A_v_min    
    
    def calculate_V_max(self) -> PlainQuantity:
        "Formula for maximum total shear capacity (V_max)"
        V_max = self.V_c + (8 * self.lambda_factor * math.sqrt(self.concrete.f_c / psi) * psi) * self._A_cv
        return V_max
    
    def calculate_V_cmax(self) -> PlainQuantity:
        "Maximum concrete shear strength"
        V_cmax = (5 * self.lambda_factor * math.sqrt(self.concrete.f_c / psi) * psi) * self._A_cv
        return V_cmax 

    def calculate_f_yt(self) -> PlainQuantity:
        """Determine the yield strength of steel based on unit system."""
        if self.concrete.unit_system == "metric":
            return min(self.steel_bar.f_yt, 400 * MPa)
        else:
            return min(self.steel_bar.f_yt, 60 * ksi)

    def _initialize_dicts_ACI_318_19(self) -> None:
        """Initialize the dictionaries used in check and design methods."""
        self._materials = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
                "Concrete density",
                "Normalweight concrete",
                "Safety factor for shear"
            ],
            "Variable": ["","fc", "fy", "γc", "λ", "Øv"],
            "Value": [self.label, round(self.concrete.f_c.to('MPa').magnitude,2), 
                      round(self.steel_bar.f_y.to('MPa').magnitude,2),self.concrete.density.magnitude,
                       self.settings.get_setting('lambda'), self.settings.get_setting('phi_v')],
            "Unit": ["", "MPa", "MPa", "kg/m³", "", ""]
        }
        self._geometry = {
            "Geometry": [
                "Section height",
                "Section width",
                "Clear cover",
                "Longitudinal tension rebar",
            ],
            "Variable": ["h", "b", "cc", "As"],
            #TODO: ver bien tema As de armadura traccionada que podria ser superior o inferior.
            "Value": [self.height.to('cm').magnitude, self.width.to('cm').magnitude, self.c_c.to('cm').magnitude,
                       round(self._A_s.to('cm**2').magnitude,2)],
            "Unit": ["cm", "cm", "cm", "cm²"]
        }
        self._forces = {
            "Design forces": [
                "Axial, positive compression",
                "Shear",
            ],
            "Variable": ["Nu", "Vu"],
            "Value": [round(self._N_u.to('kN').magnitude,2), round(self._V_u.to('kN').magnitude,2) ],
            "Unit": ["kN", "kN"]
        }
        # Example lists
        min_values = [None, None, self._A_v_min]   # Use None for items without a minimum constraint
        max_values = [self._s_max_l, self._s_max_w, None]  # Use None for items without a maximum constraint
        current_values = [self._s_l, self._s_w, self._stirrup_A_v]  # Current values to check

        # Generate check marks based on the range conditions
        checks = [
            '✔️' if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else '❌'
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._data_min_max = {
            'Check': ['Stirrup spacing along length', 'Stirrup spacing along width', 'Minimum shear reinforcement'],
            'Unit': ['cm', 'cm', 'cm²/m'],
            'Value': [round(self._s_l.to('cm').magnitude,2), round(self._s_w.to('cm').magnitude,2),
            round(self._stirrup_A_v.to('cm**2/m').magnitude,2)],
            'Min.': ["", "", round(self._A_v_min.to('cm**2/m').magnitude,2)],
            'Max.': [round(self._s_max_l.to('cm').magnitude,2), round(self._s_max_w.to('cm').magnitude,2), ""],
            'Ok?': checks
        }
        self._shear_reinforcement = {
            "Shear reinforcement strength": [
                "Number of stirrups",
                "Stirrup diameter",
                "Stirrup spacing",
                "Effective height",
                "Minimum shear reinforcing",
                "Required shear reinforcing",
                "Shear reinforcing",
                "Shear steel strength"
            ],
            "Variable": ["ns", "db", "s", "d", "Av,min","Av,req","Av", "ØVs"],
            "Value": [self._stirrup_n, self._stirrup_d_b.to('mm').magnitude, self._stirrup_s_l.to('cm').magnitude,
                    self.d.to('cm').magnitude, round(self._A_v_min.to('cm**2/m').magnitude,2),
                    round(self._A_v_req.to('cm**2/m').magnitude,2),
                    round(self._stirrup_A_v.to('cm**2/m').magnitude,2),
                    round(self._phi_V_s.to('kN').magnitude,2)],
            "Unit": ["", "mm", "cm", "cm", "cm²/m","cm²/m", "cm²/m","kN"]
        }
        check_max = '✔️' if self._max_shear_ok else '❌'
        check_FU = '✔️' if self._FUv < 1 else '❌'
        self._shear_concrete = {
            "Concrete strength": [
                "Effective shear area",
                "Longitudinal reinforcement ratio",
                "Size modification factor",
                "Axial stress",
                "Concrete effective shear stress",
                "Concrete strength",
                "Maximum shear strength",
                "Total shear strength", 
                "Max shear check",
                "Demand Capacity Ratio"
            ],
            "Variable": ["Acv", "ρw", "λs", "σNu", "kc", "ØVc", "ØVmax", "ØVn", "" ,"DCR"],
            "Value": [round(self._A_cv.to('cm**2').magnitude,2),
                      round(self._rho_w.magnitude,5),round(self._lambda_s,3),
                      round(self._sigma_u.to('MPa').magnitude,2), round(self._k_c_min.to('MPa').magnitude,2),
                      round(self._phi_V_c.to('kN').magnitude,2), round(self._phi_V_max.to('kN').magnitude,2), 
                      round(self._phi_V_n.to('kN').magnitude,2), check_max, round(self._FUv,2)],
            "Unit": ["cm²", "", "", "MPa", "MPa", "kN", "kN", "kN", "", check_FU]
        }

    # Beam results for Jupyter Notebook
    @property
    def data(self) -> None:
        markdown_content = f"Beam {self.label}, $b$={self.width.to('cm')}"\
                         f", $h$={self.height.to('cm')}, $c_{{c}}$={self.c_c.to('cm')}, \
                            Concrete {self.concrete.name}, Rebar {self.steel_bar.name}."
        # Display the combined content
        display(Markdown(markdown_content))  # type: ignore

        return None
    
    @property
    def properties(self) -> None:
        markdown_content = f"Beam {self.label}, $b$={self.width.to('cm')}"\
                         f", $h$={self.height.to('cm')}, $c_{{c}}$={self.c_c.to('cm')}, \
                            Concrete {self.concrete.name}, Rebar {self.steel_bar.name}."
        self._md_properties = markdown_content

        return None

    @property
    def shear_results(self) -> None:
        if not self._stirrup_A_v:
            warnings.warn("Shear design has not been performed yet. Call check_shear or "
                          "design_shear first.", UserWarning)
            self._md_shear_results = "Shear results are not available."
            return None

        # Create FUFormatter instance and format FU value
        formatter = Formatter()
        formatted_DCR = formatter.DCR(self._FUv )
        rebar_v = f"{int(self._stirrup_n)}eØ{self._stirrup_d_b.to('mm').magnitude}/"\
                f"{self._stirrup_s_l.to('cm').magnitude} cm"
        # Print results
        markdown_content = f"Shear reinforcing {rebar_v}, $A_v$={self._stirrup_A_v.to('cm**2/m')}"\
                         f", $V_u$={self._V_u.to('kN')}, $\\phi V_n$={self._phi_V_n.to('kN')} → {formatted_DCR}"
        self._md_shear_results = markdown_content

        return None
    
    # Beam results for Jupyter Notebook
    @property
    def results(self) -> None:
        # Ensure that both properties and shear results are available
        if not hasattr(self, '_md_properties'):
            self.properties  # This will generate _md_properties
        if not hasattr(self, '_md_shear_results'):
            self.shear_results  # This will generate _md_shear_results
        # Combine the markdown content for properties and shear results
        markdown_content = f"{self._md_properties}\n\n{self._md_shear_results}"
        
        # Display the combined content
        display(Markdown(markdown_content))  # type: ignore

        return None
    
    @property
    def shear_results_detailed(self) -> None:
        if not self._stirrup_A_v:
            warnings.warn("Shear design has not been performed yet. Call check_shear or "
                          "design_shear first.", UserWarning)
            self._md_shear_results = "Shear results are not available."
            return None
        
        # Create a TablePrinter instance and display tables
        materials_printer = TablePrinter("MATERIALS")
        materials_printer.print_table_data(self._materials, headers='keys')
        geometry_printer = TablePrinter("GEOMETRY")
        geometry_printer.print_table_data(self._geometry, headers='keys')
        forces_printer = TablePrinter("FORCES")
        forces_printer.print_table_data(self._forces, headers='keys')
        steel_printer = TablePrinter("SHEAR STRENGTH")
        steel_printer.print_table_data(self._shear_reinforcement, headers='keys')
        min_max_printer = TablePrinter("MAX AND MIN LIMIT CHECKS")
        min_max_printer.print_table_data(self._data_min_max, headers='keys')
        concrete_printer = TablePrinter("CONCRETE STRENGTH")
        concrete_printer.print_table_data(self._shear_concrete, headers='keys')

    @property
    def shear_results_detailed_doc(self) -> None:
        if not self._stirrup_A_v:
            warnings.warn("Shear design has not been performed yet. Call check_shear or "
                          "design_shear first.", UserWarning)
            self._md_shear_results = "Shear results are not available."
            return None
        
        # Convert output Dicts into DataFrames
        df_materials = pd.DataFrame(self._materials)
        df_geometry = pd.DataFrame(self._geometry)
        df_forces = pd.DataFrame(self._forces)
        df_shear_reinforcement = pd.DataFrame(self._shear_reinforcement)
        df_data_min_max = pd.DataFrame(self._data_min_max)
        df_shear_concrete = pd.DataFrame(self._shear_concrete)


        # Create a document builder instance
        doc_builder = DocumentBuilder(title='Concrete beam shear check')

        # Add first section and table
        doc_builder.add_heading('Concrete beam shear check', level=1)
        doc_builder.add_text(f'Design code: {self.concrete.design_code}')
        doc_builder.add_heading('Materials', level=2)
        doc_builder.add_table_data(df_materials)
        doc_builder.add_table_data(df_geometry)
        doc_builder.add_table_data(df_forces)

        # Add second section and another table (can use different data)
        doc_builder.add_heading('Limit checks', level=2)
        doc_builder.add_table_min_max(df_data_min_max)
        doc_builder.add_heading('Design checks', level=2)
        doc_builder.add_table_data(df_shear_reinforcement)
        doc_builder.add_table_data(df_shear_concrete)

        # Save the Word doc
        doc_builder.save(f"Concrete beam shear check {self.concrete.design_code}.docx")

def flexure() -> None:
    concrete = Concrete_ACI_318_19(name="H30",f_c=30*MPa) 
    steelBar = SteelBar(name="ADN 420", f_y=420*MPa) 
    section = RectangularConcreteBeam(
        label="B-12x24",
        concrete=concrete,
        steel_bar=steelBar,
        width=400 * mm,  
        height=500 * mm,  
    )
    debug(f"Nombre de la sección: {section.label}")
    resultados=section.design_flexure(500*kN*m)  
    debug(resultados)


def shear_ACI_metric() -> None:
    concrete= Concrete_ACI_318_19(name="C25",f_c=25*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa)
    custom_settings = {'clear_cover': 2.5*cm, 'stirrup_diameter_ini':8*mm,
                       'longitudinal_diameter_ini': 16*mm}
    section = RectangularConcreteBeam(label="V-10x16",
                                      concrete=concrete,steel_bar=steelBar,width=20*cm, height=50*cm,
                                       settings=custom_settings)
    # debug(section.settings.default_settings)
    # debug(section.get_settings)
    f = Forces(V_z=100*kN)
    section.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=20*cm) 
    results=section.check_shear(f)
    # debug(section._d, section.c_c, section._stirrup_d_b, section._long_d_b)
    # print(results)
    # results = section.design_shear(f, A_s)
    # print(results)
    # print(section.shear_design_results)
    # print(section._id)
    # print(section.shear_results)
    section.shear_results_detailed
    # section.shear_results_detailed_doc

def shear_ACI_imperial() -> None:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000*psi)  
    steelBar = SteelBar(name="ADN 420", f_y=60*ksi)  
    custom_settings = {'clear_cover': 1.5*inch, 'stirrup_diameter_ini':0.5*inch,
                       'longitudinal_diameter_ini': 1*inch} 
    section = RectangularConcreteBeam(
        label="V-10x16",
        concrete=concrete,
        steel_bar=steelBar,
        width=10*inch,  
        height=16*inch,
        settings=custom_settings  
    )
    # debug(section.settings.default_settings)
    # debug(section.get_settings)
    f = Forces(V_z=37.727*kip)
    section.set_transverse_rebar(n_stirrups=1, d_b=0.5*inch, s_l=6*inch) 
    section.check_shear(f, A_s=0.847*inch**2)
    section.design_shear(f, A_s=0.847*inch**2)
    section.shear_results_detailed  
    section.shear_results_detailed_doc

def rebar() -> None:
    concrete= Concrete_ACI_318_19(name="H30",f_c=30*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa) 
    section = RectangularConcreteBeam(
        label="V 20x50",
        concrete=concrete,
        steel_bar=steelBar,
        width=20*cm,  
        height=50*cm,  
    )
    section.settings.load_aci_318_19_settings()
    section.c_c = 30*mm
    as_req = 5 * cm**2

    beam_rebar = Rebar(section)
    long_rebar_df = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=as_req)
    print(long_rebar_df)
    best_design = beam_rebar.longitudinal_rebar_design
    debug(best_design)
    # print(long_rebar_df.iloc[0]['total_as'])
    # A_v_req = 8.045*cm**2/m
    # V_s_req = 108.602*kN
    # trans_rebar = beam_rebar.beam_transverse_rebar_ACI_318_19(A_v_req=A_v_req, V_s_req=V_s_req)
    # print(trans_rebar)

if __name__ == "__main__":
    shear_ACI_imperial()
    # shear_ACI_metric()
    # rebar()