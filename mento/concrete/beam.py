from devtools import debug
from dataclasses import dataclass
from mento.concrete.rectangular import RectangularConcreteSection
from mento.material import Concrete, SteelBar, Concrete_ACI_318_19 
from mento.settings import Settings
from mento.rebar import Rebar
from mento import MPa, ksi, psi, kip, mm, inch, kN, m, cm, Quantity
from mento.results import Formatter
from mento.forces import Forces
from IPython.display import Markdown, display
import numpy as np
import math
from typing import Optional, Dict, Any, cast
from pint.facets.plain import PlainQuantity

@dataclass
class RectangularConcreteBeam(RectangularConcreteSection):
    label: Optional[str] = None

    def __init__(self, label: Optional[str], concrete: Concrete, steel_bar: SteelBar, 
                 width: PlainQuantity, height: PlainQuantity, settings: Optional[Settings] = None):   
        super().__init__(concrete, steel_bar, width, height, settings)
        self.label = label
        self.shear_design_results: Optional[Dict[str, Any]] = None
        self.transverse_rebar: Dict[str, Any]

        if settings:
            self.settings.update(settings.settings)  # Update with any provided settings
    
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
        A_s_calc=0.85*f_c*self._width*self.d/self.steel_bar.f_y*(1-np.sqrt(1-2*R_n/(0.85*f_c)))

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
            rho=0.85*beta_1*f_c/self.steel_bar.f_y*(0.003/(self.steel_bar.epsilon_y+0.006))
            M_n_t=rho*self.steel_bar.f_y*(self.d-0.59*rho*self.steel_bar.f_y*self.d/f_c)*self.width*self.d
            M_n_prima=M_u/phi-M_n_t
            c_t=0.003*self.d/(self.steel_bar.epsilon_y+0.006)
            # HAY QUE VER DONDE ESTA DEFINIDO EL d_prima por ahora asumo
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

    def check_shear_ACI_318_19(self, Force:Forces, A_s:PlainQuantity, d_b:PlainQuantity, s:PlainQuantity, 
                               n_legs:int) -> Dict[str, Any]:
        N_u = Force.N_x
        V_u = Force.V_z
        f_c=self.concrete.f_c
        self.settings.load_aci_318_19_settings()
        phi_v = self.settings.get_setting('phi_v')
        self.lambda_factor = self.settings.get_setting('lambda')
        f_yt = self.steel_bar.f_yt

        # Minimum shear reinforcement calculation
        # 'Minimum reinforcement should be placed if the factored shear Vu 
        # is greater than half the shear capacity of the concrete,
        # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
        # Rebar needed, V_u > φ_v*V_c/2 for Imperial system
        A_v_min = max((0.75 * math.sqrt(f_c / psi) * psi/ f_yt) * self.width , (50 * psi/f_yt) * self.width)  
        
        # Shear reinforcement calculations
        A_db = (d_b ** 2) * math.pi / 4  # Area of one stirrup leg
        A_vs = n_legs * A_db  # Total area of stirrups
        A_v = A_vs / s  # Stirrup area per unit length

        V_s: PlainQuantity = A_v * f_yt * self.d  # Shear contribution of reinforcement
        phi_V_s = phi_v * V_s  # Reduced shear contribution of reinforcement

        # Effective shear area and longitudinal reinforcement ratio
        self.A_cv = self.width * self.d  # Effective shear area
        A_g = self.width * self.height  # Gross area
        rho_w = A_s / self.A_cv  # Longitudinal reinforcement ratio
        
        # Size modification factor for Imperial system
        lambda_s = math.sqrt(2 / (1 + self.d / (10*inch)))

        # Concrete shear strength calculation
        sigma_Nu = min(N_u / (6 * A_g), 0.05 * f_c)  # Axial stress influence

        # Concrete shear capacity depending on whether min rebar is present or not
        if A_v < A_v_min:
            k_c_min_rebar = 8 * lambda_s * self.lambda_factor * rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + sigma_Nu
        else:           
            k_c_min_rebar = max(2 * self.lambda_factor * math.sqrt(f_c / psi) * psi + sigma_Nu,
                            8 * self.lambda_factor * rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + sigma_Nu)

        V_cmin = 0*kip  
        V_cmax = self.calculate_V_cmax()
        
        # Calculate actual concrete shear strength
        self.V_c = min(V_cmax, max(V_cmin, k_c_min_rebar * self.A_cv))
        phi_V_c: PlainQuantity = phi_v * self.V_c  # Reduced concrete shear strength
        
        # Maximum total shear capacity
        V_max = self.calculate_V_max()
        phi_V_max: PlainQuantity = phi_v * V_max  # Reduced maximum shear capacity

        if V_u < phi_V_c/2:
            A_v_min = 0*inch
            max_shear_ok = True
        elif phi_V_c/2 < V_u < phi_V_max:
            A_v_min = self.calculate_A_v_min()
            max_shear_ok = True
        else:
            max_shear_ok = False 

        # Total shear strength
        phi_V_n: PlainQuantity = phi_v * (self.V_c + V_s)  # Total reduced shear strength (concrete + rebar)

        # Required shear reinforcing nominal strength
        V_s_req: PlainQuantity = V_u-phi_V_c
        # Required shear reinforcing
        A_v_req: PlainQuantity = max(V_s_req/(phi_v*f_yt*self.d), A_v_min)

        # Check results
        return {
            'A_v_min': A_v_min.to('cm ** 2 / m'),  # Minimum shear reinforcement area
            'A_v_req': A_v_req.to('cm ** 2 / m'), # Required shear reinforcing area
            'A_v': A_v.to('cm ** 2 / m'),  # Provided stirrup reinforcement per unit length
            'phi_V_c': phi_V_c.to('kN'),  # Concrete contribution to shear capacity
            'phi_V_s': phi_V_s.to('kN'),  # Reinforcement contribution to shear capacity
            'phi_V_n': phi_V_n.to('kN'),  # Total shear capacity
            'phi_V_max': phi_V_max.to('kN'),  # Maximum shear capacity
            'shear_ok': V_u <= phi_V_n,  # Check if applied shear is within total capacity
            'max_shear_ok': max_shear_ok,  # Check if applied shear is within max shear capacity
            "FUv" : (V_u.to('kN') / phi_V_n.to('kN')) 
        }
    
    def design_shear_ACI_318_19(self, Force:Forces, A_s:PlainQuantity) -> Dict[str, Any]:
        N_u = Force.N_x
        V_u = Force.V_z
        f_c = self.concrete.f_c
        f_yt = self.steel_bar.f_yt
        self.settings.load_aci_318_19_settings()
        phi_v = self.settings.get_setting('phi_v')
        self.lambda_factor = self.settings.get_setting('lambda')
        # Inputs into class
        self.V_u = V_u
        self.N_u = N_u
        self.A_s = A_s 
        # Effective shear area and longitudinal reinforcement ratio
        self.A_cv = self.width * self.d  # Effective shear area
        A_g = self.A_x # Gross area
        rho_w = A_s / self.A_cv  # Longitudinal reinforcement ratio
        
        # Concrete shear strength calculation
        sigma_Nu = min(N_u / (6 * A_g), 0.05 * f_c)  # Axial stress influence

        # Concrete shear capacity assuming that the beam is provided with minimum shear rebar:
        k_c_min_rebar = max(2 * self.lambda_factor * math.sqrt(f_c / psi) * psi + sigma_Nu,
                            8 * self.lambda_factor * rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + sigma_Nu)

        V_cmin = 0*kip 
        V_cmax = self.calculate_V_cmax()
        
        # Calculate actual concrete shear strength
        self.V_c = min(V_cmax, max(V_cmin, k_c_min_rebar * self.A_cv))
        phi_V_c: PlainQuantity = phi_v * self.V_c  # Reduced concrete shear strength
        
        # Maximum total shear capacity for Imperial system
        V_max = self.calculate_V_max()
        phi_V_max: PlainQuantity = phi_v * V_max  # Reduced maximum shear capacity

        # Minimum shear reinforcement calculation
        # Minimum reinforcement should be placed if the factored shear Vu 
        # is greater than half the shear capacity of the concrete, 
        # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
        if V_u < phi_V_c/2:
            A_v_min = 0*inch
            max_shear_ok = True
        elif phi_V_c/2 < V_u < phi_V_max:
            A_v_min = self.calculate_A_v_min()
            max_shear_ok = True
        else:
            max_shear_ok = False 

        # Shear reinforcement calculations
        # Required shear reinforcing nominal strength
        V_s_req: PlainQuantity = V_u-phi_V_c
        # Required shear reinforcing
        A_v_req: PlainQuantity = max(V_s_req/(phi_v*self.steel_bar.f_yt*self.d), A_v_min)

        A_v_design: dict[str,PlainQuantity] = Rebar(self).beam_transverse_rebar(A_v_req, V_s_req)
        # Store the final rebar design results as a property of the Beam instance
        self.transverse_rebar = {
            'n_stirrups': A_v_design['n_stirrups'],
            'd_b': A_v_design['d_b'].to('mm'),
            'spacing': A_v_design['spacing'].to('cm'),
            'A_v': A_v_design['A_v'].to('cm ** 2 / m')
        }

        A_v: PlainQuantity =  A_v_design['A_v']
        V_s = A_v * f_yt * self.d  # Shear contribution of reinforcement
        phi_V_s: PlainQuantity = phi_v * V_s  # Reduced shear contribution of reinforcement
        # Total shear strength
        phi_V_n: PlainQuantity = phi_v * (self.V_c + V_s)  # Total reduced shear strength (concrete + rebar)

        # Design results
        self.shear_design_results = {
            'A_v_min': A_v_min.to('cm ** 2 / m'),  # Minimum shear reinforcement area
            'A_v_req': A_v_req.to('cm ** 2 / m'), # Required shear reinforcing area
            'A_v': A_v.to('cm ** 2 / m'),  # Provided stirrup reinforcement per unit length
            'phi_V_c': phi_V_c.to('kN'),  # Concrete contribution to shear capacity
            'phi_V_s': phi_V_s.to('kN'),  # Reinforcement contribution to shear capacity
            'phi_V_n': phi_V_n.to('kN'),  # Total shear capacity
            'phi_V_max': phi_V_max.to('kN'),  # Maximum shear capacity
            'shear_ok': V_u <= phi_V_n,  # Check if applied shear is within total capacity
            'max_shear_ok': max_shear_ok,  # Check if applied shear is within max shear capacity
            "FUv" : (V_u.to('kN') / phi_V_n.to('kN')) 
        }
        return self.shear_design_results
    
    def check_shear_EN_1992(self, Force:Forces, A_s:PlainQuantity) -> None:
        return None
 
    def design_shear_EN_1992(self, Force:Forces, A_s:PlainQuantity) -> None:
        return None

    def design_shear_EHE_08(self, Force:Forces, A_s:PlainQuantity) -> None:
        return None

    def check_shear_EHE_08(self, Force:Forces, A_s:PlainQuantity) -> None:
        return None
  
    # Factory method to select the shear design method
    def design_shear(self, Force: Forces, A_s: PlainQuantity) -> Dict[str, Any]:
        if self.concrete.design_code=="ACI 318-19":
            return self.design_shear_ACI_318_19(Force, A_s)
        # elif self.concrete.design_code=="EN 1992":
        #     return self.design_shear_EN_1992(V_u, N_u, A_s)
        # elif self.concrete.design_code=="EHE-08":
        #     return self.design_shear_EHE_08(V_u, N_u, A_s)
        else:
            raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")
    
    # Factory method to select the shear check method
    def check_shear(self, Force: Forces, A_s: PlainQuantity, d_b:PlainQuantity, 
                    s:PlainQuantity, n_legs:int) -> Dict[str, PlainQuantity]:
        if self.concrete.design_code=="ACI 318-19":
            return self.check_shear_ACI_318_19(Force, A_s, d_b, s, n_legs)
        # elif self.concrete.design_code=="EN 1992":
        #     return self.check_shear_EN_1992(V_u, N_u, A_s, d_b, s, n_legs)
        # elif self.concrete.design_code=="EHE-08":
        #     return self.check_shear_EHE_08(V_u, N_u, A_s, d_b, s, n_legs)
        else:
            raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")

    def calculate_A_v_min(self) -> PlainQuantity:
        # Formula for A_v_min for Imperial system
        A_v_min = max((0.75 * math.sqrt(self.concrete.f_c / psi) * psi / self.steel_bar.f_yt) * self.width,
                      (50 * psi / self.steel_bar.f_yt) * self.width)
        return A_v_min    
    
    def calculate_V_max(self) -> PlainQuantity:
        "Formula for maximum total shear capacity (V_max)"
        V_max = self.V_c + (8 * self.lambda_factor * math.sqrt(self.concrete.f_c / psi) * psi) * self.A_cv
        return V_max
    
    def calculate_V_cmax(self) -> PlainQuantity:
        "Maximum concrete shear strength"
        V_cmax = (5 * self.lambda_factor * math.sqrt(self.concrete.f_c / psi) * psi) * self.A_cv
        return V_cmax 

    
    # Beam results forJjupyter Notebook
    @property
    def shear_results(self) -> str:
        if not self.shear_design_results:
            raise ValueError("Shear design has not been performed yet. Call design_shear first.")

        FU_v = self.shear_design_results['FUv']
        # Create FUFormatter instance and format FU value
        formatter = Formatter()
        formatted_FU = formatter.FU(FU_v)
        rebar_v = f"{self.transverse_rebar['n_stirrups']}eØ\
            {self.transverse_rebar['d_b']}/{self.transverse_rebar['spacing']}"
        # Print results
        markdown_content = f"Armadura transversal {rebar_v}, $A_v$={round(self.transverse_rebar['A_v'],2)}" \
                         f" cm²/m, $V_u$={round(self.V_u,0)},\
                              $\\phi V_n$={round(self.shear_design_results['phi_V_n'],0)} → {formatted_FU}"
        # Display the content
        display(Markdown(markdown_content)) #type: ignore

        return markdown_content


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


def shear() -> None:
    concrete= Concrete_ACI_318_19(name="C4",f_c=4000*psi) 
    steelBar= SteelBar(name="ADN 420", f_y=60*ksi) 
    section = RectangularConcreteBeam(label="V-10x16",
                                      concrete=concrete,steel_bar=steelBar,width=10*inch, height=16*inch)
    section.cc = 1.5*inch
    debug(section)
    section.stirrup_d_b = 0.5*inch
    f = Forces(V_z=37.727*kip, N_x=0*kip)
    debug(f.get_forces()) 
    A_s=0.847*inch**2 
    results=section.check_shear(f, A_s, d_b=12*mm, s=6*inch, n_legs=2) 
    debug(results)
    section.design_shear(f, A_s) 
    debug(section.shear_results)
    debug(section.shear_design_results)
    debug(section.transverse_rebar)

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
    section.cc = 30*mm
    section.stirrup_d_b = 8*mm
    as_nec = 5 * cm**2
    best_combination = Rebar(section).beam_longitudinal_rebar_ACI_318_19(A_s_req=as_nec)
    debug(best_combination)


if __name__ == "__main__":
    shear()