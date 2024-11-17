from devtools import debug
from dataclasses import dataclass
from IPython.display import Markdown, display
from typing import Optional, Dict, Any, cast, Union
from pint.facets.plain import PlainQuantity
import numpy as np
import pandas as pd
from pandas import DataFrame
import math
import warnings

from mento.rectangular import RectangularSection
from mento.material import Concrete, SteelBar, Concrete_ACI_318_19, Concrete_EHE_08
from mento.rebar import Rebar
from mento import MPa, ksi, psi, kip, mm, inch, kN, m, cm, kNm, ft
from mento.results import Formatter, TablePrinter, DocumentBuilder
from mento.forces import Forces  
from mento.node import Node

@dataclass
class RectangularBeam(RectangularSection):
    label: Optional[str] = None

    def __init__(self, label: Optional[str], concrete: Concrete, steel_bar: SteelBar, 
                 width: PlainQuantity, height: PlainQuantity, settings: Optional[Dict[str, Any]] = None):   
        super().__init__(concrete, steel_bar, width, height, settings)
        if settings:
            self.settings.update(settings)  # Update with any provided settings
        # Initialize default concrete beam attributes
        self.initialize_code_attributes()
        self.label = label
        self.shear_design_results: DataFrame = None
        self._stirrup_s_l: PlainQuantity = 0*cm
        self._stirrup_n: int = 0
        self._A_v: PlainQuantity = 0*cm**2/m
        self._A_s_req_bot: PlainQuantity = 0*cm**2
        self._A_s: PlainQuantity = 0*cm**2
        self._A_v_req: PlainQuantity = 0*cm**2/m
        self._FUv: float = 0
        self._s_l = self._stirrup_s_l
        self._s_w: PlainQuantity = 0*cm
        self._s_max_l: PlainQuantity = 0*cm
        self._s_max_w: PlainQuantity = 0*cm
        self._shear_checked = False  # Tracks if shear check or design has been done
    
    @property
    def d(self) -> PlainQuantity:
        "Effective height."
        return self._d
    
    def initialize_code_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_ACI_318_19):
            self._initialize_aci_318_attributes()
        elif isinstance(self.concrete, Concrete_EHE_08):
            self._initialize_ehe_08_attributes()

    def _initialize_aci_318_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_ACI_318_19):
            self._phi_V_n: PlainQuantity = 0*kN
            self._phi_V_s: PlainQuantity = 0*kN
            self._phi_V_c: PlainQuantity = 0*kN
            self._phi_V_max: PlainQuantity = 0*kN
            self._V_u: PlainQuantity = 0*kN
            self._N_u: PlainQuantity = 0*kN
            self._A_cv: PlainQuantity = 0*cm**2
            self._k_c_min: PlainQuantity = 0*MPa
            self._sigma_Nu: PlainQuantity = 0*MPa

    def _initialize_ehe_08_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_EHE_08):
            # Set EHE-specific attributes
            self._f_ctm = self.concrete.f_ctm
            self._f_yk = self.steel_bar.f_y
            self._f_ck = self.concrete.f_ck
            self._f_cd = self.concrete.f_cd
            self._f_ctk = self.concrete.f_ctk
            self._f_ctd = self.concrete.f_ctd
            self._f_1cd = self.concrete.f_1cd
            self._f_yda: PlainQuantity = 0*MPa
            self._V_rd_1: PlainQuantity = 0*kN
            self._V_rd_2: PlainQuantity = 0*kN
            self._N_rd: PlainQuantity = 0*kN
            self._sigma_cd: PlainQuantity = 0*MPa
            self._V_cu: PlainQuantity = 0*kN
            self._V_su: PlainQuantity = 0*kN
            self._V_u2: PlainQuantity = 0*kN
            self._K_value:float = 0
            self._rho_l: float = 0

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
        self._longitudinal_n1_b = n1
        self._longitudinal_d_b1_b = d_b1
        self._longitudinal_n2_b = n2
        self._longitudinal_d_b2_b = d_b2 if n2 > 0 else None
        self._longitudinal_n3_b = n3
        self._longitudinal_d_b3_b = d_b3 if n3 > 0 else None
        self._longitudinal_n4_b = n4
        self._longitudinal_d_b4_b = d_b4 if n4 > 0 else None
        self._total_as_b = total_as
        self._total_bars_b = total_bars
        self._clear_spacing_b = clear_spacing

# ======== ACI 318-19 methods =========

    def __maximum_flexural_reinforcement_ratio(self) -> float:
        if self.concrete.design_code=="ACI 318-19":
            concrete_aci = cast(Concrete_ACI_318_19, self.concrete)  # Cast to the specific subclass
            # Determination of maximum reinforcement ratio
            epsilon_min_rebar_ACI_318_19=self.steel_bar.epsilon_y+concrete_aci.epsilon_c # TODO: REVISAR
            rho_max=0.85*concrete_aci.beta_1*self.concrete.f_c/self.steel_bar.f_y*(concrete_aci.epsilon_c/(concrete_aci.epsilon_c+epsilon_min_rebar_ACI_318_19))
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

    def design_flexure_ACI_318_19(self, Force:Forces)-> Dict[str, Any]:
        # Si el argumento Forces no es una lista, lo convierte en una lista de un solo elemento
        if not isinstance(Force, list):
            Force = [Force]

        M_u_min=0*kNm # Negativo es el que hace que armemos arriba
        M_u_max=0*kNm # Positlivo es el que hace que armemos abajo
        for item in Force:
            if item._M_y > M_u_max:
                M_u_max=item._M_y
            if item._M_y < M_u_min:
                M_u_min=item._M_y

        self.settings.load_aci_318_19_settings()
        phi = self.settings.get_setting('phi_t')
        setting_flexural_min_reduction = self.settings.get_setting('flexural_min_reduction')
        concrete_properties=self.concrete.get_properties()
        f_c=concrete_properties['f_c']
        beta_1=concrete_properties['beta_1']
        rebar_properties=self.steel_bar.get_properties()
        f_y=rebar_properties['f_y']

        rec_mec=self.c_c + self._stirrup_d_b +1*inch # Asumption, this is the main difference between design and check.
        d_prima=self.c_c + self._stirrup_d_b +1*inch # Asumption, this is the main difference between design and check.


        #EMPEZAMOS A ITERAR HASTA QUE SEPAMOS EL VERDADERO REC MEC 
        tol=0.01*cm
        Err=2*tol
        Contador=0
        while Err>=tol:
            Contador+=1
            d=self._height-rec_mec
            b=self._width
            
            # Determination of minimum reinforcement
            A_s_min=max((3*np.sqrt(f_c / psi)*psi/self.steel_bar.f_y*self.d*self._width),
                        (200*psi/self.steel_bar.f_y*self.d*self._width))

            # Determination of maximum reinforcement
            rho_max=self.__maximum_flexural_reinforcement_ratio()
            A_s_max=rho_max*self.d*self._width

            # Determination of required reinforcement at bottom
            R_n_bot=M_u_max/(phi*b*d**2)
            A_s_calc_bot=0.85*f_c*b*d/f_y*(1-np.sqrt(1-2*R_n_bot/(0.85*f_c)))

            if A_s_calc_bot>A_s_min:
                self._A_s_calculated_bot=A_s_calc_bot
            elif  4*A_s_calc_bot/3 > A_s_min:
                self._A_s_calculated_bot=A_s_min
            else: 
                if setting_flexural_min_reduction=='True':
                    self._A_s_calculated_bot=4*A_s_calc_bot/3
                else:
                    self._A_s_calculated_bot=A_s_min
            if self._A_s_calculated_bot <= A_s_max: 
                self._A_s_comp_top=0*cm**2
            else:
                rho=0.85*beta_1*f_c/self.steel_bar.f_y*(0.003/(self.steel_bar.epsilon_y+0.006))
                M_n_t=rho*self.steel_bar.f_y*(self.d-0.59*rho*self.steel_bar.f_y*self.d/f_c)*self.width*self.d
                M_n_prima=M_u_max/phi-M_n_t
                c_t=0.003*self.d/(self.steel_bar.epsilon_y+0.006)

                f_s_prima=min(0.003*self.steel_bar.E_s*(1-d_prima/c_t),self.steel_bar.f_y)
                A_s_prima_top=M_n_prima/(f_s_prima*(self.d-d_prima))
                A_s=rho*self.width*self.d+A_s_prima_top
                self._A_s_calculated_bot=A_s
                self._A_s_comp_top=A_s_prima_top
            #Rebar bottom
            section_rebar = Rebar(self)
            self.flexure_design_results = section_rebar.longitudinal_rebar_ACI_318_19(self._A_s_calculated_bot)
            #debug(self.flexure_design_results)
            best_design=section_rebar.longitudinal_rebar_design
            self._d_b1=best_design["d_b1"] # First diameter of first layer
            self._d_b2=best_design["d_b2"] # Second diameter of first layer
            self._d_b3=best_design["d_b3"] # First diameter of second layer
            self._d_b4=best_design["d_b4"] # Second diameter of second layer
            
            self._layers_spacing = self.settings.get_setting('layers_spacing')

            if self._d_b3 is None and self._d_b4 is None:
                # Usa 0 como valor predeterminado para _d_b2 si es None
                d_b2_value = self._d_b2 if self._d_b2 is not None else 0
                heigth_of_bars = max(self._d_b1, d_b2_value) / 2
            else:
                # Usa 0 como valor predeterminado para _d_b2, _d_b3, y _d_b4 si alguno es None
                d_b2_value = self._d_b2 if self._d_b2 is not None else 0
                d_b3_value = self._d_b3 if self._d_b3 is not None else 0
                d_b4_value = self._d_b4 if self._d_b4 is not None else 0
                heigth_of_bars = (max(self._d_b1, d_b2_value) + self._layers_spacing + max(d_b3_value, d_b4_value)) / 2



            #TODO AHORA HAY QUE REHACER TODO LO MISMO PARA EL M_u_min si es negativo pero armando al reves y quedarse con los mayores en top y bottom.
            # TAMBIEN FALTA DISENAR CON REBAR EL TOP



            self.total_as_adopted=best_design["total_as"]
            self.avaiable_spaciong=best_design["clear_spacing"]

            rec_mec_calculo=self.c_c + self._stirrup_d_b + heigth_of_bars
            Err=abs(rec_mec_calculo-rec_mec)
            rec_mec=rec_mec_calculo

        results={
            'As_min_code':A_s_min.to("inch**2"),
            'As_required':self._A_s_calculated_bot.to("inch**2"),
            'As_max':A_s_max.to("inch**2"),
            'As_adopted':self.total_as_adopted_bot.to("inch**2"),
            'As_compression':self._A_s_comp_top.to("inch**2")
        }
        return pd.DataFrame([results], index=[0])

    def design_flexure_EN_1992(self, M_u: float) -> None:
        pass

    def design_flexure_EHE_08(self, M_u: float) -> None:
        pass

    def design_flexure(self,Force:Forces) -> Dict[str, Any]:
        if self.concrete.design_code=="ACI 318-19":
            return self.design_flexure_ACI_318_19(Force)
        # elif self.concrete.design_code=="EN 1992":
        #     return self.design_flexure_EN_1992(M_u)
        # elif self.concrete.design_code=="EHE-08":
        #     return self.design_flexure_EHE_08(M_u)
        else:
            raise ValueError(f"Longitudinal design method not implemented \
                    for concrete type: {type(self.concrete).__name__}")

    def check_flexure_ACI_318_19(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> None:
        pass

    def _set_initial_conditions_aci(self, Force: Forces, A_s: PlainQuantity) -> None:
        self._N_u = Force.N_x
        self._V_u = Force.V_z
        self._A_s = A_s
        self.settings.load_aci_318_19_settings()
        self.phi_v = self.settings.get_setting('phi_v')
        self.lambda_factor = self.settings.get_setting('lambda')
        self.f_yt = self._calculate_f_yt_aci()

    def _calculate_shear_reinforcement_aci(self) -> None:
        d_bs = self._stirrup_d_b
        s_l = self._stirrup_s_l
        n_legs = self._stirrup_n * 2
        A_db = (d_bs ** 2) * math.pi / 4  # Area of one stirrup leg
        A_vs = n_legs * A_db  # Total area of stirrups
        self._A_v = A_vs / s_l  # Stirrup area per unit length
        V_s = self._A_v * self.f_yt * self.d  # Shear contribution of reinforcement
        self._phi_V_s = self.phi_v * V_s  # Reduced shear contribution of reinforcement

    def _calculate_effective_shear_area_aci(self) -> None:
        self._A_cv = self.width * self.d  # Effective shear area
        self._rho_w = self._A_s.to('cm**2') / self._A_cv.to('cm**2')  # Longitudinal reinforcement ratio
        if self.concrete.unit_system == "metric":
            self._lambda_s = math.sqrt(2 / (1 + 0.004*self.d/mm))
        else:
            self._lambda_s = math.sqrt(2 / (1 + self.d / (10 * inch)))

    def _calculate_concrete_shear_strength_aci(self, A_s: PlainQuantity) -> None:
        f_c = self.concrete.f_c
        self._sigma_Nu = min(self._N_u / (6 * self.A_x), 0.05 * f_c)  # Axial stress influence
        if self.concrete.unit_system == "metric":
            V_cmin = 0 * kN
            if self._A_v < self._A_v_min:
                if A_s == 0*cm**2:
                    warnings.warn("Longitudinal rebar As cannot be zero if A_v is less than A_v_min.", UserWarning)
                self._k_c_min = 0.66 * self._lambda_s * self.lambda_factor * self._rho_w ** (1/3)\
                    * math.sqrt(f_c / MPa)* MPa + self._sigma_Nu
            else:
                self._k_c_min = max(0.17 * self.lambda_factor * math.sqrt(f_c / MPa) * MPa + self._sigma_Nu,
                                0.66 * self.lambda_factor * self._rho_w ** (1/3) * math.sqrt(f_c / MPa) * MPa\
                                      + self._sigma_Nu)
        else:
            V_cmin = 0 * kip
            if self._A_v < self._A_v_min:
                self._k_c_min = 8 * self._lambda_s * self.lambda_factor * self._rho_w ** (1/3) * math.sqrt(f_c / psi)\
                      * psi + self._sigma_Nu
            else:
                self._k_c_min = max(2 * self.lambda_factor * math.sqrt(f_c / psi) * psi + self._sigma_Nu,
                                8 * self.lambda_factor * self._rho_w ** (1/3) * math.sqrt(f_c / psi) * psi +\
                                      self._sigma_Nu)
        # Maximum concrete shear strength
        if self.concrete.unit_system == "metric":
            V_cmax = (0.42 * self.lambda_factor * math.sqrt(self.concrete.f_c / MPa) * MPa) * self._A_cv
        else:
            V_cmax = (5 * self.lambda_factor * math.sqrt(self.concrete.f_c / psi) * psi) * self._A_cv
        self.V_c = min(V_cmax, max(V_cmin, self._k_c_min * self._A_cv))
        self._phi_V_c = self.phi_v * self.V_c

    def _calculate_max_shear_capacity_aci(self) -> None:
        "Formula for maximum total shear capacity (V_max)"
        if self.concrete.unit_system == "metric":
            V_max = self.V_c + (0.66 * self.lambda_factor * math.sqrt(self.concrete.f_c / MPa) * MPa) * self._A_cv
        else:
            V_max = self.V_c + (8 * self.lambda_factor * math.sqrt(self.concrete.f_c / psi) * psi) * self._A_cv
        self._phi_V_max = self.phi_v * V_max
        self._max_shear_ok = self._V_u < self._phi_V_max  

    def _check_minimum_reinforcement_requirement_aci(self) -> None:
        if self.concrete.unit_system == "metric":
            if self._V_u < 0.083*self.phi_v*self.lambda_factor*math.sqrt(self.concrete.f_c / MPa) * MPa*self._A_cv:
                self._A_v_req = 0 * cm**2/m
                self._A_v_min = 0 * cm**2/m
                self._max_shear_ok = True
            elif 0.083*self.phi_v*self.lambda_factor*math.sqrt(self.concrete.f_c / MPa) * MPa*self._A_cv\
                < self._V_u < self._phi_V_max:
                self._calculate_A_v_min_ACI(self.concrete.f_c)
                self._max_shear_ok = True
            else:
                self._calculate_A_v_min_ACI(self.concrete.f_c)
                self._max_shear_ok = False
        else:
            if self._V_u < self.phi_v*self.lambda_factor*math.sqrt(self.concrete.f_c / psi) * psi*self._A_cv:
                self._A_v_req = 0 * inch**2/ft
                self._A_v_min = 0 * inch**2/ft
                self._max_shear_ok = True
                # if self._V_u != 0*kN:
                #     self._stirrup_d_b = 0*inch
            elif self.phi_v*self.lambda_factor*math.sqrt(self.concrete.f_c / psi) * psi*self._A_cv\
                < self._V_u < self._phi_V_max:
                self._calculate_A_v_min_ACI(self.concrete.f_c)
                self._max_shear_ok = True
            else:
                self._calculate_A_v_min_ACI(self.concrete.f_c)
                self._max_shear_ok = False

    def _calculate_V_s_req(self) -> None:
        self._V_s_req = self._V_u - self._phi_V_c
        self._A_v_req = max(self._V_s_req / (self.phi_v * self.f_yt * self.d), self._A_v_min)

    def _calculate_total_shear_strength_aci(self) -> None:
        self._phi_V_n = self.phi_v * (self.V_c + self._A_v * self.f_yt * self.d)
        self._FUv = (self._V_u.to('kN') / self._phi_V_n.to('kN'))

    def _calculate_rebar_spacing_aci(self) -> None:
        section_rebar = Rebar(self)
        n_legs_actual = self._stirrup_n * 2  # Ensure legs are even
        self._s_l = self._stirrup_s_l
        self._s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs_actual - 1)
        self._s_max_l, self._s_max_w = section_rebar.calculate_max_spacing_ACI(self._V_u - self._phi_V_c, self._A_cv)
        self._s_l = max(self._s_l, 0 * inch)
        self._s_w = max(self._s_w, 0 * inch)

    def _compile_results_aci_shear(self, force: Forces) -> Dict[str, Any]:
        return {
            'Section Label': self.label,
            'Load Combo': force.label,
            'Av,min': self._A_v_min.to('cm ** 2 / m'),
            'Av,req': self._A_v_req.to('cm ** 2 / m'),
            'Av': self._A_v.to('cm ** 2 / m'),
            'Vu': self._V_u.to('kN'),
            'Nu': self._N_u.to('kN'),
            'ØVc': self._phi_V_c.to('kN'),
            'ØVs': self._phi_V_s.to('kN'),
            'ØVn': self._phi_V_n.to('kN'),
            'ØVmax': self._phi_V_max.to('kN'),
            'Vu<ØVmax': self._max_shear_ok,
            'Vu<ØVn': self._V_u <= self._phi_V_n,
            'DCR': self._FUv
        }

    def _calculate_A_v_min_ACI(self, f_c: PlainQuantity) -> None:
        """Calculate the minimum shear reinforcement based on unit system."""
        # 'Minimum reinforcement should be placed if the factored shear Vu 
        # is greater than half the shear capacity of the concrete,
        # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
        # Rebar needed, V_u > φ_v*V_c/2 for Imperial system
        f_yt = self._calculate_f_yt_aci()
        
        if self.concrete.unit_system == "metric":
            self._A_v_min = max(
                (0.062 * math.sqrt(f_c / MPa) * MPa / f_yt) * self.width,
                (0.35 * MPa / f_yt) * self.width
            )
        else:
            self._A_v_min = max(
                (0.75 * math.sqrt(f_c / psi) * psi / f_yt) * self.width,
                (50 * psi / f_yt) * self.width
            )   
       
    def _calculate_f_yt_aci(self) -> PlainQuantity:
        """Determine the yield strength of steel based on unit system."""
        if self.concrete.unit_system == "metric":
            return min(self.steel_bar.f_y, 400 * MPa)
        else:
            return min(self.steel_bar.f_y, 60 * ksi)

    def check_shear_ACI_318_19(self, force: Forces, A_s: PlainQuantity = 0 * cm ** 2) -> pd.DataFrame:
        # Set the initial variables
        self._set_initial_conditions_aci(force, A_s)

        # Minimum shear reinforcement calculation
        self._calculate_A_v_min_ACI(self.concrete.f_c)
        if self._stirrup_n > 0:
            # Shear reinforcement calculations
            self._calculate_shear_reinforcement_aci()

        # Effective shear area and longitudinal reinforcement ratio
        self._calculate_effective_shear_area_aci()

        # Check if minimum reinforcement is required
        self._check_minimum_reinforcement_requirement_aci()

        # Concrete shear strength calculation
        self._calculate_concrete_shear_strength_aci(A_s)

        # Maximum total shear capacity
        self._calculate_max_shear_capacity_aci()

        # Calculate required shear reinforcement
        self._calculate_V_s_req()

        # Total shear strength
        self._calculate_total_shear_strength_aci()

        # Rebar spacing checks
        self._calculate_rebar_spacing_aci()

        # Check results and return DataFrame
        results = self._compile_results_aci_shear(force)
        self._initialize_dicts_ACI_318_19()
        return pd.DataFrame([results], index=[0])

    def design_shear_ACI_318_19(self, force: Forces, A_s: PlainQuantity = 0 * cm ** 2) -> pd.DataFrame:
        # Set the initial variables
        self._set_initial_conditions_aci(force, A_s)

        # Minimum shear reinforcement calculation
        self._calculate_A_v_min_ACI(self.concrete.f_c)
        # Consider that the beam has minimum reinforcement
        self._A_v = self._A_v_min

        # Effective shear area and longitudinal reinforcement ratio
        self._calculate_effective_shear_area_aci()

        # Concrete shear strength calculation
        self._calculate_concrete_shear_strength_aci(A_s)

        # Maximum total shear capacity
        self._calculate_max_shear_capacity_aci()

        # Check if minimum reinforcement is required
        self._check_minimum_reinforcement_requirement_aci()

        # Calculate required shear reinforcement
        self._calculate_V_s_req()

        # Shear reinforcement calculations
        section_rebar = Rebar(self)
        section_rebar.transverse_rebar(self._A_v_req, self._V_s_req)
        best_design = section_rebar.transverse_rebar_design
        self._s_l = best_design['s_l']
        self._s_w = best_design['s_w']
        self._s_max_l = best_design['s_max_l']
        self._s_max_w = best_design['s_max_w']
        self.set_transverse_rebar(best_design['n_stir'], best_design['d_b'], best_design['s_l'])
        self._A_v = best_design['A_v']
        self._A_v = best_design['A_v']
        V_s = self._A_v * self.f_yt * self.d  # Shear contribution of reinforcement
        self._phi_V_s = self.phi_v * V_s  # Reduced shear contribution of reinforcement

        # Total shear strength
        self._calculate_total_shear_strength_aci()
        self._FUv = (self._V_u.to('kN') / self._phi_V_n.to('kN'))

        # Design results and return DataFrame
        results = self._compile_results_aci_shear(force)
        self._initialize_dicts_ACI_318_19()
        return pd.DataFrame([results], index=[0])

# ======== EHE-08 methods =========

    def _initialize_variables_ehe(self, Force: Forces, A_s: PlainQuantity) -> None:
        if isinstance(self.concrete, Concrete_EHE_08):
            # Set the initial variables
            self._N_rd = Force.N_x
            self._V_rd_1 = Force.V_z  # Consider the same shear at the edge of support and in d
            self._V_rd_2 = Force.V_z  # Consider the same shear at the edge of support and in d
            self._M_rd = Force.M_y
            self._A_s = A_s

            # Load settings for gamma factors
            self.settings.load_ehe_08_settings()
            self._gamma_c = self.settings.get_setting('gamma_c')
            self._gamma_s = self.settings.get_setting('gamma_s')
            self._f_yd = self._f_yk / self._gamma_s
            self._f_yda = min(400 * MPa, self._f_yd)

            # Minimum shear reinforcement calculation
            self._A_v_min = self._f_ctm * self.width / (7.5 * self._f_yda)

            # Compression stress, positive
            self._A_p = 0*cm**2 # No prestressing for now
            self._rho_l = min((A_s + self._A_p) / (self.width * self.d), 0.02)

            # Shear calculation for sections without rebar
            self._xi = min(1 + math.sqrt(200 * mm / self.d), 2)
            self._f_cv = self._f_ck

    def _calculate_V_u1(self) -> PlainQuantity:
        self._alpha = math.radians(90)
        self._theta = math.radians(45)
        self._cot_theta = 1 / math.tan(self._theta)
        self._cot_alpha = 1 / math.tan(self._alpha)
        self._sigma_cd =  self._N_rd / self.A_x # Without compression reinforcement considered 
        self._K_value = self._calculate_axial_coefficient_ehe(self._sigma_cd, self._f_cd)
        return self._K_value * self._f_1cd * self.width * self.d\
              * (self._cot_theta + self._cot_alpha) / (1 + self._cot_theta ** 2)

    def _shear_without_rebar_EHE(self) -> PlainQuantity:
        W_y = self.width * self.height ** 2 / 6
        # Positive of compression
        sigma_t_min = -self._M_rd / W_y + self._N_rd / self.A_x
        self._stirrup_d_b = 0*mm
        self._A_v_min = 0*cm**2/m
        debug(sigma_t_min, -self._f_ctd)

        if sigma_t_min > -self._f_ctd:
            # Total shear capacity without rebar (case without cracking)
            alpha_l = 1
            S_y = self.width * self.height ** 2 / 8
            # Stress must be prestressing force really and not N_rd compression
            self._sigma_cd =  0*MPa if self._A_p == 0*cm**2 else self._N_rd / self.A_x
            return (self.I_y * self.width / S_y) *\
                  math.sqrt(self._f_ctd.to('MPa').magnitude ** 2\
                             + alpha_l * self._sigma_cd.to('MPa').magnitude * self._f_ctd.to('MPa').magnitude)*MPa
        else:
            # Total shear capacity without rebar (case with cracking)
            xi = min(1 + math.sqrt(200 * mm / self.d), 2)
            self._sigma_cd = min(self._N_rd / self.A_x, 0.3 * self._f_cd)
            V_u2_min = (0.075 / self._gamma_c * xi ** (3 / 2) * (self._f_ck / MPa) ** (1 / 2)\
                         + 0.15 * self._sigma_cd / MPa) * MPa * self.width * self.d
            V_u2 = (0.18 / self._gamma_c * xi * (100 * self._rho_l * self._f_ck / MPa) ** (1 / 3)\
                     + 0.15 * self._sigma_cd/ MPa) * MPa * self.width * self.d
            return max(V_u2_min, V_u2)
    
    def _calculate_axial_coefficient_ehe(self, sigma_cd: PlainQuantity, f_cd: PlainQuantity) -> float:
        if sigma_cd == 0:
            return 1
        elif sigma_cd <= 0.25 * f_cd:
            return 1 + sigma_cd / f_cd
        elif sigma_cd < 0.5 * f_cd:
            return 1.25
        else:
            return 2.5 * (1 - sigma_cd / f_cd)
        
    def check_shear_EHE_08(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> DataFrame:
        if isinstance(self.concrete, Concrete_EHE_08):
            # Initialize all the code related variables
            self._initialize_variables_ehe(Force,A_s)
            
            # Maximum shear strength check
            self._V_u1 = self._calculate_V_u1()
            self._max_shear_ok = self._V_rd_1 < self._V_u1
            
            if self._stirrup_n == 0:
                # Calculate V_u2 for no rebar scenario
                self._V_u2 = self._shear_without_rebar_EHE()
                self._V_cu = self._V_u2

            else:
                # Shear reinforcement calculations
                d_bs = self._stirrup_d_b
                s_l = self._stirrup_s_l
                n_legs = self._stirrup_n*2
                self._A_s = A_s
                A_db = (d_bs ** 2) * math.pi / 4  # Area of one stirrup leg
                A_vs = n_legs * A_db  # Total area of stirrups
                self._A_v = A_vs / s_l  # Stirrup area per unit length
                # Total shear strength with rebar (case with rebar and cracking)
                theta_e = self._theta # Cracks angle (assumed 45 degrees)
                cot_theta_e = 1 / math.tan(theta_e)

                if 0.5 <= self._cot_theta < cot_theta_e:
                    beta = (2 * self._cot_theta - 1) / (2 * cot_theta_e - 1)
                elif cot_theta_e <= self._cot_theta <= 2:
                    beta = (self._cot_theta - 2) / (cot_theta_e - 2)
                else:
                    beta = 1  # Default value if condition is not met

                V_cu = (0.15 / self._gamma_c * self._xi * (100 * self._rho_l * self._f_cv / MPa) ** (1 / 3)\
                         + 0.15 * self._sigma_cd / MPa)\
                    * MPa * beta * self.width * self.d
                V_u2_min = (0.075 / self._gamma_c * self._xi ** (3 / 2) * (self._f_cv / MPa) ** (1 / 2)\
                             + 0.15 * self._sigma_cd / MPa)\
                    * MPa * self.width * self.d
                self._V_cu = max(V_cu, V_u2_min)
                
                z = 0.9 * self.d
                self._V_su = z * math.sin(self._alpha) * (self._cot_alpha + self._cot_theta) * self._A_v * self._f_yda
                self._V_u2 = self._V_cu + self._V_su

                # Required shear reinforcing strength
                V_s_req = self._V_rd_2 - self._V_cu

                # Required shear reinforcing area
                self._A_v_req = max(V_s_req / (z * math.sin(self._alpha)\
                                                * (self._cot_alpha + self._cot_theta) * self._f_yda), self._A_v_min)

                # Rebar spacing checks
                section_rebar = Rebar(self)
                n_legs_actual = self._stirrup_n * 2      # Ensure legs are even
                self._s_l = self._stirrup_s_l
                self._s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs_actual - 1) 
                self._s_max_l, self._s_max_w =\
                      section_rebar.calculate_max_spacing_EHE(self._V_rd_2, self._V_u1, self._alpha)

            self._FUv = (self._V_rd_2.to('kN') / self._V_u2.to('kN'))
            # Design results
            results = {
                'Label': self.label, #Beam label
                'Av,min': self._A_v_min.to('cm ** 2 / m'),  # Minimum shear reinforcement area
                'Av,req': self._A_v_req.to('cm ** 2 / m'), # Required shear reinforcing area
                'Av': self._A_v.to('cm ** 2 / m'),  # Provided stirrup reinforcement per unit length
                'Vrd,1': self._V_rd_1.to('kN'), # Max Vu for the design at the support
                'Vrd,2': self._V_rd_2.to('kN'), # Max Vu for the design at d from the support
                'Vcu': self._V_cu.to('kN'),  # Concrete contribution to shear capacity
                'Vsu': self._V_su.to('kN'),  # Reinforcement contribution to shear capacity
                'Vu2': self._V_u2.to('kN'),  # Total shear capacity
                'Vu1': self._V_u1.to('kN'),  # Maximum shear capacity
                'Vrd,1<Vu1': self._max_shear_ok,  # Check if applied shear is within max shear capacity
                'Vrd,2<Vu2': self._V_rd_2 <= self._V_u2,  # Check if applied shear is within total capacity
                "DCR" :  self._FUv
            }
            self._initialize_dicts_EHE_08()
            return pd.DataFrame([results], index=[0])
        else:
            raise ValueError("Concrete type is not compatible with EHE-08 shear check.")
  
    def design_shear_EHE_08(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> None:
        return None
# ======== EN 1992 methods =========
    def check_shear_EN_1992(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> None:
        return None
 
    def design_shear_EN_1992(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> None:
        return None

    # Factory method to select the shear design method
    def design_shear(self, A_s: PlainQuantity = 0*cm**2) -> DataFrame:
        if not self.node or not self.node.forces:
            raise ValueError("No Node or forces list associated with this beam.")
        
        self._shear_results_list = []  # Store individual results for each force
        self._shear_results_detailed_list = {}  # Store detailed results by force ID
        max_A_v_req = 0*cm # Track the maximum A_v_req to identify the limiting case

        for force in self.node.forces:
            if self.concrete.design_code=="ACI 318-19":
                result =  self.design_shear_ACI_318_19(force, A_s)
            # elif self.concrete.design_code=="EN 1992":
            #     result =  self.design_shear_EN_1992(V_u, N_u, A_s)
            # elif self.concrete.design_code=="EHE-08":
            #     result =  self.design_shear_EHE_08(V_u, N_u, A_s)
            else:
                raise ValueError(f"Shear design method not implemented for concrete type:"\
                    f"{type(self.concrete).__name__}")
            self._shear_results_list.append(result)
            self._shear_results_detailed_list[force.id] = {
                'forces': self._forces.copy(),
                'shear_reinforcement': self._shear_reinforcement.copy(),
                'min_max': self._data_min_max.copy(),
                'shear_concrete': self._shear_concrete.copy(),
            }
            # Check if this result is the limiting case
            current_A_v_req = result['Av,req'][0]
            if current_A_v_req > max_A_v_req:
                max_A_v_req = current_A_v_req
                self._limiting_case_shear = result
                self._limiting_case_shear_details = self._shear_results_detailed_list[force.id]

                # Update shear design results for the worst case
                section_rebar = Rebar(self)
                self.shear_design_results = section_rebar.transverse_rebar(self._A_v_req, self._V_s_req)
                self._best_rebar_design = section_rebar.transverse_rebar_design

        # Compile all results into a single DataFrame
        all_results = pd.concat(self._shear_results_list, ignore_index=True)

        # Identify the most limiting case by Av,required
        self.limiting_case_shear = all_results.loc[all_results['Av,req'].idxmax()]  # Select row with highest Av,req
        # Mark shear as checked
        self._shear_checked = True  
        return all_results
    
    # Factory method to select the shear check method
    def check_shear(self, A_s: PlainQuantity = 0*cm**2) -> DataFrame:
        if not self.node or not self.node.forces:
            raise ValueError("No Node or forces list associated with this beam.")
        
        self._shear_results_list = []  # Store individual results for each force
        self._shear_results_detailed_list = {}  # Store detailed results by force ID
        max_dcr = 0  # Track the maximum DCR to identify the limiting case

        for force in self.node.forces:
            # Select the method based on design code
            if self.concrete.design_code=="ACI 318-19":
                result = self.check_shear_ACI_318_19(force, A_s)
            # elif self.concrete.design_code=="EN 1992":
            #     result =  self.check_shear_EN_1992(V_u, N_u, A_s, d_b, s, n_legs)
            elif self.concrete.design_code=="EHE-08":
                result = self.check_shear_EHE_08(force, A_s)
            else:
                raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")  # noqa: E501
            self._shear_results_list.append(result)
            self._shear_results_detailed_list[force.id] = {
                'forces': self._forces.copy(),
                'shear_reinforcement': self._shear_reinforcement.copy(),
                'min_max': self._data_min_max.copy(),
                'shear_concrete': self._shear_concrete.copy(),
            }

            # Check if this result is the limiting case
            current_dcr = result['DCR'][0]
            if current_dcr > max_dcr:
                max_dcr = current_dcr
                self._limiting_case_shear = result
                self._limiting_case_shear_details = self._shear_results_detailed_list[force.id]

        # Compile all results into a single DataFrame
        all_results = pd.concat(self._shear_results_list, ignore_index=True)

        # Identify the most limiting case by Demand-to-Capacity Ratio (DCR) or other criteria
        self.limiting_case_shear = all_results.loc[all_results['DCR'].idxmax()]  # Select row with highest DCR

        self._shear_checked = True  # Mark shear as checked
        return all_results

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
                      round(self.steel_bar.f_y.to('MPa').magnitude,2),round(self.concrete.density.to('kg/m**3').magnitude,1),
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
        # Min max lists
        min_values = [None, None, self._A_v_min]   # Use None for items without a minimum constraint
        max_values = [self._s_max_l, self._s_max_w, None]  # Use None for items without a maximum constraint
        current_values = [self._s_l, self._s_w, self._A_v]  # Current values to check

        # Generate check marks based on the range conditions
        checks = [
            '✔️' if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else '❌'
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._data_min_max = {
            'Check': ['Stirrup spacing along length', 'Stirrup spacing along width', 'Minimum shear reinforcement'],
            'Unit': ['cm', 'cm', 'cm²/m'],
            'Value': [round(self._s_l.to('cm').magnitude,2), round(self._s_w.to('cm').magnitude,2),
            round(self._A_v.to('cm**2/m').magnitude,2)],
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
                "Defined shear reinforcing",
                "Shear steel strength"
            ],
            "Variable": ["ns", "db", "s", "d", "Av,min","Av,req","Av", "ØVs"],
            "Value": [self._stirrup_n, self._stirrup_d_b.to('mm').magnitude, self._stirrup_s_l.to('cm').magnitude,
                    self.d.to('cm').magnitude, round(self._A_v_min.to('cm**2/m').magnitude,2),
                    round(self._A_v_req.to('cm**2/m').magnitude,2),
                    round(self._A_v.to('cm**2/m').magnitude,2),
                    round(self._phi_V_s.to('kN').magnitude,2)],
            "Unit": ["", "mm", "cm", "cm", "cm²/m","cm²/m", "cm²/m","kN"]
        }
        check_max = '✔️' if self._max_shear_ok else '❌'
        check_FU = '✔️' if self._FUv < 1 else '❌'
        self._shear_concrete: dict[str, list[Union[str, float, None]]] = {
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
                      round(self._sigma_Nu.to('MPa').magnitude,2), round(self._k_c_min.to('MPa').magnitude,2),
                      round(self._phi_V_c.to('kN').magnitude,2), round(self._phi_V_max.to('kN').magnitude,2), 
                      round(self._phi_V_n.to('kN').magnitude,2), check_max, round(self._FUv,2)],
            "Unit": ["cm²", "", "", "MPa", "MPa", "kN", "kN", "kN", "", check_FU]
        }
    
    def _initialize_dicts_EHE_08(self) -> None:
        if isinstance(self.concrete, Concrete_EHE_08):
            """Initialize the dictionaries used in check and design methods."""
            self._materials = {
                "Materiales": [
                    "Marca de la sección",
                    "Resistencia característica del hormigón",
                    "Resistencia de diseño del hormigón",
                    "Resistencia a tracción media del hormigón",
                    "Resistencia a tracción de diseño del hormigón",
                    "Resistencia a compresión de la biela",
                    "Tipo de control del hormigón",
                    "Resistencia característica del acero",
                    "Resistencia de diseño del tirante",
                ],
                "Variable": ["","fck", "fcd", "fctm", "fctd", "f1cd", "", "fyk", "fydα"],
                "Valor": [self.label, round(self.concrete.f_ck.to('MPa').magnitude,2),
                        round(self.concrete.f_cd.to('MPa').magnitude,2),
                        round(self.concrete.f_ctm.to('MPa').magnitude,2),
                        round(self.concrete.f_ctd.to('MPa').magnitude,2),
                        round(self.concrete.f_1cd.to('MPa').magnitude,2),
                        "Directo",
                        round(self.steel_bar.f_y.to('MPa').magnitude,2),
                        round(self._f_yda.to('MPa').magnitude,2),
                        ],
                "Unidad": ["", "MPa", "MPa", "MPa","MPa","MPa", "","MPa","MPa"]
            }
            self._geometry = {
                "Geometría": [
                    "Altura de la sección",
                    "Ancho de la sección",
                    "Recubrimiento geométrico",
                    "Armadura longitudinal traccionada",
                ],
                "Variable": ["h", "b", "rgeom", "As"],
                #TODO: ver bien tema As de armadura traccionada que podria ser superior o inferior.
                "Valor": [self.height.to('cm').magnitude, self.width.to('cm').magnitude, self.c_c.to('cm').magnitude,
                        round(self._A_s.to('cm**2').magnitude,2)],
                "Unidad": ["cm", "cm", "cm", "cm²"]
            }
            self._forces = {
                "Fuerzas de diseño": [
                    "Axial, positivo de compresión",
                    "Cortante",
                    "Momento flector"
                ],
                "Variable": ["Nrd", "Vrd,2", "Mrd"],
                "Value": [round(self._N_rd.to('kN').magnitude,2), round(self._V_rd_2.to('kN').magnitude,2), 
                          round(self._M_rd.to('kN*m').magnitude,2)],
                "Unit": ["kN", "kN", "kNm"]
            }
            # Min max lists
            min_values = [None, None, self._A_v_min]   # Use None for items without a minimum constraint
            max_values = [self._s_max_l, self._s_max_w, None]  # Use None for items without a maximum constraint
            current_values = [self._s_l, self._s_w, self._A_v]  # Current values to check

            # Generate check marks based on the range conditions
            checks = [
                '✔️' if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else '❌'
                for curr, min_val, max_val in zip(current_values, min_values, max_values)
            ]
            self._data_min_max = {
                'Check': ['Separación de estribos longitudinal', 
                          'Separación de estribos transversal', 'Armadura transversal mínima'],
                'Unidad': ['cm', 'cm', 'cm²/m'],
                'Valor': [round(self._s_l.to('cm').magnitude,2), round(self._s_w.to('cm').magnitude,2),
                round(self._A_v.to('cm**2/m').magnitude,2)],
                'Min.': ["", "", round(self._A_v_min.to('cm**2/m').magnitude,2)],
                'Max.': [round(self._s_max_l.to('cm').magnitude,2), round(self._s_max_w.to('cm').magnitude,2), ""],
                'Ok?': checks
            }
            self._shear_reinforcement = {
                "Capacidad de la armadura": [
                    "Cantidad de estribos",
                    "Diámetro de estribo",
                    "Separación de estribos",
                    "Altura efectiva",
                    "Armadura transversal mínima",
                    "Armadura transversal requerida",
                    "Armadura transversal dispuesta",
                    "Capacidad de la armadura a corte"
                ],
                "Variable": ["ns", "db", "s", "d", "Av,min","Av,req","Av", "Vus"],
                "Valor": [self._stirrup_n, self._stirrup_d_b.to('mm').magnitude, self._stirrup_s_l.to('cm').magnitude,
                        self.d.to('cm').magnitude, round(self._A_v_min.to('cm**2/m').magnitude,2),
                        round(self._A_v_req.to('cm**2/m').magnitude,2),
                        round(self._A_v.to('cm**2/m').magnitude,2),
                        round(self._V_su.to('kN').magnitude,2)],
                "Unidad": ["", "mm", "cm", "cm", "cm²/m","cm²/m", "cm²/m","kN"]
            }
            check_max = '✔️' if self._max_shear_ok else '❌'
            check_FU = '✔️' if self._FUv < 1 else '❌'
            self._shear_concrete = {
                "Capacidad del hormigón": [
                    "Cuantía de armadura longitudinal",
                    "Tensión de compresión",
                    "Factor que depende de la compresión",
                    "Capacidad del hormigón",
                    "Capacidad máxima de la sección",
                    "Capacidad total de la sección", 
                    "Cortante máximo check",
                    "Factor de Utilización"
                ],
                "Variable": ["ρl", "σcd", "K", "Vcu", "Vu1", "Vu2", "" ,"FU"],
                "Valor": [round(self._rho_l,5),
                        round(self._sigma_cd.to('MPa').magnitude,2), round(self._K_value,2),
                        round(self._V_cu.to('kN').magnitude,2), round(self._V_u1.to('kN').magnitude,2), 
                        round(self._V_u2.to('kN').magnitude,2), check_max, round(self._FUv,2)],
                "Unidad": ["", "MPa", "", "kN", "kN", "kN", "", check_FU]
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
        if not self._shear_checked:
            warnings.warn("Shear design has not been performed yet. Call check_shear or "
                          "design_shear first.", UserWarning)
            self._md_shear_results = "Shear results are not available."
            return None
        # Use limiting case results
        result_data = self._limiting_case_shear_details
        limiting_reinforcement = result_data['shear_reinforcement']
        limiting_forces = result_data['forces']
        limiting_shear_concrete = result_data['shear_concrete']

        # Create FUFormatter instance and format FU value
        formatter = Formatter()
        formatted_DCR = formatter.DCR(limiting_shear_concrete['Value'][-1])
        # rebar_v = f"{int(limiting_reinforcement['ns'])}eØ{limiting_reinforcement['db']} mm/"\
        # f"{limiting_reinforcement['s']} cm"
        rebar_v = f"{int(limiting_reinforcement['Value'][0])}eØ{limiting_reinforcement['Value'][1]}/"\
        f"{limiting_reinforcement['Value'][2]} cm"

        # Print results for the limiting case
        # markdown_content = f"Shear reinforcing {rebar_v}, $A_v$={limiting_reinforcement['Av']} cm²/m"\
        #                  f", $V_u$={limiting_forces['Vu']} kN, $\\phi V_n$={limiting_forces['ØVn']} kN → {formatted_DCR}"  # noqa: E501
        markdown_content = f"Shear reinforcing {rebar_v}, $A_v$={limiting_reinforcement['Value'][6]} cm²/m"\
                         f", $V_u$={limiting_forces['Value'][1]} kN, $\\phi V_n$={limiting_shear_concrete['Value'][7]} kN → {formatted_DCR}"  # noqa: E501

        self._md_shear_results = markdown_content

        # # Create FUFormatter instance and format FU value
        # formatter = Formatter()
        # formatted_DCR = formatter.DCR(self._FUv )
        # rebar_v = f"{int(self._stirrup_n)}eØ{self._stirrup_d_b.to('mm').magnitude}/"\
        #         f"{self._stirrup_s_l.to('cm').magnitude} cm"
        # # Print results
        # markdown_content = f"Shear reinforcing {rebar_v}, $A_v$={self._A_v.to('cm**2/m')}"\
        #                  f", $V_u$={self._V_u.to('kN')}, $\\phi V_n$={self._phi_V_n.to('kN')} → {formatted_DCR}"
        # self._md_shear_results = markdown_content

        return None
    
    # Beam results for Jupyter Notebook
    @property
    def results(self) -> None:
        # Ensure that both properties and shear results are available
        if not hasattr(self, '_md_properties'):
            self.properties  # This will generate _md_properties
        if not hasattr(self, '_md_shear_results'):
            self.shear_results  # This will generate _md_shear_results
        #TODO: Flexure results missing
        # Combine the markdown content for properties and shear results
        markdown_content = f"{self._md_properties}\n\n{self._md_shear_results}"
        
        # Display the combined content
        display(Markdown(markdown_content))  # type: ignore

        return None
    
    def shear_results_detailed(self, force: Optional[Forces] = None) -> None:
        """
        Displays detailed shear results.

        Parameters
        ----------
        forces : Forces, optional
            The specific Forces object to display results for. If None, displays results for the limiting case.
        """
        if not self._shear_checked:
            warnings.warn("Shear check has not been performed yet. Call check_shear or " 
                          "design_shear first.", UserWarning)
            self._md_shear_results = "Shear results are not available."
            return None
         # Determine which results to display (limiting case by default)
        if force:
            force_id = force.id
            if force_id not in self._shear_results_detailed_list:
                raise ValueError(f"No results found for Forces object with ID {force_id}.")
            result_data = self._shear_results_detailed_list[force_id]
        else:
            # Default to limiting case
            result_data = self._limiting_case_shear_details

        # Create a TablePrinter instance and display tables
        materials_printer = TablePrinter("MATERIALS (LIMITING CASE)")
        materials_printer.print_table_data(self._materials, headers='keys')
        geometry_printer = TablePrinter("GEOMETRY (LIMITING CASE)")
        geometry_printer.print_table_data(self._geometry, headers='keys')
        forces_printer = TablePrinter("FORCES")
        forces_printer.print_table_data(result_data['forces'], headers='keys')
        steel_printer = TablePrinter("SHEAR STRENGTH")
        steel_printer.print_table_data(result_data['shear_reinforcement'], headers='keys')
        min_max_printer = TablePrinter("MAX AND MIN LIMIT CHECKS")
        min_max_printer.print_table_data(result_data['min_max'], headers='keys')
        concrete_printer = TablePrinter("CONCRETE STRENGTH")
        concrete_printer.print_table_data(result_data['shear_concrete'], headers='keys')

    def shear_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        """
        Prints detailed shear results in Word.

        Parameters
        ----------
        forces : Forces, optional
            The specific Forces object to display results for. If None, displays results for the limiting case.
        """
        if not self._shear_checked:
            warnings.warn("Shear check has not been performed yet. Call check_shear or " 
                          "design_shear first.", UserWarning)
            self._md_shear_results = "Shear results are not available."
            return None
         # Determine which results to display (limiting case by default)
        if force:
            force_id = force.id
            if force_id not in self._shear_results_detailed_list:
                raise ValueError(f"No results found for Forces object with ID {force_id}.")
            result_data = self._shear_results_detailed_list[force_id]
        else:
            # Default to limiting case
            result_data = self._limiting_case_shear_details
        
        # Convert output Dicts into DataFrames
        df_materials = pd.DataFrame(self._materials)
        df_geometry = pd.DataFrame(self._geometry)
        df_forces = pd.DataFrame(result_data['forces'])
        df_shear_reinforcement = pd.DataFrame(result_data['shear_reinforcement'])
        df_data_min_max = pd.DataFrame(result_data['min_max'])
        df_shear_concrete = pd.DataFrame(result_data['shear_concrete'])


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
    #Example 6.6 CRSI GUIDE
    concrete = Concrete_ACI_318_19(name="fc4000",f_c=4000*psi) 
    steelBar = SteelBar(name="G60", f_y=60000*psi) 
    section = RectangularBeam(
        label="B-12x24",
        concrete=concrete,
        steel_bar=steelBar,
        width=12 * inch,  
        height=24 * inch,   
    )
    debug(f"Nombre de la sección: {section.label}")
    f = Forces(M_y=258.3*kip*ft)
    resultados=section.design_flexure(f)  
    debug(resultados)


def shear_ACI_metric() -> None:
    concrete= Concrete_ACI_318_19(name="C25",f_c=25*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa)
    custom_settings = {'clear_cover': 2.5*cm, 'stirrup_diameter_ini':8*mm,
                       'longitudinal_diameter_ini': 16*mm}
    beam = RectangularBeam(label="101",
                                      concrete=concrete,steel_bar=steelBar,width=20*cm, height=50*cm,
                                       settings=custom_settings)
    f1 = Forces(label='1.4D', V_z=50*kN)
    f2 = Forces(label='1.2D+1.6L', V_z=155*kN)
    f3 = Forces(label='W', V_z=220*kN)
    f4 = Forces(label='S', V_z=80*kN)
    f5 = Forces(label='E', V_z=10*kN)
    Node(section=beam, forces=[f1, f2, f3, f4, f5])
    # beam.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=20*cm) 
    # results, limiting_case = beam.check_shear(A_s=5*cm**2)
    results = beam.design_shear(A_s=5*cm**2)
    print(results)
    beam.shear_results_detailed()
    print(beam.shear_design_results)
    print(beam.results)
    # beam.shear_results_detailed_doc(f2)

def shear_ACI_imperial() -> None:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000*psi)  
    steelBar = SteelBar(name="ADN 420", f_y=60*ksi)  
    custom_settings = {'clear_cover': 1.5*inch, 'stirrup_diameter_ini':0.5*inch,
                       'longitudinal_diameter_ini': 1*inch} 
    beam = RectangularBeam(
        label="102",
        concrete=concrete,
        steel_bar=steelBar,
        width=10*inch,  
        height=16*inch,
        settings=custom_settings  
    )

    f1 = Forces(label='D', V_z=37.727*kip, N_x=20*kip)
    f2 = Forces(label='L', V_z=6*kip) # No shear reinforcing
    Node(section=beam, forces=[f1, f2])
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.5*inch, s_l=6*inch) 
    results = beam.check_shear(A_s=0.847*inch**2)
    print(results)
    # section.design_shear(f, A_s=0.847*inch**2)
    # section.shear_results_detailed  
    # section.shear_results_detailed_doc

def rebar() -> None:
    concrete= Concrete_ACI_318_19(name="H30",f_c=30*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa) 
    section = RectangularBeam(
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

def shear_EHE_08() -> None:
    concrete= Concrete_EHE_08(name="C25",f_ck=25*MPa) 
    steelBar= SteelBar(name="B500S", f_y=500*MPa)
    custom_settings = {'clear_cover': 2.6*cm, 'stirrup_diameter_ini':8*mm,
                       'longitudinal_diameter_ini': 16*mm}
    beam = RectangularBeam(label="V-20x60",
                                      concrete=concrete,steel_bar=steelBar,width=20*cm, height=60*cm,
                                       settings=custom_settings)
    # f = Forces(V_z=100*kN, M_y=100*kNm)
    f = Forces(V_z=150*kN, M_y=50*kNm)
    A_s = 8.04*cm**2
    beam.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=20*cm)
    Node(beam, f)
    beam.check_shear(A_s)
    # section.shear_results_detailed
    # section.shear_results_detailed_doc

if __name__ == "__main__":
    # flexure()
    # shear_ACI_imperial()
    # shear_EHE_08()
    shear_ACI_metric()
    # rebar()
