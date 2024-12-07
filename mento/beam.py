from devtools import debug
from dataclasses import dataclass
from IPython.display import Markdown, display
from typing import Optional, Dict, Any, cast, Union, List
from pint.facets.plain import PlainQuantity
import numpy as np
import pandas as pd
from pandas import DataFrame
import math
import warnings

from mento.rectangular import RectangularSection
from mento.material import Concrete, SteelBar, Concrete_ACI_318_19, Concrete_EN_1992_2004
from mento.rebar import Rebar
from mento.units import MPa, ksi, psi, kip, mm, inch, kN, m, cm, kNm, ft, deg, dimensionless
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
        self._DCRv: float = 0
        self._s_l = self._stirrup_s_l
        self._s_w: PlainQuantity = 0*cm
        self._s_max_l: PlainQuantity = 0*cm
        self._s_max_w: PlainQuantity = 0*cm
        self._shear_checked = False  # Tracks if shear check or design has been done
        self.layers_spacing = self.settings.get_setting('layers_spacing')
        self._n1_bot = 0
        self._d_b1_b = 0*mm
        self._n2_b = 0
        self._d_b2_b = 0*mm
        self._n3_b = 0
        self._d_b3_b = 0*mm
        self._n4_b = 0
        self._d_b4_b = 0*mm
        self._n1_t = 0
        self._d_b1_t = 0*mm
        self._n2_t = 0
        self._d_b2_t = 0*mm
        self._n3_t = 0
        self._d_b3_t = 0*mm
        self._n4_t = 0
        self._d_b4_t = 0*mm
        self._A_s_bot=0*mm**2
        self._A_s_top=0*mm**2
        self._bot_rebar_centroid=0*mm
        self._top_rebar_centroid=0*mm

    @property
    def d(self) -> PlainQuantity:
        "Effective height."
        return self._d
    
    def initialize_code_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_ACI_318_19):
            self._initialize_aci_318_attributes()
        elif isinstance(self.concrete, Concrete_EN_1992_2004):
            self._initialize_en_1992_2004_attributes()

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

    def _initialize_en_1992_2004_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_EN_1992_2004):
            self._f_yk = self.steel_bar.f_y
            self._f_ck = self.concrete.f_ck
            self._V_Ed_1: PlainQuantity = 0*kN
            self._V_Ed_2: PlainQuantity = 0*kN
            self._N_Ed: PlainQuantity = 0*kN
            self._M_Ed: PlainQuantity = 0*kNm
            self._sigma_cd: PlainQuantity = 0*MPa
            self._V_Rd_c: PlainQuantity = 0*kN
            self._V_Rd_s: PlainQuantity = 0*kN
            self._V_Rd_max: PlainQuantity = 0*kN
            self._k_value: float = 0
            self._rho_l: PlainQuantity = 0 * dimensionless

    def set_transverse_rebar(self, n_stirrups: int = 0, d_b:PlainQuantity = 0*mm, s_l:PlainQuantity = 0*cm) -> None:
        """Sets the transverse rebar in the object."""
        self._stirrup_n = n_stirrups
        self._stirrup_d_b = d_b
        self._stirrup_s_l = s_l
        # Update effective height d with new values
        self._d = self._height -(self.c_c+self._stirrup_d_b+self._long_d_b/2) # Initial value
        n_legs = n_stirrups * 2
        A_db = (d_b ** 2) * math.pi / 4  # Area of one stirrup leg
        A_vs = n_legs * A_db  # Total area of stirrups
        self._A_v = A_vs / s_l  # Stirrup area per unit length

    def set_longitudinal_rebar_bot(self, n1: int, d_b1: PlainQuantity, n2: int = 0, d_b2: PlainQuantity=0*mm, 
                                n3: int=0, d_b3: PlainQuantity=0*mm, n4: int=0, d_b4: PlainQuantity=0*mm, 
                                ) -> None:
        """Sets the longitudinal rebar in the object."""
        self._n1_bot = n1
        self._d_b1_b = d_b1
        self._n2_b = n2
        self._d_b2_b = d_b2 if n2 > 0 else None
        self._n3_b = n3
        self._d_b3_b = d_b3 if n3 > 0 else None
        self._n4_b = n4
        self._d_b4_b = d_b4 if n4 > 0 else None
        self._A_s_bot = n1*d_b1**2*np.pi/4 + n2*d_b2**2*np.pi/4 + n3*d_b3**2*np.pi/4 + n4*d_b4**2*np.pi/4
        self._bot_rebar_centroid=self.__calculate_long_rebar_centroid(n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4)

    def set_longitudinal_rebar_top(self, n1: int, d_b1: PlainQuantity, n2: int=0, d_b2: PlainQuantity=0*mm, 
                                n3: int=0, d_b3: PlainQuantity = 0*mm, n4: int=0, d_b4: PlainQuantity=0*mm 
                                ) -> None:
        """Sets the longitudinal rebar in the object."""
        self._n1_t = n1
        self._d_b1_t = d_b1
        self._n2_t = n2
        self._d_b2_t = d_b2 if n2 > 0 else None
        self._n3_t = n3
        self._d_b3_t = d_b3 if n3 > 0 else None
        self._n4_t = n4
        self._d_b4_t = d_b4 if n4 > 0 else None
        self._A_s_top = n1*d_b1**2*np.pi/4 + n2*d_b2**2*np.pi/4 + n3*d_b3**2*np.pi/4 + n4*d_b4**2*np.pi/4
        self._top_rebar_centroid=self.__calculate_long_rebar_centroid(n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4)

    def __calculate_long_rebar_centroid(
        self, 
        n1: int, 
        d_b1: PlainQuantity, 
        n2: int = 0, 
        d_b2: PlainQuantity = 0 * mm, 
        n3: int = 0, 
        d_b3: PlainQuantity = 0 * mm, 
        n4: int = 0, 
        d_b4: PlainQuantity = 0 * mm) -> PlainQuantity:
        """
        Calculates the centroid (baricenter) of a group of rebars based on their diameters, quantities, 
        and layer spacing.

        Returns:
            float: The calculated centroid height of the rebar group.
        """
        # Default to 0 if any diameter or quantity is None
        d_b1 = d_b1 if d_b1 is not None else 0*mm
        d_b2 = d_b2 if d_b2 is not None else 0*mm
        d_b3 = d_b3 if d_b3 is not None else 0*mm
        d_b4 = d_b4 if d_b4 is not None else 0*mm

        n1 = n1 if n1 is not None else 0
        n2 = n2 if n2 is not None else 0
        n3 = n3 if n3 is not None else 0
        n4 = n4 if n4 is not None else 0

        # Calculate the vertical positions of the bar layers
        y1 = d_b1 / 2
        y2 = d_b2 / 2
        y3 = max(d_b1, d_b2) + self.layers_spacing + d_b3 / 2
        y4 = max(d_b1, d_b2) + self.layers_spacing + d_b4 / 2

        # Calculate the total area of each layer
        A1 = n1 * d_b1  # Area proportional to number of bars and their diameter
        A2 = n2 * d_b2
        A3 = n3 * d_b3
        A4 = n4 * d_b4

        # Calculate the centroid as a weighted average
        total_area = A1 + A2 + A3 + A4
        if total_area == 0:
            return 0*mm  # Avoid division by zero if no bars are present

        centroid = (A1 * y1 + A2 * y2 + A3 * y3 + A4 * y4) / total_area
        return centroid

# ======== ACI 318-19 methods =========

    def __calculate_phi_ACI_318_19(self, epsilon_most_strained: float) -> float:
        """
        Calculates the strength reduction factor (φ) for flexural design 
        based on ACI 318-19.
        It is used for columns; for beams, it is not required since beams 
        are always designed to be tension-controlled, with φ=0.9.

        Parameters:
            epsilon_most_strained (float): Strain in the most strained steel fiber.

        Returns:
            float: The strength reduction factor (φ), ranging from 0.65 to 0.9.

        Description:
            - φ = 0.65 if ε_most_strained ≤ ε_y (yield strain).
            - φ transitions linearly from 0.65 to 0.9 if ε_y < ε_most_strained ≤ ε_y + ε_c 
            (concrete crushing strain).
            - φ = 0.9 if ε_most_strained > ε_y + ε_c.
        """
        # Retrieve concrete crushing strain (ε_c)
        epsilon_c = self.concrete.get_properties()["epsilon_c"]

        # Calculate φ based on ε_most_strained
        if epsilon_most_strained <= self.steel_bar.epsilon_y:
            return 0.65
        elif epsilon_most_strained <= self.steel_bar.epsilon_y + epsilon_c:
            return (
                (0.9 - 0.65) 
                * (epsilon_most_strained - self.steel_bar.epsilon_y) / epsilon_c
                + 0.65
            )
        else:
            return 0.9

    def __maximum_flexural_reinforcement_ratio_ACI_318_19(self) -> float:
        """
        Calculates the maximum flexural reinforcement ratio (ρ_max) according to the 
        ACI 318-19 design code.

        Returns:
            float: The maximum reinforcement ratio (ρ_max) for the section, 
            or 0 if the design code is not ACI 318-19.

        Description:
            This function determines the maximum reinforcement ratio (ρ_max) 
            allowed by the ACI 318-19 design code to ensure ductile behavior 
            of reinforced concrete sections. The calculation depends on the 
            properties of the concrete (β1 and ε_c) and steel (ε_y).

            Note:
            - The parameter `epsilon_min_rebar_ACI_318_19` is calculated as the 
            sum of the steel's yield strain (ε_y) and the concrete's strain 
            at crushing (ε_c).
            - This function only works if the `design_code` of the concrete is 
            set to "ACI 318-19". For other codes, it returns 0.

            TODO:
            - Verify the calculation of `epsilon_min_rebar_ACI_318_19` to ensure 
            it complies with the ACI 318-19 provisions.
        """
        if self.concrete.design_code == "ACI 318-19":
            # Cast the concrete object to the specific ACI subclass
            concrete_aci = cast(Concrete_ACI_318_19, self.concrete)
            
            # Calculate minimum steel strain for ductility
            epsilon_min_rebar_ACI_318_19 = (
                self.steel_bar.epsilon_y + concrete_aci.epsilon_c
            )
            
            # Calculate maximum reinforcement ratio (ρ_max)
            rho_max = (
                0.85 * concrete_aci.beta_1 * self.concrete.f_c / self.steel_bar.f_y
                * (concrete_aci.epsilon_c /
                (concrete_aci.epsilon_c + epsilon_min_rebar_ACI_318_19))
            )
            return rho_max
        else:
            return 0

    def __minimum_flexural_reinforcement_ratio_ACI_318_19(self)->PlainQuantity:
        concrete_properties = self.concrete.get_properties()
        f_c = concrete_properties['f_c']
        rebar_properties = self.steel_bar.get_properties()
        f_y = rebar_properties['f_y']
        minimum_ratio = max((3 * np.sqrt(f_c / psi) * psi / f_y), 
                    (200 * psi / f_y))
        return minimum_ratio

    def __calculate_flexural_reinforcement_ACI_318_19(self, M_u: float, d: float, d_prima: float)-> tuple[PlainQuantity, PlainQuantity, PlainQuantity, PlainQuantity]:
        """
        Calculates the flexural reinforcement (bottom tensile and top compression, if required) 
        for a given factored moment.

        Parameters:
            M_u (float): Factored moment (must always be passed as a positive value).
            d (float): Effective depth of the tensile reinforcement.
            d_prima (float): Effective depth of the compression reinforcement.

        Returns:
            Tuple:
                - A_s_min (float): Minimum reinforcement area required by the code.
                - A_s_max (float): Maximum reinforcement area allowed by the code.
                - A_s_final (float): Final reinforcement area adopted for the tensile zone.
                - A_s_comp (float): Compression reinforcement area (if required).

        Description:
            This is a private function that calculates the reinforcement for a given factored moment.
            The moment \( M_u \) must always be passed as a **positive value**, regardless of whether 
            it corresponds to a positive or negative moment. To study a positive moment, invoke this 
            function with \( d \) as the tensile reinforcement depth and \( d_prima \) as the 
            compression reinforcement depth. To study a negative moment, pass \( d_prima \) as the 
            tensile depth and \( d \) as the compression depth.

            This function is meant to be used internally by higher-level methods and is not intended 
            for direct use by the user.

        Example:
            # For a positive moment:
            A_s_min, A_s_max, A_s_final, A_s_comp = self.__calculate_flexural_reinforcement(M_u_max, d, d_prima)
                
            # For a negative moment:
            A_s_min, A_s_max, A_s_final, A_s_comp = self.__calculate_flexural_reinforcement(abs(M_u_min), d_prima, d)

        """
        # Extract relevant properties from self
        phi = self.settings.get_setting('phi_t')
        setting_flexural_min_reduction = self.settings.get_setting('flexural_min_reduction')
        concrete_properties = self.concrete.get_properties()
        f_c = concrete_properties['f_c']
        beta_1 = concrete_properties['beta_1']
        rebar_properties = self.steel_bar.get_properties()
        f_y = rebar_properties['f_y']
        b = self._width

        # Determine minimum and maximum reinforcement
        rho_min=self.__minimum_flexural_reinforcement_ratio_ACI_318_19()
        A_s_min = rho_min * d * b
        rho_max = self.__maximum_flexural_reinforcement_ratio_ACI_318_19()
        A_s_max = rho_max * d * b

        # Calculate required reinforcement
        R_n = M_u / (phi * b * d**2)
        A_s_calc = 0.85 * f_c * b * d / f_y * (1 - np.sqrt(1 - 2 * R_n / (0.85 * f_c)))
        
        self._c=A_s_calc*f_y/(0.85*f_c*b*beta_1) # C=T -> 0.85*f_c*c*beta_1*b = A_s *f_y -> a = A_s *f_y / (0.85*f_c*b)


        # Adjust required reinforcement to limits
        if A_s_calc > A_s_min:
            A_s_final = A_s_calc
        elif 4 * A_s_calc / 3 > A_s_min:
            A_s_final = A_s_min
        else:
            if setting_flexural_min_reduction == 'True':
                A_s_final = 4 * A_s_calc / 3
            else:
                A_s_final = A_s_min

        # Check if compression reinforcement is required
        if A_s_final <= A_s_max:
            A_s_comp = 0 * cm**2
        else:
            rho = 0.85 * beta_1 * f_c / f_y * (0.003 / (self.steel_bar.epsilon_y + 0.006))
            M_n_t = rho * f_y * (d - 0.59 * rho * f_y * d / f_c) * b * d
            M_n_prima = M_u / phi - M_n_t
            c_t = 0.003 * d / (self.steel_bar.epsilon_y + 0.006)
            f_s_prima = min(0.003 * self.steel_bar.E_s * (1 - d_prima / c_t), f_y)
            A_s_comp = M_n_prima / (f_s_prima * (d - d_prima))
            A_s_final = rho * b * d + A_s_comp

        return A_s_min, A_s_max, A_s_final, A_s_comp


    def _compile_results_aci_flexure_metric(self, force: Forces) -> List[Dict[str, Any]]:
        # Create dictionaries for bottom and top rows
        if force._M_y>=0:
            result = {
                'Section Label': self.label,
                'Load Combo': force.label,
                'Position': 'Bottom',
                'As,min': self._A_s_min_bot.to('cm ** 2'),
                'As,req bot': self._A_s_req_bot.to('cm ** 2'),
                'As,req top': self._A_s_req_top.to('cm ** 2'),
                'As bot': self._A_s_bot.to('cm ** 2'),
                'c/d': self._c.to('mm')/self.d.to('mm'),
                'Mu': force._M_y.to('kN*m'),
                'ØMn': self._Phi_M_n_positive.to('kN*m'),
                'Mu<ØMn': force._M_y <= self._Phi_M_n_positive,
                'DCR': force._M_y.to('kN*m') / self._Phi_M_n_positive.to('kN*m')
            }
        else:
            result = {
                'Section Label': self.label,
                'Load Combo': force.label,
                'Position': 'Top',
                'As,min': self._A_s_min_top.to('cm ** 2'),
                'As,req bot': self._A_s_req_bot.to('cm ** 2'),
                'As,req top': self._A_s_req_top.to('cm ** 2'),
                'As top': self._A_s_top.to('cm ** 2'),
                'c/d': self._c.to('mm')/(self.height-self._d_prime).to('mm'),
                'Mu': force._M_y.to('kN*m'),
                'ØMn': self._Phi_M_n_negative.to('kN*m'),
                'Mu<ØMn': force._M_y <= self._Phi_M_n_negative,
                'DCR': force._M_y.to('kN*m') / self._Phi_M_n_negative.to('kN*m')
            }

        # Return both rows as a list of dictionaries
        return result

    def _determine_nominal_moment_simple_reinf_ACI_318_19(self, A_s: PlainQuantity, d: PlainQuantity) -> PlainQuantity:
        #TODO REDACTAR BIEN ESTO
        '''
            Se usa esta formula SOLO cuando el acero dispuesto es menor o igual que A_s_max
            
        Equilibrium: C=T
        0.85*f_c*a*b = A_s *f_y -> a = A_s *f_y / (0.85*f_c*b)
        '''
        a = A_s * self.steel_bar.f_y /(0.85*self.concrete.f_c*self._width)
        M_n=A_s* self.steel_bar.f_y *(d-a/2)
        return M_n
    
    def _determine_nominal_moment_double_reinf_ACI_318_19(self, A_s: PlainQuantity, d: PlainQuantity, c_mec: PlainQuantity, d_prime: PlainQuantity, A_s_prime: PlainQuantity) -> PlainQuantity:
        '''
            Se usa solo cuando la viga tiene mas refuerzo que el maximo, y tiene armadura de compresion.
        '''
        f_c = self.concrete.f_c
        if isinstance(self.concrete, Concrete_ACI_318_19):
            beta_1 = self.concrete._beta_1
            epsilon_c = self.concrete._epsilon_c

        f_y= self.steel_bar.f_y
        E_s = self.steel_bar._E_s
        epsilon_y=self.steel_bar._epsilon_y
        b = self._width

        '''
        Equilibrium: C=T --> Cc + Cs = T
        Sin embargo la forma de calcuar Cs depende de si el acero comprimido está en fluencia o no.
        Se asume que si, y se chequea
        0.85*f_c*a*b + A_s_prime * f_y = A_s * f_y
        '''

        # Assuming f_s_prime=f_y:
        c_assumed = (A_s * f_y - A_s_prime * f_y) /(0.85*f_c*b*beta_1)
        epsilon_s=(c_assumed-d_prime)/c_assumed*epsilon_c

        #Now we check our assumption
        if(epsilon_s>=epsilon_y):
            # Our assumption was right so:
            a_assumed=c_assumed*beta_1
            M_n= 0.85*f_c*a_assumed*b * (d-a_assumed/2) +  A_s_prime * f_y * (d-d_prime)
            return M_n
        else:
            # Our assumption was wrong so we have to determine the real c value from quadratic equation:
            # A_s * f_y = 0.85 * f_c *b * β1 * c + A_s_prime * (c-d_prime)/c * εc * E_s
            A = 0.85*f_c*b*beta_1
            B = A_s_prime * epsilon_c * E_s - A_s * f_y
            C = -d_prime*A_s_prime*epsilon_c*E_s
            c=(-B+np.sqrt(B**2-4*A*C))/(2*A)
            a=c*beta_1
            epsilon_s_prime=(c-d_prime)/c*epsilon_c
            f_s_prime=epsilon_s_prime*E_s
            A_s_2=A_s_prime*f_s_prime/f_y
            A_s_1=A_s-A_s_2
            M_n_1=A_s_1*f_y*(d-a/2)
            M_n_2=A_s_prime*f_s_prime*(d-d_prime)
            M_n= M_n_1+M_n_2
            return M_n

    def _determine_nominal_moment_ACI_318_19(self) -> tuple[PlainQuantity, PlainQuantity]:
        '''
            Para una determinada seccion con sus dos armados, top y bottom, encuentra el momento nominal para 
            momentos positivos y negativos
        '''

        concrete_properties = self.concrete.get_properties()
        f_c = concrete_properties['f_c']
        beta_1 = concrete_properties['beta_1']
        rebar_properties = self.steel_bar.get_properties()
        f_y = rebar_properties['f_y']
        E_s = rebar_properties['E_s']
        epsilon_y=rebar_properties['epsilon_y']
        b = self._width
        h=self.height
        # Retrieve concrete crushing strain (ε_c)
        epsilon_c = self.concrete.get_properties()["epsilon_c"]

        # Load design settings for ACI 318-19
        self.settings.load_aci_318_19_settings()

        # Informacion de la seccion:
        A_s_bot=self._A_s_bot
        A_s_top=self._A_s_top
        c_mec = self.c_c + self._stirrup_d_b + self._bot_rebar_centroid
        self._d_prime = self.c_c + self._stirrup_d_b + self._top_rebar_centroid
        d = self._height - c_mec
        self._rho_min = self.__minimum_flexural_reinforcement_ratio_ACI_318_19()
        self._rho_max = self.__maximum_flexural_reinforcement_ratio_ACI_318_19()

        # DETERMINACION DE LA CAPACIDAD PARA MOMENTO POSITIVO (TRACCIONA ABAJO):
        self._A_s_min_bot = self._rho_min * d * b
        self._A_s_max_bot = self._rho_max * d * b
        if (A_s_bot<=self._A_s_max_bot):
            M_n_positive = self._determine_nominal_moment_simple_reinf_ACI_318_19(self._A_s_bot, d)
        elif (A_s_top==0*inch**2):
            M_n_positive = self._determine_nominal_moment_simple_reinf_ACI_318_19(self._A_s_max_bot, d)
        else:
            M_n_positive=self._determine_nominal_moment_double_reinf_ACI_318_19(A_s_bot, d, c_mec, self._d_prime, A_s_top)

        # DETERMINACION DE LA CAPACIDAD PARA MOMENTO NEGATIVO (TRACCIONA ARRIBA):
        self._A_s_min_top = self._rho_min * (h-self._d_prime) * b
        self._A_s_max_top = self._rho_max *  (h-self._d_prime) * b
        if (A_s_top==0):
            M_n_negative=0
        elif (A_s_top<=self._A_s_max_top):
            M_n_negative=-self._determine_nominal_moment_simple_reinf_ACI_318_19(A_s_top, (h-self._d_prime), self._d_prime)
        elif (A_s_bot==0):
            M_n_negative=-self._determine_nominal_moment_simple_reinf_ACI_318_19(self._A_s_max_top, (h-self._d_prime), self._d_prime)
        else:
            M_n_negative=-self._determine_nominal_moment_double_reinf_ACI_318_19(A_s_top, (h-self._d_prime), self._d_prime, c_mec, A_s_bot)

        return M_n_positive, M_n_negative

    def check_flexure_ACI_318_19(self, Force: Forces) -> pd.DataFrame:
        '''
            flexure_design_results_top: es el diccionario de armado de rebar que tiene las 4 posiciones
        '''

        # Load design settings for ACI 318-19
        self.settings.load_aci_318_19_settings()
        phi = self.settings.get_setting('phi_t')
        # Initial assumptions for mechanical cover and compression depth

        c_mec = self.c_c + self._stirrup_d_b + self._bot_rebar_centroid
        d_prime = self.c_c + self._stirrup_d_b + self._top_rebar_centroid
        d = self._height - c_mec

        self._M_n_positive, self._M_n_negative=self._determine_nominal_moment_ACI_318_19()
        self._Phi_M_n_positive = phi * self._M_n_positive
        self._Phi_M_n_negative = phi * self._M_n_negative

        if Force._M_y>=0:
            A_s_min, A_s_max, self._A_s_req_bot, self._A_s_req_top=self.__calculate_flexural_reinforcement_ACI_318_19(Force._M_y, d, d_prime)
        else:
            A_s_min, A_s_max, self._A_s_req_top, self._A_s_req_bot=self.__calculate_flexural_reinforcement_ACI_318_19(-Force._M_y, self._height-d_prime, c_mec)

        results=self._compile_results_aci_flexure_metric(Force)
        # self._initialize_dicts_ACI_318_19_flexure()
        return pd.DataFrame([results], index=[0])

    def design_flexure_ACI_318_19(self, Force: Forces) -> Dict[str, Any]:
        """
        Designs the flexural reinforcement for a beam cross-section 
        according to ACI 318-19.

        Parameters:
            Force (Forces): A list of moment forces acting on the section. If a 
            single moment is passed, it will be converted into a list.

        Returns:
            Dict[str, Any]: A DataFrame containing the adopted bottom and top 
            reinforcement areas and their respective bar spacings.

        Description:
            This function determines the flexural reinforcement for positive 
            (tension at the bottom) and negative moments (tension at the top) 
            iteratively. The process adjusts the effective depths (`d` and 
            `d_prima`) and the mechanical cover (`rec_mec`) until convergence.

            The reinforcement design starts with the positive moment (M_u_max), 
            determining the required bottom reinforcement and, if needed, 
            compression reinforcement at the top. For sections experiencing 
            negative moments (M_u_min), the function recalculates the top 
            reinforcement requirements by swapping the positions of `d` and 
            `d_prima`.

            The function outputs a summary of the final adopted reinforcement 
            areas and their spacing.
        """

        # Ensure Force is a list
        if not isinstance(Force, list):
            Force = [Force]

        # Initialize the maximum and minimum moments
        M_u_min = 0 * kNm  # Negative moment (tension at the top)
        M_u_max = 0 * kNm  # Positive moment (tension at the bottom)
        for item in Force:
            if item._M_y > M_u_max:
                M_u_max = item._M_y
            if item._M_y < M_u_min:
                M_u_min = item._M_y

        # Load design settings for ACI 318-19
        self.settings.load_aci_318_19_settings()

        # Initial assumptions for mechanical cover and compression depth
        rec_mec = self.c_c + self._stirrup_d_b + 1 * inch
        d_prima = self.c_c + self._stirrup_d_b + 1 * inch

        # Start the iterative process
        tol = 0.01 * cm  # Tolerance for convergence
        Err = 2 * tol
        iteration_count = 0

        while Err >= tol:
            iteration_count += 1
            # Update the effective depth for bottom tension reinforcement
            d = self._height - rec_mec

            # Calculate reinforcement for the positive moment
            (
                self._A_s_min_bot,
                self._A_s_max_bot,
                self._A_s_final_bot_Positive_M,
                self._A_s_comp_top,
            ) = self.__calculate_flexural_reinforcement_ACI_318_19(M_u_max, d, d_prima)

            # Initialize bottom and top reinforcement
            self._A_s_bottom = self._A_s_final_bot_Positive_M
            self._A_s_top = self._A_s_comp_top

            # If there is a negative moment, calculate the top reinforcement
            if M_u_min < 0:
                (
                    self._A_s_min_top,
                    self._A_s_max_top,
                    self._A_s_final_top_Negative_M,
                    self._A_s_comp_bot,
                ) = self.__calculate_flexural_reinforcement_ACI_318_19(abs(M_u_min), d_prima, d)

                # Adjust reinforcement areas based on positive and negative moments
                self._A_s_bottom = max(self._A_s_final_bot_Positive_M, self._A_s_comp_bot)
                self._A_s_top = max(self._A_s_comp_top, self._A_s_final_top_Negative_M)

            # Design bottom reinforcement
            section_rebar_bot = Rebar(self)
            self.flexure_design_results_bot = section_rebar_bot.longitudinal_rebar_ACI_318_19(self._A_s_bottom)
            best_design = section_rebar_bot.longitudinal_rebar_design

            # Extract bar information
            d_b1_bot=best_design["d_b1"]
            d_b2_bot=best_design["d_b2"]
            d_b3_bot=best_design["d_b3"]
            d_b4_bot=best_design["d_b4"]
            n_1_bot=best_design["n_1"]
            n_2_bot=best_design["n_2"]
            n_3_bot=best_design["n_3"]
            n_4_bot=best_design["n_4"]
            total_as_adopted_bot = best_design["total_as"]
            total_bars_bot=n_1_bot+n_2_bot+n_3_bot+n_4_bot
            clearing_spacing_bot = best_design["clear_spacing"]


            # Set rebar information to section
            self.set_longitudinal_rebar_bot(n_1_bot, d_b1_bot, n_2_bot, d_b2_bot, 
                                n_3_bot, d_b3_bot, n_4_bot, d_b4_bot, 
                                total_as_adopted_bot, total_bars_bot, clearing_spacing_bot)

            c_mec_calc = self.c_c + self._stirrup_d_b + self._bot_rebar_centroid

            # Design top reinforcement
            if self._A_s_top > 0:
                section_rebar_top = Rebar(self)
                self.flexure_design_results_top = section_rebar_top.longitudinal_rebar_ACI_318_19(self._A_s_top)
                best_design_top = section_rebar_top.longitudinal_rebar_design

                # Extract bar information for top reinforcement
                d_b1_top=best_design_top["d_b1"]
                d_b2_top=best_design_top["d_b2"]
                d_b3_top=best_design_top["d_b3"]
                d_b4_top=best_design_top["d_b4"]
                n_1_top=best_design_top["n_1"]
                n_2_top=best_design_top["n_2"]
                n_3_top=best_design_top["n_3"]
                n_4_top=best_design_top["n_4"]
                total_as_adopted_top=n_1_top+n_2_top+n_3_top+n_4_top
                total_bars_top = best_design_top["total_as"]
                clearing_spacing_top = best_design_top["clear_spacing"]

                # Set rebar information to section
                self.set_longitudinal_rebar_top(n_1_top, d_b1_top, n_2_top, d_b2_top, 
                                    n_3_top, d_b3_top, n_4_top, d_b4_top, 
                                    total_as_adopted_top, total_bars_top, clearing_spacing_top)

                d_prima_calculo = self.c_c + self._stirrup_d_b + self._top_rebar_centroid
            else:
                # If no top reinforcement is required
                d_prima_calculo = d_prima
                total_as_adopted_top = 0 * inch**2
                clearing_spacing_top = self._width - 2 * self.c_c - 2 * self._stirrup_d_b
                # Set rebar information to section
                self.set_longitudinal_rebar_top(0, 0*mm, 0, 0*mm, 
                                    0, 0*mm, 0, 0*mm,total_as_adopted_top, 0, clearing_spacing_top)

            # Update error for iteration
            Err = max(abs(c_mec_calc - rec_mec), abs(d_prima_calculo - d_prima))
            rec_mec = c_mec_calc
            d_prima = d_prima_calculo

        # Return results as a DataFrame
        results = {
            "Bottom_As_adopted": self._total_as_b.to("inch**2"),
            "Bottom separation of bars": self._clear_spacing_b.to("inch"),
            "As_compression_adopted": self._total_as_t.to("inch**2"),
            "Top separation of bars": self._clear_spacing_t.to("inch"),
        }
        return pd.DataFrame([results], index=[0])

    def design_flexure_EN_1992(self, M_u: float) -> None:
        pass

    def design_flexure(self) -> pd.DataFrame:
        """
        Designs flexural reinforcement for the beam using the provided forces and design code.
        Identifies the limiting cases for top and bottom reinforcement, designs for those cases, 
        and then checks flexural capacity for all forces.
        
        Returns
        -------
        DataFrame
            A DataFrame summarizing the flexural design results for all forces.
        """
        if not self.node or not self.node.forces:
            raise ValueError("No Node or forces list associated with this beam.")

        self._flexure_results_list = []  # Store results for each force
        self._flexure_results_detailed_list = {}  # Store detailed results by force ID

        # Initialize limiting cases
        max_M_y_top = 0 * kN * m  # For negative M_y (top reinforcement design)
        max_M_y_bot = 0 * kN * m  # For positive M_y (bottom reinforcement design)

        # Identify the limiting cases
        for force in self.node.forces:
            # For top reinforcement, consider the minimum (most negative) moment
            if force.M_y < max_M_y_top:
                max_M_y_top = force.M_y
            # For bottom reinforcement, consider the maximum positive moment
            if force.M_y > max_M_y_bot:
                max_M_y_bot = force.M_y

        # Design flexural reinforcement for the limiting cases
        if self.concrete.design_code == "ACI 318-19":
            # Only design for top if there are negative moments
            top_result = None
            if max_M_y_top < 0 * kN * m:
                top_result = self.design_flexure_ACI_318_19(max_M_y_top)

            bot_result = self.design_flexure_ACI_318_19(max_M_y_bot)
        else:
            raise ValueError(f"Longitudinal design method not implemented "
                            f"for concrete type: {type(self.concrete).__name__}")

        # Store the limiting case designs
        self._limiting_case_top = top_result
        self._limiting_case_bottom = bot_result

        # Assign the designed reinforcement to the section (only if top reinforcement is required)
        self.set_flexure_rebar(top_result['As'] if top_result else None, bot_result['As'])

        # Check flexural capacity for all forces with the assigned reinforcement
        for force in self.node.forces:
            result = self.check_flexure_ACI_318_19(force)  # Assuming ACI 318-19 for simplicity
            self._flexure_results_list.append(result)

            # Store detailed results for each force
            self._flexure_results_detailed_list[force.id] = {
                'forces': self._forces_flexure.copy(),
                'min_max': self._data_min_max.copy(),
                'flexure_capacity_top': self._flexure_capacity_top.copy(),
                'flexure_capacity_bottom': self._flexure_capacity_bot.copy(),
            }

        # Compile all results into a single DataFrame
        all_results = pd.DataFrame(self._flexure_results_list)

        # Mark flexure as checked
        self._flexure_checked = True

        return all_results




    def _initialize_dicts_ACI_318_19_flexure(self) -> None:
        """Initialize the dictionaries used in check and design methods."""
        self._materials_flexure = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
                "Concrete density",
                "Normalweight concrete",
                "Safety factor for bending"
            ],
            "Variable": ["","fc", "fy", "γc", "λ", "Øt"],
            "Value": [self.label, round(self.concrete.f_c.to('MPa').magnitude,2), 
                      round(self.steel_bar.f_y.to('MPa').magnitude,2),round(self.concrete.density.to('kg/m**3').magnitude,1),
                       self.settings.get_setting('lambda'), self.settings.get_setting('phi_t')],
            "Unit": ["", "MPa", "MPa", "kg/m³", "", ""]
        }
        self._geometry_flexure = {
            "Geometry": [
                "Section height",
                "Section width",
                "Clear cover",
                "Mechanical top cover",
                "Mechanical bottom cover",
            ],
            "Variable": ["h", "b", "cc", "cm,top", "cm,bot"],
            #TODO: ver bien tema As de armadura traccionada que podria ser superior o inferior.
            "Value": [self.height.to('cm').magnitude, self.width.to('cm').magnitude, self.c_c.to('cm').magnitude,
                       self._c_mec_top.to('cm').magnitude, self._c_mec_bot.to('cm').magnitude],
            "Unit": ["cm", "cm", "cm", "cm", "cm"]
        }
        self._forces_flexure = {
            "Design forces": [
                "Top max moment",
                "Bottom max moment",
            ],
            "Variable": ["Mu,top", "Mu,bot"],
            "Value": [round(self._M_u_top.to('kN*m').magnitude,2), round(self._M_u_bot.to('kN*m').magnitude,2) ],
            "Unit": ["kNm", "kNm"]
        }
        # Min max lists
        min_spacing_top = max(self.settings.get_setting('clear_spacing'), self.settings.get_setting('vibrator size'),
                              self._d_b_max_top)
        min_spacing_bot = max(self.settings.get_setting('clear_spacing'), self._d_b_max_bot)
        min_values = [self._A_s_min_top, None, min_spacing_top, self._A_s_min_bot,
                      None,min_spacing_bot]   # Use None for items without a minimum constraint
        max_values = [None, self._A_s_max_top, None, None,
                      self._A_s_max_bot,None] # Use None for items without a maximum constraint
        current_values = [self._A_s_top, self._A_s_top, self._avaialable_s_top, self._A_s_bot, self._A_s_bot,
                          self._available_s_bot]  # Current values to check

        # Generate check marks based on the range conditions
        checks = [
            '✔️' if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else '❌'
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._data_min_max_flexure = {
            'Check': ['Minimum rebar top', 'Maximum rebar top', 'Minimum spacing top','Minimum rebar bottom', 
                      'Maximum rebar bottom', 'Minimum spacing bottom'],
            'Unit': ['cm²', 'cm²','cm²','cm²','mm', 'mm'],
            'Value': [self._A_s_top.to('cm**2').magnitude, self._A_s_top.to('cm**2').magnitude,
            self._avaialable_s_top.to('mmm').magnitude, self._A_s_bot.to('cm**2').magnitude, 
            self._A_s_bot.to('cm**2').magnitude,self._avaialable_s_bot.to('mmm').magnitude,],
            'Min.': [round(self._A_s_min_top.to('cm**2').magnitude,2), "", min_spacing_top.to('mm').magnitude,
                     round(self._A_s_min_bot.to('cm**2').magnitude,2), "", min_spacing_bot.to('mm').magnitude],
            'Max.': ["", round(self._A_s_max_top.to('cm**2').magnitude,2), "", "",
                     round(self._A_s_max_bot.to('cm**2').magnitude,2), ""],
            'Ok?': checks
        }
        a_d_top = self._a_top/self._d_top
        a_d_bot = self._a_bot/self._d_bot
        check_DCR_top = '✔️' if self._DCR_top < 1 else '❌'
        check_DCR_bot = '✔️' if self._DCR_bot < 1 else '❌'
        self._flexure_capacity_top = {
            "Top reinforcement check": [
                "First layer bars",
                "Second layer bars",
                "Effective height",
                "Depth of equivalent strength block ratio",
                "Minimum rebar reinforcing",
                "Required rebar reinforcing",
                "Defined rebar reinforcing",
                "Longitudinal reinforcement ratio",
                "Total flexural strength", 
                "Demand Capacity Ratio"
            ],
            "Variable": ["n1+n2", "n3+n4", "d", "a/d", "As,min","As,req","As", "ρl", "ØMn", "DCR"],
            "Value": ["", "", self._d_top.to('cm').magnitude,
                    a_d_top, round(self._A_s_min_top.to('cm**2').magnitude,2),
                    round(self._A_s_req_top.to('cm**2').magnitude,2),
                    round(self._A_s_top.to('cm**2').magnitude,2), round(self._rho_l_top.magnitude,5),
                    round(self._phi_M_n_top.to('kN*m').magnitude,2), round(self._DCR_top,2)],
            "Unit": ["", "", "cm", "", "cm²","cm²", "cm²","","kNm", check_DCR_top]
        }
        self._flexure_capacity_bot = {
            "Bottom reinforcement check": [
                "First layer bars",
                "Second layer bars",
                "Effective height",
                "Depth of equivalent strength block ratio",
                "Minimum rebar reinforcing",
                "Required rebar reinforcing",
                "Defined rebar reinforcing",
                "Longitudinal reinforcement ratio",
                "Total flexural strength", 
                "Demand Capacity Ratio"
            ],
            "Variable": ["n1+n2", "n3+n4", "d", "a/d", "As,min","As,req","As", "ρl", "ØMn", "DCR"],
            "Value": ["", "", self._d_bot.to('cm').magnitude,
                    a_d_bot, round(self._A_s_min_bot.to('cm**2').magnitude,2),
                    round(self._A_s_req_bot.to('cm**2').magnitude,2),
                    round(self._A_s_bot.to('cm**2').magnitude,2), round(self._rho_l_bot.magnitude,5),
                    round(self._phi_M_n_bot.to('kN*m').magnitude,2), round(self._DCR_bot,2)],
            "Unit": ["", "", "cm", "", "cm²","cm²", "cm²","","kNm", check_DCR_bot]
        }

    def _set_initial_conditions_aci_shear(self, Force: Forces, A_s: PlainQuantity) -> None:
        self._N_u = Force.N_x
        self._V_u = Force.V_z
        #TODO: Cambiar por bot o top según signo de momento de Force
        self._A_s = A_s
        self.settings.load_aci_318_19_settings()
        self.phi_v = self.settings.get_setting('phi_v')
        self.lambda_factor = self.settings.get_setting('lambda')
        self.f_yt = self._calculate_f_yt_aci()

    def _calculate_shear_reinforcement_aci(self) -> None:
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
        self._set_initial_conditions_aci_shear(force, A_s)

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
        self._initialize_dicts_ACI_318_19_shear()
        return pd.DataFrame([results], index=[0])

    def design_shear_ACI_318_19(self, force: Forces, A_s: PlainQuantity = 0 * cm ** 2) -> pd.DataFrame:
        # Set the initial variables
        self._set_initial_conditions_aci_shear(force, A_s)

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
        self._initialize_dicts_ACI_318_19_shear()
        return pd.DataFrame([results], index=[0])

# ======== EN-1992-2004 methods =========

    def _initialize_variables_EN_1992_2004(self, Force: Forces, A_s: PlainQuantity) -> None:
        if isinstance(self.concrete, Concrete_EN_1992_2004):
            # Set the initial variables
            self._N_Ed = Force.N_x
            self._V_Ed_1 = Force.V_z  # Consider the same shear at the edge of support and in d
            self._V_Ed_2 = Force.V_z  # Consider the same shear at the edge of support and in d
            self._A_s = A_s

            # Load settings for gamma factors
            self.settings.load_en_1992_2004_settings()
            self.settings.load_en_1992_2004_settings()
            self._alpha_cc = self.settings.get_setting('alpha_cc')
            self._gamma_c = self.settings.get_setting('gamma_c')
            self._gamma_s = self.settings.get_setting('gamma_s')
            self._f_ywk = self._f_yk
            self._f_ywd = self._f_ywk/self._gamma_s
            self._f_ywk = self._f_yk
            self._f_ywd = self._f_ywk/self._gamma_s
            self._f_cd = self._alpha_cc*self._f_ck/self._gamma_c

            # Minimum shear reinforcement calculation
            self._alpha = math.radians(90)
            rho_min = 0.08*math.sqrt(self._f_ck.to('MPa').magnitude) / (self._f_ywk)*MPa
            self._A_v_min = rho_min * self.width * math.sin(self._alpha)
            # Compression stress, positive
            self._A_p = 0*cm**2 # No prestressing for now
            self._rho_l = min((self._A_s + self._A_p) / (self.width * self.d), 0.02)
            # Shear calculation for sections without rebar
            self._k_value = min(1 + math.sqrt(200 * mm / self.d), 2)
            # Positive of compression
            self._sigma_cp = min(self._N_Ed / self.A_x, 0.2*self._f_cd)

    def _shear_without_rebar_EN_1992_2004(self) -> PlainQuantity:
        self._stirrup_d_b = 0*mm

        self._theta = 0
        # Total shear capacity without rebar
        C_rdc = 0.18/self._gamma_c
        v_min = 0.035*self._k_value**(3/2)*math.sqrt(self._f_ck.to('MPa').magnitude)
        k_1 = 0.15
        V_Rd_c_min = ((v_min+k_1*self._sigma_cp.to('MPa').magnitude)* self.width * self.d * MPa).to('kN')
        V_Rd_c = ((C_rdc*self._k_value*(100*self._rho_l*self._f_ck.to('MPa').magnitude)**(1/3)*MPa\
                  +k_1*self._sigma_cp.to('MPa'))* self.width * self.d).to('kN')        
        return max(V_Rd_c_min, V_Rd_c)
        
    def check_shear_EN_1992_2004(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> DataFrame:
        if isinstance(self.concrete, Concrete_EN_1992_2004):
            # Initialize all the code related variables
            self._initialize_variables_EN_1992_2004(Force,A_s)            
            
            if self._stirrup_n == 0:
                # Calculate V_Rd_c
                self._V_Rd_c = self._shear_without_rebar_EN_1992_2004()
                # According to EN1992-1-1 §6.2.1(4) minimum shear reinforcement should nevertheless be provided
                # according to EN1992-1-1 §9.2.2. The minimum shear reinforcement may be omitted in members where
                # transverse redistribution of loads is possible (such as slabs) and members of minor importance
                # which do not contribute significantly to the overall resistance and stability of the structure.
                self._A_v_req = self._A_v_min
                #Maximum shear capacity is the same as the concrete capacity
                self._V_Rd = self._V_Rd_c
                self._V_Rd_max = self._V_Rd
                self._max_shear_ok = self._V_Ed_1 <= self._V_Rd_max

            else:
                # Shear reinforcement calculations
                d_bs = self._stirrup_d_b
                s_l = self._stirrup_s_l
                n_legs = self._stirrup_n*2
                self._A_s = A_s
                A_db = (d_bs ** 2) * math.pi / 4  # Area of one stirrup leg
                A_vs = n_legs * A_db  # Total area of stirrups
                self._A_v = A_vs / s_l  # Stirrup area per unit length
                # Total shear strength with rebar
                alpha_cw = 1  # Non-prestressed members or members subject to tensile stress due to axial force
                v_1 = 0.6 * (1 - self._f_ck.to('MPa').magnitude / 250)  # Strength reduction factor for concrete struts
                z = 0.9 * self.d  # Lever arm

                # The θ angle is lmited between 21,8° ≤ θ ≤ 45°(1 ≤ cot(θ) ≤ 2.5)
                # Check the minimum strut angle θ = 21.8° (cot(θ) = 2.5)
                theta_min = math.radians(21.8)
                cot_theta_min = 1 / math.tan(theta_min)

                V_Rd_max_min_angle = (alpha_cw * self.width * z * v_1 * self._f_cd / (cot_theta_min +
                                                                                       math.tan(theta_min))).to('kN')
                debug(V_Rd_max_min_angle, cot_theta_min, math.tan(theta_min))

                if self._V_Ed_1 <= V_Rd_max_min_angle:
                    # If within the minimum angle
                    self._theta = theta_min
                    self._cot_theta = cot_theta_min
                    self._V_Rd_max = V_Rd_max_min_angle
                    self._max_shear_ok = True
                else:
                    # Check the maximum strut angle θ = 45° (cot(θ) = 1.0)
                    theta_max = math.radians(45)
                    cot_theta_max = 1 / math.tan(theta_max)
                    V_Rd_max_max_angle: PlainQuantity = (alpha_cw * self.width * z * v_1 * self._f_cd / (cot_theta_max +
                                                                                           math.tan(theta_max))).to('kN')

                    if self._V_Ed_1 > V_Rd_max_max_angle:
                        self._theta = theta_max
                        self._cot_theta = 1 / math.tan(self._theta)
                        self._V_Rd_max = V_Rd_max_max_angle
                        self._max_shear_ok = False
                    else:
                        self._max_shear_ok = True
                        # Determine the angle θ of the strut based on the shear force
                        self._theta = 0.5 * math.asin((self._V_Ed_1 / V_Rd_max_max_angle))
                        debug(self._theta)
                        self._cot_theta = 1 / math.tan(self._theta)
                        self._V_Rd_max = (alpha_cw * self.width * z * v_1 * self._f_cd / (self._cot_theta +
                                                                                           math.tan(self._theta))).to('kN')
                # Required shear reinforcing area
                self._A_v_req = max((self._V_Ed_2 / (z * self._f_ywd * self._cot_theta)), self._A_v_min)
                self._V_Rd_s = (self._A_v * z * self._f_ywd * self._cot_theta)
                #Maximum shear capacity is the same as the steel capacity
                self._V_Rd = self._V_Rd_s

                # Rebar spacing checks
                section_rebar = Rebar(self)
                n_legs_actual = self._stirrup_n * 2      # Ensure legs are even
                self._s_l = self._stirrup_s_l
                self._s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs_actual - 1) 
                self._s_max_l, self._s_max_w =\
                      section_rebar.calculate_max_spacing_EN_1992_2004(self._alpha)

            self._DCRv = (self._V_Ed_2.to('kN') / self._V_Rd.to('kN'))
            # Design results
            results = {
                'Label': self.label, #Beam label
                'Av,min': self._A_v_min.to('cm ** 2 / m'),  # Minimum shear reinforcement area
                'Av,req': self._A_v_req.to('cm ** 2 / m'), # Required shear reinforcing area
                'Av': self._A_v.to('cm ** 2 / m'),  # Provided stirrup reinforcement per unit length
                'VEd,1': self._V_Ed_1.to('kN'), # Max Vu for the design at the support
                'VEd,2': self._V_Ed_2.to('kN'), # Max Vu for the design at d from the support
                'VRd,c': self._V_Rd_c.to('kN'),  # Concrete contribution to shear capacity
                'VRd,s': self._V_Rd_s.to('kN'),  # Reinforcement contribution to shear capacity
                'VRd': self._V_Rd.to('kN'),  # Total shear capacity
                'VRd,max': self._V_Rd_max.to('kN'),  # Maximum shear capacity
                'VEd,1<VRd,max': self._max_shear_ok,  # Check if applied shear is within max shear capacity
                'VEd,2<VRd': self._V_Ed_2 <= self._V_Rd,  # Check if applied shear is within total capacity
                "DCR" :  self._DCRv
            }
            self._initialize_dicts_EN_1992_2004()
            return pd.DataFrame([results], index=[0])
        else:
            raise ValueError("Concrete type is not compatible with EHE-08 shear check.")
  
    def design_shear_EN_1992_2004(self, Force:Forces, A_s:PlainQuantity = 0*cm**2) -> None:
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
            elif self.concrete.design_code=="EN 1992-2004":
                result =  self.design_shear_EN_1992_2004(force, A_s)
            else:
                raise ValueError(f"Shear design method not implemented for concrete type:"\
                    f"{type(self.concrete).__name__}")
            self._shear_results_list.append(result)
            self._shear_results_detailed_list[force.id] = {
                'forces': self._forces_shear.copy(),
                'shear_reinforcement': self._shear_reinforcement.copy(),
                'min_max': self._data_min_max_shear.copy(),
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
            elif self.concrete.design_code=="EN 1992-2004":
                result =  self.check_shear_EN_1992_2004(force, A_s)
            else:
                raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")  # noqa: E501
            self._shear_results_list.append(result)
            self._shear_results_detailed_list[force.id] = {
                'forces': self._forces_shear.copy(),
                'shear_reinforcement': self._shear_reinforcement.copy(),
                'min_max': self._data_min_max_shear.copy(),
                'checks_pass': self._all_shear_checks_passed,
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

    def _initialize_dicts_ACI_318_19_shear(self) -> None:
        """Initialize the dictionaries used in check and design methods."""
        self._materials_shear = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
                "Concrete density",
                "Normalweight concrete",
                "Safety factor for shear"
            ],
            "Variable": ["","fc", "fy", "wc", "λ", "Øv"],
            "Value": [self.label, round(self.concrete.f_c.to('MPa').magnitude,2), 
                      round(self.steel_bar.f_y.to('MPa').magnitude,2),round(self.concrete.density.to('kg/m**3').magnitude,1),
                       self.settings.get_setting('lambda'), self.settings.get_setting('phi_v')],
            "Unit": ["", "MPa", "MPa", "kg/m³", "", ""]
        }
        self._geometry_shear = {
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
        self._forces_shear = {
            "Design forces": [
                "Axial, positive for compression",
                "Shear",
            ],
            "Variable": ["Nu", "Vu"],
            "Value": [round(self._N_u.to('kN').magnitude,2), round(self._V_u.to('kN').magnitude,2) ],
            "Unit": ["kN", "kN"]
        }
        # Min max lists
        db_min = 10 * mm if self.concrete.unit_system == "metric" else 3 / 8 * inch
        min_values = [None, None, self._A_v_min, db_min]   # Use None for items without a minimum constraint
        max_values = [self._s_max_l, self._s_max_w, None, None]  # Use None for items without a maximum constraint
        current_values = [self._s_l, self._s_w, self._A_v, self._stirrup_d_b]  # Current values to check

        # Generate check marks based on the range conditions
        checks = [
            '✔️' if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else '❌'
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_shear_checks_passed = all(check == '✔️' for check in checks)
        self._data_min_max_shear = {
            'Check': ['Stirrup spacing along length', 'Stirrup spacing along width', 'Minimum shear reinforcement',
                       'Minimum rebar diameter'],
            'Unit': ['cm', 'cm', 'cm²/m', 'mm'],
            'Value': [round(self._s_l.to('cm').magnitude,2), round(self._s_w.to('cm').magnitude,2),
            round(self._A_v.to('cm**2/m').magnitude,2), round(self._stirrup_d_b.magnitude,0)],
            'Min.': ["", "", round(self._A_v_min.to('cm**2/m').magnitude,2), round(db_min.magnitude,0)],
            'Max.': [round(self._s_max_l.to('cm').magnitude,2), round(self._s_max_w.to('cm').magnitude,2), "", ""],
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
                "Shear rebar strength"
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
            "Shear strength": [
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
    
    def _initialize_dicts_EN_1992_2004(self) -> None:
        if isinstance(self.concrete, Concrete_EN_1992_2004):
            """Initialize the dictionaries used in check and design methods."""
            self._materials_shear = {
                "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
                "Safety factor for concrete",
                "Safety factor for steel",
                "Coefficient for long term effects and loading effects"
                ],
            "Variable": ["","fck", "fywk", "γc", "γs", "αcc"],
                "Value": [self.label, round(self.concrete.f_ck.to('MPa').magnitude,2),
                        round(self.steel_bar.f_y.to('MPa').magnitude,2),
                        self._gamma_c, self._gamma_s, self._alpha_cc
                        ],
                "Unit": ["", "MPa", "MPa", "","",""]
            }
            self._geometry_shear = {
                    "Geometry": [
                    "Section height",
                    "Section width",
                    "Clear cover",
                    "Longitudinal tension rebar",
                ],
                "Variable": ["h", "b", "rgeom", "As"],
                #TODO: ver bien tema As de armadura traccionada que podria ser superior o inferior.
                "Value": [self.height.to('cm').magnitude, self.width.to('cm').magnitude, self.c_c.to('cm').magnitude,
                        round(self._A_s.to('cm**2').magnitude,2)],
                "Unit": ["cm", "cm", "cm", "cm²"]
            }
            self._forces_shear = {
                "Design forces": [
                    "Axial, positive for compression",
                    "Shear",
                ],
                "Variable": ["NEd", "VEd,2"],
                "Value": [round(self._N_Ed.to('kN').magnitude,2), round(self._V_Ed_2.to('kN').magnitude,2)],
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
            self._all_shear_checks_passed = all(check == '✔️' for check in checks)
            self._data_min_max_shear = {
                'Check': ['Stirrup spacing along length', 'Stirrup spacing along width', 'Minimum shear reinforcement'],
                'Unit': ['cm', 'cm', 'cm²/m'],
                'Valor': [round(self._s_l.to('cm').magnitude,2), round(self._s_w.to('cm').magnitude,2),
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
                    "Shear rebar strength"
                ],
                "Variable": ["ns", "db", "s", "d", "Asw,min","Asw,req","Asw", "VRd,s"],
                "Value": [self._stirrup_n, self._stirrup_d_b.to('mm').magnitude, self._stirrup_s_l.to('cm').magnitude,
                        self.d.to('cm').magnitude, round(self._A_v_min.to('cm**2/m').magnitude,2),
                        round(self._A_v_req.to('cm**2/m').magnitude,2),
                        round(self._A_v.to('cm**2/m').magnitude,2),
                        round(self._V_Rd_s.to('kN').magnitude,2)],
                "Unit": ["", "mm", "cm", "cm", "cm²/m","cm²/m", "cm²/m","kN"]
            }
            check_max = '✔️' if self._max_shear_ok else '❌'
            check_DCR = '✔️' if self._DCRv < 1 else '❌'
            self._shear_concrete = {
                    "Shear strength": [
                    "Longitudinal reinforcement ratio",
                    'k value',
                    "Axial stress",
                    "Concrete strut angle",
                    "Concrete strength",
                    "Maximum shear strength",
                    "Total shear strength", 
                    "Max shear check",
                    "Demand Capacity Ratio"
                ],
                "Variable": ["ρl", "k", "σcd","Θ","VRd,c", "VRd,max", "VRd","" ,"DCR"],
                "Value": [round(self._rho_l.magnitude,4),
                        round(self._k_value,2),round(self._sigma_cp.to('MPa').magnitude,2),
                        round(math.degrees(self._theta),1),
                        round(self._V_Rd_c.to('kN').magnitude,2), round(self._V_Rd_max.to('kN').magnitude,2), 
                        round(self._V_Rd.to('kN').magnitude,2), check_max, round(self._DCRv,2)],
                "Unit": ["", "", "MPa", "deg", "kN", "kN", "kN", "", check_DCR]
            }
    # Beam results for Jupyter Notebook
    @property
    def data(self) -> None:
        markdown_content = f"Beam {self.label}, $b$={self.width.to('cm')}"\
                         f", $h$={self.height.to('cm')}, $c_{{c}}$={self.c_c.to('cm')}, \
                            Concrete {self.concrete.name}, Rebar {self.steel_bar.name}."
        self._md_data = markdown_content
        # Display the combined content
        display(Markdown(markdown_content))  # type: ignore

        return None
    
    @property
    def flexure_results(self) -> None:
        if not self._flexure_checked:
            warnings.warn("Flexural design has not been performed yet. Call check_flexure or "
                        "design_flexure first.", UserWarning)
            self._md_flexure_results = "Flexural results are not available."
            return None

        # Use limiting case results
        result_data = self._limiting_case_flexure_details
        top_reinforcement = result_data['top_reinforcement']
        bottom_reinforcement = result_data['bottom_reinforcement']
        limiting_forces = result_data['forces']
        flexure_capacity = result_data['flexure_capacity']

        # Create FUFormatter instance and format FU value
        formatter = Formatter()
        formatted_DCR_top = formatter.DCR(flexure_capacity['Value'][0])
        formatted_DCR_bottom = formatter.DCR(flexure_capacity['Value'][1])

        # Format top and bottom reinforcement details
        rebar_top = f"{int(top_reinforcement['Value'][0])}Ø{top_reinforcement['Value'][1]}"
        area_top = top_reinforcement['Value'][2]
        moment_top = limiting_forces['Value'][0]
        capacity_top = flexure_capacity['Value'][0]

        rebar_bottom = f"{int(bottom_reinforcement['Value'][0])}Ø{bottom_reinforcement['Value'][1]}"
        area_bottom = bottom_reinforcement['Value'][2]
        moment_bottom = limiting_forces['Value'][1]
        capacity_bottom = flexure_capacity['Value'][1]

        # Construct markdown output
        markdown_content = (
            f"Bottom longitudinal rebar  {rebar_bottom}, $A_s.bot$ = {area_bottom} cm², "
            f"$M_u$ = {moment_bottom} kNm, $\\phi M_n$ = {capacity_bottom} kNm, "
            f"DCR = {formatted_DCR_bottom}\n"
            f"Top longitudinal rebar  {rebar_top}, $A_s.top$ = {area_top} cm², "
            f"$M_u$ = {moment_top} kNm, $\\phi M_n$ = {capacity_top} kNm, "
            f"DCR = {formatted_DCR_top}"
        )

        self._md_flexure_results = markdown_content
        display(Markdown(markdown_content))

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
        if self._A_v == 0*cm:
            rebar_v = "not assigned"
        else:
            rebar_v = f"{int(limiting_reinforcement['Value'][0])}eØ{limiting_reinforcement['Value'][1]}/"\
                    f"{limiting_reinforcement['Value'][2]} cm"
        # Limitng cases checks 
        warning = "⚠️ Some checks failed, see detailed results." if not self._all_shear_checks_passed else "" 
        markdown_content = f"Shear reinforcing {rebar_v}, $A_v$={limiting_reinforcement['Value'][6]} cm²/m"\
                         f", $V_u$={limiting_forces['Value'][1]} kN, $\\phi V_n$={limiting_shear_concrete['Value'][7]} kN → {formatted_DCR} {warning}"  # noqa: E501

        self._md_shear_results = markdown_content
        display(Markdown(markdown_content))

        return None
    
    # Beam results for Jupyter Notebook
    @property
    def results(self) -> None:
        # Ensure that both properties and shear results are available
        if not hasattr(self, '_md_properties'):
            self.data  # This will generate _md_properties
        if not hasattr(self, '_md_flexure_results'):
            self.flexure_results  # This will generate _md_flexure_results
        if not hasattr(self, '_md_shear_results'):
            self.shear_results  # This will generate _md_shear_results
        # Combine the markdown content for properties and shear results
        markdown_content = f"{self._md_data}\n{self._md_flexure_results}\n{self._md_shear_results}"
        
        # Display the combined content
        display(Markdown(markdown_content))  # type: ignore

        return None

    def flexure_results_detailed(self, force: Optional[Forces] = None) -> None:
        """
        Displays detailed flexure results.

        Parameters
        ----------
        forces : Forces, optional
            The specific Forces object to display results for. If None, displays results for the limiting case.
        Returns
        -------
        None
        """
        if not self._flexure_checked:
            warnings.warn("Flexural check has not been performed yet. Call check_flexure or "
                        "design_flexure first.", UserWarning)
            self._md_flexure_results = "Flexure results are not available."
            return None

        # Determine which results to display (limiting case by default)
        if force:
            force_id = force.id
            if force_id not in self._flexure_results_detailed_list:
                raise ValueError(f"No results found for Forces object with ID {force_id}.")
            result_data = self._flexure_results_detailed_list[force_id]
        else:
            # Default to limiting case
            result_data = self._limiting_case_flexure_details

        # Create TablePrinter instances for detailed display
        materials_printer = TablePrinter("MATERIALS")
        materials_printer.print_table_data(self._materials_flexure, headers='keys')

        geometry_printer = TablePrinter("GEOMETRY")
        geometry_printer.print_table_data(self._geometry_flexure, headers='keys')

        forces_printer = TablePrinter("FORCES")
        forces_printer.print_table_data(result_data['forces'], headers='keys')

        min_max_printer = TablePrinter("MAX AND MIN LIMIT CHECKS")
        min_max_printer.print_table_data(result_data['min_max'], headers='keys')

        capacity_printer = TablePrinter("FLEXURAL CAPACITY - TOP")
        capacity_printer.print_table_data(result_data['flexure_capacity_top'], headers='keys')
        capacity_printer = TablePrinter("FLEXURAL CAPACITY - BOTTOM")
        capacity_printer.print_table_data(result_data['flexure_capacity_bot'], headers='keys')


    def flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        """
        Prints detailed flexure results in Word.

        Parameters
        ----------
        forces : Forces, optional
            The specific Forces object to display results for. If None, displays results for the limiting case.
        """
        if not self._flexure_checked:
            warnings.warn("Flexural check has not been performed yet. Call check_flexure or "
                        "design_flexure first.", UserWarning)
            self._md_flexure_results = "Flexural results are not available."
            return None

        # Determine which results to display (limiting case by default)
        if force:
            force_id = force.id
            if force_id not in self._flexure_results_detailed_list:
                raise ValueError(f"No results found for Forces object with ID {force_id}.")
            result_data = self._flexure_results_detailed_list[force_id]
        else:
            # Default to limiting case
            result_data = self._limiting_case_flexure_details

        # Convert output Dicts into DataFrames
        df_materials = pd.DataFrame(self._materials_flexure)
        df_geometry = pd.DataFrame(self._geometry_flexure)
        df_forces = pd.DataFrame(result_data['forces'])
        df_data_min_max = pd.DataFrame(result_data['min_max'])
        df_flexure_capacity_top = pd.DataFrame(result_data['flexure_capacity_top'])
        df_flexure_capacity_bottom = pd.DataFrame(result_data['flexure_capacity_bottom'])
        
        # Create a document builder instance
        doc_builder = DocumentBuilder(title='Concrete beam flexure check')

        # Add first section and table
        doc_builder.add_heading('Concrete beam flexure check', level=1)
        doc_builder.add_text(f'Design code: {self.concrete.design_code}')
        doc_builder.add_heading('Materials', level=2)
        doc_builder.add_table_data(df_materials)
        doc_builder.add_table_data(df_geometry)
        doc_builder.add_table_data(df_forces)

        # Add third section for limit checks
        doc_builder.add_heading('Limit Checks', level=2)
        doc_builder.add_table_data(df_data_min_max)

        # Add second section for flexural checks
        doc_builder.add_heading('Flexural Capacity Top', level=2)
        doc_builder.add_table_data(df_flexure_capacity_top)
        doc_builder.add_heading('Flexural Capacity Bottom', level=2)
        doc_builder.add_table_data(df_flexure_capacity_bottom)

        # Save the Word doc
        doc_builder.save(f"Concrete beam flexure check {self.concrete.design_code}.docx")

    
    def shear_results_detailed(self, force: Optional[Forces] = None) -> None:
        """
        Displays detailed shear results.

        Parameters
        ----------
        forces : Forces, optional
            The specific Forces object to display results for. If None, displays results for the limiting case.
        Returns
        -------
        None
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
        materials_printer = TablePrinter("MATERIALS")
        materials_printer.print_table_data(self._materials_shear, headers='keys')
        geometry_printer = TablePrinter("GEOMETRY")
        geometry_printer.print_table_data(self._geometry_shear, headers='keys')
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
        df_materials = pd.DataFrame(self._materials_shear)
        df_geometry = pd.DataFrame(self._geometry_shear)
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

def flexure_design_test() -> None:
    #Example 6.6 CRSI GUIDE
    concrete = Concrete_ACI_318_19(name="C4",f_c=4000*psi) 
    steelBar = SteelBar(name="G60", f_y=60000*psi) 
    beam = RectangularBeam(
        label="101",
        concrete=concrete,
        steel_bar=steelBar,
        width=12 * inch,  
        height=24 * inch,   
    )
    f1 = Forces(label='D', M_y=258.3*kip*ft)
    f2 = Forces(label='L', M_y=58.3*kip*ft)
    Node(section=beam, forces=[f1, f2])
    results=beam.design_flexure() 
    print(results)

def flexure_check_test() -> None:
    concrete = Concrete_ACI_318_19(name="H-25",f_c=25*MPa) 
    steelBar = SteelBar(name="420", f_y=420*MPa) 

    beam = RectangularBeam(
        label="101",
        concrete=concrete,
        steel_bar=steelBar,
        width=20*cm,  
        height=60*cm,   
    )

    beam.set_longitudinal_rebar_bot(n1=2,d_b1=25*mm)
    f1 = Forces(label='D', M_y=100*kNm)
    f2 = Forces(label='D', M_y=120*kNm)
    results=beam.check_flexure_ACI_318_19(f1)
    print(results)
    results=beam.check_flexure_ACI_318_19(f2)
    print(results)


def flexure_Mn() -> None:
    # MOMENTO NOMINAL SIMPLEMENTE ARMADO
    #Example from https://www.google.com/search?q=ACI+318+19+nominal+moment+of+section&oq=ACI+318+19+nominal+moment+of+section&gs_lcrp=EgZjaHJvbWUyBggAEEUYOTIHCAEQIRigAdIBCTIyNzkzajBqN6gCALACAA&sourceid=chrome&ie=UTF-8#fpstate=ive&vld=cid:dab5f5e7,vid:m4H0QbGDYIg,st:0
    concrete = Concrete_ACI_318_19(name="C6",f_c=6000*psi) 
    steelBar = SteelBar(name="G80", f_y=80000*psi) 
    beam = RectangularBeam(
        label="B20x30",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * inch,  
        height=30 * inch,   
    )
    A_s=10.92 * inch**2
    d=27*inch
    c_mec=3*inch
    result=beam._determine_nominal_moment_simple_reinf_ACI_318_19(A_s,d,c_mec)
    print(result.to(kip*ft))


    # MOMENTO NOMINAL DOBLEMENTE ARMADO
    concrete2 = Concrete_ACI_318_19(name="C4",f_c=4000*psi) 
    steelBar2 = SteelBar(name="G60", f_y=60000*psi) 
    beam2 = RectangularBeam(
        label="B14x27",
        concrete=concrete2,
        steel_bar=steelBar2,
        width=14 * inch,  
        height=27 * inch,   
    )
    A_s=6 * inch**2
    d=24*inch
    c_mec=3*inch
    d_prime=2.5*inch
    A_s_prime=1.8*inch**2
    result=beam2._determine_nominal_moment_double_reinf_ACI_318_19(A_s, d, c_mec, d_prime, A_s_prime)
    print(result.to(kip*ft))

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
    beam.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=20*cm) 
    results = beam.check_shear(A_s=8.04*cm**2)
    # results = beam.design_shear(A_s=5*cm**2)
    print(results)
    # beam.shear_results_detailed()
    # print(beam.shear_design_results)
    # print(beam.results)
    beam.shear_results_detailed_doc()

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

def shear_EN_1992() -> None:
    concrete= Concrete_EN_1992_2004(name="C25",f_ck=25*MPa) 
    steelBar= SteelBar(name="B500S", f_y=500*MPa)
    custom_settings = {'clear_cover': 2.6*cm, 'stirrup_diameter_ini':8*mm,
                       'longitudinal_diameter_ini': 16*mm}
    beam = RectangularBeam(label="101",
                                      concrete=concrete,steel_bar=steelBar,width=20*cm, height=60*cm,
                                       settings=custom_settings)
    # f = Forces(V_z=100*kN, M_y=100*kNm)
    f = Forces(V_z=500*kN, N_x=0*kN)
    A_s = 8.04*cm**2
    beam.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=25*cm)
    Node(beam, f)
    beam.check_shear(A_s)
    beam.shear_results_detailed()
    beam.shear_results_detailed_doc()

if __name__ == "__main__":
    flexure_check_test()
    # flexure_Mn()
    # shear_ACI_imperial()
    shear_EN_1992()
    # shear_ACI_metric()
    # rebar()