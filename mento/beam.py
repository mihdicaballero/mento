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
        self.layers_spacing = self.settings.get_setting('layers_spacing')
                    
    
    @property
    def d(self) -> PlainQuantity:
        "Effective height."
        return self._d
    
    def initialize_code_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_ACI_318_19):
            self._initialize_aci_318_attributes()
        elif isinstance(self.concrete, Concrete_EN_1992_2004):
            self._initialize_en_1992_attributes()

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


    def _initialize_en_1992_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_EN_1992_2004):
            self._f_yk = self.steel_bar.f_y
            self._f_ck = self.concrete.f_ck
            self._f_cd = self.concrete.f_cd
            self._V_Ed_1: PlainQuantity = 0*kN
            self._V_Ed_2: PlainQuantity = 0*kN
            self._N_Ed: PlainQuantity = 0*kN
            self._M_Ed: PlainQuantity = 0*kNm
            self._sigma_cd: PlainQuantity = 0*MPa
            self._V_Rd_c: PlainQuantity = 0*kN
            self._V_Rd_s: PlainQuantity = 0*kN
            self._V_Rd_max: PlainQuantity = 0*kN
            self._k_value:float = 0
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
        self._bot_rebar_centroid=self.__calculate_long_rebar_centroid(n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4)

    def set_longitudinal_rebar_top(self, n1: int, d_b1: PlainQuantity, n2: int, d_b2: PlainQuantity, 
                                n3: int, d_b3: PlainQuantity, n4: int, d_b4: PlainQuantity, 
                                total_as: PlainQuantity, total_bars: int, clear_spacing: PlainQuantity) -> None:
        """Sets the longitudinal rebar in the object."""
        self._longitudinal_n1_t = n1
        self._longitudinal_d_b1_t = d_b1
        self._longitudinal_n2_t = n2
        self._longitudinal_d_b2_t = d_b2 if n2 > 0 else None
        self._longitudinal_n3_t = n3
        self._longitudinal_d_b3_t = d_b3 if n3 > 0 else None
        self._longitudinal_n4_t = n4
        self._longitudinal_d_b4_t = d_b4 if n4 > 0 else None
        self._total_as_t = total_as
        self._total_bars_t = total_bars
        self._clear_spacing_t = clear_spacing
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

    def __maximum_flexural_reinforcement_ratio(self) -> float:
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
            )  # TODO: REVISAR
            
            # Calculate maximum reinforcement ratio (ρ_max)
            rho_max = (
                0.85 * concrete_aci.beta_1 * self.concrete.f_c / self.steel_bar.f_y
                * (concrete_aci.epsilon_c /
                (concrete_aci.epsilon_c + epsilon_min_rebar_ACI_318_19))
            )
            return rho_max
        else:
            return 0

    def __minimum_flexural_reinforcement_ratio(self)->PlainQuantity:
        concrete_properties = self.concrete.get_properties()
        f_c = concrete_properties['f_c']
        rebar_properties = self.steel_bar.get_properties()
        f_y = rebar_properties['f_y']
        minimum_ratio = max((3 * np.sqrt(f_c / psi) * psi / f_y), 
                    (200 * psi / f_y))
        return minimum_ratio
    

    def __calculate_phi_ACI_318_19(self, epsilon_most_strained: float) -> float:
        """
        Calculates the strength reduction factor (φ) for flexural design 
        based on ACI 318-19.

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
        rho_min=self.__minimum_flexural_reinforcement_ratio()
        A_s_min = rho_min * d * b
        rho_max = self.__maximum_flexural_reinforcement_ratio()
        A_s_max = rho_max * d * b

        # Calculate required reinforcement
        R_n = M_u / (phi * b * d**2)
        A_s_calc = 0.85 * f_c * b * d / f_y * (1 - np.sqrt(1 - 2 * R_n / (0.85 * f_c)))

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


    def _compile_results_aci_flexure_metric(self, force: Forces, A_s_min_bot: PlainQuantity, A_s_req_bot: PlainQuantity, A_s_bot: PlainQuantity,
                                        c_d_bot: float, M_u_bot: PlainQuantity, phi_M_n_bot: PlainQuantity,
                                        A_s_min_top: PlainQuantity, A_s_req_top: PlainQuantity, A_s_top: PlainQuantity,
                                        c_d_top: float, M_u_top: PlainQuantity, phi_M_n_top: PlainQuantity) -> List[Dict[str, Any]]:
        # Create dictionaries for bottom and top rows
        bottom_result = {
            'Section Label': self.label,
            'Load Combo': force.label,
            'Position': 'Bottom',
            'As,min': A_s_min_bot.to('cm ** 2'),
            'As,req': A_s_req_bot.to('cm ** 2'),
            'As': A_s_bot.to('cm ** 2'),
            'c/d': c_d_bot,
            'Mu': M_u_bot.to('kN*m'),
            'ØMn': phi_M_n_bot.to('kN*m'),
            'Mu<ØMn': M_u_bot <= phi_M_n_bot,
            'DCR': M_u_bot / phi_M_n_bot
        }
        
        top_result = {
            'Section Label': self.label,
            'Load Combo': force.label,
            'Position': 'Top',
            'As,min': A_s_min_top.to('cm ** 2'),
            'As,req': A_s_req_top.to('cm ** 2'),
            'As': A_s_top.to('cm ** 2'),
            'c/d': c_d_top,
            'Mu': M_u_top.to('kN*m'),
            'ØMn': phi_M_n_top.to('kN*m'),
            'Mu<ØMn': M_u_top <= phi_M_n_top,
            'DCR': M_u_top / phi_M_n_top
        }
    
        # Return both rows as a list of dictionaries
        return [bottom_result, top_result]


    def check_flexure_ACI_318_19(self, Force: Forces) -> Dict[str, Any]:
        '''
            flexure_design_results_top: es el diccionario de armado de rebar que tiene las 4 posiciones
        '''

        # Load design settings for ACI 318-19
        self.settings.load_aci_318_19_settings()

        # Initial assumptions for mechanical cover and compression depth

        rec_mec = self.c_c + self._stirrup_d_b + self._bot_rebar_centroid
        d_prima = self.c_c + self._stirrup_d_b + self._top_rebar_centroid
        d = self._height - rec_mec


        for force in Force:
            # RECORREMOS TODAS LAS FUERZAS ASI SACAMOS LOS DCR RATIO PARA CADA UNA
            if force._M_y > 0:
                (self._A_s_min_bot, self._A_s_max_bot, self._A_s_final_bot_Positive_M,
                self._A_s_comp_top) = self.__calculate_flexural_reinforcement(force._M_y, d, d_prima)

            else:
                (self._A_s_min_top, self._A_s_max_top, self._A_s_final_top_Negative_M,
                self._A_s_comp_bot) = self.__calculate_flexural_reinforcement(abs(force._M_y), d_prima, d)

        results=self._compile_results_aci_flexure_metric()

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

            rec_mec_calculo = self.c_c + self._stirrup_d_b + self._bot_rebar_centroid

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
            Err = max(abs(rec_mec_calculo - rec_mec), abs(d_prima_calculo - d_prima))
            rec_mec = rec_mec_calculo
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

    def design_flexure_EHE_08(self, M_u: float) -> None:
        pass

    def design_flexure(self) ->  DataFrame:
        if not self.node or not self.node.forces:
            raise ValueError("No Node or forces list associated with this beam.")
        self._flexure_results_list = []  # Store individual results for each force
        self._flexure_results_detailed_list = {}  # Store detailed results by force ID
        max_A_s_bottom_req = 0*cm**2 # Track the maximum A_s_bottom_req to identify the limiting case
        max_A_s_top_req = 0*cm**2 # Track the maximum A_s_top_req to identify the limiting case

        for force in self.node.forces:
            if self.concrete.design_code=="ACI 318-19":
                result = self.design_flexure_ACI_318_19(force)
            # elif self.concrete.design_code=="EN 1992":
            #     result = self.design_flexure_EN_1992(force)
            else:
                raise ValueError(f"Longitudinal design method not implemented \
                        for concrete type: {type(self.concrete).__name__}")
            self._flexure_results_list.append(result)
            #TODO: Esto falta ajustar a flexión
            # self._flexure_results_detailed_list[force.id] = {
            #     'forces': self._forces.copy(),
            #     'shear_reinforcement': self._shear_reinforcement.copy(),
            #     'min_max': self._data_min_max.copy(),
            #     'shear_concrete': self._shear_concrete.copy(),
            # }
            # Check if this result is the limiting case
            # current_A_v_req = result['Av,req'][0]
            # if current_A_v_req > max_A_v_req:
            #     max_A_v_req = current_A_v_req
            #     self._limiting_case_shear = result
            #     self._limiting_case_shear_details = self._shear_results_detailed_list[force.id]

            #     # Update shear design results for the worst case
            #     section_rebar = Rebar(self)
            #     self.shear_design_results = section_rebar.transverse_rebar(self._A_v_req, self._V_s_req)
            #     self._best_rebar_design = section_rebar.transverse_rebar_design

        # Compile all results into a single DataFrame
        all_results = pd.concat(self._flexure_results_list, ignore_index=True)

        # Identify the most limiting case by Av,required
        # self.limiting_case_shear = all_results.loc[all_results['Av,req'].idxmax()]  # Select row with highest Av,req
        # Mark shear as checked
        self._flexure_checked = True  
        return all_results

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
        self._initialize_dicts_ACI_318_19_shear()
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
            self.settings.load_en_1992_settings()
            self._gamma_c = self.settings.get_setting('gamma_c')
            self._gamma_s = self.settings.get_setting('gamma_s')
            self._f_ywk = self._f_yk
            self._f_ywd = self._f_ywk/self._gamma_s

            # Minimum shear reinforcement calculation
            self._A_v_min = 0.08*math.sqrt(self._f_ck.to('MPa').magnitude) / (self._f_ywk)*MPa

            # Compression stress, positive
            self._A_p = 0*cm**2 # No prestressing for now
            self._rho_l = min((A_s + self._A_p) / (self.width * self.d), 0.02)

            # Shear calculation for sections without rebar
            self._k_value = min(1 + math.sqrt(200 * mm / self.d), 2)

    def _calculate_V_u1(self) -> PlainQuantity:
        self._alpha = math.radians(90)
        self._theta = math.radians(45)
        self._cot_theta = 1 / math.tan(self._theta)
        self._cot_alpha = 1 / math.tan(self._alpha)
        self._sigma_cd =  self._N_rd / self.A_x # Without compression reinforcement considered 
        self._K_value = self._calculate_axial_coefficient_ehe(self._sigma_cd, self._f_cd)
        return self._K_value * self._f_1cd * self.width * self.d\
              * (self._cot_theta + self._cot_alpha) / (1 + self._cot_theta ** 2)

    def _shear_without_rebar_EN_1992_2004(self) -> PlainQuantity:
        
        self._stirrup_d_b = 0*mm
        self._A_v_min = 0*cm**2/m
        # Positive of compression
        self._sigm_cp = min(self._N_Ed / self.A_x, 0.2*self._f_cd)

        # Total shear capacity without rebar
        C_rdc = 0.18/self._gamma_c
        v_min = 0.035*self._k_value**(3/2)*math.sqrt(self._f_ck.to('MPa').magnitude)
        k_1 = 0.15
        V_Rd_c_min = (v_min+k_1*self._sigm_cp.to('MPa').magnitude)* self.width * self.d * MPa
        V_Rd_c = (C_rdc*self._k_value*(100*self._rho_l*self._f_ck.to('MPa').magnitude)**(1/3)*MPa\
                  +k_1*self._sigm_cp.to('MPa').magnitude)* self.width * self.d * MPa
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
                self._max_shear_ok = self._V_Ed_1 < self._V_Rd_c

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
                alpha_cw = 1 #For non-prestressed members or members subject to tensile stress due to axial force
                v_1 = 0.6*(1 - self._f_ck.to('MPa').magnitude/250)
                # The θ angle is lmited between 21,8° ≤ θ ≤ 45°(1 ≤ cot(θ) ≤ 2.5)
                self._theta = 21.8*deg # Cracks angle (assumed 45 degrees)
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
                V_s_req = self._V_Ed_2 - self._V_cu

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
                'Vrd,1': self._V_Ed_1.to('kN'), # Max Vu for the design at the support
                'Vrd,2': self._V_Ed_2.to('kN'), # Max Vu for the design at d from the support
                'Vcu': self._V_cu.to('kN'),  # Concrete contribution to shear capacity
                'Vsu': self._V_su.to('kN'),  # Reinforcement contribution to shear capacity
                'Vu2': self._V_u2.to('kN'),  # Total shear capacity
                'Vu1': self._V_u1.to('kN'),  # Maximum shear capacity
                'Vrd,1<Vu1': self._max_shear_ok,  # Check if applied shear is within max shear capacity
                'Vrd,2<Vu2': self._V_Ed_2 <= self._V_u2,  # Check if applied shear is within total capacity
                "DCR" :  self._FUv
            }
            self._initialize_dicts_EHE_08()
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

    def _initialize_dicts_ACI_318_19_shear(self) -> None:
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
        rebar_v = f"{int(limiting_reinforcement['Value'][0])}eØ{limiting_reinforcement['Value'][1]}/"\
        f"{limiting_reinforcement['Value'][2]} cm"

        markdown_content = f"Shear reinforcing {rebar_v}, $A_v$={limiting_reinforcement['Value'][6]} cm²/m"\
                         f", $V_u$={limiting_forces['Value'][1]} kN, $\\phi V_n$={limiting_shear_concrete['Value'][7]} kN → {formatted_DCR}"  # noqa: E501

        self._md_shear_results = markdown_content

        return None
    
    # Beam results for Jupyter Notebook
    @property
    def results(self) -> None:
        # Ensure that both properties and shear results are available
        if not hasattr(self, '_md_properties'):
            self.properties  # This will generate _md_properties
        if not hasattr(self, '_md_flexure_results'):
            self.flexure_results  # This will generate _md_flexure_results
        if not hasattr(self, '_md_shear_results'):
            self.shear_results  # This will generate _md_shear_results
        # Combine the markdown content for properties and shear results
        markdown_content = f"{self._md_properties}\n{self._md_flexure_results}\n{self._md_shear_results}"
        
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
        materials_printer.print_table_data(self._materials, headers='keys')

        geometry_printer = TablePrinter("GEOMETRY")
        geometry_printer.print_table_data(self._geometry, headers='keys')

        forces_printer = TablePrinter("FORCES")
        forces_printer.print_table_data(result_data['forces'], headers='keys')

        reinforcement_printer = TablePrinter("FLEXURAL REINFORCEMENT")
        reinforcement_printer.print_table_data(result_data['flexure_reinforcement'], headers='keys')

        capacity_printer = TablePrinter("FLEXURAL CAPACITY")
        capacity_printer.print_table_data(result_data['flexure_capacity'], headers='keys')

        min_max_printer = TablePrinter("MAX AND MIN LIMIT CHECKS")
        min_max_printer.print_table_data(result_data['min_max'], headers='keys')

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
        df_materials = pd.DataFrame(self._materials)
        df_geometry = pd.DataFrame(self._geometry)
        df_forces = pd.DataFrame(result_data['forces'])
        df_flexure_reinforcement = pd.DataFrame(result_data['flexure_reinforcement'])
        df_flexure_capacity = pd.DataFrame(result_data['flexure_capacity'])
        df_data_min_max = pd.DataFrame(result_data['min_max'])

        # Create a document builder instance
        doc_builder = DocumentBuilder(title='Concrete beam flexure check')

        # Add first section and table
        doc_builder.add_heading('Concrete beam flexure check', level=1)
        doc_builder.add_text(f'Design code: {self.concrete.design_code}')
        doc_builder.add_heading('Materials', level=2)
        doc_builder.add_table_data(df_materials)
        doc_builder.add_table_data(df_geometry)
        doc_builder.add_table_data(df_forces)

        # Add second section for flexural checks
        doc_builder.add_heading('Flexural Reinforcement', level=2)
        doc_builder.add_table_data(df_flexure_reinforcement)
        doc_builder.add_heading('Flexural Capacity', level=2)
        doc_builder.add_table_data(df_flexure_capacity)

        # Add third section for limit checks
        doc_builder.add_heading('Limit Checks', level=2)
        doc_builder.add_table_data(df_data_min_max)

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
        materials_printer.print_table_data(self._materials, headers='keys')
        geometry_printer = TablePrinter("GEOMETRY")
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

def shear_EN_1992() -> None:
    concrete= Concrete_EN_1992_2004(name="C25",f_ck=25*MPa) 
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
    flexure()
    # shear_ACI_imperial()
    # shear_EHE_08()
    #shear_ACI_metric()
    # rebar()