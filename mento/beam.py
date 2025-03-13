import os # Cleaning console
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
from mento.material import Concrete, SteelBar, Concrete_ACI_318_19, Concrete_EN_1992_2004, Concrete_CIRSOC_201_25
from mento.rebar import Rebar
from mento.units import MPa, psi, mm, inch, kN, m, cm, kNm, dimensionless
from mento.results import Formatter, TablePrinter, DocumentBuilder
from mento.forces import Forces  

from mento.codes.EN_1992_2004_beam import _check_shear_EN_1992_2004
from mento.codes.ACI_318_19_beam import _check_shear_ACI_318_19, _design_shear_ACI_318_19


@dataclass
class RectangularBeam(RectangularSection):
    label: Optional[str] = None

    def __init__(self, label: Optional[str], concrete: Concrete, steel_bar: SteelBar, 
                 width: PlainQuantity, height: PlainQuantity, settings: Optional[Dict[str, Any]] = None):   
        super().__init__(concrete, steel_bar, width, height, settings)
        if settings:
            self.settings.update(settings)  # Update with any provided settings
        self.label = label
        self.layers_spacing = self.settings.get_setting('layers_spacing')
        
        # Centralized attribute initialization
        self._initialize_attributes()
        
        # self.shear_design_results: DataFrame = None

    def _initialize_attributes(self) -> None:
        """Initialize all attributes of the beam."""
        # Stirrups and shear attributes
        self._stirrup_s_l: PlainQuantity = 0 * cm
        self._stirrup_s_w: PlainQuantity = 0 * cm
        self._stirrup_s_max_l: PlainQuantity = 0 * cm
        self._stirrup_s_max_w: PlainQuantity = 0 * cm
        self._stirrup_n: int = 0
        self._A_v_min: PlainQuantity = 0 * cm**2/m
        self._A_v: PlainQuantity = 0 * cm**2 / m
        self._A_s_req_bot: PlainQuantity = 0 * cm**2
        self._A_s_req_top: PlainQuantity = 0 * cm**2
        self._A_v_req: PlainQuantity = 0 * cm**2 / m
        self._DCRv: float = 0
        self._DCRb_top: float = 0
        self._DCRb_bot: float = 0

        # Design checks and effective heights
        self._rho_l_bot: PlainQuantity = 0 * dimensionless
        self._rho_l_top: PlainQuantity = 0 * dimensionless
        self._bot_rebar_centroid = 0 * mm
        self._top_rebar_centroid = 0 * mm
        self._c_d_top: float = 0
        self._c_d_bot: float = 0
        self._shear_checked = False  # Tracks if shear check or design has been done
        self._flexure_checked = False  # Tracks if shear check or design has been done

        # Initialize default concrete beam attributes
        self._initialize_code_attributes()
        # Longitudinal rebar attributes
        self._initialize_longitudinal_rebar_attributes()

        # Results attributes
        self._materials_shear: Dict = {}
        self._geometry_shear: Dict = {}
        self._forces_shear: Dict = {}
        self._all_shear_checks_passed: bool = False
        self._data_min_max_shear: Dict = {} 
        self._shear_reinforcement: Dict = {}
        self._shear_concrete: Dict = {}
        self._shear_all_checks: bool = False

    def _initialize_code_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_ACI_318_19):
            self._initialize_ACI_318_attributes()
        elif isinstance(self.concrete, Concrete_EN_1992_2004):
            self._initialize_en_1992_2004_attributes()

    def _initialize_longitudinal_rebar_attributes(self) -> None:
        """Initialize all rebar-related attributes with default values."""
        # Bottom rebar defaults
        self._n2_b, self._d_b2_b = 0, 0 * mm
        self._n3_b, self._d_b3_b = 0, 0 * mm
        self._n4_b, self._d_b4_b = 0, 0 * mm

        # Top rebar defaults
        self._n2_t, self._d_b2_t = 0, 0 * mm
        self._n3_t, self._d_b3_t = 0, 0 * mm
        self._n4_t, self._d_b4_t = 0, 0 * mm

        # Unit system default minimum rebar
        if self.concrete.unit_system == "metric":
            self._n1_b, self._d_b1_b = 2, 8 * mm
            self._n1_t, self._d_b1_t = 2, 8 * mm
        else:
            self._n1_b, self._d_b1_b = 2, 3/8 * inch
            self._n1_t, self._d_b1_t = 2, 3/8 * inch

        # Update dependent attributes
        self._update_longitudinal_rebar_attributes()

    def _initialize_ACI_318_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_ACI_318_19):
            self._phi_V_n: PlainQuantity = 0*kN
            self._phi_V_s: PlainQuantity = 0*kN
            self._phi_V_c: PlainQuantity = 0*kN
            self._phi_V_max: PlainQuantity = 0*kN
            self._V_u: PlainQuantity = 0*kN
            self._M_u: PlainQuantity = 0*kNm
            self._M_u_bot: PlainQuantity = 0*kNm
            self._M_u_top: PlainQuantity = 0*kNm
            self._N_u: PlainQuantity = 0*kN
            self._A_cv: PlainQuantity = 0*cm**2
            self._k_c_min: PlainQuantity = 0*MPa
            self._sigma_Nu: PlainQuantity = 0*MPa
            self.V_c: PlainQuantity = 0*kN
            self._V_s_req: PlainQuantity = 0*kN
            self.phi_v: float = 0
            self.phi_t: float = 0
            self.lambda_factor = 0
            self._rho_w: PlainQuantity = 0*dimensionless
            self._lambda_s: float = 0
            self.f_yt: PlainQuantity = 0*MPa
            self._A_s_tension: PlainQuantity = 0*cm**2
            self._max_shear_ok: bool =False

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
            self._V_Rd: PlainQuantity = 0*kN
            self._k_value: float = 0
            self._alpha_cc:float  = 0
            self._gamma_c: float = 0
            self._gamma_s:float  = 0
            self._f_ywk: PlainQuantity = 0*MPa
            self._f_ywd: PlainQuantity = 0*MPa
            self._f_cd: PlainQuantity = 0*MPa
            self._alpha: float = 0
            self._A_p = 0*cm**2 # No prestressing for now
            self._sigma_cp: PlainQuantity = 0*MPa
            self._theta: float = 0
            self._cot_theta: float = 0

    def _update_longitudinal_rebar_attributes(self) -> None:
        """Recalculate attributes dependent on rebar configuration for both top and bottom reinforcing."""
        self._calculate_longitudinal_rebar_area()
        self._calculate_long_rebar_centroid()
        self._calculate_min_clear_spacing()
        self._update_effective_heights()

    def _update_effective_heights(self) -> None:
        """Update effective heights and depths for moment and shear calculations."""
        self._c_mec_bot = self.c_c + self._stirrup_d_b + self._bot_rebar_centroid
        self._c_mec_top = self.c_c + self._stirrup_d_b + self._top_rebar_centroid
        self._d_bot = self._height - self._c_mec_bot
        self._d_top = self._height - self._c_mec_top
        # Use bottom or top effective height
        self._d_shear = min(self._d_bot, self._d_top)

    def set_transverse_rebar(self, n_stirrups: int = 0, d_b:PlainQuantity = 0*mm, s_l:PlainQuantity = 0*cm) -> None:
        """Sets the transverse rebar in the beam section."""
        self._stirrup_n = n_stirrups
        self._stirrup_d_b = d_b
        self._stirrup_s_l = s_l
        n_legs = n_stirrups * 2
        A_db = (d_b ** 2) * math.pi / 4  # Area of one stirrup leg
        A_vs = n_legs * A_db  # Total area of stirrups
        self._A_v = A_vs / s_l  # Stirrup area per unit length

        # Update effective heights
        self._update_effective_heights()

    def set_longitudinal_rebar_bot(self, n1: int=0, d_b1: PlainQuantity=0*mm, n2: int = 0, d_b2: PlainQuantity=0*mm, 
                                n3: int=0, d_b3: PlainQuantity=0*mm, n4: int=0, d_b4: PlainQuantity=0*mm, 
                                ) -> None:
        """Update the bottom rebar configuration and recalculate attributes."""
        self._n1_b = n1 or self._n1_b
        self._d_b1_b = d_b1 or self._d_b1_b
        self._n2_b = n2 or self._n2_b
        self._d_b2_b = d_b2 or self._d_b2_b
        self._n3_b = n3 or self._n3_b
        self._d_b3_b = d_b3 or self._d_b3_b
        self._n4_b = n4 or self._n4_b
        self._d_b4_b = d_b4 or self._d_b4_b
        self._update_longitudinal_rebar_attributes()

    def set_longitudinal_rebar_top(self, n1: int, d_b1: PlainQuantity, n2: int=0, d_b2: PlainQuantity=0*mm, 
                                n3: int=0, d_b3: PlainQuantity = 0*mm, n4: int=0, d_b4: PlainQuantity=0*mm 
                                ) -> None:
        """Update the top rebar configuration and recalculate attributes."""
        self._n1_t = n1 or self._n1_t
        self._d_b1_t = d_b1 or self._d_b1_t
        self._n2_t = n2 or self._n2_t
        self._d_b2_t = d_b2 or self._d_b2_t
        self._n3_t = n3 or self._n3_t
        self._d_b3_t = d_b3 or self._d_b3_t
        self._n4_t = n4 or self._n4_t
        self._d_b4_t = d_b4 or self._d_b4_t
        self._update_longitudinal_rebar_attributes()
    
    def _calculate_longitudinal_rebar_area(self) -> None:
        """Calculate the total rebar area for a given configuration."""

        # AREA OF BOTTOM BARS
        self._A_s_bot = (self._n1_b * self._d_b1_b**2 * np.pi / 4 +
                self._n2_b * self._d_b2_b**2 * np.pi / 4 +
                self._n3_b * self._d_b3_b**2 * np.pi / 4 +
                self._n4_b * self._d_b4_b**2 * np.pi / 4)
    
        # AREA OF TOP BARS
        self._A_s_top = (self._n1_t * self._d_b1_t**2 * np.pi / 4 +
                self._n2_t * self._d_b2_t**2 * np.pi / 4 +
                self._n3_t * self._d_b3_t**2 * np.pi / 4 +
                self._n4_t * self._d_b4_t**2 * np.pi / 4)
        
    def _calculate_min_clear_spacing(self) -> None:
            """
            Calculates the maximum clear spacing between bars for the bottom rebar layers.

            Returns:
                PlainQuantity: The maximum clear spacing between bars in either the first or second layer.
            """

            def layer_clear_spacing(n_a: int, d_a: PlainQuantity, n_b: int, d_b: PlainQuantity) -> PlainQuantity:
                """
                Helper function to calculate clear spacing for a given layer.

                Parameters:
                    n_a (int): Number of bars in the first group of the layer.
                    d_a (PlainQuantity): Diameter of bars in the first group of the layer.
                    n_b (int): Number of bars in the second group of the layer.
                    d_b (PlainQuantity): Diameter of bars in the second group of the layer.

                Returns:
                    PlainQuantity: Clear spacing for the given layer.
                """
                effective_width = self.width - 2 * (self.c_c + self._stirrup_d_b)
                total_bars = n_a + n_b
                if total_bars <= 1:
                    return effective_width - max(d_a, d_b)  # Clear space for one bar
                total_bar_width = n_a * d_a + n_b * d_b
                return (effective_width - total_bar_width) / (total_bars - 1)

            # AVAIABLE CLEAR SPACING FOR BOTTOM BARS
            # Calculate clear spacing for each layer
            spacing_layer1_b = layer_clear_spacing(self._n1_b, self._d_b1_b, self._n2_b, self._d_b2_b)
            spacing_layer2_b = layer_clear_spacing(self._n3_b, self._d_b3_b, self._n4_b, self._d_b4_b)

            # Return the maximum clear spacing between the two layers
            self._available_s_bot=min(spacing_layer1_b, spacing_layer2_b)

            # AVAIABLE CLEAR SPACING FOR TOP BARS
            # Calculate clear spacing for each layer
            spacing_layer1_t = layer_clear_spacing(self._n1_t, self._d_b1_t, self._n2_t, self._d_b2_t)
            spacing_layer2_t = layer_clear_spacing(self._n3_t, self._d_b3_t, self._n4_t, self._d_b4_t)

            # Return the maximum clear spacing between the two layers
            self._available_s_top=min(spacing_layer1_t, spacing_layer2_t)

    def _calculate_long_rebar_centroid(self) -> None:
        """
        Calculates the centroid (baricenter) of a group of rebars based on their diameters, quantities, 
        and layer spacing.

        Returns:
            float: The calculated centroid height of the rebar group.
        """
        # BOTTOM BARS CENTROID
        # Calculate the vertical positions of the bar layers
        y1_b = self._d_b1_b / 2
        y2_b = self._d_b2_b / 2
        y3_b = max(self._d_b1_b, self._d_b2_b) + self.layers_spacing + self._d_b3_b / 2
        y4_b = max(self._d_b1_b, self._d_b2_b) + self.layers_spacing + self._d_b4_b / 2
        # Calculate the total area of each layer
        area_1_b = self._n1_b * self._d_b1_b**2*np.pi/4  # Area proportional to number of bars and their diameter
        area_2_b = self._n2_b * self._d_b2_b**2*np.pi/4
        area_3_b = self._n3_b * self._d_b3_b**2*np.pi/4
        area_4_b = self._n4_b * self._d_b4_b**2*np.pi/4

        # Calculate the centroid as a weighted average
        total_area_b = area_1_b + area_2_b + area_3_b + area_4_b
        if total_area_b == 0:
            return 0*mm  # Avoid division by zero if no bars are present

        self._bot_rebar_centroid = (area_1_b * y1_b + area_2_b * y2_b + area_3_b * y3_b + area_4_b * y4_b) / total_area_b
        
        #TOP BARS CENTROID
        # Calculate the vertical positions of the bar layers
        y1_t = self._d_b1_t / 2
        y2_t = self._d_b2_t / 2
        y3_t = max(self._d_b1_t, self._d_b2_t) + self.layers_spacing + self._d_b3_t / 2
        y4_t = max(self._d_b1_t, self._d_b2_t) + self.layers_spacing + self._d_b4_t / 2

        # Calculate the total area of each layer
        area_1_t = self._n1_t * self._d_b1_t**2*np.pi/4  # Area proportional to number of bars and their diameter
        area_2_t = self._n2_t * self._d_b2_t**2*np.pi/4
        area_3_t = self._n3_t * self._d_b3_t**2*np.pi/4
        area_4_t = self._n4_t * self._d_b4_t**2*np.pi/4

        # Calculate the centroid as a weighted average
        total_area_t = area_1_t + area_2_t + area_3_t + area_4_t
        if total_area_t == 0:
            return 0*mm  # Avoid division by zero if no bars are present

        self._top_rebar_centroid = (area_1_t * y1_t + area_2_t * y2_t + area_3_t * y3_t + area_4_t * y4_t) / total_area_t
        
    def _calculate_phi_ACI_318_19(self, epsilon_most_strained: float) -> float:
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

    def _maximum_flexural_reinforcement_ratio_ACI_318_19(self) -> float:
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

        """
        # Cast the concrete object to the specific ACI subclass
        concrete_aci = cast(Concrete_ACI_318_19, self.concrete)
        
        # Calculate minimum steel strain for ductility (tension controled according 9.3.3.1 and 21.2.2)
        epsilon_min_rebar_ACI_318_19 = (
            self.steel_bar.epsilon_y + concrete_aci._epsilon_c
        )
        
        # Calculate maximum reinforcement ratio (ρ_max)
        rho_max = (
            0.85 * concrete_aci.beta_1 * self.concrete.f_c / self.steel_bar.f_y
            * (concrete_aci._epsilon_c /
            (concrete_aci._epsilon_c + epsilon_min_rebar_ACI_318_19))
        )

        return rho_max

    def _minimum_flexural_reinforcement_ratio_ACI_318_19(self, M_u:PlainQuantity) -> PlainQuantity:
        """
        Calculates the minimum flexural reinforcement ratio according to ACI 318-19 
        provisions based on the factored moment, M_u.

        This method determines the minimum amount of tensile reinforcement 
        (in terms of a reinforcement ratio) that should be provided in a 
        reinforced concrete section according to ACI 318-19. If the factored 
        moment M_u is zero, it means there is no flexural demand, and hence 
        no minimum flexural reinforcement is required (the ratio is zero). 
        If M_u is not zero, the method checks the unit system (metric or imperial) 
        and computes the required minimum ratio accordingly. These calculations 
        depend on the compressive strength of the concrete (f_c) and the yield 
        strength of the reinforcing steel (f_y).

        Parameters
        ----------
        M_u : PlainQuantity
            The factored moment for the section where the minimum flexural 
            reinforcement ratio is required. The unit should be consistent with 
            the chosen system (e.g., kNm in metric).

        Returns
        -------
        PlainQuantity
            The minimum flexural reinforcement ratio (dimensionless).

        Notes
        -----
        - If M_u = 0, it indicates no moment demand, thus no minimum flexural 
        reinforcement is required (resulting in a zero ratio).
        - For M_u > 0, the minimum ratio is determined using formulas involving 
        the square root of f_c and the value of f_y, in accordance with 
        ACI 318-19.
        - The result is a dimensionless ratio representing the minimum area 
        of steel to the area of the concrete section.

        References
        ----------
        ACI Committee 318. "Building Code Requirements for Structural Concrete 
        (ACI 318-19) and Commentary", American Concrete Institute, 2019.
        """
        if M_u==0*kNm:
            minimum_ratio = 0*dimensionless
        else:
            if self.concrete.unit_system == "metric":
                minimum_ratio = max((0.25 * np.sqrt(self.concrete.f_c / MPa) * MPa / self.steel_bar.f_y), 
                    (1.4 * MPa / self.steel_bar.f_y))
            else:
                minimum_ratio = max((3 * np.sqrt(self.concrete.f_c / psi) * psi / self.steel_bar.f_y), 
                    (200 * psi / self.steel_bar.f_y))
        return minimum_ratio

    def _calculate_flexural_reinforcement_ACI_318_19(self, M_u: PlainQuantity, d: float, d_prima: float) -> tuple[PlainQuantity, PlainQuantity, PlainQuantity, PlainQuantity, float]:
        """
        Calculates the flexural reinforcement for a given factored moment according to ACI 318-19.

        This function computes the required reinforcement areas (minimum, maximum, and final) and
        the compression reinforcement (if required) for a given factored moment. The moment M_u must 
        always be provided as a positive value. For a positive moment, pass 'd' as the effective depth 
        of the tensile reinforcement and 'd_prima' as the effective depth of the compression reinforcement.
        For a negative moment, reverse the roles of 'd' and 'd_prima'.

        Parameters:
            M_u (PlainQuantity): The factored moment (always a positive value).
            d (float): Effective depth of the tensile reinforcement.
            d_prima (float): Effective depth of the compression reinforcement.

        Returns:
            tuple: A tuple containing:
                - A_s_min (PlainQuantity): Minimum reinforcement area required by the code.
                - A_s_max (PlainQuantity): Maximum reinforcement area allowed by the code.
                - A_s_final (PlainQuantity): Final reinforcement area adopted for the tensile zone.
                - A_s_comp (PlainQuantity): Compression reinforcement area (if required).
                - c_d (float): Ratio of the calculated neutral axis depth to the effective depth (c/d).
        """

        # Extract relevant properties and settings
        setting_flexural_min_reduction = self.settings.get_setting('flexural_min_reduction')
        beta_1 = self.concrete.get_properties()["beta_1"]
        b = self._width

        self.settings.load_ACI_318_19_settings()
        self.phi_t = self.settings.get_setting('phi_t')

        # Determine minimum and maximum reinforcement areas
        rho_min = self._minimum_flexural_reinforcement_ratio_ACI_318_19(M_u)
        A_s_min = rho_min * d * b
        rho_max = self._maximum_flexural_reinforcement_ratio_ACI_318_19()
        A_s_max = rho_max * d * b

        # Calculate required reinforcement based on the nominal moment capacity
        R_n = M_u / (self.phi_t * b * d**2)
        # TODO: REVIEW WHAT TO DO WHEN THE VALUE INSIDE THE SQUARE ROOT IS NEGATIVE.
        # Verify if the value under the square root is negative
        sqrt_value = 1 - 2 * R_n / (0.85 * self.concrete.f_c)
        if sqrt_value < 0:
            # Raise exception if the square root is negative;
            # Here we assign A_s_max so that the calculation does not break,
            # resulting in a DCR greater than 1.
            A_s_calc = A_s_max
        else:
            A_s_calc = 0.85 * self.concrete.f_c * b * d / self.steel_bar.f_y * (1 - np.sqrt(sqrt_value))
        
        # Calculate the neutral axis depth based on equilibrium: 0.85 * f_c * c * beta_1 * b = A_s * f_y
        c = A_s_calc * self.steel_bar.f_y / (0.85 * self.concrete.f_c * b * beta_1)

        # Helper function to clean near-zero values
        def clean_zero(value: float, tolerance: float = 1e-6) -> float:
            return 0.0 if abs(value) < tolerance else value

        c_d = clean_zero(c / d)

        # Adjust the required reinforcement to meet the limits
        if A_s_calc > A_s_min:
            A_s_final = A_s_calc
        elif 4 * A_s_calc / 3 > A_s_min:
            A_s_final = A_s_min
        else:
            A_s_final = clean_zero(4 * A_s_calc.to('cm**2').magnitude / 3) * cm**2 if setting_flexural_min_reduction == 'True' else A_s_min

        # Determine if compression reinforcement is required
        if A_s_final <= A_s_max:
            A_s_comp = 0 * cm**2
        else:
            rho = 0.85 * beta_1 * self.concrete.f_c.to(psi).magnitude / self.steel_bar.f_y.to(psi).magnitude * (0.003 / (self.steel_bar.epsilon_y + 0.006))
            M_n_t = rho * self.steel_bar.f_y * (d - 0.59 * rho * self.steel_bar.f_y * d / self.concrete.f_c) * b * d
            M_n_prima = M_u / self.phi_t - M_n_t
            c_t = 0.003 * d / (self.steel_bar.epsilon_y + 0.006)
            c_d = clean_zero(c_t / d)
            f_s_prima = min(0.003 * self.steel_bar.E_s * (1 - d_prima / c_t), self.steel_bar.f_y)
            A_s_comp = M_n_prima / (f_s_prima * (d - d_prima))
            A_s_final = rho * b * d + A_s_comp

        if sqrt_value < 0:
            # In the case where the quadratic equation fails (negative square root),
            # use an approximate value to provide an estimate.
            A_s_final = M_u / (self.phi_t * 0.9 * d * self.steel_bar.f_y)

        return A_s_min, A_s_max, A_s_final, A_s_comp, c_d

    def _determine_nominal_moment_simple_reinf_ACI_318_19(self, A_s: PlainQuantity, d: PlainQuantity) -> PlainQuantity:
        """
        Determines the nominal moment for a simply reinforced section according to ACI 318-19.

        This formula is used ONLY when the provided reinforcement area (A_s) is less than or equal to A_s_max.

        The equilibrium of forces is assumed (compression equals tension):
            0.85 * f_c * a * b = A_s * f_y
        which implies:
            a = (A_s * f_y) / (0.85 * f_c * b)

        Parameters:
            A_s (PlainQuantity): The area of reinforcement.
            d (PlainQuantity): The effective depth of the section.

        Returns:
            PlainQuantity: The nominal moment (M_n) calculated as A_s * f_y * (d - a/2).
        """
        # Calculate the depth of the equivalent rectangular stress block (a)
        a = A_s * self.steel_bar.f_y / (0.85 * self.concrete.f_c * self._width)
        # Calculate the nominal moment (M_n)
        M_n = A_s * self.steel_bar.f_y * (d - a / 2)
        return M_n
    
    def _determine_nominal_moment_double_reinf_ACI_318_19(self, A_s: PlainQuantity, d: PlainQuantity, d_prime: PlainQuantity, A_s_prime: PlainQuantity) -> PlainQuantity:
        """
        Determines the nominal moment for a doubly reinforced beam section according to ACI 318-19.

        This method is used only when the beam has reinforcement exceeding the maximum limit 
        and includes compression reinforcement.

        Equilibrium is assumed (Compression = Tension):
            C_total = Concrete compression (Cc) + Compression reinforcement (Cs)
            Tension (T) = A_s * f_y

        Initially, it is assumed that the compression reinforcement yields (i.e., f_s_prime = f_y),
        so the equilibrium equation becomes:
            0.85 * f_c * a * b + A_s_prime * f_y = A_s * f_y

        Parameters:
            A_s (PlainQuantity): Area of the tensile reinforcement.
            d (PlainQuantity): Effective depth of the beam section.
            d_prime (PlainQuantity): Effective depth (cover) to the compression reinforcement.
            A_s_prime (PlainQuantity): Area of the compression reinforcement.

        Returns:
            PlainQuantity: The nominal moment (M_n) of the doubly reinforced section.
        """
        f_c = self.concrete.f_c
        if isinstance(self.concrete, Concrete_ACI_318_19):
            beta_1 = self.concrete._beta_1
            epsilon_c = self.concrete._epsilon_c

        f_y = self.steel_bar.f_y
        E_s = self.steel_bar._E_s
        epsilon_y = self.steel_bar._epsilon_y
        b = self._width

        # -------------------------------------------------------------------------
        # Step 1: Assume that the compression steel is yielding (f_s_prime = f_y)
        # Based on equilibrium:
        #   0.85 * f_c * a * b + A_s_prime * f_y = A_s * f_y
        # Solving for the neutral axis depth 'c_assumed':
        c_assumed = (A_s * f_y - A_s_prime * f_y) / (0.85 * f_c * b * beta_1)
        
        # Compute the strain in the compression reinforcement with the assumed c value:
        epsilon_s = (c_assumed - d_prime) / c_assumed * epsilon_c if c_assumed > 0 else 0
        # 0 is an arbitrary value, for example, when As top and Bottom are equal, c_assumed is 0, 
        # and compression steel is not yielding, so epsilon_s must be calculated.

        # -------------------------------------------------------------------------
        # Step 2: Check if the assumed compression steel strain exceeds the yield strain.
        if epsilon_s >= epsilon_y:
            # The assumption is valid (compression reinforcement yields).
            a_assumed = c_assumed * beta_1
            M_n = 0.85 * f_c * a_assumed * b * (d - a_assumed / 2) + A_s_prime * f_y * (d - d_prime)
            return M_n
        else:
            # The assumption is invalid, so determine the actual neutral axis depth 'c' using the quadratic equation:
            # Based on equilibrium:
            #   A_s * f_y = 0.85 * f_c * b * beta_1 * c + A_s_prime * ( (c - d_prime) / c * epsilon_c * E_s )
            # Rearranging into a quadratic form: A*c^2 + B*c + C = 0, where:

            A = 0.85 * f_c * b * beta_1
            B = A_s_prime * epsilon_c * E_s - A_s * f_y
            C = -d_prime * A_s_prime * epsilon_c * E_s
            
            # Solve for c using the quadratic formula:
            c = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)
            
            # Compute the corresponding depth of the equivalent rectangular stress block:
            a = c * beta_1
            
            # Recompute the strain and stress in the compression reinforcement:
            epsilon_s_prime = (c - d_prime) / c * epsilon_c
            f_s_prime = epsilon_s_prime * E_s

            
            # Determine the adjusted areas for tension reinforcement:
            # A portion of the tensile reinforcement is balanced by the compression reinforcement.
            A_s_2 = A_s_prime * f_s_prime / f_y
            A_s_1 = A_s - A_s_2
            
            # Calculate the nominal moment contributions:
            M_n_1 = A_s_1 * f_y * (d - a / 2)
            M_n_2 = A_s_prime * f_s_prime * (d - d_prime)
            M_n = M_n_1 + M_n_2
            
            return M_n

    def _determine_nominal_moment_ACI_318_19(self, force: Forces) -> None:
        """
        Determines the nominal moment for a given section with both top and bottom reinforcement,
        calculating the nominal moment for both positive and negative moment scenarios.

        For positive moments, the tension is in the bottom reinforcement.
        For negative moments, the tension is in the top reinforcement.

        Parameters:
            force (Forces): An object containing the forces acting on the section, including the moment M_y.

        Returns:
            None
        """

        # Load design settings for ACI 318-19
        self.settings.load_ACI_318_19_settings()
        self.phi_t = self.settings.get_setting('phi_t')

        # Calculate minimum and maximum reinforcement ratios
        rho_min = self._minimum_flexural_reinforcement_ratio_ACI_318_19(force._M_y)
        rho_max = self._maximum_flexural_reinforcement_ratio_ACI_318_19()
        
        # For positive moments (tension in the bottom), set minimum reinforcement accordingly.
        if force._M_y > 0 * kNm:
            rho_min_top = 0 * dimensionless
            rho_min_bot = rho_min
        else:
            rho_min_top = rho_min
            rho_min_bot = 0 * dimensionless

        # Calculate minimum and maximum bottom reinforcement areas
        self._A_s_min_bot = rho_min_bot * self._d_bot * self._width
        self._A_s_max_bot = rho_max * self._d_bot * self._width

        # Determine the nominal moment for positive moments
        if (self._A_s_bot <= self._A_s_max_bot):
            M_n_positive = self._determine_nominal_moment_simple_reinf_ACI_318_19(self._A_s_bot, self._d_bot)
        elif (self._A_s_top == 0 * cm**2):
            M_n_positive = self._determine_nominal_moment_simple_reinf_ACI_318_19(self._A_s_max_bot, self._d_bot)
        else:
            M_n_positive = self._determine_nominal_moment_double_reinf_ACI_318_19(
                self._A_s_bot, self._d_bot, self._c_mec_top, self._A_s_top
            )

        # Determine capacity for negative moment (tension at the top)
        self._A_s_min_top = rho_min_top * self._d_top * self._width
        self._A_s_max_top = rho_max * self._d_top * self._width
        
        if (self._A_s_top == 0 * cm**2):
            M_n_negative = 0 * kNm
        elif (self._A_s_top <= self._A_s_max_top):
            M_n_negative = self._determine_nominal_moment_simple_reinf_ACI_318_19(self._A_s_top, self._d_top)
        elif (self._A_s_bot == 0 * cm**2):
            M_n_negative = self._determine_nominal_moment_simple_reinf_ACI_318_19(self._A_s_max_top, self._d_top)
        else:
            M_n_negative = self._determine_nominal_moment_double_reinf_ACI_318_19(
                self._A_s_top, self._d_top, self._c_mec_bot, self._A_s_bot
            )

        # Calculate the design moment capacities for both bottom and top reinforcement
        self._phi_M_n_bot: PlainQuantity = self.phi_t * M_n_positive
        self._phi_M_n_top: PlainQuantity = self.phi_t * M_n_negative

        return None

    def _check_flexure_ACI_318_19(self, force: Forces) -> pd.DataFrame:
        """
        Checks the flexural capacity of the section according to ACI 318-19 guidelines.

        This function accepts a single force and performs the flexural check of the section 
        following the ACI 318-19 requirements. It initializes the design variables, computes the 
        nominal moments for both top and bottom reinforcement, determines the required reinforcement 
        areas, and calculates the design capacity ratios. Finally, the results are compiled into a 
        Pandas DataFrame.

        Parameters:
            force (Forces): The force acting on the section, which must include a single moment value.

        Returns:
            pd.DataFrame: A DataFrame containing the flexural design metrics and results.
        """

        # Initialize the design variables according to ACI 318-19 requirements using the provided force.
        self._initialize_variables_ACI_318_19(force)

        # Calculate the nominal moments for both top and bottom reinforcement.
        self._determine_nominal_moment_ACI_318_19(force)

        if self._M_u >= 0:
            # For positive moments, calculate the reinforcement requirements for the bottom tension side.
            (self._A_s_min_bot, self._A_s_max_bot, self._A_s_req_bot, self._A_s_req_top, 
            self._c_d_bot) = self._calculate_flexural_reinforcement_ACI_318_19(self._M_u_bot, self._d_bot, self._c_mec_top)
            self._c_d_top = 0
            # Calculate the design capacity ratio for the bottom side.
            self._DCRb_bot = round(self._M_u_bot.to('kN*m').magnitude / self._phi_M_n_bot.to('kN*m').magnitude, 3)
            self._DCRb_top = 0
        else:
            # For negative moments, calculate the reinforcement requirements for the top tension side.
            (self._A_s_min_top, self._A_s_max_top, self._A_s_req_top, self._A_s_req_bot, 
            self._c_d_top) = self._calculate_flexural_reinforcement_ACI_318_19(abs(self._M_u_top), self._d_top, self._c_mec_bot)
            self._c_d_bot = 0
            # Calculate the design capacity ratio for the top side.
            self._DCRb_top = round(-self._M_u_top.to('kN*m').magnitude / self._phi_M_n_top.to('kN*m').magnitude, 3)
            self._DCRb_bot = 0

        # Determine the maximum detailing cover dimensions for top and bottom.
        self._d_b_max_top = max(self._d_b1_t, self._d_b2_t, self._d_b3_t, self._d_b4_t)
        self._d_b_max_bot = max(self._d_b1_b, self._d_b2_b, self._d_b3_b, self._d_b4_b)
        
        # Calculate the longitudinal reinforcement ratios for both sides.
        self._rho_l_bot = self._A_s_bot / (self._d_bot * self._width)
        self._rho_l_top = self._A_s_bot / (self._d_top * self._width)
        
        # Compile the design results into a dictionary.
        results = self._compile_results_ACI_flexure_metric(force)
        
        # Initialize any additional dictionaries required for ACI 318-19 flexural checks.
        self._initialize_dicts_ACI_318_19_flexure()
        
        # Return the results as a Pandas DataFrame.
        return pd.DataFrame([results], index=[0])

    def _design_flexure_ACI_318_19(self, max_M_y_bot: PlainQuantity, max_M_y_top: PlainQuantity) -> Dict[str, Any]:
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

        # Load design settings for ACI 318-19
        self.settings.load_ACI_318_19_settings()
        self.phi_t = self.settings.get_setting('phi_t')
        # Initial assumptions for mechanical cover and compression depth
        rec_mec = self.c_c + self._stirrup_d_b + 1*cm
        d_prima = self.c_c + self._stirrup_d_b + 1*cm
        # Start the iterative process
        tol = 0.01 * cm  # Tolerance for convergence
        Err = 2 * tol
        iteration_count = 0

        while Err >= tol:
            iteration_count += 1
            # Update the effective depth for bottom tension reinforcement
            d = self._height - rec_mec

            # Calculate reinforcement for the positive moment, even if it is 0
            (
                self._A_s_min_bot,
                self._A_s_max_bot,
                A_s_final_bot_Positive_M,
                A_s_comp_top,
                self._c_d_bot
            ) = self._calculate_flexural_reinforcement_ACI_318_19(max_M_y_bot, d, d_prima)
            # Initialize bottom and top reinforcement
            self._A_s_bot = A_s_final_bot_Positive_M
            self._A_s_top = A_s_comp_top
            # If there is a negative moment, calculate the top reinforcement
            if max_M_y_top < 0:
                (
                    self._A_s_min_top,
                    self._A_s_max_top,
                    A_s_final_top_Negative_M,
                    A_s_comp_bot,
                    self._c_d_top
                ) = self._calculate_flexural_reinforcement_ACI_318_19(abs(max_M_y_top), self.height-d_prima, rec_mec)
                
                # Adjust reinforcement areas based on positive and negative moments
                self._A_s_bot = max(A_s_final_bot_Positive_M, A_s_comp_bot)                 
                self._A_s_top = max(A_s_comp_top,A_s_final_top_Negative_M)
                
            # Design bottom reinforcement
            section_rebar_bot = Rebar(self)
            self.flexure_design_results_bot = section_rebar_bot.longitudinal_rebar_ACI_318_19(self._A_s_bot)
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

            # Set rebar information to section
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

                # Set rebar information to section
                self.set_longitudinal_rebar_bot(n_1_bot, d_b1_bot, n_2_bot, d_b2_bot, 
                                n_3_bot, d_b3_bot, n_4_bot, d_b4_bot)
                self.set_longitudinal_rebar_top(n_1_top, d_b1_top, n_2_top, d_b2_top, 
                                    n_3_top, d_b3_top, n_4_top, d_b4_top)

                d_prima_calculo = self.c_c + self._stirrup_d_b + self._top_rebar_centroid
            else:
                # If no top reinforcement is required
                d_prima_calculo = d_prima
                # Set rebar information to section
                self.set_longitudinal_rebar_bot(n_1_bot, d_b1_bot, n_2_bot, d_b2_bot, 
                                n_3_bot, d_b3_bot, n_4_bot, d_b4_bot)
                self.set_longitudinal_rebar_top(0, 0*mm, 0, 0*mm, 0, 0*mm, 0, 0*mm)

            # Update error for iteration
            Err = max(abs(c_mec_calc - rec_mec), abs(d_prima_calculo - d_prima))
            rec_mec = c_mec_calc
            d_prima = d_prima_calculo

        # Return results as a DataFrame
        results = {
            "Bottom_As_adopted": self._A_s_bot.to("inch**2"),
            "Bottom separation of bars": self._available_s_bot.to("inch"),
            "As_compression_adopted": self._A_s_top.to("inch**2"),
            "Top separation of bars": self._available_s_top.to("inch"),
        }
        return pd.DataFrame([results], index=[0])

    def _design_flexure_EN_1992(self, force: Forces) -> None:
        pass

    def design_flexure(self, forces: list[Forces]) -> pd.DataFrame:
        """
        Designs flexural reinforcement for the beam using the provided forces and design code.
        Identifies the limiting cases for top and bottom reinforcement, designs for those cases, 
        and then checks flexural capacity for all forces.
        
        Returns
        -------
        DataFrame
            A DataFrame summarizing the flexural design results for all forces.
        """

        # Initialize limiting cases
        max_M_y_top = 0 * kN * m  # For negative M_y (top reinforcement design)
        max_M_y_bot = 0 * kN * m  # For positive M_y (bottom reinforcement design)

        # Identify the limiting cases
        for force in forces:
            # For top reinforcement, consider the minimum (most negative) moment
            if force._M_y <= max_M_y_top:
                max_M_y_top = force._M_y
                self._limiting_case_bot=force

            # For bottom reinforcement, consider the maximum positive moment
            if force._M_y > max_M_y_bot:
                max_M_y_bot = force._M_y
                self._limiting_case_top=force
        # Design flexural reinforcement for the limiting cases
        if self.concrete.design_code=="ACI 318-19" or self.concrete.design_code=="CIRSOC 201-25":
            self._design_flexure_ACI_318_19(max_M_y_bot, max_M_y_top)
        else:
            raise ValueError(f"Longitudinal design method not implemented "
                            f"for concrete type: {type(self.concrete).__name__}")
        
        # Check flexural capacity for all forces with the assigned reinforcement
        all_results = self._check_flexure(forces)
        return all_results

    def check_flexure(self, forces: list[Forces]) -> DataFrame:
        
        # Initialize variables to track limiting cases
        max_dcr_top: float = 0
        max_dcr_bot: float = 0
        limiting_case_top = None
        limiting_case_bot = None
        limiting_case_top_details = None
        limiting_case_bot_details = None

        # To compile results for all forces
        self._flexure_results_list = []  # Store individual results for each force
        self._flexure_results_detailed_list = {}  # Store detailed results by force ID

        for force in forces:
            # Select the method based on design code
            if self.concrete.design_code=="ACI 318-19" or self.concrete.design_code=="CIRSOC 201-25":
                result = self._check_flexure_ACI_318_19(force)
            # elif self.concrete.design_code=="EN 1992-2004":
            #     result =  self._check_shear_EN_1992_2004(force)
            else:
                raise ValueError(f"Flexure design method not implemented for concrete type: {type(self.concrete).__name__}")  # noqa: E501
            self._flexure_results_list.append(result)

            # Store detailed results for this force
            self._flexure_results_detailed_list[force.id] = {
                'forces': self._forces_flexure.copy(),
                'min_max': self._data_min_max_flexure.copy(),
                'flexure_capacity_top': self._flexure_capacity_top.copy(),
                'flexure_capacity_bot': self._flexure_capacity_bot.copy(),
                'checks_pass': self._all_flexure_checks_passed
            }

            # Extract the DCR values for top and bottom from the results
            current_dcr_top = self._DCRb_top
            current_dcr_bot = self._DCRb_bot

            # Update top limiting case
            if current_dcr_top > max_dcr_top:
                max_dcr_top = current_dcr_top
                limiting_case_top = result
                limiting_case_top_details = self._flexure_results_detailed_list[force.id]

            # Update bottom limiting case
            if current_dcr_bot > max_dcr_bot:
                max_dcr_bot = current_dcr_bot
                limiting_case_bot = result
                limiting_case_bot_details = self._flexure_results_detailed_list[force.id]

        # Compile all results into a single DataFrame
        all_results = pd.concat(self._flexure_results_list, ignore_index=True)

        # Store limiting cases
        self._limiting_case_flexure_top = limiting_case_top
        self._limiting_case_flexure_bot = limiting_case_bot
        self._limiting_case_flexure_top_details = limiting_case_top_details
        self._limiting_case_flexure_bot_details = limiting_case_bot_details

        # Store maximum DCRs for easy access
        self._max_dcr_top = max_dcr_top
        self._max_dcr_bot = max_dcr_bot
        
        # Mark shear as checked
        self._flexure_checked = True
        return all_results
    
    def _compile_results_ACI_flexure_metric(self, force: Forces) -> Dict[str, Any]:
        # Create dictionaries for bottom and top rows
        if self._M_u>=0:
            result = {
                'Section Label': self.label,
                'Load Combo': force.label,
                'Position': 'Bottom',
                'As,min': self._A_s_min_bot.to('cm ** 2'),
                'As,req top': self._A_s_req_top.to('cm ** 2'),
                'As,req bot': self._A_s_req_bot.to('cm ** 2'),
                'As': self._A_s_bot.to('cm ** 2'),
                # 'c/d': self._c_d_bot,
                'Mu': self._M_u_bot.to('kN*m'),
                'ØMn': self._phi_M_n_bot.to('kN*m'),
                'Mu<ØMn': self._M_u_bot <= self._phi_M_n_bot,
                'DCR': self._DCRb_bot
            }
        else:
            result = {
                'Section Label': self.label,
                'Load Combo': force.label,
                'Position': 'Top',
                'As,min': self._A_s_min_top.to('cm ** 2'),
                'As,req top': self._A_s_req_top.to('cm ** 2'),
                'As,req bot': self._A_s_req_bot.to('cm ** 2'),
                'As': self._A_s_top.to('cm ** 2'),
                # 'c/d': self._c_d_top,
                'Mu': self._M_u_top.to('kN*m'),
                'ØMn': self._phi_M_n_top.to('kN*m'),
                'Mu<ØMn': -self._M_u_top <= self._phi_M_n_top,
                'DCR': self._DCRb_top
            }
        return result
    
    # Factory method to select the shear design method
    def design_shear(self, forces: list[Forces]) -> DataFrame:
        # Track the maximum A_v_req to identify the limiting case
        max_A_v_req = 0*cm**2/m 

        # Step 1: Identify the worst-case force
        for force in forces:
            if self.concrete.design_code=="ACI 318-19" or self.concrete.design_code=="CIRSOC 201-25":
                _design_shear_ACI_318_19(self, force)
            # elif self.concrete.design_code=="EN 1992-2004":
                # self._design_shear_EN_1992_2004(force)
            else:
                raise ValueError(f"Shear design method not implemented for concrete type:"\
                    f"{type(self.concrete).__name__}")
            # Check if this result is the limiting case
            current_A_v_req = self._A_v_req
            if current_A_v_req >= max_A_v_req:
                max_A_v_req = current_A_v_req
                max_V_s_req = self._V_s_req

        # Step 2: Perform rebar design for the worst-case force
        section_rebar = Rebar(self)
        self.shear_design_results = section_rebar.transverse_rebar(max_A_v_req, max_V_s_req)
        self._best_rebar_design = section_rebar.transverse_rebar_design
        self._stirrup_s_l = self._best_rebar_design['s_l']
        self._stirrup_s_w = self._best_rebar_design['s_w']
        self._stirrup_s_max_l = self._best_rebar_design['s_max_l']
        self._stirrup_s_max_w = self._best_rebar_design['s_max_w']
        self.set_transverse_rebar(self._best_rebar_design['n_stir'],
                                  self._best_rebar_design['d_b'],
                                  self._best_rebar_design['s_l'])

        # Step 3: Check shear adequacy for all forces using the designed rebar
        all_results = self.check_shear(forces)  
        return all_results
    
    # Factory method to select the shear check method
    def check_shear(self, forces: list[Forces]) -> DataFrame:

        self._shear_results_list = []  # Store individual results for each force
        self._shear_results_detailed_list = {}  # Store detailed results by force ID
        max_dcr = 0  # Track the maximum DCR to identify the limiting case
        self._limiting_case_shear_details = None

        for force in forces:
            # Select the method based on design code
            if self.concrete.design_code=="ACI 318-19" or self.concrete.design_code=="CIRSOC 201-25":
                result = _check_shear_ACI_318_19(self, force)
            elif self.concrete.design_code=="EN 1992-2004":
                result =  _check_shear_EN_1992_2004(self, force)
            else:
                raise ValueError(f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}")  # noqa: E501
            self._shear_results_list.append(result)
            self._shear_results_detailed_list[force.id] = {
                'forces': self._forces_shear.copy(),
                'shear_reinforcement': self._shear_reinforcement.copy(),
                'min_max': self._data_min_max_shear.copy(),
                'checks_pass': self._all_shear_checks_passed,
                'shear_concrete': self._shear_concrete.copy()
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

    def _format_longitudinal_rebar_string(self, n1: int, d_b1: PlainQuantity, n2: int = 0, d_b2: PlainQuantity = 0*mm) -> str:
        """
        Returns a formatted string representing the rebars and their diameters.
        
        Example output:
            "2Ø10" (if only n1 and d_b1 are defined)
            "2Ø10+2Ø8" (if both n1/d_b1 and n2/d_b2 are defined)
        """
        if n1 == 0:
            rebar_string ="-"
        else:
            rebar_string = f"{n1}Ø{int(d_b1.magnitude)}"
            if n2 > 0 and d_b2.magnitude > 0:  # Check if n2 and d_b2 are defined
                rebar_string += f"+{n2}Ø{int(d_b2.magnitude)}"
        return rebar_string


##########################################################
# RESULTS
##########################################################


    def _initialize_dicts_ACI_318_19_flexure(self) -> None:
        """Initialize the dictionaries used in check and design methods."""
        self._materials_flexure = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
            ],
            "Variable": ["","fc", "fy"],
            "Value": [self.label, round(self.concrete.f_c.to('MPa').magnitude,2), 
                      round(self.steel_bar.f_y.to('MPa').magnitude,2)],
            "Unit": ["", "MPa", "MPa"]
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
                       round(self._c_mec_top.to('cm').magnitude,2), round(self._c_mec_bot.to('cm').magnitude,2)],
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
        min_spacing_top: PlainQuantity = max(self.settings.get_setting('clear_spacing'), self.settings.get_setting('vibrator_size'),
                              self._d_b_max_top)
        min_spacing_bot: PlainQuantity = max(self.settings.get_setting('clear_spacing'), self._d_b_max_bot)
        min_values = [self._A_s_min_top, min_spacing_top, self._A_s_min_bot,
                      min_spacing_bot]   # Use None for items without a minimum constraint
        max_values = [self._A_s_max_top, None, 
                      self._A_s_max_bot,None] # Use None for items without a maximum constraint
        current_values = [self._A_s_top, self._available_s_top, self._A_s_bot, 
                          self._available_s_bot]  # Current values to check

        # Generate check marks based on the range conditions
        checks = [
            '✔️' if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else '❌'
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_flexure_checks_passed = all(check == '✔️' for check in checks)
        self._data_min_max_flexure = {
            'Check': ['Min/Max As rebar top', 'Minimum spacing top','Min/Max As rebar bottom', 
                      'Minimum spacing bottom'],
            'Unit': ['cm²', 'mm', 'cm²', 'mm'],
            'Value': [round(self._A_s_top.to('cm**2').magnitude,2), 
            round(self._available_s_top.to('mm').magnitude,2), round(self._A_s_bot.to('cm**2').magnitude,2), 
            round(self._available_s_bot.to('mm').magnitude,2),],
            'Min.': [round(self._A_s_min_top.to('cm**2').magnitude,2), round(min_spacing_top.to('mm').magnitude,2),
                     round(self._A_s_min_bot.to('cm**2').magnitude,2),  round(min_spacing_bot.to('mm').magnitude,2)],
            'Max.': [round(self._A_s_max_top.to('cm**2').magnitude,2), "", 
                     round(self._A_s_max_bot.to('cm**2').magnitude,2), ""],
            'Ok?': checks
        }
        check_DCR_top = '✔️' if self._DCRb_top < 1 else '❌'
        check_DCR_bot = '✔️' if self._DCRb_bot < 1 else '❌'
        self._flexure_capacity_top = {
            "Top reinforcement check": [
                "First layer bars",
                "Second layer bars",
                "Effective height",
                "Depth of equivalent strength block ratio",
                "Minimum rebar reinforcing",
                "Required rebar reinforcing top",
                "Required rebar reinforcing bottom",
                "Defined rebar reinforcing top",
                "Longitudinal reinforcement ratio",
                "Total flexural strength", 
                "Demand Capacity Ratio"
            ],
            "Variable": ["n1+n2", "n3+n4", "d", "c/d", "As,min","As,req top","As,req bot","As", "ρl", "ØMn", "DCR"],
            "Value": [self._format_longitudinal_rebar_string(self._n1_t, self._d_b1_t, self._n2_t, self._d_b2_t),
                    self._format_longitudinal_rebar_string(self._n3_t, self._d_b3_t, self._n4_t, self._d_b4_t),
                    round(self._d_top.to('cm').magnitude,2),
                    self._c_d_top, round(self._A_s_min_top.to('cm**2').magnitude,2),
                    round(self._A_s_req_top.to('cm**2').magnitude,2),
                    round(self._A_s_req_bot.to('cm**2').magnitude,2),
                    round(self._A_s_top.to('cm**2').magnitude,2), round(self._rho_l_top.magnitude,5),
                    round(self._phi_M_n_top.to('kN*m').magnitude,2), round(self._DCRb_top,2)],
            "Unit": ["", "", "cm", "", "cm²","cm²", "cm²","cm²","","kNm", check_DCR_top]
        }
        self._flexure_capacity_bot = {
            "Bottom reinforcement check": [
                "First layer bars",
                "Second layer bars",
                "Effective height",
                "Depth of equivalent strength block ratio",
                "Minimum rebar reinforcing",
                "Required rebar reinforcing bottom",
                "Required rebar reinforcing top",
                "Defined rebar reinforcing bottom",
                "Longitudinal reinforcement ratio",
                "Total flexural strength", 
                "Demand Capacity Ratio"
            ],
            "Variable": ["n1+n2", "n3+n4", "d", "c/d", "As,min","As,req top", "As,req bot","As", "ρl", "ØMn", "DCR"],
            "Value": [self._format_longitudinal_rebar_string(self._n1_b, self._d_b1_b, self._n2_b, self._d_b2_b),
                    self._format_longitudinal_rebar_string(self._n3_b, self._d_b3_b, self._n4_b, self._d_b4_b),
                    round(self._d_bot.to('cm').magnitude,2),
                    self._c_d_bot, round(self._A_s_min_bot.to('cm**2').magnitude,2),
                    round(self._A_s_req_top.to('cm**2').magnitude,2),
                    round(self._A_s_req_bot.to('cm**2').magnitude,2),
                    round(self._A_s_bot.to('cm**2').magnitude,2), round(self._rho_l_bot.magnitude,5),
                    round(self._phi_M_n_bot.to('kN*m').magnitude,2), round(self._DCRb_bot,2)],
            "Unit": ["", "", "cm", "", "cm²","cm²", "cm²","cm²","","kNm", check_DCR_bot]
        }
        self._flexure_all_checks = self._all_flexure_checks_passed and check_DCR_bot and check_DCR_top
    
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
            warnings.warn("Flexural design has not been performed yet. Call _check_flexure() or "
                        "design_flexure() first.", UserWarning)
            self._md_flexure_results = "Flexural results are not available."
            return None
        # Check if limiting case details exist
        top_details = self._limiting_case_flexure_top_details or {}
        bot_details = self._limiting_case_flexure_bot_details or {}
        # Use limiting case results
        top_result_data = top_details.get('flexure_capacity_top')
        bot_result_data = bot_details.get('flexure_capacity_bot')

        checks_pass_top = top_details.get('checks_pass')
        checks_pass_bot = bot_details.get('checks_pass')
        warning_top = "⚠️ Some checks failed, see detailed results." if not checks_pass_top else ""
        warning_bot = "⚠️ Some checks failed, see detailed results." if not checks_pass_bot else ""

        # Pending for approval
        # all_checks_top = top_details.get('flexure_check')
        # all_checks_bot = bot_details.get('flexure_check')
        # all_checks = all_checks_top and all_checks_bot

        # Formatter instance for DCR formatting
        markdown_content = ""
        formatter = Formatter()
        
        # Handle top result data
        if top_result_data:
            top_rebar_1 = top_result_data['Value'][0]
            top_rebar_2 = top_result_data['Value'][1]
            area_top = top_result_data['Value'][7]
            Mu_top = self._limiting_case_flexure_top_details['forces']['Value'][0]
            Mn_top = top_result_data['Value'][9]
            DCR_top = top_result_data['Value'][10]

            rebar_top = f"{top_rebar_1}" + (f" ++ {top_rebar_2}" if top_rebar_2 != '-' else "")
            formatted_DCR_top = formatter.DCR(DCR_top)

            markdown_content += (
                f"Top longitudinal rebar: {rebar_top}, $A_{{s,top}}$ = {area_top} cm², $M_u$ = {Mu_top} kNm, "
                f"$\\phi M_n$ = {Mn_top} kNm → {formatted_DCR_top} {warning_top}\n\n"
            )
        else:
            markdown_content += "No top moment to check.\n\n"

        # Handle bottom result data
        if bot_result_data:
            bot_rebar_1 = bot_result_data['Value'][0]
            bot_rebar_2 = bot_result_data['Value'][1]
            area_bot = bot_result_data['Value'][7]
            Mu_bot = self._limiting_case_flexure_bot_details['forces']['Value'][1]
            Mn_bot = bot_result_data['Value'][9]
            DCR_bot = bot_result_data['Value'][10]

            rebar_bot = f"{bot_rebar_1}" + (f" ++ {bot_rebar_2}" if bot_rebar_2 != '-' else "")
            formatted_DCR_bot = formatter.DCR(DCR_bot)

            markdown_content += (
                f"Bottom longitudinal rebar: {rebar_bot}, $A_{{s,bot}}$ = {area_bot} cm², $M_u$ = {Mu_bot} kNm, "
                f"$\\phi M_n$ = {Mn_bot} kNm → {formatted_DCR_bot} {warning_bot}"
            )
        else:
            markdown_content += "No bottom moment to check."

        # markdown_content += 'Beam flexure checks PASS ✔️' if all_checks else "Beam flexure checks FAIL ❌"

        self._md_flexure_results = markdown_content
        display(Markdown(markdown_content))

    @property
    def shear_results(self) -> None:
        if not self._shear_checked:
            warnings.warn("Shear design has not been performed yet. Call check_shear() or "
                          "design_shear() first.", UserWarning)
            self._md_shear_results = "Shear results are not available."
            return None
        

        # Check if limiting case details exist
        shear_details = self._limiting_case_shear_details or {}
        # Use limiting case results
        limiting_reinforcement = shear_details.get('shear_reinforcement')
        limiting_forces = shear_details.get('forces')
        limiting_shear_concrete = shear_details.get('shear_concrete')
        checks_pass = shear_details.get('checks_pass')
        markdown_content = ""
        if shear_details:
            # Create FUFormatter instance and format FU value
            formatter = Formatter()
            formatted_DCR = formatter.DCR(limiting_shear_concrete['Value'][-1])
            if self._A_v == 0*cm:
                rebar_v = "not assigned"
            else:
                rebar_v = f"{int(limiting_reinforcement['Value'][0])}eØ{limiting_reinforcement['Value'][1]}/"\
                        f"{limiting_reinforcement['Value'][2]} cm"
            # Limitng cases checks 
            warning = "⚠️ Some checks failed, see detailed results." if not checks_pass else "" 
            markdown_content = f"Shear reinforcing {rebar_v}, $A_v$={limiting_reinforcement['Value'][6]} cm²/m"\
                            f", $V_u$={limiting_forces['Value'][1]} kN, $\\phi V_n$={limiting_shear_concrete['Value'][7]} kN → {formatted_DCR} {warning}"  # noqa: E501
        else:
            markdown_content += "No shear to check."
        self._md_shear_results = markdown_content
        display(Markdown(markdown_content))

        return None
    
    # Beam results for Jupyter Notebook
    @property
    def results(self) -> None:
        """
        Ensure that properties, flexure results, and shear results are available and display them.
        Handles cases where flexure or shear results are not yet available.
        """
        if not hasattr(self, '_md_properties'):
            self.data  # This will generate _md_properties
        if self._flexure_checked:
            self.flexure_results  # This will generate _md_flexure_results
        if self._shear_checked:
            self.shear_results  # This will generate _md_shear_results
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
            warnings.warn("Flexural check has not been performed yet. Call _check_flexure or "
                        "design_flexure first.", UserWarning)
            self._md_flexure_results = "Flexure results are not available."
            return None

        # Determine which results to display (limiting case by default)
        if force:
            if force.id not in self._flexure_results_detailed_list:
                raise ValueError(f"No results found for Forces object with ID {force.id}.")
            result_data = self._flexure_results_detailed_list[force.id]
            top_result_data = result_data['flexure_capacity_top']
            bot_result_data = result_data['flexure_capacity_bot']
            forces_result = result_data['forces']
            min_max_result = result_data['min_max']
        else:
            if self._limiting_case_flexure_top_details is None:
                raise ValueError("Top limiting case details are not available.")
            if self._limiting_case_flexure_bot_details is None:
                raise ValueError("Bottom limiting case details are not available.")
            # Use the worst-case top and bottom scenarios
            top_result_data = self._limiting_case_flexure_top_details['flexure_capacity_top']
            bot_result_data = self._limiting_case_flexure_bot_details['flexure_capacity_bot']
            forces_result = {
                "Design forces": [
                    "Top max moment",
                    "Bottom max moment",
                ],
                "Variable": ["Mu,top", "Mu,bot"],
                "Value": [
                    round(self._limiting_case_flexure_top_details['forces']['Value'][0], 2),
                    round(self._limiting_case_flexure_bot_details['forces']['Value'][1], 2),
                ],
                "Unit": ["kNm", "kNm"]
            }
            min_max_result = {
            'Check': [
                'Min/Max As rebar top',
                'Minimum spacing top',
                'Min/Max As rebar bottom',
                'Minimum spacing bottom'
            ],
            'Unit': ['cm²', 'mm', 'cm²', 'mm'],
            'Value': [
                round(self._limiting_case_flexure_top_details['min_max']['Value'][0], 2),  # Top limiting case As
                round(self._limiting_case_flexure_top_details['min_max']['Value'][1], 2),  # Top limiting case spacing
                round(self._limiting_case_flexure_bot_details['min_max']['Value'][2], 2),  # Bottom limiting case As
                round(self._limiting_case_flexure_bot_details['min_max']['Value'][3], 2)   # Bottom limiting case spacing
            ],
            'Min.': [
                round(self._limiting_case_flexure_top_details['min_max']['Min.'][0], 2),  # Top limiting case As_min
                self._limiting_case_flexure_top_details['min_max']['Min.'][1],           # Top limiting case spacing min
                round(self._limiting_case_flexure_bot_details['min_max']['Min.'][2], 2),  # Bottom limiting case As_min
                self._limiting_case_flexure_bot_details['min_max']['Min.'][3]            # Bottom limiting case spacing min
            ],
            'Max.': [
                round(self._limiting_case_flexure_top_details['min_max']['Max.'][0], 2),  # Top limiting case As_max
                "",  # No max constraint for spacing
                round(self._limiting_case_flexure_bot_details['min_max']['Max.'][2], 2),  # Bottom limiting case As_max
                ""   # No max constraint for spacing
            ],
            'Ok?': [
                self._limiting_case_flexure_top_details['min_max']['Ok?'][0],  # Top As check
                self._limiting_case_flexure_top_details['min_max']['Ok?'][1],  # Top spacing check
                self._limiting_case_flexure_bot_details['min_max']['Ok?'][2],  # Bottom As check
                self._limiting_case_flexure_bot_details['min_max']['Ok?'][3]   # Bottom spacing check
            ]
            }

        # Create TablePrinter instances for detailed display
        print('===== BEAM FLEXURE DETAILED RESULTS =====')
        materials_printer = TablePrinter("MATERIALS")
        materials_printer.print_table_data(self._materials_flexure, headers='keys')

        geometry_printer = TablePrinter("GEOMETRY")
        geometry_printer.print_table_data(self._geometry_flexure, headers='keys')

        forces_printer = TablePrinter("FORCES")
        forces_printer.print_table_data(forces_result, headers='keys')

        min_max_printer = TablePrinter("MAX AND MIN LIMIT CHECKS")
        min_max_printer.print_table_data(min_max_result, headers='keys')

        capacity_printer = TablePrinter("FLEXURAL CAPACITY - TOP")
        capacity_printer.print_table_data(top_result_data, headers='keys')
        capacity_printer = TablePrinter("FLEXURAL CAPACITY - BOTTOM")
        capacity_printer.print_table_data(bot_result_data, headers='keys')

    def flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        """
        Prints detailed flexure results in Word.

        Parameters
        ----------
        forces : Forces, optional
            The specific Forces object to display results for. If None, displays results for the limiting case.
        """
        if not self._flexure_checked:
            warnings.warn("Flexural check has not been performed yet. Call _check_flexure or "
                        "design_flexure first.", UserWarning)
            self._md_flexure_results = "Flexural results are not available."
            return None

        # Determine which results to display (limiting case by default)
        if force:
            if force.id not in self._flexure_results_detailed_list:
                raise ValueError(f"No results found for Forces object with ID {force.id}.")
            result_data = self._flexure_results_detailed_list[force.id]
            top_result_data = result_data['flexure_capacity_top']
            bot_result_data = result_data['flexure_capacity_bot']
            forces_result = result_data['forces']
            min_max_result = result_data['min_max']
        else:
            if self._limiting_case_flexure_top_details is None:
                raise ValueError("Top limiting case details are not available.")
            if self._limiting_case_flexure_bot_details is None:
                raise ValueError("Bottom limiting case details are not available.")
            # Use the worst-case top and bottom scenarios
            top_result_data = self._limiting_case_flexure_top_details['flexure_capacity_top']
            bot_result_data = self._limiting_case_flexure_bot_details['flexure_capacity_bot']
            forces_result = {
                "Design forces": [
                    "Top max moment",
                    "Bottom max moment",
                ],
                "Variable": ["Mu,top", "Mu,bot"],
                "Value": [
                    round(self._limiting_case_flexure_top_details['forces']['Value'][0], 2),
                    round(self._limiting_case_flexure_bot_details['forces']['Value'][1], 2),
                ],
                "Unit": ["kNm", "kNm"]
            }
            min_max_result = {
            'Check': [
                'Min/Max As rebar top',
                'Minimum spacing top',
                'Min/Max As rebar bottom',
                'Minimum spacing bottom'
            ],
            'Unit': ['cm²', 'mm', 'cm²', 'mm'],
            'Value': [
                round(self._limiting_case_flexure_top_details['min_max']['Value'][0], 2),  # Top limiting case As
                round(self._limiting_case_flexure_top_details['min_max']['Value'][1], 2),  # Top limiting case spacing
                round(self._limiting_case_flexure_bot_details['min_max']['Value'][2], 2),  # Bottom limiting case As
                round(self._limiting_case_flexure_bot_details['min_max']['Value'][3], 2)   # Bottom limiting case spacing
            ],
            'Min.': [
                round(self._limiting_case_flexure_top_details['min_max']['Min.'][0], 2),  # Top limiting case As_min
                self._limiting_case_flexure_top_details['min_max']['Min.'][1],           # Top limiting case spacing min
                round(self._limiting_case_flexure_bot_details['min_max']['Min.'][2], 2),  # Bottom limiting case As_min
                self._limiting_case_flexure_bot_details['min_max']['Min.'][3]            # Bottom limiting case spacing min
            ],
            'Max.': [
                round(self._limiting_case_flexure_top_details['min_max']['Max.'][0], 2),  # Top limiting case As_max
                "",  # No max constraint for spacing
                round(self._limiting_case_flexure_bot_details['min_max']['Max.'][2], 2),  # Bottom limiting case As_max
                ""   # No max constraint for spacing
            ],
            'Ok?': [
                self._limiting_case_flexure_top_details['min_max']['Ok?'][0],  # Top As check
                self._limiting_case_flexure_top_details['min_max']['Ok?'][1],  # Top spacing check
                self._limiting_case_flexure_bot_details['min_max']['Ok?'][2],  # Bottom As check
                self._limiting_case_flexure_bot_details['min_max']['Ok?'][3]   # Bottom spacing check
            ]
            }

        # Convert output Dicts into DataFrames
        df_materials = pd.DataFrame(self._materials_flexure)
        df_geometry = pd.DataFrame(self._geometry_flexure)
        df_forces = pd.DataFrame(forces_result)
        df_data_min_max = pd.DataFrame(min_max_result)
        df_flexure_capacity_top = pd.DataFrame(top_result_data)
        df_flexure_capacity_bottom = pd.DataFrame(bot_result_data)
        
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
        doc_builder.add_table_dcr(df_flexure_capacity_top)
        doc_builder.add_heading('Flexural Capacity Bottom', level=2)
        doc_builder.add_table_dcr(df_flexure_capacity_bottom)

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
        print('===== BEAM SHEAR DETAILED RESULTS =====')
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
        doc_builder.add_table_dcr(df_shear_concrete)

        # Save the Word doc
        doc_builder.save(f"Concrete beam shear check {self.concrete.design_code}.docx")


##########################################################################################################


def clear_console() -> None:
    """
    Clears the console based on the operating system.
    """
    if os.name == 'nt':  # For Windows
        os.system('cls')
    else:  # For macOS and Linux
        os.system('clear')



def test_on_determine_nominal_moment_ACI_318_19()->None:
    concrete = Concrete_ACI_318_19(name="fc 4000", f_c=4000*psi)  
    steelBar = SteelBar(name="fy 60000", f_y=60*ksi)  
    custom_settings = {'clear_cover': 1.5*inch} 
    beam = RectangularBeam(
        label="B-12x24",
        concrete=concrete,
        steel_bar=steelBar,
        width=12*inch,  
        height=24*inch,
        settings=custom_settings  
    )
    
    beam.set_longitudinal_rebar_bot(n1=2,d_b1=1.128*inch, n2=1, d_b2=1.128*inch,  n3=2,d_b3=1*inch, n4=1, d_b4=1*inch)
    beam.set_longitudinal_rebar_top(n1=2,d_b1=1.128*inch, n2=1, d_b2=1.128*inch,  n3=2,d_b3=1*inch, n4=1, d_b4=1*inch)

    f = Forces(label='Test_01', M_y=400*kip*ft)
    beam._determine_nominal_moment_ACI_318_19(f)





def flexure_design_test() -> None:
    # clear_console()
    concrete= Concrete_ACI_318_19(name="C25",f_c=25*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa)
    custom_settings = {'clear_cover': 2.5*cm}
    beam = RectangularBeam(label="101",
                                      concrete=concrete,steel_bar=steelBar,width=15*cm, height=50*cm,
                                       settings=custom_settings)
    f1 = Forces(label='C1', M_y=20*kNm)
    f2 = Forces(label='C2', M_y=-20*kNm)
    forces=[f1, f2]
    results=beam._design_flexure(forces)
    # print(beam.flexure_design_results_bot,'\n', beam.flexure_design_results_top)
    print(results)
    # beam.flexure_results_detailed()

def flexure_design_test_calcpad_example() -> None:
    concrete = Concrete_ACI_318_19(name="fc 4000", f_c=4000*psi)  
    steelBar = SteelBar(name="fy 60000", f_y=60*ksi)  
    custom_settings = {'clear_cover': 1.5*inch} 
    beam = RectangularBeam(
        label="B-12x24",
        concrete=concrete,
        steel_bar=steelBar,
        width=12*inch,  
        height=24*inch,
        settings=custom_settings  
    )
    
    #beam.set_longitudinal_rebar_bot(n1=2,d_b1=1.375*inch, n3=2,d_b3=1.27*inch)
    #beam.set_longitudinal_rebar_top(n1=2,d_b1=1.375*inch, n3=2,d_b3=1.27*inch)

    f = Forces(label='Test_01', V_z = 40*kip, M_y=400*kip*ft)
    f2 = Forces(label='Test_01', V_z = 100*kip, M_y=-400*kip*ft)
    forces=[f,f2]

    flexure_results = beam._design_flexure(forces)
    print(flexure_results)
 


def flexure_check_test() -> None:
    clear_console()
    concrete = Concrete_ACI_318_19(name="H-25",f_c=25*MPa) 
    steelBar = SteelBar(name="420", f_y=420*MPa) 

    beam = RectangularBeam(
        label="101",
        concrete=concrete,
        steel_bar=steelBar,
        width=20*cm,  
        height=60*cm,   
    )

    # beam.set_longitudinal_rebar_bot(n1=2,d_b1=20*mm)
    beam.set_longitudinal_rebar_bot(n1=2,d_b1=12*mm, n2=1, d_b2=12*mm, n3=2,d_b3=12*mm, n4=1, d_b4=10*mm)
    beam.set_longitudinal_rebar_top(n1=2,d_b1=16*mm)
    f1 = Forces(label='D', M_y=0*kNm, V_z=50*kN)
    f2 = Forces(label='L', M_y=-100*kNm)
    f3 = Forces(label='W', M_y=-50*kNm)
    f4 = Forces(label='S', M_y=110*kNm)
    forces=[f1,f2,f3,f4]
    print(beam._check_flexure(forces))
    
    # beam.check_shear()
    beam.flexure_results_detailed()
    # beam.flexure_results_detailed_doc()
    # beam.shear_results_detailed_doc()

def flexure_Mn() -> None:
    clear_console()
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
    result=beam._determine_nominal_moment_simple_reinf_ACI_318_19(A_s,d)
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
    d_prime=2.5*inch
    A_s_prime=1.8*inch**2
    result=beam2._determine_nominal_moment_double_reinf_ACI_318_19(A_s, d, d_prime, A_s_prime)
    print(result.to(kip*ft))




def shear_ACI_metric() -> None:
    concrete= Concrete_ACI_318_19(name="C30",f_c=30*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa)
    custom_settings = {'clear_cover': 30*mm}
    beam = RectangularBeam(label="101",
                                      concrete=concrete,steel_bar=steelBar,width=20*cm, height=50*cm,
                                       settings=custom_settings)
    f1 = Forces(label='1.4D', V_z=100*kN)
    f2 = Forces(label='1.2D+1.6L', V_z=1.55*kN)
    f3 = Forces(label='W', V_z=2.20*kN)
    f4 = Forces(label='S', V_z=8.0*kN)
    f5 = Forces(label='E', V_z=1.0*kN)
    forces=[f1, f2, f3, f4, f5]
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)
    # beam.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=20*cm) 
    # results = beam.check_shear()
    results = beam._design_shear(forces)
    print(results)
    print(beam.shear_design_results)
    # beam.shear_results_detailed()
    # print(beam.shear_design_results)
    # print(beam.results)
    # beam.shear_results_detailed_doc()

def shear_ACI_imperial() -> None:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000*psi)  
    steelBar = SteelBar(name="ADN 420", f_y=60*ksi)  
    custom_settings = {'clear_cover': 1.5*inch} 
    beam = RectangularBeam(
        label="102",
        concrete=concrete,
        steel_bar=steelBar,
        width=10*inch,  
        height=16*inch,
        settings=custom_settings  
    )

    # f1 = Forces(label='D', V_z=37.727*kip, N_x=20*kip)
    f1 = Forces(label='D', V_z=37.727*kip)
    # f1 = Forces(label='D', V_z=8*kip)
    # f2 = Forces(label='L', V_z=6*kip) # No shear reinforcing
    forces=[f1]
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.5*inch, s_l=6*inch)

    # beam.set_longitudinal_rebar_bot(n1=2, d_b1=0.625*inch)
    print(beam._A_v)
    results = beam._check_shear(forces)
    # results = beam.design_shear()
    print(results)
    # section.design_shear(f, A_s=0.847*inch**2)
    beam.shear_results_detailed()  
    # section.shear_results_detailed_doc()

def rebar() -> None:
    concrete= Concrete_ACI_318_19(name="H30",f_c=30*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa)
    custom_settings = {'clear_cover': 30*mm, 'stirrup_diameter_ini':8*mm}
    section = RectangularBeam(
        label="V 20x50",
        concrete=concrete,
        steel_bar=steelBar,
        width=20*cm,  
        height=50*cm,
        settings=custom_settings  
    )
    as_req = 7*cm**2

    beam_rebar = Rebar(section)
    long_rebar_df = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=as_req)
    best_design = beam_rebar.longitudinal_rebar_design
    print(long_rebar_df)

def rebar_df() -> None:
    concrete= Concrete_ACI_318_19(name="H30",f_c=30*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa)
    custom_settings = {'clear_cover': 30*mm, 'stirrup_diameter_ini':8*mm}
    section = RectangularBeam(
        label="V 20x50",
        concrete=concrete,
        steel_bar=steelBar,
        width=20*cm,  
        height=50*cm,
        settings=custom_settings  
    )
    beam_rebar = Rebar(section)
    # Create a list of required steel areas from 0.5 to 10 with a step of 0.5 cm²
    as_req_list = np.arange(0.5, 10.5, 0.5)*cm**2
    # Initialize an empty DataFrame to store the results
    results_df = pd.DataFrame()

    # Loop through each required steel area
    for as_req in as_req_list:
        # Run the longitudinal_rebar_ACI_19 method
        long_rebar_df = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=as_req)
        
        # Extract the first row of the resulting DataFrame
        first_row = long_rebar_df.iloc[0:1].copy()
        
        # Add a column to store the required steel area
        first_row['as_req'] = as_req.magnitude  # Store the magnitude (value without units)
        
        # Append the first row to the results DataFrame
        results_df = pd.concat([results_df, first_row], ignore_index=True)

    # Display the results DataFrame
    print(results_df)

def shear_EN_1992() -> None:
    concrete= Concrete_EN_1992_2004(name="C25",f_ck=25*MPa) 
    steelBar= SteelBar(name="B500S", f_y=500*MPa)
    custom_settings = {'clear_cover': 2.6*cm}
    beam = RectangularBeam(label="101",
                                      concrete=concrete,steel_bar=steelBar,width=20*cm, height=60*cm,
                                       settings=custom_settings)
    # f = Forces(V_z=100*kN, M_y=100*kNm)
    f = Forces(V_z=30*kN, N_x=0*kN)
    forces=[f]
    beam.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=25*cm)
    beam.set_longitudinal_rebar_bot(n1=4,d_b1=16*mm)
    beam.check_shear(forces)
    beam.shear_results_detailed()
    # beam.shear_results_detailed_doc()

if __name__ == "__main__":
    shear_EN_1992()

