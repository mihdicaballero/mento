from dataclasses import dataclass, field
from IPython.display import Markdown, display
from typing import Optional, Dict
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch
from pint import Quantity
import numpy as np
import pandas as pd
from pandas import DataFrame
import math
import warnings
from importlib.metadata import version

from mento.rectangular import RectangularSection
from mento.material import (
    Concrete_ACI_318_19,
    Concrete_EN_1992_2004,
)
from mento.rebar import Rebar
from mento.units import MPa, mm, inch, kN, m, cm, kNm, dimensionless
from mento.results import Formatter, TablePrinter, DocumentBuilder, CUSTOM_COLORS
from mento.forces import Forces
from mento.settings import BeamSettings

from mento.codes.EN_1992_2004_beam import (
    _check_shear_EN_1992_2004,
    _design_shear_EN_1992_2004,
)
from mento.codes.ACI_318_19_beam import (
    _check_shear_ACI_318_19,
    _design_shear_ACI_318_19,
    _check_flexure_ACI_318_19,
    _design_flexure_ACI_318_19,
)

MENTO_VERSION = version("mento")


@dataclass
class RectangularBeam(RectangularSection):
    """
    Represents a reinforced concrete rectangular beam section with methods for design, checking,
    and visualization of longitudinal and transverse reinforcement according to various design codes.

    Attributes:
        settings (BeamSettings): Access to global design rules and settings.

    Methods:
        set_transverse_rebar(n_stirrups, d_b, s_l):
            Sets the transverse (stirrup) rebar configuration for the beam.
        set_longitudinal_rebar_bot(n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4):
            Sets the bottom longitudinal rebar configuration.
        set_longitudinal_rebar_top(n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4):
            Sets the top longitudinal rebar configuration.
        design_flexure(forces):
            Designs the flexural reinforcement for the beam based on provided forces and design code.
        check_flexure(forces):
            Checks the flexural capacity for all provided forces and stores results.
        design_shear(forces):
            Designs the shear reinforcement for the beam based on provided forces and design code.
        check_shear(forces):
            Checks the shear capacity for all provided forces and stores results.
        data:
            Property. Displays basic beam data in Markdown format.
        flexure_results:
            Property. Displays summary of flexural design/check results in Markdown format.
        shear_results:
            Property. Displays summary of shear design/check results in Markdown format.
        results:
            Property. Displays all available results (properties, flexure, shear) in Markdown format.
        flexure_results_detailed(force=None):
            Displays detailed flexure results for a specific force or the limiting case.
        flexure_results_detailed_doc(force=None):
            Exports detailed flexure results to a Word document.
        shear_results_detailed(force=None):
            Displays detailed shear results for a specific force or the limiting case.
        shear_results_detailed_doc(force=None):
            Exports detailed shear results to a Word document.
        plot():
            Plots the beam section with its rebar.

    Usage:
        - Instantiate with required section and material properties.
        - Use set_longitudinal_rebar_* and set_transverse_rebar to configure reinforcement.
        - Call design_flexure and design_shear to perform design.
        - Use flexure_results, shear_results, and results properties for summary output.
        - Use flexure_results_detailed and shear_results_detailed for detailed output.
        - Use plot() to visualize the rebar arrangement.

    Note:
        This class is intended for use in reinforced concrete beam design workflows,
        supporting multiple design codes and detailed reporting.
    """

    settings: BeamSettings = field(
        default_factory=BeamSettings
    )  # allow user to pass settings

    def __post_init__(self) -> None:
        super().__post_init__()  # Call parent attributes
        if not self.settings:
            # create defaults based on concrete's unit system
            self.settings = BeamSettings(unit_system=self.concrete.unit_system)
        else:
            # Fill missing fields with defaults from the given unit system
            # (ensures your partial BeamSettings still gets the rest of defaults)
            self.settings.unit_system = self.concrete.unit_system
            self.settings.__post_init__()  # recompute missing fields
        self._initialize_attributes()

    ##########################################################
    # INITIALIZE ATTRIBUTES
    ##########################################################

    def _initialize_attributes(self) -> None:
        """Initialize all attributes of the beam."""
        # Stirrups and shear attributes
        self._stirrup_d_b = self.settings.stirrup_diameter_ini
        self._stirrup_s_l: Quantity = 0 * cm
        self._stirrup_s_w: Quantity = 0 * cm
        self._stirrup_s_max_l: Quantity = 0 * cm
        self._stirrup_s_max_w: Quantity = 0 * cm
        self._stirrup_n: int = 0
        self._A_v_min: Quantity = 0 * cm**2 / m
        self._A_v: Quantity = 0 * cm**2 / m
        self._A_s_req_bot: Quantity = 0 * cm**2
        self._A_s_req_top: Quantity = 0 * cm**2
        self._A_v_req: Quantity = 0 * cm**2 / m
        self._A_s_tension: Quantity = 0 * cm**2
        self._DCRv: float = 0
        self._DCRb_top: float = 0
        self._DCRb_bot: float = 0
        self._alpha: float = math.radians(90)
        self._V_s_req: Quantity = 0 * kN

        # Design checks and effective heights
        self._rho_l_bot: Quantity = 0 * dimensionless
        self._rho_l_top: Quantity = 0 * dimensionless
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

        self._materials_flexure: Dict = {}
        self._geometry_flexure: Dict = {}
        self._forces_flexure: Dict = {}
        self._all_flexure_checks_passed: bool = False
        self._data_min_max_flexure: Dict = {}
        self._flexure_capacity_bot: Dict = {}
        self._flexure_capacity_top: Dict = {}
        self._flexure_all_checks: bool = False

    def _initialize_code_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_ACI_318_19):
            self._initialize_ACI_318_attributes()
        elif isinstance(self.concrete, Concrete_EN_1992_2004):
            self._initialize_EN_1992_2004_attributes()

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

        # Unit system default starter rebar for initial effective height
        if self.concrete.unit_system == "metric":
            self._n1_b, self._d_b1_b = 2, 8 * mm
            self._n1_t, self._d_b1_t = 2, 8 * mm
        else:
            self._n1_b, self._d_b1_b = 2, 3 / 8 * inch
            self._n1_t, self._d_b1_t = 2, 3 / 8 * inch

        # Update dependent attributes
        self._update_longitudinal_rebar_attributes()

    def _initialize_ACI_318_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_ACI_318_19):
            self._phi_V_n: Quantity = 0 * kN
            self._phi_V_s: Quantity = 0 * kN
            self._phi_V_c: Quantity = 0 * kN
            self._phi_V_max: Quantity = 0 * kN
            self._V_u: Quantity = 0 * kN
            self._M_u: Quantity = 0 * kNm
            self._M_u_bot: Quantity = 0 * kNm
            self._M_u_top: Quantity = 0 * kNm
            self._N_u: Quantity = 0 * kN
            self._A_cv: Quantity = 0 * cm**2
            self._k_c_min: Quantity = 0 * MPa
            self._sigma_Nu: Quantity = 0 * MPa
            self.V_c: Quantity = 0 * kN
            self._rho_w: Quantity = 0 * dimensionless
            self._lambda_s: float = 0
            self.f_yt: Quantity = 0 * MPa
            self._max_shear_ok: bool = False
            self._A_s_min_bot: Quantity = 0 * cm**2
            self._A_s_min_top: Quantity = 0 * cm**2
            self._A_s_max_bot: Quantity = 0 * cm**2
            self._A_s_max_top: Quantity = 0 * cm**2
            self._phi_M_n_bot: Quantity = 0 * kNm
            self._phi_M_n_top: Quantity = 0 * kNm
            self._d_b_max_bot: Quantity = 0 * mm
            self._d_b_max_top: Quantity = 0 * mm
            self.flexure_design_results_bot: DataFrame = None
            self.flexure_design_results_top: DataFrame = None

    def _initialize_EN_1992_2004_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_EN_1992_2004):
            # self._f_yk = self.steel_bar.f_y
            # self._f_ck = self.concrete.f_ck
            # self._f_ctm = self.concrete.f_ctm
            # self._epsilon_cu3 = self.concrete._epsilon_cu3
            # self._E_s = self.steel_bar.E_s
            self._V_Ed_1: Quantity = 0 * kN
            self._V_Ed_2: Quantity = 0 * kN
            self._N_Ed: Quantity = 0 * kN
            self._M_Ed: Quantity = 0 * kNm
            self._sigma_cd: Quantity = 0 * MPa
            self._V_Rd_c: Quantity = 0 * kN
            self._V_Rd_s: Quantity = 0 * kN
            self._V_Rd_max: Quantity = 0 * kN
            self._V_Rd: Quantity = 0 * kN
            self._k_value: float = 0
            self._f_ywk: Quantity = 0 * MPa
            self._f_ywd: Quantity = 0 * MPa
            self._f_yd: Quantity = 0 * MPa
            self._f_cd: Quantity = 0 * MPa
            self._A_p = 0 * cm**2  # No prestressed for now
            self._sigma_cp: Quantity = 0 * MPa
            self._theta: float = 0
            self._cot_theta: float = 0
            self._z: Quantity = 0 * cm

    ##########################################################
    # SET LONGITUDINAL AND TRANSVERSE REBAR AND UPDATE ATTRIBUTES
    ##########################################################

    def set_transverse_rebar(
        self,
        n_stirrups: int = 0,
        d_b: Quantity = 0 * mm,
        s_l: Quantity = 0 * cm,
    ) -> None:
        """Sets the transverse rebar in the beam section."""
        self._stirrup_n = n_stirrups
        self._stirrup_d_b = d_b
        self._stirrup_s_l = s_l
        n_legs = n_stirrups * 2
        A_db = (d_b**2) * math.pi / 4  # Area of one stirrup leg
        A_vs = n_legs * A_db  # Total area of stirrups
        self._A_v = A_vs / s_l  # Stirrup area per unit length

        # Update effective heights
        self._update_effective_heights()

    def set_longitudinal_rebar_bot(
        self,
        n1: int = 0,
        d_b1: Quantity = 0 * mm,
        n2: int = 0,
        d_b2: Quantity = 0 * mm,
        n3: int = 0,
        d_b3: Quantity = 0 * mm,
        n4: int = 0,
        d_b4: Quantity = 0 * mm,
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

    def set_longitudinal_rebar_top(
        self,
        n1: int,
        d_b1: Quantity,
        n2: int = 0,
        d_b2: Quantity = 0 * mm,
        n3: int = 0,
        d_b3: Quantity = 0 * mm,
        n4: int = 0,
        d_b4: Quantity = 0 * mm,
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
        self._A_s_bot = (
            self._n1_b * self._d_b1_b**2 * np.pi / 4
            + self._n2_b * self._d_b2_b**2 * np.pi / 4
            + self._n3_b * self._d_b3_b**2 * np.pi / 4
            + self._n4_b * self._d_b4_b**2 * np.pi / 4
        )

        # AREA OF TOP BARS
        self._A_s_top = (
            self._n1_t * self._d_b1_t**2 * np.pi / 4
            + self._n2_t * self._d_b2_t**2 * np.pi / 4
            + self._n3_t * self._d_b3_t**2 * np.pi / 4
            + self._n4_t * self._d_b4_t**2 * np.pi / 4
        )

    def _calculate_min_clear_spacing(self) -> None:
        """
        Calculates the maximum clear spacing between bars for the bottom rebar layers.

        Returns:
            Quantity: The maximum clear spacing between bars in either the first or second layer.
        """

        def layer_clear_spacing(
            n_a: int, d_a: Quantity, n_b: int, d_b: Quantity
        ) -> Quantity:
            """
            Helper function to calculate clear spacing for a given layer.

            Parameters:
                n_a (int): Number of bars in the first group of the layer.
                d_a (Quantity): Diameter of bars in the first group of the layer.
                n_b (int): Number of bars in the second group of the layer.
                d_b (Quantity): Diameter of bars in the second group of the layer.

            Returns:
                Quantity: Clear spacing for the given layer.
            """
            effective_width = self.width - 2 * (self.c_c + self._stirrup_d_b)
            total_bars = n_a + n_b
            if total_bars <= 1:
                return effective_width - max(d_a, d_b)  # Clear space for one bar
            total_bar_width = n_a * d_a + n_b * d_b
            return (effective_width - total_bar_width) / (total_bars - 1)

        # AVAIABLE CLEAR SPACING FOR BOTTOM BARS
        # Calculate clear spacing for each layer
        spacing_layer1_b = layer_clear_spacing(
            self._n1_b, self._d_b1_b, self._n2_b, self._d_b2_b
        )
        spacing_layer2_b = layer_clear_spacing(
            self._n3_b, self._d_b3_b, self._n4_b, self._d_b4_b
        )

        # Return the maximum clear spacing between the two layers
        self._available_s_bot = min(spacing_layer1_b, spacing_layer2_b)

        # AVAIABLE CLEAR SPACING FOR TOP BARS
        # Calculate clear spacing for each layer
        spacing_layer1_t = layer_clear_spacing(
            self._n1_t, self._d_b1_t, self._n2_t, self._d_b2_t
        )
        spacing_layer2_t = layer_clear_spacing(
            self._n3_t, self._d_b3_t, self._n4_t, self._d_b4_t
        )

        # Return the maximum clear spacing between the two layers
        self._available_s_top = min(spacing_layer1_t, spacing_layer2_t)

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
        y3_b = (
            max(self._d_b1_b, self._d_b2_b)
            + self.settings.layers_spacing
            + self._d_b3_b / 2
        )
        y4_b = (
            max(self._d_b1_b, self._d_b2_b)
            + self.settings.layers_spacing
            + self._d_b4_b / 2
        )
        # Calculate the total area of each layer
        area_1_b = (
            self._n1_b * self._d_b1_b**2 * np.pi / 4
        )  # Area proportional to number of bars and their diameter
        area_2_b = self._n2_b * self._d_b2_b**2 * np.pi / 4
        area_3_b = self._n3_b * self._d_b3_b**2 * np.pi / 4
        area_4_b = self._n4_b * self._d_b4_b**2 * np.pi / 4

        # Calculate the centroid as a weighted average
        total_area_b = area_1_b + area_2_b + area_3_b + area_4_b
        if total_area_b == 0:
            return 0 * mm  # Avoid division by zero if no bars are present

        self._bot_rebar_centroid = (
            area_1_b * y1_b + area_2_b * y2_b + area_3_b * y3_b + area_4_b * y4_b
        ) / total_area_b

        # TOP BARS CENTROID
        # Calculate the vertical positions of the bar layers
        y1_t = self._d_b1_t / 2
        y2_t = self._d_b2_t / 2
        y3_t = (
            max(self._d_b1_t, self._d_b2_t)
            + self.settings.layers_spacing
            + self._d_b3_t / 2
        )
        y4_t = (
            max(self._d_b1_t, self._d_b2_t)
            + self.settings.layers_spacing
            + self._d_b4_t / 2
        )

        # Calculate the total area of each layer
        area_1_t = (
            self._n1_t * self._d_b1_t**2 * np.pi / 4
        )  # Area proportional to number of bars and their diameter
        area_2_t = self._n2_t * self._d_b2_t**2 * np.pi / 4
        area_3_t = self._n3_t * self._d_b3_t**2 * np.pi / 4
        area_4_t = self._n4_t * self._d_b4_t**2 * np.pi / 4

        # Calculate the centroid as a weighted average
        total_area_t = area_1_t + area_2_t + area_3_t + area_4_t
        if total_area_t == 0:
            return 0 * mm  # Avoid division by zero if no bars are present

        self._top_rebar_centroid = (
            area_1_t * y1_t + area_2_t * y2_t + area_3_t * y3_t + area_4_t * y4_t
        ) / total_area_t

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
        self._d_bot = self.height - self._c_mec_bot
        self._d_top = self.height - self._c_mec_top
        # Use bottom or top effective height
        self._d_shear = min(self._d_bot, self._d_top)

    ##########################################################
    # CHECK & DESIGN FLEXURE
    ##########################################################

    def design_flexure(self, forces: list[Forces]) -> DataFrame:
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
                self._limiting_case_bot = force
            # For bottom reinforcement, consider the maximum positive moment
            if force._M_y > max_M_y_bot:
                max_M_y_bot = force._M_y
                self._limiting_case_top = force
        # Design flexural reinforcement for the limiting cases
        if (
            self.concrete.design_code == "ACI 318-19"
            or self.concrete.design_code == "CIRSOC 201-25"
        ):
            _design_flexure_ACI_318_19(self, max_M_y_bot, max_M_y_top)
        else:
            raise ValueError(
                f"Longitudinal design method not implemented "
                f"for concrete type: {type(self.concrete).__name__}"
            )

        # Check flexural capacity for all forces with the assigned reinforcement
        all_results = self.check_flexure(forces)
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
            if (
                self.concrete.design_code == "ACI 318-19"
                or self.concrete.design_code == "CIRSOC 201-25"
            ):
                result = _check_flexure_ACI_318_19(self, force)
            # elif self.concrete.design_code=="EN 1992-2004":
            #     result =  _check_shear_EN_1992_2004(self, force)
            else:
                raise ValueError(
                    f"Flexure design method not implemented for concrete type: {type(self.concrete).__name__}"
                )  # noqa: E501
            self._flexure_results_list.append(result)

            # Store detailed results for this force
            self._flexure_results_detailed_list[force.id] = {
                "forces": self._forces_flexure.copy(),
                "min_max": self._data_min_max_flexure.copy(),
                "flexure_capacity_top": self._flexure_capacity_top.copy(),
                "flexure_capacity_bot": self._flexure_capacity_bot.copy(),
                "checks_pass": self._all_flexure_checks_passed,
            }

            # Extract the DCR values for top and bottom from the results
            current_dcr_top = self._DCRb_top
            current_dcr_bot = self._DCRb_bot

            # Update top limiting case
            if current_dcr_top >= max_dcr_top:
                max_dcr_top = current_dcr_top
                limiting_case_top = result
                limiting_case_top_details = self._flexure_results_detailed_list[
                    force.id
                ]

            # Update bottom limiting case
            if current_dcr_bot >= max_dcr_bot:
                max_dcr_bot = current_dcr_bot
                limiting_case_bot = result
                limiting_case_bot_details = self._flexure_results_detailed_list[
                    force.id
                ]

        # Compile all results into a single DataFrame
        all_data = pd.concat(self._flexure_results_list, ignore_index=True)
        units_row = self._get_units_row_flexure()
        all_results = pd.concat([units_row, all_data], ignore_index=True)

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

    ##########################################################
    # CHECK & DESIGN SHEAR
    ##########################################################

    # Factory method to select the shear design method
    def design_shear(self, forces: list[Forces]) -> DataFrame:
        # Track the maximum A_v_req to identify the limiting case
        max_A_v_req = 0 * cm**2 / m

        # Step 1: Identify the worst-case force
        for force in forces:
            if (
                self.concrete.design_code == "ACI 318-19"
                or self.concrete.design_code == "CIRSOC 201-25"
            ):
                _design_shear_ACI_318_19(self, force)
            elif self.concrete.design_code == "EN 1992-2004":
                _design_shear_EN_1992_2004(self, force)
            else:
                raise ValueError(
                    f"Shear design method not implemented for concrete type:"
                    f"{type(self.concrete).__name__}"
                )
            # Check if this result is the limiting case
            current_A_v_req = self._A_v_req
            if current_A_v_req >= max_A_v_req:
                max_A_v_req = current_A_v_req
                max_V_s_req = self._V_s_req

        # Step 2: Perform rebar design for the worst-case force
        section_rebar = Rebar(self)
        self.shear_design_results = section_rebar.transverse_rebar(
            max_A_v_req, max_V_s_req, self._alpha
        )
        self._best_rebar_design = section_rebar.transverse_rebar_design
        self._stirrup_s_l = self._best_rebar_design["s_l"]
        self._stirrup_s_w = self._best_rebar_design["s_w"]
        self._stirrup_s_max_l = self._best_rebar_design["s_max_l"]
        self._stirrup_s_max_w = self._best_rebar_design["s_max_w"]
        self.set_transverse_rebar(
            self._best_rebar_design["n_stir"],
            self._best_rebar_design["d_b"],
            self._best_rebar_design["s_l"],
        )

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
            if (
                self.concrete.design_code == "ACI 318-19"
                or self.concrete.design_code == "CIRSOC 201-25"
            ):
                result = _check_shear_ACI_318_19(self, force)
            elif self.concrete.design_code == "EN 1992-2004":
                result = _check_shear_EN_1992_2004(self, force)
            else:
                raise ValueError(
                    f"Shear design method not implemented for concrete type: {type(self.concrete).__name__}"
                )  # noqa: E501

            self._shear_results_list.append(result)

            self._shear_results_detailed_list[force.id] = {
                "forces": self._forces_shear.copy(),
                "shear_reinforcement": self._shear_reinforcement.copy(),
                "min_max": self._data_min_max_shear.copy(),
                "checks_pass": self._all_shear_checks_passed,
                "shear_concrete": self._shear_concrete.copy(),
            }

            # Check if this result is the limiting case
            current_dcr = result["DCR"][0]
            if current_dcr > max_dcr:
                max_dcr = current_dcr
                self._limiting_case_shear = result
                self._limiting_case_shear_details = self._shear_results_detailed_list[
                    force.id
                ]

        # Compile all results into a single DataFrame
        all_data = pd.concat(self._shear_results_list, ignore_index=True)
        units_row = self._get_units_row_shear()
        all_results = pd.concat([units_row, all_data], ignore_index=True)

        # Identify the most limiting case by Demand-to-Capacity Ratio (DCR) or other criteria
        self.limiting_case_shear = all_data.loc[
            all_data["DCR"].idxmax()
        ]  # Select row with highest DCR

        self._shear_checked = True  # Mark shear as checked
        return all_results

    def _get_units_row_shear(self) -> pd.DataFrame:
        if isinstance(self.concrete, Concrete_EN_1992_2004):
            # Orden exacto de columnas para EN 1992
            return pd.DataFrame(
                [
                    {
                        "Label": "",
                        "Comb.": "",
                        "Av,min": "cm²/m",
                        "Av,req": "cm²/m",
                        "Av": "cm²/m",
                        "VEd,1": "kN",
                        "VEd,2": "kN",
                        "VRd,c": "kN",
                        "VRd,s": "kN",
                        "VRd": "kN",
                        "VRd,max": "kN",
                        "VEd,1≤VRd,max": "",
                        "VEd,2≤VRd": "",
                        "DCR": "",
                    }
                ]
            )
        else:
            return pd.DataFrame(
                [
                    {
                        "Label": "",
                        "Comb.": "",
                        "Av,min": "cm²/m",
                        "Av,req": "cm²/m",
                        "Av": "cm²/m",
                        "Vu": "kN",
                        "Nu": "kN",
                        "ØVc": "kN",
                        "ØVs": "kN",
                        "ØVn": "kN",
                        "ØVmax": "kN",
                        "Vu≤ØVmax": "",
                        "Vu≤ØVn": "",
                        "DCR": "",
                    }
                ]
            )

    def _get_units_row_flexure(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "Label": "",
                    "Comb.": "",
                    "Position": "",
                    "As,min": "cm²",
                    "As,req top": "cm²",
                    "As,req bot": "cm²",
                    "As": "cm²",
                    # "c/d": "",  # Uncomment if you include this field later
                    "Mu": "kNm",
                    "ØMn": "kNm",
                    "Mu≤ØMn": "",
                    "DCR": "",
                }
            ]
        )

    ##########################################################
    # RESULTS
    ##########################################################

    # Beam results for Jupyter Notebook
    @property
    def data(self) -> None:
        markdown_content = (
            f"Beam {self.label}, $b$={self.width.to('cm')}"
            f", $h$={self.height.to('cm')}, $c_{{c}}$={self.c_c.to('cm')}, \
                            Concrete {self.concrete.name}, Rebar {self.steel_bar.name}."
        )
        self._md_data = markdown_content
        # Display the combined content
        display(Markdown(markdown_content))  # type: ignore

        return None

    @property
    def flexure_results(self) -> None:
        if not self._flexure_checked:
            warnings.warn(
                "Flexural design has not been performed yet. Call _check_flexure() or "
                "design_flexure() first.",
                UserWarning,
            )
            self._md_flexure_results = "Flexural results are not available."
            return None
        # Check if limiting case details exist
        top_details = self._limiting_case_flexure_top_details or {}
        bot_details = self._limiting_case_flexure_bot_details or {}
        # Use limiting case results
        top_result_data = top_details.get("flexure_capacity_top")
        bot_result_data = bot_details.get("flexure_capacity_bot")

        checks_pass_top = top_details.get("checks_pass")
        checks_pass_bot = bot_details.get("checks_pass")
        warning_top = (
            "⚠️ Some checks failed, see detailed results."
            if not checks_pass_top
            else ""
        )
        warning_bot = (
            "⚠️ Some checks failed, see detailed results."
            if not checks_pass_bot
            else ""
        )

        # Pending for approval
        # all_checks_top = top_details.get('flexure_check')
        # all_checks_bot = bot_details.get('flexure_check')
        # all_checks = all_checks_top and all_checks_bot

        # Formatter instance for DCR formatting
        markdown_content = ""
        formatter = Formatter()

        # Handle top result data
        if top_result_data:
            top_rebar_1 = top_result_data["Value"][0]
            top_rebar_2 = top_result_data["Value"][1]
            area_top = top_result_data["Value"][7]
            Mu_top = self._limiting_case_flexure_top_details["forces"]["Value"][0]
            Mn_top = top_result_data["Value"][9]
            DCR_top = top_result_data["Value"][10]

            rebar_top = f"{top_rebar_1}" + (
                f" ++ {top_rebar_2}" if top_rebar_2 != "-" else ""
            )
            formatted_DCR_top = formatter.DCR(DCR_top)

            markdown_content += (
                f"Top longitudinal rebar: {rebar_top}, $A_{{s,top}}$ = {area_top} cm², $M_u$ = {Mu_top} kNm, "
                f"$\\phi M_n$ = {Mn_top} kNm → {formatted_DCR_top} {warning_top}\n\n"
            )
        else:
            markdown_content += "No top moment to check.\n\n"

        # Handle bottom result data
        if bot_result_data:
            bot_rebar_1 = bot_result_data["Value"][0]
            bot_rebar_2 = bot_result_data["Value"][1]
            area_bot = bot_result_data["Value"][7]
            Mu_bot = self._limiting_case_flexure_bot_details["forces"]["Value"][1]
            Mn_bot = bot_result_data["Value"][9]
            DCR_bot = bot_result_data["Value"][10]

            rebar_bot = f"{bot_rebar_1}" + (
                f" ++ {bot_rebar_2}" if bot_rebar_2 != "-" else ""
            )
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
            warnings.warn(
                "Shear design has not been performed yet. Call check_shear() or "
                "design_shear() first.",
                UserWarning,
            )
            self._md_shear_results = "Shear results are not available."
            return None

        # Check if limiting case details exist
        shear_details = self._limiting_case_shear_details or {}
        # Use limiting case results
        limiting_reinforcement = shear_details.get("shear_reinforcement")
        limiting_forces = shear_details.get("forces")
        limiting_shear_concrete = shear_details.get("shear_concrete")
        checks_pass = shear_details.get("checks_pass")
        markdown_content = ""
        if shear_details:
            # Create FUFormatter instance and format FU value
            formatter = Formatter()
            formatted_DCR = formatter.DCR(limiting_shear_concrete["Value"][-1])
            if self._A_v == 0 * cm:
                rebar_v = "not assigned"
            else:
                rebar_v = (
                    f"{int(limiting_reinforcement['Value'][0])}eØ{limiting_reinforcement['Value'][1]}/"
                    f"{limiting_reinforcement['Value'][2]} cm"
                )
            # Limitng cases checks
            warning = (
                "⚠️ Some checks failed, see detailed results."
                if not checks_pass
                else ""
            )
            if (
                self.concrete.design_code == "ACI 318-19"
                or self.concrete.design_code == "CIRSOC 201-25"
            ):
                markdown_content = (
                    f"Shear reinforcing {rebar_v}, $A_v$={limiting_reinforcement['Value'][6]} cm²/m"
                    f", $V_u$={limiting_forces['Value'][1]} kN, $\\phi V_n$={limiting_shear_concrete['Value'][7]} kN → {formatted_DCR} {warning}"
                )  # noqa: E501
            else:  # self.concrete.design_code == "EN 1992-2004"
                markdown_content = (
                    f"Shear reinforcing {rebar_v}, $A_{{sw}}$={limiting_reinforcement['Value'][6]} cm²/m"
                    f", $V_{{Ed,2}}$={limiting_forces['Value'][1]} kN, $V_{{Rd}}$={limiting_shear_concrete['Value'][6]} kN → {formatted_DCR} {warning}"
                )  # noqa: E501
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
        if not hasattr(self, "_md_properties"):
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
            warnings.warn(
                "Flexural check has not been performed yet. Call _check_flexure or "
                "design_flexure first.",
                UserWarning,
            )
            self._md_flexure_results = "Flexure results are not available."
            return None

        # Determine which results to display (limiting case by default)
        if force:
            if force.id not in self._flexure_results_detailed_list:
                raise ValueError(
                    f"No results found for Forces object with ID {force.id}."
                )
            result_data = self._flexure_results_detailed_list[force.id]
            top_result_data = result_data["flexure_capacity_top"]
            bot_result_data = result_data["flexure_capacity_bot"]
            forces_result = result_data["forces"]
            min_max_result = result_data["min_max"]
        else:
            if self._limiting_case_flexure_top_details is None:
                raise ValueError("Top limiting case details are not available.")
            if self._limiting_case_flexure_bot_details is None:
                raise ValueError("Bottom limiting case details are not available.")
            # Use the worst-case top and bottom scenarios
            top_result_data = self._limiting_case_flexure_top_details[
                "flexure_capacity_top"
            ]
            bot_result_data = self._limiting_case_flexure_bot_details[
                "flexure_capacity_bot"
            ]
            forces_result = {
                "Design_forces": [
                    "Top max moment",
                    "Bottom max moment",
                ],
                "Variable": ["Mu,top", "Mu,bot"],
                "Value": [
                    round(
                        self._limiting_case_flexure_top_details["forces"]["Value"][0], 2
                    ),
                    round(
                        self._limiting_case_flexure_bot_details["forces"]["Value"][1], 2
                    ),
                ],
                "Unit": ["kNm", "kNm"],
            }
            min_max_result = {
                "Check": [
                    "Min/Max As rebar top",
                    "Minimum spacing top",
                    "Min/Max As rebar bottom",
                    "Minimum spacing bottom",
                ],
                "Unit": ["cm²", "mm", "cm²", "mm"],
                "Value": [
                    round(
                        self._limiting_case_flexure_top_details["min_max"]["Value"][0],
                        2,
                    ),  # Top limiting case As
                    round(
                        self._limiting_case_flexure_top_details["min_max"]["Value"][1],
                        2,
                    ),  # Top limiting case spacing
                    round(
                        self._limiting_case_flexure_bot_details["min_max"]["Value"][2],
                        2,
                    ),  # Bottom limiting case As
                    round(
                        self._limiting_case_flexure_bot_details["min_max"]["Value"][3],
                        2,
                    ),  # Bottom limiting case spacing
                ],
                "Min.": [
                    round(
                        self._limiting_case_flexure_top_details["min_max"]["Min."][0], 2
                    ),  # Top limiting case As_min
                    self._limiting_case_flexure_top_details["min_max"]["Min."][
                        1
                    ],  # Top limiting case spacing min
                    round(
                        self._limiting_case_flexure_bot_details["min_max"]["Min."][2], 2
                    ),  # Bottom limiting case As_min
                    self._limiting_case_flexure_bot_details["min_max"]["Min."][
                        3
                    ],  # Bottom limiting case spacing min
                ],
                "Max.": [
                    round(
                        self._limiting_case_flexure_top_details["min_max"]["Max."][0], 2
                    ),  # Top limiting case As_max
                    "",  # No max constraint for spacing
                    round(
                        self._limiting_case_flexure_bot_details["min_max"]["Max."][2], 2
                    ),  # Bottom limiting case As_max
                    "",  # No max constraint for spacing
                ],
                "Ok?": [
                    self._limiting_case_flexure_top_details["min_max"]["Ok?"][
                        0
                    ],  # Top As check
                    self._limiting_case_flexure_top_details["min_max"]["Ok?"][
                        1
                    ],  # Top spacing check
                    self._limiting_case_flexure_bot_details["min_max"]["Ok?"][
                        2
                    ],  # Bottom As check
                    self._limiting_case_flexure_bot_details["min_max"]["Ok?"][
                        3
                    ],  # Bottom spacing check
                ],
            }

        # Create TablePrinter instances for detailed display
        print("===== BEAM FLEXURE DETAILED RESULTS =====")
        materials_printer = TablePrinter("MATERIALS")
        materials_printer.print_table_data(self._materials_flexure, headers="keys")

        geometry_printer = TablePrinter("GEOMETRY")
        geometry_printer.print_table_data(self._geometry_flexure, headers="keys")

        forces_printer = TablePrinter("FORCES")
        forces_printer.print_table_data(forces_result, headers="keys")

        min_max_printer = TablePrinter("MAX AND MIN LIMIT CHECKS")
        min_max_printer.print_table_data(min_max_result, headers="keys")

        capacity_printer = TablePrinter("FLEXURAL CAPACITY - TOP")
        capacity_printer.print_table_data(top_result_data, headers="keys")
        capacity_printer = TablePrinter("FLEXURAL CAPACITY - BOTTOM")
        capacity_printer.print_table_data(bot_result_data, headers="keys")

    def flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        """
        Prints detailed flexure results in Word.

        Parameters
        ----------
        forces : Forces, optional
            The specific Forces object to display results for. If None, displays results for the limiting case.
        """
        if not self._flexure_checked:
            warnings.warn(
                "Flexural check has not been performed yet. Call _check_flexure or "
                "design_flexure first.",
                UserWarning,
            )
            self._md_flexure_results = "Flexural results are not available."
            return None

        # Determine which results to display (limiting case by default)
        if force:
            if force.id not in self._flexure_results_detailed_list:
                raise ValueError(
                    f"No results found for Forces object with ID {force.id}."
                )
            result_data = self._flexure_results_detailed_list[force.id]
            top_result_data = result_data["flexure_capacity_top"]
            bot_result_data = result_data["flexure_capacity_bot"]
            forces_result = result_data["forces"]
            min_max_result = result_data["min_max"]
        else:
            if self._limiting_case_flexure_top_details is None:
                raise ValueError("Top limiting case details are not available.")
            if self._limiting_case_flexure_bot_details is None:
                raise ValueError("Bottom limiting case details are not available.")
            # Use the worst-case top and bottom scenarios
            top_result_data = self._limiting_case_flexure_top_details[
                "flexure_capacity_top"
            ]
            bot_result_data = self._limiting_case_flexure_bot_details[
                "flexure_capacity_bot"
            ]
            forces_result = {
                "Design forces": [
                    "Top max moment",
                    "Bottom max moment",
                ],
                "Variable": ["Mu,top", "Mu,bot"],
                "Value": [
                    round(
                        self._limiting_case_flexure_top_details["forces"]["Value"][0], 2
                    ),
                    round(
                        self._limiting_case_flexure_bot_details["forces"]["Value"][1], 2
                    ),
                ],
                "Unit": ["kNm", "kNm"],
            }
            min_max_result = {
                "Check": [
                    "Min/Max As rebar top",
                    "Minimum spacing top",
                    "Min/Max As rebar bottom",
                    "Minimum spacing bottom",
                ],
                "Unit": ["cm²", "mm", "cm²", "mm"],
                "Value": [
                    round(
                        self._limiting_case_flexure_top_details["min_max"]["Value"][0],
                        2,
                    ),  # Top limiting case As
                    round(
                        self._limiting_case_flexure_top_details["min_max"]["Value"][1],
                        2,
                    ),  # Top limiting case spacing
                    round(
                        self._limiting_case_flexure_bot_details["min_max"]["Value"][2],
                        2,
                    ),  # Bottom limiting case As
                    round(
                        self._limiting_case_flexure_bot_details["min_max"]["Value"][3],
                        2,
                    ),  # Bottom limiting case spacing
                ],
                "Min.": [
                    round(
                        self._limiting_case_flexure_top_details["min_max"]["Min."][0], 2
                    ),  # Top limiting case As_min
                    self._limiting_case_flexure_top_details["min_max"]["Min."][
                        1
                    ],  # Top limiting case spacing min
                    round(
                        self._limiting_case_flexure_bot_details["min_max"]["Min."][2], 2
                    ),  # Bottom limiting case As_min
                    self._limiting_case_flexure_bot_details["min_max"]["Min."][
                        3
                    ],  # Bottom limiting case spacing min
                ],
                "Max.": [
                    round(
                        self._limiting_case_flexure_top_details["min_max"]["Max."][0], 2
                    ),  # Top limiting case As_max
                    "",  # No max constraint for spacing
                    round(
                        self._limiting_case_flexure_bot_details["min_max"]["Max."][2], 2
                    ),  # Bottom limiting case As_max
                    "",  # No max constraint for spacing
                ],
                "Ok?": [
                    self._limiting_case_flexure_top_details["min_max"]["Ok?"][
                        0
                    ],  # Top As check
                    self._limiting_case_flexure_top_details["min_max"]["Ok?"][
                        1
                    ],  # Top spacing check
                    self._limiting_case_flexure_bot_details["min_max"]["Ok?"][
                        2
                    ],  # Bottom As check
                    self._limiting_case_flexure_bot_details["min_max"]["Ok?"][
                        3
                    ],  # Bottom spacing check
                ],
            }

        # Convert output Dicts into DataFrames
        df_materials = pd.DataFrame(self._materials_flexure)
        df_geometry = pd.DataFrame(self._geometry_flexure)
        df_forces = pd.DataFrame(forces_result)
        df_data_min_max = pd.DataFrame(min_max_result)
        df_flexure_capacity_top = pd.DataFrame(top_result_data)
        df_flexure_capacity_bottom = pd.DataFrame(bot_result_data)

        # Create a document builder instance
        doc_builder = DocumentBuilder(title="Concrete beam flexure check")

        # Add first section and table
        doc_builder.add_heading(f"Beam {self.label} flexure check", level=1)
        doc_builder.add_text(
            f"Made with mento {MENTO_VERSION}. Design code: {self.concrete.design_code}"
        )
        doc_builder.add_heading("Materials", level=2)
        doc_builder.add_table_data(df_materials)
        doc_builder.add_table_data(df_geometry)
        doc_builder.add_table_data(df_forces)

        # Add third section for limit checks
        doc_builder.add_heading("Limit Checks", level=2)
        doc_builder.add_table_data(df_data_min_max)

        # Add second section for flexural checks
        doc_builder.add_heading("Flexural Capacity Top", level=2)
        doc_builder.add_table_dcr(df_flexure_capacity_top)
        doc_builder.add_heading("Flexural Capacity Bottom", level=2)
        doc_builder.add_table_dcr(df_flexure_capacity_bottom)

        # Save the Word doc
        doc_builder.save(
            f"Beam {self.label} flexure check {self.concrete.design_code}.docx"
        )

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
            warnings.warn(
                "Shear check has not been performed yet. Call check_shear or "
                "design_shear first.",
                UserWarning,
            )
            self._md_shear_results = "Shear results are not available."
            return None
        # Determine which results to display (limiting case by default)
        if force:
            force_id = force.id
            if force_id not in self._shear_results_detailed_list:
                raise ValueError(
                    f"No results found for Forces object with ID {force_id}."
                )
            result_data = self._shear_results_detailed_list[force_id]
        else:
            # Default to limiting case
            result_data = self._limiting_case_shear_details

        # Create a TablePrinter instance and display tables
        print("===== BEAM SHEAR DETAILED RESULTS =====")
        materials_printer = TablePrinter("MATERIALS")
        materials_printer.print_table_data(self._materials_shear, headers="keys")
        geometry_printer = TablePrinter("GEOMETRY")
        geometry_printer.print_table_data(self._geometry_shear, headers="keys")
        forces_printer = TablePrinter("FORCES")
        forces_printer.print_table_data(result_data["forces"], headers="keys")
        steel_printer = TablePrinter("SHEAR STRENGTH")
        steel_printer.print_table_data(
            result_data["shear_reinforcement"], headers="keys"
        )
        min_max_printer = TablePrinter("MAX AND MIN LIMIT CHECKS")
        min_max_printer.print_table_data(result_data["min_max"], headers="keys")
        concrete_printer = TablePrinter("CONCRETE STRENGTH")
        concrete_printer.print_table_data(result_data["shear_concrete"], headers="keys")

    def shear_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        """
        Prints detailed shear results in Word.

        Parameters
        ----------
        forces : Forces, optional
            The specific Forces object to display results for. If None, displays results for the limiting case.
        """
        if not self._shear_checked:
            warnings.warn(
                "Shear check has not been performed yet. Call check_shear or "
                "design_shear first.",
                UserWarning,
            )
            self._md_shear_results = "Shear results are not available."
            return None
        # Determine which results to display (limiting case by default)
        if force:
            force_id = force.id
            if force_id not in self._shear_results_detailed_list:
                raise ValueError(
                    f"No results found for Forces object with ID {force_id}."
                )
            result_data = self._shear_results_detailed_list[force_id]
        else:
            # Default to limiting case
            result_data = self._limiting_case_shear_details

        # Convert output Dicts into DataFrames
        df_materials = pd.DataFrame(self._materials_shear)
        df_geometry = pd.DataFrame(self._geometry_shear)
        df_forces = pd.DataFrame(result_data["forces"])
        df_shear_reinforcement = pd.DataFrame(result_data["shear_reinforcement"])
        df_data_min_max = pd.DataFrame(result_data["min_max"])
        df_shear_concrete = pd.DataFrame(result_data["shear_concrete"])

        # Create a document builder instance
        doc_builder = DocumentBuilder(title="Concrete beam shear check")

        # Add first section and table
        doc_builder.add_heading(f"Beam {self.label} shear check", level=1)
        doc_builder.add_text(
            f"Made with mento {MENTO_VERSION}. Design code: {self.concrete.design_code}"
        )
        doc_builder.add_heading("Materials", level=2)
        doc_builder.add_table_data(df_materials)
        doc_builder.add_table_data(df_geometry)
        doc_builder.add_table_data(df_forces)

        # Add second section and another table (can use different data)
        doc_builder.add_heading("Limit checks", level=2)
        doc_builder.add_table_min_max(df_data_min_max)
        doc_builder.add_heading("Design checks", level=2)
        doc_builder.add_table_data(df_shear_reinforcement)
        doc_builder.add_table_dcr(df_shear_concrete)

        # Save the Word doc
        doc_builder.save(
            f"Beam {self.label} shear check {self.concrete.design_code}.docx"
        )

    def _format_longitudinal_rebar_string(
        self, n1: int, d_b1: Quantity, n2: int = 0, d_b2: Quantity = 0 * mm
    ) -> str:
        """
        Returns a formatted string representing the rebars and their diameters.

        Example output:
            "2Ø10" (if only n1 and d_b1 are defined)
            "2Ø10+2Ø8" (if both n1/d_b1 and n2/d_b2 are defined)
        """
        if n1 == 0:
            rebar_string = "-"
        else:
            rebar_string = f"{n1}Ø{int(d_b1.magnitude)}"
            if n2 > 0 and d_b2.magnitude > 0:  # Check if n2 and d_b2 are defined
                rebar_string += f"+{n2}Ø{int(d_b2.magnitude)}"
        return rebar_string

    ##########################################################
    # PLOT BEAM SECTION WITH REBAR
    ##########################################################

    def _plot_rebar_layer(
        self,
        width_cm: float,
        height_cm: float,
        c_c_cm: float,
        stirrup_d_b_cm: float,
        layers_spacing_cm: float,
        n1: int,
        d_b1: Quantity,
        n2: int,
        d_b2: Quantity,
        max_db: Quantity,
        is_bottom: bool = True,
        is_second_layer: bool = False,
    ) -> None:
        """
        Helper method to plot a single layer of rebars.
        """
        # Calculate y-position based on layer and bottom/top
        y_base = (
            c_c_cm + stirrup_d_b_cm
            if is_bottom
            else height_cm - c_c_cm - stirrup_d_b_cm
        )

        if is_second_layer:
            y_base += (
                layers_spacing_cm + max_db.to("cm").magnitude
                if is_bottom
                else -layers_spacing_cm - max_db.to("cm").magnitude
            )

        # Plot side bars (position 1 or 3)
        if n1 > 0:
            diameter_cm = d_b1.to("cm").magnitude
            y_center = (
                y_base + diameter_cm / 2 if is_bottom else y_base - diameter_cm / 2
            )
            for i in range(n1):
                x = (
                    c_c_cm
                    + stirrup_d_b_cm
                    + diameter_cm / 2
                    + (
                        i
                        * (width_cm - 2 * (c_c_cm + stirrup_d_b_cm + diameter_cm / 2))
                        / (n1 - 1)
                    )
                )
                circle = Circle(
                    (x, y_center),
                    diameter_cm / 2,
                    color=CUSTOM_COLORS["dark_gray"],
                    fill=True,
                )
                self._ax.add_patch(circle)

        # Plot intermediate bars (position 2 or 4)
        if n2 > 0:
            diameter_cm = d_b2.to("cm").magnitude
            y_center = (
                y_base + diameter_cm / 2 if is_bottom else y_base - diameter_cm / 2
            )
            for i in range(n2):
                x = (
                    c_c_cm
                    + stirrup_d_b_cm
                    + diameter_cm / 2
                    + (
                        (i + 1)
                        * (width_cm - 2 * (c_c_cm + stirrup_d_b_cm + diameter_cm / 2))
                        / (n2 + 1)
                    )
                )
                circle = Circle(
                    (x, y_center),
                    diameter_cm / 2,
                    color=CUSTOM_COLORS["dark_gray"],
                    fill=True,
                )
                self._ax.add_patch(circle)

    def plot(self) -> None:
        """
        Plots the rectangular section with a dark gray border, light gray hatch, and dimensions.
        Also plots the stirrup with rounded corners and thickness.
        """

        # Convert dimensions to consistent units (cm)
        width_cm: float = self.width.to("cm").magnitude
        height_cm: float = self.height.to("cm").magnitude
        c_c_cm: float = self.c_c.to("cm").magnitude
        stirrup_d_b_cm: float = self._stirrup_d_b.to("cm").magnitude
        layers_spacing_cm: float = self.settings.layers_spacing.to("cm").magnitude

        # Create figure and axis
        fig, self._ax = plt.subplots()

        # Create a rectangle patch for the section
        rect = Rectangle(
            (0, 0),
            width_cm,
            height_cm,
            linewidth=1.3,
            edgecolor=CUSTOM_COLORS["dark_gray"],
            facecolor=CUSTOM_COLORS["light_gray"],
        )
        self._ax.add_patch(rect)

        # Calculate stirrup dimensions
        c_c = self.c_c.to("cm").magnitude
        stirrup_width = self.width.to("cm").magnitude - 2 * c_c
        stirrup_height = self.height.to("cm").magnitude - 2 * c_c
        stirrup_thickness = self._stirrup_d_b.to("cm").magnitude

        # Create rounded corners for the stirrup
        inner_radius = stirrup_thickness * 2
        outer_radius = stirrup_thickness * 3

        # Create the outer rounded rectangle for the stirrup
        outer_rounded_rect = FancyBboxPatch(
            (c_c, c_c),  # Bottom-left corner
            stirrup_width,  # Width
            stirrup_height,  # Height
            boxstyle=f"Round, pad=0, rounding_size={outer_radius}",  # Rounded corners
            edgecolor=CUSTOM_COLORS["dark_blue"],
            facecolor="white",
            linewidth=1,
        )
        self._ax.add_patch(outer_rounded_rect)

        # Create the inner rounded rectangle for the stirrup (offset by thickness)
        inner_rounded_rect = FancyBboxPatch(
            (c_c + stirrup_thickness, c_c + stirrup_thickness),  # Bottom-left corner
            stirrup_width - 2 * stirrup_thickness,  # Width
            stirrup_height - 2 * stirrup_thickness,  # Height
            boxstyle=f"Round, pad=0, rounding_size={inner_radius}",  # Rounded corners
            edgecolor=CUSTOM_COLORS["dark_blue"],
            facecolor=CUSTOM_COLORS["light_gray"],
            linewidth=1,
        )
        self._ax.add_patch(inner_rounded_rect)

        # Set plot limits with some padding
        padding = max(width_cm, height_cm) * 0.2
        self._ax.set_xlim(-padding, width_cm + padding)
        self._ax.set_ylim(-padding, height_cm + padding)

        # Text and dimension offsets
        dim_offset = 2.5
        text_offset = dim_offset + 2
        # Add width dimension
        self._ax.annotate(
            "",  # No text here, text is added separately
            xy=(0, -dim_offset),  # Start of arrow (left side)
            xytext=(width_cm, -dim_offset),  # End of arrow (right side)
            arrowprops={
                "arrowstyle": "<->",
                "lw": 1,
                "color": CUSTOM_COLORS["dark_blue"],
            },
        )
        if self.concrete.unit_system == "imperial":
            # Example: format to 2 decimal places, then use pint's compact (~P) format
            width = "{:.0f~P}".format(self.width.to("inch"))
            height = "{:.0f~P}".format(self.height.to("inch"))
        else:
            width = "{:.0f~P}".format(self.width.to("cm"))
            height = "{:.0f~P}".format(self.height.to("cm"))
        # Add width dimension text below the arrow
        self._ax.text(
            width_cm / 2,  # Center of the arrow
            -text_offset,  # Slightly below the arrow
            width,
            ha="center",
            va="top",
            color=CUSTOM_COLORS["dark_gray"],
        )

        # Add height dimension
        self._ax.annotate(
            "",  # No text here, text is added separately
            xy=(-dim_offset, 0),  # Start of arrow (bottom)
            xytext=(-dim_offset, height_cm),  # End of arrow (top)
            arrowprops={
                "arrowstyle": "<->",
                "lw": 1,
                "color": CUSTOM_COLORS["dark_blue"],
            },
        )
        # Add height dimension text to the left of the arrow
        self._ax.text(
            -text_offset,  # Slightly to the left of the arrow
            height_cm / 2,  # Center of the arrow
            height,
            ha="right",
            va="center",
            color=CUSTOM_COLORS["dark_gray"],
            rotation=90,  # Rotate text vertically
        )

        # Set aspect of the plot to be equal
        self._ax.set_aspect("equal")
        # Remove axes for better visualization
        self._ax.axis("off")

        # Calculate rebar positions
        # Bottom rebars
        self._plot_rebar_layer(
            width_cm,
            height_cm,
            c_c_cm,
            stirrup_d_b_cm,
            layers_spacing_cm,
            self._n1_b,
            self._d_b1_b,
            self._n2_b,
            self._d_b2_b,
            max_db=self._d_b1_b,
            is_bottom=True,
        )
        self._plot_rebar_layer(
            width_cm,
            height_cm,
            c_c_cm,
            stirrup_d_b_cm,
            layers_spacing_cm,
            self._n3_b,
            self._d_b3_b,
            self._n4_b,
            self._d_b4_b,
            max_db=self._d_b1_b,
            is_bottom=True,
            is_second_layer=True,
        )

        # Top rebars
        self._plot_rebar_layer(
            width_cm,
            height_cm,
            c_c_cm,
            stirrup_d_b_cm,
            layers_spacing_cm,
            self._n1_t,
            self._d_b1_t,
            self._n2_t,
            self._d_b2_t,
            max_db=self._d_b1_t,
            is_bottom=False,
        )
        self._plot_rebar_layer(
            width_cm,
            height_cm,
            c_c_cm,
            stirrup_d_b_cm,
            layers_spacing_cm,
            self._n3_t,
            self._d_b3_t,
            self._n4_t,
            self._d_b4_t,
            max_db=self._d_b1_t,
            is_bottom=False,
            is_second_layer=True,
        )

        # Show plot
        plt.show()
