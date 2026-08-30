from dataclasses import dataclass
from IPython.display import Markdown, display
from typing import TYPE_CHECKING, Any, Optional, Dict, Tuple

if TYPE_CHECKING:
    from matplotlib.figure import Figure
from pint import Quantity
import numpy as np
import pandas as pd
from pandas import DataFrame
import math
import warnings
from numbers import Integral
# from devtools import debug

from mento.rectangular import RectangularSection
from mento.material import (
    Concrete_ACI_318_19,
    Concrete_EN_1992_2004,
)
from mento.rebar import Rebar
from mento.units import MPa, mm, inch, kN, m, cm, kNm, dimensionless
from mento.results import Formatter, TablePrinter
from mento.i18n import get_language, translate
from mento.forces import Forces
from mento.settings import BeamSettings
from mento.documents import flexure_report_doc, shear_report_doc
from mento.plots import plot_beam_section
from mento.report_tables import build_flexure_report, build_shear_report
from mento.design_results import (
    FlexureCheck,
    FlexureDesign,
    SectionReinforcement,
    ShearCheck,
    ShearDesign,
    build_flexure_design,
    build_reinforcement,
    build_shear_design,
    capture_flexure_check,
    capture_shear_check,
)

from mento.codes.EN_1992_2004_beam import (
    _check_shear_EN_1992_2004,
    _check_flexure_EN_1992_2004,
    _design_shear_EN_1992_2004,
    _design_flexure_EN_1992_2004,
)
from mento.codes.ACI_318_19_beam import (
    _check_shear_ACI_318_19,
    _design_shear_ACI_318_19,
    _check_flexure_ACI_318_19,
    _design_flexure_ACI_318_19,
)


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

    settings: Optional[BeamSettings] = None

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
        self.mode = "beam"
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
        # One immutable result per load combination, appended by check_flexure and
        # check_shear. The attributes above only hold the combination that ran last;
        # enveloping is a pure operation over these lists.
        self._flexure_checks: list[FlexureCheck] = []
        self._shear_checks: list[ShearCheck] = []
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
        self._doubly_reinforced = False  # Tracks if doubly reinforced section is used

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

    def _create_rebar_designer(self) -> Rebar:
        """Factory for the longitudinal reinforcement optimizer."""
        return Rebar(self)

    def _apply_longitudinal_design_bot(self, design: Any) -> None:
        """Apply a discrete design to the bottom reinforcement.

        ``design`` is any mapping with the rebar-designer keys (n_1..n_4,
        d_b1..d_b4); the designer hands back a pandas row, tests pass plain
        dicts.
        """
        self.set_longitudinal_rebar_bot(
            int(design.get("n_1", 0)),
            design.get("d_b1"),
            int(design.get("n_2", 0)),
            design.get("d_b2"),
            int(design.get("n_3", 0)),
            design.get("d_b3"),
            int(design.get("n_4", 0)),
            design.get("d_b4"),
        )

    def _apply_longitudinal_design_top(self, design: Any) -> None:
        """Apply a discrete design to the top reinforcement.

        See :meth:`_apply_longitudinal_design_bot` for the accepted shape.
        """
        self.set_longitudinal_rebar_top(
            int(design.get("n_1", 0)),
            design.get("d_b1"),
            int(design.get("n_2", 0)),
            design.get("d_b2"),
            int(design.get("n_3", 0)),
            design.get("d_b3"),
            int(design.get("n_4", 0)),
            design.get("d_b4"),
        )

    def _clear_top_longitudinal(self) -> None:
        """Reset the top reinforcement to the default placeholder bars."""
        if self.concrete.unit_system == "metric":
            self.set_longitudinal_rebar_top(0, 0 * mm)
        else:
            self.set_longitudinal_rebar_top(0, 0 * inch)

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
            self.flexure_design_results_bot: Any = None
            self.flexure_design_results_top: Any = None
            self._A_s_bool_bot: bool = False
            self._A_s_bool_top: bool = False

    def _initialize_EN_1992_2004_attributes(self) -> None:
        if isinstance(self.concrete, Concrete_EN_1992_2004):
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
            self._f_ywk = self.steel_bar.f_y
            self._f_ywd: Quantity = 0 * MPa
            self._f_cd: Quantity = 0 * MPa
            self._f_cd_shear: Quantity = 0 * MPa
            self._A_p = 0 * cm**2  # No prestressed for now
            self._sigma_cp: Quantity = 0 * MPa
            self._theta: float = 0
            self._cot_theta: float = 0
            self._z: Quantity = 0 * cm
            self._M_Rd_bot: Quantity = 0 * kNm
            self._M_Rd_top: Quantity = 0 * kNm
            self._M_Ed_bot: Quantity = 0 * kNm
            self._M_Ed_top: Quantity = 0 * kNm
            # Shared with the ACI branch: the flexural design driver and the
            # results tables read these before any check has run.
            self._A_s_min_bot: Quantity = 0 * cm**2
            self._A_s_min_top: Quantity = 0 * cm**2
            self._A_s_max_bot: Quantity = 0 * cm**2
            self._A_s_max_top: Quantity = 0 * cm**2
            self.flexure_design_results_bot: Any = None
            self.flexure_design_results_top: Any = None

    ##########################################################
    # SET LONGITUDINAL AND TRANSVERSE REBAR AND UPDATE ATTRIBUTES
    ##########################################################

    def set_transverse_rebar(
        self,
        n_stirrups: int = 0,
        d_b: Quantity = 0 * mm,
        s_l: Quantity = 0 * cm,
    ) -> None:
        """Set transverse reinforcement or clear it with an all-zero input."""

        # Reject booleans and non-integer stirrup counts.
        if isinstance(n_stirrups, bool) or not isinstance(n_stirrups, Integral):
            raise TypeError("n_stirrups must be an integer.")
        n_stirrups = int(n_stirrups)

        # Diameter and spacing must be physical lengths.
        if not isinstance(d_b, Quantity) or not d_b.check("[length]"):
            raise TypeError("d_b must be a length Quantity.")
        if not isinstance(s_l, Quantity) or not s_l.check("[length]"):
            raise TypeError("s_l must be a length Quantity.")

        diameter_mm = d_b.to("mm").magnitude
        spacing_mm = s_l.to("mm").magnitude
        if not math.isfinite(diameter_mm):
            raise ValueError("d_b must be finite.")
        if not math.isfinite(spacing_mm):
            raise ValueError("s_l must be finite.")

        # An all-zero input explicitly removes the transverse reinforcement.
        if n_stirrups == 0 and diameter_mm == 0 and spacing_mm == 0:
            self._stirrup_n = 0
            self._stirrup_d_b = d_b
            self._stirrup_s_l = s_l
            self._A_v = 0 * cm**2 / m
            self._update_effective_heights()
            return

        # Every non-empty reinforcement configuration must be strictly positive.
        if n_stirrups <= 0:
            raise ValueError("n_stirrups must be greater than zero.")
        if diameter_mm <= 0:
            raise ValueError("d_b must be greater than zero.")
        if spacing_mm <= 0:
            raise ValueError("s_l must be greater than zero.")

        # Store the inputs only after all validations pass.
        self._stirrup_n = n_stirrups
        self._stirrup_d_b = d_b
        self._stirrup_s_l = s_l

        # A closed stirrup contributes two vertical legs.
        n_legs = n_stirrups * 2
        A_db = d_b**2 * math.pi / 4
        A_vs = n_legs * A_db

        # Calculate the transverse reinforcement area per unit length.
        self._A_v = A_vs / s_l

        # Recalculate effective depths using the new stirrup diameter.
        self._update_effective_heights()

    def _len_unit(self) -> Quantity:
        # Base length unit for this section
        return 1 * mm if getattr(self.concrete, "unit_system", "metric") == "metric" else 1 * inch

    def _area_unit(self) -> Quantity:
        L = self._len_unit()
        return L**2

    def set_longitudinal_rebar_bot(
        self,
        n1: int,
        d_b1: Quantity | None,
        n2: int = 0,
        d_b2: Quantity | None = None,
        n3: int = 0,
        d_b3: Quantity | None = None,
        n4: int = 0,
        d_b4: Quantity | None = None,
    ) -> None:
        L = self._len_unit()
        self._n1_b = n1
        self._d_b1_b = d_b1 if d_b1 is not None else 0 * L
        self._n2_b = n2
        self._d_b2_b = d_b2 if d_b2 is not None else 0 * L
        self._n3_b = n3
        self._d_b3_b = d_b3 if d_b3 is not None else 0 * L
        self._n4_b = n4
        self._d_b4_b = d_b4 if d_b4 is not None else 0 * L
        self._update_longitudinal_rebar_attributes()

    def set_longitudinal_rebar_top(
        self,
        n1: int,
        d_b1: Quantity | None,
        n2: int = 0,
        d_b2: Quantity | None = None,
        n3: int = 0,
        d_b3: Quantity | None = None,
        n4: int = 0,
        d_b4: Quantity | None = None,
    ) -> None:
        L = self._len_unit()
        self._n1_t = n1
        self._d_b1_t = d_b1 if d_b1 is not None else 0 * L
        self._n2_t = n2
        self._d_b2_t = d_b2 if d_b2 is not None else 0 * L
        self._n3_t = n3
        self._d_b3_t = d_b3 if d_b3 is not None else 0 * L
        self._n4_t = n4
        self._d_b4_t = d_b4 if d_b4 is not None else 0 * L
        self._update_longitudinal_rebar_attributes()

    def _calculate_longitudinal_rebar_area(self) -> None:
        """Calculate total rebar area (robust to None or zero diameters)."""
        A0 = 0 * self._area_unit()

        def area(n: int, d: Quantity | None) -> Quantity:
            # Treat None or zero diameter as zero area; works with pint
            if d is None:
                return A0
            mag = getattr(d, "magnitude", None)
            if mag is not None and mag == 0:
                return A0
            return n * (d**2) * np.pi / 4

        # Bottom
        self._A_s_bot = (
            area(self._n1_b, self._d_b1_b)
            + area(self._n2_b, self._d_b2_b)
            + area(self._n3_b, self._d_b3_b)
            + area(self._n4_b, self._d_b4_b)
        )

        # Top
        self._A_s_top = (
            area(self._n1_t, self._d_b1_t)
            + area(self._n2_t, self._d_b2_t)
            + area(self._n3_t, self._d_b3_t)
            + area(self._n4_t, self._d_b4_t)
        )

    def _calculate_min_clear_spacing(self) -> None:
        """
        Calculates the maximum clear spacing between bars for the bottom rebar layers.

        Returns:
            Quantity: The maximum clear spacing between bars in either the first or second layer.
        """

        def layer_clear_spacing(n_a: int, d_a: Quantity, n_b: int, d_b: Quantity) -> Quantity:
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
        spacing_layer1_b = layer_clear_spacing(self._n1_b, self._d_b1_b, self._n2_b, self._d_b2_b)
        spacing_layer2_b = layer_clear_spacing(self._n3_b, self._d_b3_b, self._n4_b, self._d_b4_b)

        # Return the maximum clear spacing between the two layers
        self._available_s_bot = min(spacing_layer1_b, spacing_layer2_b)

        # AVAIABLE CLEAR SPACING FOR TOP BARS
        # Calculate clear spacing for each layer
        spacing_layer1_t = layer_clear_spacing(self._n1_t, self._d_b1_t, self._n2_t, self._d_b2_t)
        spacing_layer2_t = layer_clear_spacing(self._n3_t, self._d_b3_t, self._n4_t, self._d_b4_t)

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
        y3_b = max(self._d_b1_b, self._d_b2_b) + self.settings.layers_spacing + self._d_b3_b / 2
        y4_b = max(self._d_b1_b, self._d_b2_b) + self.settings.layers_spacing + self._d_b4_b / 2
        # Calculate the total area of each layer
        area_1_b = self._n1_b * self._d_b1_b**2 * np.pi / 4.0
        area_2_b = self._n2_b * self._d_b2_b**2 * np.pi / 4.0
        area_3_b = self._n3_b * self._d_b3_b**2 * np.pi / 4.0
        area_4_b = self._n4_b * self._d_b4_b**2 * np.pi / 4.0

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
        y3_t = max(self._d_b1_t, self._d_b2_t) + self.settings.layers_spacing + self._d_b3_t / 2
        y4_t = max(self._d_b1_t, self._d_b2_t) + self.settings.layers_spacing + self._d_b4_t / 2

        # Calculate the total area of each layer
        area_1_t = self._n1_t * self._d_b1_t**2 * np.pi / 4  # Area proportional to number of bars and their diameter
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
        if self.concrete.design_code == "ACI 318-19" or self.concrete.design_code == "CIRSOC 201-25":
            _design_flexure_ACI_318_19(self, max_M_y_bot, max_M_y_top)
        elif self.concrete.design_code == "EN 1992-2004":
            _design_flexure_EN_1992_2004(self, max_M_y_bot, max_M_y_top)

        # Check flexural capacity for all forces with the assigned reinforcement
        all_results = self.check_flexure(forces)
        return all_results

    @property
    def flexure_checks(self) -> Tuple[FlexureCheck, ...]:
        """One immutable flexure result per combination of the last check.

        The beam's own attributes describe whichever combination ran last;
        these describe all of them, so enveloping is a plain ``max`` over the
        tuple rather than something the check has to accumulate as it goes.
        """
        return tuple(self._flexure_checks)

    @property
    def shear_checks(self) -> Tuple[ShearCheck, ...]:
        """One immutable shear result per combination of the last check."""
        return tuple(self._shear_checks)

    def flexure_check_results(self, forces: list[Forces]) -> Tuple[FlexureCheck, ...]:
        """Check flexure and return one result per combination, building no report.

        The same numbers as :meth:`check_flexure`, without the report tables or
        the result DataFrame — which is most of what a check costs. This is the
        entry point for looping over many sections::

            worst = max(r.bottom.DCR for r in beam.flexure_check_results(combos))

        The section is still updated in place while the check runs, so a loop
        over many sections needs one section object per station, not one shared.
        """
        self._flexure_checks = []
        for force in forces:
            self._run_flexure_check(force, report=False)
            self._flexure_checks.append(capture_flexure_check(self, force.label))
        self._flexure_checked = True
        return tuple(self._flexure_checks)

    def shear_check_results(self, forces: list[Forces]) -> Tuple[ShearCheck, ...]:
        """Check shear and return one result per combination, building no report.

        The counterpart of :meth:`flexure_check_results`; see it for the
        trade-off and the caveat about sharing a section across stations.
        """
        self._shear_checks = []
        for force in forces:
            self._run_shear_check(force, report=False)
            self._shear_checks.append(capture_shear_check(self, force.label))
        self._shear_checked = True
        return tuple(self._shear_checks)

    def _run_flexure_check(self, force: Forces, *, report: bool) -> Optional[DataFrame]:
        """Compute one flexure combination, and build its report only if asked.

        The design code does the calculation; assembling the row and the detail
        tables is presentation and lives in :mod:`mento.report_tables`.
        """
        if self.concrete.design_code in ("ACI 318-19", "CIRSOC 201-25"):
            _check_flexure_ACI_318_19(self, force)
        else:
            _check_flexure_EN_1992_2004(self, force)
        return build_flexure_report(self, force) if report else None

    def _run_shear_check(self, force: Forces, *, report: bool) -> Optional[DataFrame]:
        """Compute one shear combination, and build its report only if asked."""
        if self.concrete.design_code in ("ACI 318-19", "CIRSOC 201-25"):
            _check_shear_ACI_318_19(self, force)
        else:
            _check_shear_EN_1992_2004(self, force)
        return build_shear_report(self, force) if report else None

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
        self._flexure_results_detailed_list: Dict[Any, Dict[str, Any]] = {}  # Store detailed results by force ID
        self._flexure_checks = []

        for force in forces:
            result = self._run_flexure_check(force, report=True)
            self._flexure_results_list.append(result)

            # Store detailed results for this force
            self._flexure_results_detailed_list[force.id] = {
                "forces": self._forces_flexure.copy(),
                "min_max": self._data_min_max_flexure.copy(),
                "flexure_capacity_top": self._flexure_capacity_top.copy(),
                "flexure_capacity_bot": self._flexure_capacity_bot.copy(),
                "checks_pass": self._all_flexure_checks_passed,
            }

            # Capture this combination's result as a value, before the next one
            # overwrites the attributes it just left on the beam.
            self._flexure_checks.append(capture_flexure_check(self, force.label))

            # Extract the DCR values for top and bottom from the results
            current_dcr_top = self._DCRb_top
            current_dcr_bot = self._DCRb_bot

            # Update top limiting case
            if current_dcr_top >= max_dcr_top:
                max_dcr_top = current_dcr_top
                limiting_case_top = result
                limiting_case_top_details = self._flexure_results_detailed_list[force.id]

            # Update bottom limiting case
            if current_dcr_bot >= max_dcr_bot:
                max_dcr_bot = current_dcr_bot
                limiting_case_bot = result
                limiting_case_bot_details = self._flexure_results_detailed_list[force.id]

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
            if self.concrete.design_code == "ACI 318-19" or self.concrete.design_code == "CIRSOC 201-25":
                _design_shear_ACI_318_19(self, force)
            elif self.concrete.design_code == "EN 1992-2004":
                _design_shear_EN_1992_2004(self, force)
            # Check if this result is the limiting case
            current_A_v_req = self._A_v_req
            if current_A_v_req >= max_A_v_req:
                max_A_v_req = current_A_v_req
                max_V_s_req = self._V_s_req

        # Step 2: Perform rebar design for the worst-case force
        section_rebar = Rebar(self)
        self.shear_design_results = section_rebar.transverse_rebar(max_A_v_req, max_V_s_req, self._alpha)
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

        # Update longitudinal rebar attributes
        self._update_longitudinal_rebar_attributes()

        # Step 3: Check shear adequacy for all forces using the designed rebar
        all_results = self.check_shear(forces)
        return all_results

    # Factory method to select the shear check method
    def check_shear(self, forces: list[Forces]) -> DataFrame:
        self._shear_results_list = []  # Store individual results for each force
        self._shear_results_detailed_list: Dict[Any, Dict[str, Any]] = {}  # Store detailed results by force ID
        max_dcr = 0  # Track the maximum DCR to identify the limiting case
        self._limiting_case_shear_details = None
        self._shear_checks = []

        for force in forces:
            result = self._run_shear_check(force, report=True)
            self._shear_results_list.append(result)

            self._shear_results_detailed_list[force.id] = {
                "forces": self._forces_shear.copy(),
                "shear_reinforcement": self._shear_reinforcement.copy(),
                "min_max": self._data_min_max_shear.copy(),
                "checks_pass": self._all_shear_checks_passed,
                "shear_concrete": self._shear_concrete.copy(),
            }

            # Capture this combination's result as a value, before the next one
            # overwrites the attributes it just left on the beam.
            self._shear_checks.append(capture_shear_check(self, force.label))

            # Check if this result is the limiting case
            current_dcr = result["DCR"][0]
            if current_dcr >= max_dcr:
                max_dcr = current_dcr
                self._limiting_case_shear = result
                self._limiting_case_shear_details = self._shear_results_detailed_list[force.id]

        # Compile all results into a single DataFrame
        all_data = pd.concat(self._shear_results_list, ignore_index=True)
        units_row = self._get_units_row_shear()
        all_results = pd.concat([units_row, all_data], ignore_index=True)

        # Identify the most limiting case by Demand-to-Capacity Ratio (DCR) or other criteria
        self.limiting_case_shear = all_data.loc[all_data["DCR"].idxmax()]  # Select row with highest DCR

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
                        "NEd": "kN",
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
        # TODO: Add imperial units row output
        if self.concrete.design_code == "ACI 318-19" or self.concrete.design_code == "CIRSOC 201-25":
            units_row = pd.DataFrame(
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
        elif self.concrete.design_code == "EN 1992-2004":
            units_row = pd.DataFrame(
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
                        "MEd": "kNm",
                        "MRd": "kNm",
                        "MEd≤MRd": "",
                        "DCR": "",
                    }
                ]
            )
        else:
            raise ValueError(f"Flexure design method not implemented for concrete type: {type(self.concrete).__name__}")  # noqa: E501

        return units_row

    ##########################################################
    # CHECK & DESIGN ALL
    ##########################################################

    def design(self, forces: list[Forces]) -> None:
        """
        Complete design: flexure + shear + flexure check.
        """
        self.design_flexure(forces)
        self.design_shear(forces)
        self.check_flexure(forces)

    def check(self, forces: list[Forces]) -> None:
        """
        Complete check: flexure + shear.
        """
        self.check_flexure(forces)
        self.check_shear(forces)

    ##########################################################
    # RESULTS
    ##########################################################

    @property
    def reinforcement(self) -> SectionReinforcement:
        """The reinforcement this beam currently carries, as plain data.

        Readable at any time — it describes the section, not a result, so it
        does not need a check to have run::

            beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
            beam.reinforcement.bottom.A_s
            beam.reinforcement.transverse.n_legs

        For what a check *demanded* of that reinforcement (``A_s_req``,
        ``DCR``), use :attr:`flexure_design` and :attr:`shear_design`.
        """
        return build_reinforcement(self)

    @property
    def flexure_design(self) -> FlexureDesign:
        """Longitudinal reinforcement of this beam, as plain data.

        Available after a flexure check or design has been run. Unlike
        :attr:`flexure_results`, which prints a table for a notebook, this
        returns an object you can read values from::

            beam.flexure_design.bottom.A_s
            beam.flexure_design.bottom.layers[0].n

        Raises:
            DesignNotRunError: if no flexure check or design has been run.
        """
        return build_flexure_design(self)

    @property
    def shear_design(self) -> ShearDesign:
        """Transverse reinforcement of this beam, as plain data.

        Available after a shear check or design has been run::

            beam.shear_design.s_l
            beam.shear_design.A_v

        Raises:
            DesignNotRunError: if no shear check or design has been run.
        """
        return build_shear_design(self)

    # Beam results for Jupyter Notebook
    @property
    def data(self) -> None:
        type = self.mode.capitalize()
        markdown_content = (
            f"{type} {self.label}, $b$={self.width.to('cm')}"
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
                "Flexural design has not been performed yet. Call _check_flexure() or design_flexure() first.",
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
        warning_top = "⚠️ Some checks failed, see detailed results." if not checks_pass_top else ""
        warning_bot = "⚠️ Some checks failed, see detailed results." if not checks_pass_bot else ""

        # Pending for approval
        # all_checks_top = top_details.get('flexure_check')
        # all_checks_bot = bot_details.get('flexure_check')
        # all_checks = all_checks_top and all_checks_bot

        # Formatter instance for DCR formatting
        markdown_content = ""
        formatter = Formatter()
        symbols = self._flexure_symbols
        md_demand = symbols["md_demand"]
        md_capacity = symbols["md_capacity"]

        # Handle top result data
        if top_result_data:
            top_rebar_1 = top_result_data["Value"][0]
            top_rebar_2 = top_result_data["Value"][1]
            area_top = top_result_data["Value"][7]
            Mu_top = self._limiting_case_flexure_top_details["forces"]["Value"][0]
            Mn_top = top_result_data["Value"][9]
            DCR_top = top_result_data["Value"][10]

            rebar_top = f"{top_rebar_1}" + (f" ++ {top_rebar_2}" if top_rebar_2 != "-" else "")
            formatted_DCR_top = formatter.DCR(DCR_top)

            markdown_content += (
                f"Top longitudinal rebar: {rebar_top}, $A_{{s,top}}$ = {area_top} cm², "
                f"${md_demand}$ = {Mu_top} kNm, "
                f"${md_capacity}$ = {Mn_top} kNm → {formatted_DCR_top} {warning_top}\n\n"
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

            rebar_bot = f"{bot_rebar_1}" + (f" ++ {bot_rebar_2}" if bot_rebar_2 != "-" else "")
            formatted_DCR_bot = formatter.DCR(DCR_bot)

            markdown_content += (
                f"Bottom longitudinal rebar: {rebar_bot}, $A_{{s,bot}}$ = {area_bot} cm², "
                f"${md_demand}$ = {Mu_bot} kNm, "
                f"${md_capacity}$ = {Mn_bot} kNm → {formatted_DCR_bot} {warning_bot}"
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
                "Shear design has not been performed yet. Call check_shear() or design_shear() first.",
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
            warning = "⚠️ Some checks failed, see detailed results." if not checks_pass else ""
            if self.concrete.design_code == "ACI 318-19" or self.concrete.design_code == "CIRSOC 201-25":
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

    @property
    def _flexure_symbols(self) -> Dict[str, str]:
        """Flexural symbols of the active design code.

        ACI 318-19 and CIRSOC 201-25 work with a factored demand :math:`M_u` and a
        reduced capacity :math:`\\phi M_n`; EN 1992-2004 works with :math:`M_{Ed}`
        and :math:`M_{Rd}`. The per-code result dictionaries already carry the right
        wording, but the limiting-case summaries assembled here have to pick it.

        ``md_*`` entries are LaTeX for the Jupyter markdown line; the others are the
        plain variable names used in the detailed tables and Word reports.
        """
        if self.concrete.design_code == "EN 1992-2004":
            return {
                "demand": "MEd",
                "demand_top": "MEd,top",
                "demand_bot": "MEd,bot",
                "md_demand": "M_{Ed}",
                "md_capacity": "M_{Rd}",
            }
        return {
            "demand": "Mu",
            "demand_top": "Mu,top",
            "demand_bot": "Mu,bot",
            "md_demand": "M_u",
            "md_capacity": "\\phi M_n",
        }

    @property
    def _report_text(self) -> Dict[str, str]:
        """Report wording of this element, so a slab is not reported as a beam.

        ``ShearWall`` overrides the detailed report methods with its own wording,
        so only the beam and the slab are covered here. The strings double as the
        keys of the language catalog in :mod:`mento.i18n`.
        """
        if self.mode == "slab":
            return {
                "flexure_banner": "===== SLAB FLEXURE DETAILED RESULTS =====",
                "shear_banner": "===== SLAB SHEAR DETAILED RESULTS =====",
                "flexure_doc_title": "Concrete slab flexure check",
                "shear_doc_title": "Concrete slab shear check",
                "flexure_heading": "Slab {label} flexure check",
                "shear_heading": "Slab {label} shear check",
            }
        return {
            "flexure_banner": "===== BEAM FLEXURE DETAILED RESULTS =====",
            "shear_banner": "===== BEAM SHEAR DETAILED RESULTS =====",
            "flexure_doc_title": "Concrete beam flexure check",
            "shear_doc_title": "Concrete beam shear check",
            "flexure_heading": "Beam {label} flexure check",
            "shear_heading": "Beam {label} shear check",
        }

    def _report_file_name(self, heading_key: str) -> str:
        """Name of the Word file of a report.

        Built from the English heading whatever the report language is, so a
        project keeps one naming scheme.
        """
        heading = self._report_text[heading_key].format(label=self.label)
        return f"{heading} {self.concrete.design_code}.docx"

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
                "Flexural check has not been performed yet. Call _check_flexure or design_flexure first.",
                UserWarning,
            )
            self._md_flexure_results = "Flexure results are not available."
            return None

        # Determine which results to display (limiting case by default)
        if force:
            if force.id not in self._flexure_results_detailed_list:
                raise ValueError(f"No results found for Forces object with ID {force.id}.")
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
            top_result_data = self._limiting_case_flexure_top_details["flexure_capacity_top"]
            bot_result_data = self._limiting_case_flexure_bot_details["flexure_capacity_bot"]
            forces_result = {
                "Design forces": [
                    "Top max moment",
                    "Bottom max moment",
                ],
                "Variable": [
                    self._flexure_symbols["demand_top"],
                    self._flexure_symbols["demand_bot"],
                ],
                "Value": [
                    round(self._limiting_case_flexure_top_details["forces"]["Value"][0], 2),
                    round(self._limiting_case_flexure_bot_details["forces"]["Value"][1], 2),
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
                    round(self._limiting_case_flexure_top_details["min_max"]["Min."][0], 2),  # Top limiting case As_min
                    self._limiting_case_flexure_top_details["min_max"]["Min."][1],  # Top limiting case spacing min
                    round(
                        self._limiting_case_flexure_bot_details["min_max"]["Min."][2], 2
                    ),  # Bottom limiting case As_min
                    self._limiting_case_flexure_bot_details["min_max"]["Min."][3],  # Bottom limiting case spacing min
                ],
                "Max.": [
                    round(self._limiting_case_flexure_top_details["min_max"]["Max."][0], 2),  # Top limiting case As_max
                    "",  # No max constraint for spacing
                    round(
                        self._limiting_case_flexure_bot_details["min_max"]["Max."][2], 2
                    ),  # Bottom limiting case As_max
                    "",  # No max constraint for spacing
                ],
                "Ok?": [
                    self._limiting_case_flexure_top_details["min_max"]["Ok?"][0],  # Top As check
                    self._limiting_case_flexure_top_details["min_max"]["Ok?"][1],  # Top spacing check
                    self._limiting_case_flexure_bot_details["min_max"]["Ok?"][2],  # Bottom As check
                    self._limiting_case_flexure_bot_details["min_max"]["Ok?"][3],  # Bottom spacing check
                ],
            }

        # Create TablePrinter instances for detailed display
        language = get_language()
        print(translate(self._report_text["flexure_banner"], language))
        materials_printer = TablePrinter("MATERIALS", language)
        materials_printer.print_table_data(self._materials_flexure, headers="keys")

        geometry_printer = TablePrinter("GEOMETRY", language)
        geometry_printer.print_table_data(self._geometry_flexure, headers="keys")

        forces_printer = TablePrinter("FORCES", language)
        forces_printer.print_table_data(forces_result, headers="keys")

        min_max_printer = TablePrinter("MAX AND MIN LIMIT CHECKS", language)
        min_max_printer.print_table_data(min_max_result, headers="keys")

        capacity_printer = TablePrinter("FLEXURAL CAPACITY - TOP", language)
        capacity_printer.print_table_data(top_result_data, headers="keys")
        capacity_printer = TablePrinter("FLEXURAL CAPACITY - BOTTOM", language)
        capacity_printer.print_table_data(bot_result_data, headers="keys")

    def flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        """Write the detailed flexure results to a Word document.

        The assembly lives in :mod:`mento.documents`.
        """
        flexure_report_doc(self, force)

    def shear_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        """Write the detailed shear results to a Word document.

        The assembly lives in :mod:`mento.documents`.
        """
        shear_report_doc(self, force)

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
                "Shear check has not been performed yet. Call check_shear or design_shear first.",
                UserWarning,
            )
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
        language = get_language()
        print(translate(self._report_text["shear_banner"], language))
        materials_printer = TablePrinter("MATERIALS", language)
        materials_printer.print_table_data(self._materials_shear, headers="keys")
        geometry_printer = TablePrinter("GEOMETRY", language)
        geometry_printer.print_table_data(self._geometry_shear, headers="keys")
        forces_printer = TablePrinter("FORCES", language)
        forces_printer.print_table_data(result_data["forces"], headers="keys")
        steel_printer = TablePrinter("SHEAR STRENGTH", language)
        steel_printer.print_table_data(result_data["shear_reinforcement"], headers="keys")
        min_max_printer = TablePrinter("MAX AND MIN LIMIT CHECKS", language)
        min_max_printer.print_table_data(result_data["min_max"], headers="keys")
        concrete_printer = TablePrinter("CONCRETE STRENGTH", language)
        concrete_printer.print_table_data(result_data["shear_concrete"], headers="keys")

    def _format_longitudinal_rebar_string(self, n1: int, d_b1: Quantity, n2: int = 0, d_b2: Quantity = 0 * mm) -> str:
        """
        Returns a formatted string representing the rebars and their diameters.

        Rules:
        - If n1 and n2 have the same diameter → combine (e.g., 2Ø16 + 1Ø16 → 3Ø16)
        - If they differ → show both groups (e.g., 2Ø16+2Ø10)
        - If no bars exist → "-"
        """
        # Convert diameters safely
        phi1 = int(d_b1.to("mm").magnitude) if (d_b1 is not None and d_b1.magnitude > 0) else 0
        phi2 = int(d_b2.to("mm").magnitude) if (d_b2 is not None and d_b2.magnitude > 0) else 0

        # No bars at all
        if n1 == 0 and n2 == 0:
            return "-"

        # Only one group
        if n2 == 0 or phi2 == 0:
            return f"{n1}Ø{phi1}"
        if n1 == 0 or phi1 == 0:
            return f"{n2}Ø{phi2}"

        # Same diameter → combine quantities
        if phi1 == phi2:
            return f"{n1 + n2}Ø{phi1}"

        # Different diameters → write both
        return f"{n1}Ø{phi1}+{n2}Ø{phi2}"

    ##########################################################
    # PLOT BEAM SECTION WITH REBAR
    ##########################################################

    def plot(self, show: bool = False) -> "Figure":
        """Draw the cross-section with its reinforcement.

        The drawing itself lives in :mod:`mento.plots`; keeping it out of this
        class is what stops the element from being part matplotlib.
        """
        return plot_beam_section(self, show=show)
