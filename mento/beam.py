from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, Optional, Dict, Tuple

if TYPE_CHECKING:
    from matplotlib.figure import Figure
from pint import Quantity
import numpy as np
import pandas as pd
from pandas import DataFrame
import math
from numbers import Integral
# from devtools import debug

from mento.rectangular import RectangularSection
from mento.codes.registry import design_code
from mento.precompute import refresh_section_floats
from mento.rebar import Rebar
from mento.units import mm, inch, kN, m, cm, dimensionless
from mento.forces import Forces
from mento.settings import BeamSettings
from mento.reports import views
from mento.reports.documents import flexure_report_doc, shear_report_doc
from mento.plots.sections import plot_beam_section
from mento.reports.tables import build_flexure_report, build_shear_report
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


class _DesignCodeAttributes:
    """Attributes the design codes write onto a section, declared for type checkers.

    Kept off ``RectangularBeam`` itself: that class is a dataclass, so an
    annotation in its body -- even one guarded by ``TYPE_CHECKING`` -- counts as a
    field to a type checker, and every private name below showed up in the editor
    as a required constructor argument. A plain base class contributes no fields,
    so the synthesized signature stays the handful the caller actually passes.
    """

    # ------------------------------------------------------------------
    # The ADR-0001 compatibility layer: declared, never assigned here.
    #
    # A check no longer writes its results to the section. These exist because
    # the report tables still read them off the element, and each design code
    # zeroes its own set through ``DesignCode.initialize_attributes``, so a
    # table asked for before any check has run finds a number instead of an
    # AttributeError.
    #
    # Declared under TYPE_CHECKING so the type checker knows they exist and
    # what they hold, without adding anything at runtime.
    #
    # They leave with the compatibility layer. A new design code should not add
    # to this list -- its report tables should read the check state instead.
    # ------------------------------------------------------------------
    if TYPE_CHECKING:
        # ACI 318-19 and CIRSOC 201-25
        _phi_V_n: Quantity
        _phi_V_s: Quantity
        _phi_V_c: Quantity
        _phi_V_max: Quantity
        _V_u: Quantity
        _M_u: Quantity
        _M_u_bot: Quantity
        _M_u_top: Quantity
        _N_u: Quantity
        _A_cv: Quantity
        _k_c_min: Quantity
        _sigma_Nu: Quantity
        V_c: Quantity
        _rho_w: Quantity
        f_yt: Quantity
        _phi_M_n_bot: Quantity
        _phi_M_n_top: Quantity
        _d_b_max_bot: Quantity
        _d_b_max_top: Quantity
        _lambda_s: float
        _max_shear_ok: bool
        _A_s_bool_bot: bool
        _A_s_bool_top: bool
        # EN 1992-2004
        _V_Ed_1: Quantity
        _V_Ed_2: Quantity
        _N_Ed: Quantity
        _M_Ed: Quantity
        _M_Ed_bot: Quantity
        _M_Ed_top: Quantity
        _M_Rd_bot: Quantity
        _M_Rd_top: Quantity
        _sigma_cd: Quantity
        _sigma_cp: Quantity
        _V_Rd_c: Quantity
        _V_Rd_s: Quantity
        _V_Rd_max: Quantity
        _V_Rd: Quantity
        _f_ywk: Quantity
        _f_ywd: Quantity
        _f_cd: Quantity
        _f_cd_shear: Quantity
        _A_p: Quantity
        _z: Quantity
        _k_value: float
        _theta: float
        _cot_theta: float
        # Both
        _A_s_min_bot: Quantity
        _A_s_min_top: Quantity
        _A_s_max_bot: Quantity
        _A_s_max_top: Quantity
        flexure_design_results_bot: Any
        flexure_design_results_top: Any


@dataclass
class RectangularBeam(RectangularSection, _DesignCodeAttributes):
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

        # Rendered Markdown of the last view asked for, kept so a caller can
        # read the text back instead of only seeing it displayed.
        self._md_data: str = ""
        self._md_flexure_results: str = ""
        self._md_shear_results: str = ""

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
        design_code(self.concrete).initialize_attributes(self)

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

    def _finalize_longitudinal_design(
        self,
        A_req_bot: Quantity,
        A_req_top: Quantity,
        resists: Callable[[], bool],
    ) -> None:
        """Last word on the layout, once both faces have been designed.

        Nothing to do for a beam: its two faces are detailed independently, so
        the layout the design arrived at is the layout it keeps. It exists for
        the elements whose faces are tied to each other -- a footing is placed
        as one mat, and chooses the module and both bars here.

        Args:
            A_req_bot: Steel the bottom face has to carry, minimum included.
            A_req_top: The same for the top.
            resists: Whether the layout now on the section resists both design
                moments. An element that re-selects bars has to ask, since a
                thicker bar sits deeper and lowers the effective depth it was
                chosen against.
        """

    def _clear_top_longitudinal(self) -> None:
        """Reset the top reinforcement to the default placeholder bars."""
        if self.concrete.unit_system == "metric":
            self.set_longitudinal_rebar_top(0, 0 * mm)
        else:
            self.set_longitudinal_rebar_top(0, 0 * inch)

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

    def _leg_spacing_across_width(self) -> Quantity:
        """Distance between adjacent stirrup legs, measured across the section.

        A beam's legs are the sides of its stirrup cage, so their spacing is not
        detailed directly: it falls out of how many legs are spread between the
        centres of the outermost pair. A slab's is a detailing decision of its
        own, which is why the two sections answer this differently and the shear
        checks ask the section instead of assuming the cage.

        Zero when the section carries no stirrups, so that a section without
        shear reinforcement reports no spacing rather than a negative one.
        """
        n_legs = self._stirrup_n * 2
        if n_legs < 2:
            return 0 * self._len_unit()
        return max(
            (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs - 1),
            0 * self._len_unit(),
        )

    def _apply_transverse_design(self, design: Any) -> None:
        """Store the row the transverse rebar search chose.

        A hook rather than a call inlined in :meth:`design_shear`, because a
        section is detailed the way it is reinforced: a beam by a stirrup count
        and a spacing, a slab by a spacing in each direction. Only the section
        knows which of the two the row means.
        """
        self.set_transverse_rebar(design["n_stir"], design["d_b"], design["s_l"])

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
            self._bot_rebar_centroid = 0 * mm
        else:
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
            self._top_rebar_centroid = 0 * mm
        else:
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
        # Every geometry and reinforcement change funnels through here, so this
        # is the one place the float view can go stale. Rebuilt eagerly so that
        # no check ever has to create it -- see mento.precompute.
        refresh_section_floats(self)

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
        design_code(self.concrete).design_flexure(self, max_M_y_bot, max_M_y_top)

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

        Nothing is written to the section on the ACI and CIRSOC path: the check
        returns its result and only the reporting path copies it back.
        """
        self._flexure_checks = []
        for force in forces:
            state = self._run_flexure_check(force, report=False)
            self._flexure_checks.append(capture_flexure_check(self, force.label, state))
        self._flexure_checked = True
        return tuple(self._flexure_checks)

    def shear_check_results(self, forces: list[Forces]) -> Tuple[ShearCheck, ...]:
        """Check shear and return one result per combination, building no report.

        Nothing is written to the section: the ACI check returns its result as a
        value and only the reporting path copies it back. That is what makes a
        loop over many stations safe — see :meth:`flexure_check_results` for the
        flexure counterpart, which still writes.
        """
        self._shear_checks = []
        for force in forces:
            state = self._run_shear_check(force, report=False)
            self._shear_checks.append(capture_shear_check(self, force.label, state))
        self._shear_checked = True
        return tuple(self._shear_checks)

    def _run_flexure_check(self, force: Forces, *, report: bool) -> Optional[Any]:
        """Compute one flexure combination, and build its report only if asked.

        The design code does the calculation and returns it; assembling the row
        and the detail tables is presentation and lives in
        :mod:`mento.reports.tables`, which still reads the element — so the
        state is copied back only when a report is wanted.
        """
        code = design_code(self.concrete)
        state: Any = code.check_flexure(self, force)
        if report:
            code.apply_flexure_state(self, state)
        if report:
            self._flexure_report_row = build_flexure_report(self, force)
        return state

    def _run_shear_check(self, force: Forces, *, report: bool) -> Optional[Any]:
        """Compute one shear combination, and build its report only if asked.

        Returns the ACI state so a values-only caller never needs the section to
        have been written to. Building a report does need that, because the
        tables still read the element — the compatibility layer of ADR-0001.
        """
        code = design_code(self.concrete)
        state = code.check_shear(self, force)
        if report:
            code.apply_shear_state(self, state)
        if report:
            self._shear_report_row = build_shear_report(self, force)
        return state

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
            state = self._run_flexure_check(force, report=True)
            result = self._flexure_report_row
            self._flexure_results_list.append(result)

            # Store detailed results for this force
            self._flexure_results_detailed_list[force.id] = {
                "forces": self._forces_flexure.copy(),
                "min_max": self._data_min_max_flexure.copy(),
                "flexure_capacity_top": self._flexure_capacity_top.copy(),
                "flexure_capacity_bot": self._flexure_capacity_bot.copy(),
                "checks_pass": self._all_flexure_checks_passed,
            }

            # The result is a value of the check itself, not a reading of the
            # attributes it left on the beam -- those describe the last
            # combination only, and are on their way out with them.
            self._flexure_checks.append(capture_flexure_check(self, force.label, state))

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

    def _clear_transverse_rebar_design(self) -> None:
        """Record a shear design that places no stirrups.

        Mirrors what the stirrup designer leaves behind for a real cage, so the
        report tables and ``BeamSummary`` read the same keys either way -- they
        just read zeros.
        """
        length_unit = "cm" if self.concrete.unit_system == "metric" else "inch"
        area_unit = "cm**2/m" if self.concrete.unit_system == "metric" else "inch**2/ft"
        empty = {
            "n_stir": 0,
            "d_b": (0 * mm).to(length_unit),
            "s_l": (0 * mm).to(length_unit),
            "s_w": (0 * mm).to(length_unit),
            "A_v": (0 * cm**2 / m).to(area_unit),
            "s_max_l": (0 * mm).to(length_unit),
            "s_max_w": (0 * mm).to(length_unit),
        }
        self.shear_design_results = pd.DataFrame([empty])
        self._best_rebar_design = self.shear_design_results.iloc[0]
        self._stirrup_s_w = empty["s_w"]
        self._stirrup_s_max_l = empty["s_max_l"]
        self._stirrup_s_max_w = empty["s_max_w"]
        self.set_transverse_rebar(0, 0 * mm, 0 * cm)

    @property
    def _stirrups_optional(self) -> bool:
        """Whether this element may be built with no transverse reinforcement at all.

        A beam never may. Stirrups run the whole length of a beam whatever the
        shear diagram says, so designing one always places at least the minimum
        area -- ACI 318-19 Table 9.6.3.4, EN 1992-1-1 Eq. (9.5N) -- even where
        the concrete alone would carry Vu. ACI 9.6.3.1 waives that minimum below
        a threshold and the check reports the waiver, but the design does not
        follow it down.

        A one-way slab is the exception both codes name: ACI 318-19 7.6.3.1 only
        calls for shear reinforcement where Vu exceeds phi*Vc, and EN 1992-1-1
        6.2.1(4) drops the minimum in members where the load can redistribute
        transversally. So a slab the concrete alone carries is designed with no
        stirrups at all, rather than given a cage it does not need.
        """
        return self.mode == "slab"

    # Factory method to select the shear design method
    def design_shear(self, forces: list[Forces]) -> DataFrame:
        # Track the maximum A_v_req to identify the limiting case
        max_A_v_req = 0 * cm**2 / m
        max_V_s_req = 0 * kN

        # Step 1: Identify the worst-case force
        for force in forces:
            design_code(self.concrete).design_shear(self, force)
            # Check if this result is the limiting case
            current_A_v_req = self._A_v_req
            if current_A_v_req >= max_A_v_req:
                max_A_v_req = current_A_v_req
                max_V_s_req = self._V_s_req

        # Step 2: Perform rebar design for the worst-case force
        if self._stirrups_optional and max_A_v_req <= 0 * cm**2 / m:
            # No combination asks for shear reinforcement and the code does not
            # impose a minimum here, so the section is built without stirrups.
            self._clear_transverse_rebar_design()
            self._update_longitudinal_rebar_attributes()
            return self.check_shear(forces)

        section_rebar = Rebar(self)
        self.shear_design_results = section_rebar.transverse_rebar(max_A_v_req, max_V_s_req, self._alpha)
        self._best_rebar_design = section_rebar.transverse_rebar_design
        self._stirrup_s_l = self._best_rebar_design["s_l"]
        self._stirrup_s_w = self._best_rebar_design["s_w"]
        self._stirrup_s_max_l = self._best_rebar_design["s_max_l"]
        self._stirrup_s_max_w = self._best_rebar_design["s_max_w"]
        self._apply_transverse_design(self._best_rebar_design)

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
            state = self._run_shear_check(force, report=True)
            result = self._shear_report_row
            self._shear_results_list.append(result)

            self._shear_results_detailed_list[force.id] = {
                "forces": self._forces_shear.copy(),
                "shear_reinforcement": self._shear_reinforcement.copy(),
                "min_max": self._data_min_max_shear.copy(),
                "checks_pass": self._all_shear_checks_passed,
                "shear_concrete": self._shear_concrete.copy(),
            }

            # As in check_flexure: the value of the check, not the beam's
            # attributes afterwards.
            self._shear_checks.append(capture_shear_check(self, force.label, state))

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
        """The unit row of the shear summary, in the active code's own names."""
        return pd.DataFrame([design_code(self.concrete).units_row_shear])

    def _get_units_row_flexure(self) -> pd.DataFrame:
        """The unit row of the flexure summary, in the active code's own names."""
        # TODO: Add imperial units row output
        return pd.DataFrame([design_code(self.concrete).units_row_flexure])

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

    ##########################################################
    # PRESENTATION — every one of these delegates
    ##########################################################
    # Rendering a result is a presentation choice, not part of the beam, so the
    # notebook views live in `mento.reports.views`, the Word reports in
    # `mento.reports.documents` and the drawing in `mento.plots.sections`.

    @property
    def data(self) -> None:
        """Show the section's basic data as Markdown."""
        return views.data(self)

    @property
    def flexure_results(self) -> None:
        """Show a summary of the flexure results as Markdown."""
        return views.flexure_results(self)

    @property
    def shear_results(self) -> None:
        """Show a summary of the shear results as Markdown."""
        return views.shear_results(self)

    @property
    def results(self) -> None:
        """Show every available result as Markdown."""
        return views.results(self)

    def flexure_results_detailed(self, force: Optional[Forces] = None) -> None:
        """Print the detailed flexure tables for one combination."""
        return views.flexure_results_detailed(self, force)

    def shear_results_detailed(self, force: Optional[Forces] = None) -> None:
        """Print the detailed shear tables for one combination."""
        return views.shear_results_detailed(self, force)

    @property
    def _report_text(self) -> Dict[str, str]:
        """Translated headings shared by the views and the Word reports."""
        return views._report_text(self)

    def _report_file_name(self, heading_key: str) -> str:
        return views._report_file_name(self, heading_key)

    # Beam results for Jupyter Notebook

    @property
    def _flexure_symbols(self) -> Dict[str, str]:
        r"""Flexural symbols of the active design code.

        ACI 318-19 and CIRSOC 201-25 work with a factored demand :math:`M_u` and a
        reduced capacity :math:`\phi M_n`; EN 1992-2004 works with :math:`M_{Ed}`
        and :math:`M_{Rd}`. Each code declares its own in its registry entry.
        """
        return design_code(self.concrete).flexure_symbols

    def flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        """Write the detailed flexure results to a Word document.

        The assembly lives in :mod:`mento.reports.documents`.
        """
        flexure_report_doc(self, force)

    def shear_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        """Write the detailed shear results to a Word document.

        The assembly lives in :mod:`mento.reports.documents`.
        """
        shear_report_doc(self, force)

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

        The drawing itself lives in :mod:`mento.plots.sections`; keeping it out of this
        class is what stops the element from being part matplotlib.
        """
        return plot_beam_section(self, show=show)
