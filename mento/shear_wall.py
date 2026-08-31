from __future__ import annotations

import math
from typing import TYPE_CHECKING, Dict, Optional

if TYPE_CHECKING:
    from matplotlib.figure import Figure

import pandas as pd
from pandas import DataFrame
from pint import Quantity

from mento.beam import RectangularBeam
from mento.forces import Forces
from mento.material import Concrete, SteelBar
from mento.settings import BeamSettings
from mento.units import cm, dimensionless, kN, mm

from mento.codes.registry import design_code
from mento.plots.walls import plot_wall_elevation
from mento.reports import walls as wall_reports


class ShearWall(RectangularBeam):
    """
    Reinforced concrete structural wall — shear check and design.

    The design code is whatever the concrete declares; only codes whose
    registry entry supplies the wall hooks can check one.

    Geometry:
        thickness — wall thickness  (t)         [maps to parent's ``width``]
        length    — wall in-plane length  (lw)  [maps to parent's ``height``]
        height    — wall story height  (hw)     [exposed via property; replaces ``hw``]

    Reinforcement:
        Horizontal distributed bars resist in-plane shear (ρt).
        Vertical distributed bars provide minimum vertical steel (ρl).
        Use set_horizontal_rebar() and set_vertical_rebar() instead of stirrups.
    """

    def __init__(
        self,
        *,
        concrete: Concrete,
        steel_bar: SteelBar,
        c_c: Quantity,
        thickness: Quantity,
        length: Quantity,
        height: Quantity,
        label: Optional[str] = None,
        level: Optional[str] = None,
        settings: Optional[BeamSettings] = None,
    ) -> None:
        self.level: Optional[str] = level

        # Pre-seed `_length` and `_wall_height` so parent's `__post_init__` can read
        # `self.length`/`self.height` (our properties) during cross-section computation.
        # `_wall_height` is bootstrapped to `length` so parent's `_A_x = width * height`
        # produces the correct in-plane cross-section area; we replace it after super().
        self._length: Quantity = length
        self._wall_height: Quantity = length

        super().__init__(
            concrete=concrete,
            steel_bar=steel_bar,
            c_c=c_c,
            width=thickness,
            height=length,
            label=label,
            settings=settings,
        )

        # Replace bootstrap with the actual wall story height.
        self._wall_height = height
        # _initialize_wall_attributes was already called via __post_init__ above.

    # ------------------------------------------------------------------
    # Wall-friendly dimension properties
    # ------------------------------------------------------------------

    @property
    def thickness(self) -> Quantity:
        """Wall thickness — alias for the parent's ``width`` attribute."""
        return self.width

    @property
    def length(self) -> Quantity:
        """Wall in-plane length (the dimension that resists in-plane shear)."""
        return self._length

    @property
    def height(self) -> Quantity:  # type: ignore[override]
        """Wall story height — replaces the legacy ``hw`` field."""
        return self._wall_height

    @height.setter
    def height(self, value: Quantity) -> None:
        # Parent's dataclass-generated ``__init__`` assigns ``self.height = length``;
        # we accept that without complaint because ``_wall_height`` is bootstrapped to
        # the same value by our ``__init__``. Subsequent user assignments update the
        # wall story height directly.
        self._wall_height = value

    def __post_init__(self) -> None:
        super().__post_init__()
        self._initialize_wall_attributes()

    # ------------------------------------------------------------------
    # Wall-specific initialization
    # ------------------------------------------------------------------

    def _initialize_wall_attributes(self) -> None:
        self.mode = "shear_wall"

        # Distributed mesh is placed on BOTH faces of the wall (each face = E.F.).
        # The reinforcement ratio counts the bars from every curtain:
        #   ρ = n_curtains · A_b / (t · s)
        self._n_curtains: int = 2

        # Horizontal (transverse) distributed rebar
        self._d_b_h: Quantity = 0 * mm
        self._s_h: Quantity = 0 * mm
        self._rho_t: Quantity = 0 * dimensionless

        # Vertical (longitudinal) distributed rebar
        self._d_b_v: Quantity = 0 * mm
        self._s_v: Quantity = 0 * mm
        self._rho_l: Quantity = 0 * dimensionless

        # Wall shear result quantities
        self._Acv: Quantity = 0 * cm**2
        self._alpha_c: float = 0.0
        self._hw_lw: float = 0.0
        self._f_yt_wall: Quantity = 0 * mm / mm * kN / kN  # typed as dimensionless placeholder; set on first check
        self._V_u: Quantity = 0 * kN
        self._N_u: Quantity = 0 * kN
        self._V_c_wall: Quantity = 0 * kN
        self._V_s_wall: Quantity = 0 * kN
        self._V_n_wall: Quantity = 0 * kN
        self._V_n_max: Quantity = 0 * kN
        self._phi_V_n_wall: Quantity = 0 * kN
        self._phi_V_n_max_wall: Quantity = 0 * kN
        self._DCRv_wall: float = 0.0

        # Minimum ratios and spacing limits
        self._rho_t_min: Quantity = 0.0025 * dimensionless
        self._rho_l_min: Quantity = 0.0025 * dimensionless
        self._rho_t_req: Quantity = 0 * dimensionless
        self._s_h_max: Quantity = 0 * mm
        self._s_v_max: Quantity = 0 * mm

        # Status flags
        self._shear_wall_checked: bool = False
        self._all_wall_shear_checks_passed: bool = False

        # Detail dicts (mirrors beam pattern)
        self._materials_shear_wall: Dict = {}
        self._geometry_shear_wall: Dict = {}
        self._forces_shear_wall: Dict = {}
        self._shear_capacity_wall: Dict = {}
        self._data_min_max_wall: Dict = {}

    # ------------------------------------------------------------------
    # Rebar setters
    # ------------------------------------------------------------------

    def set_horizontal_rebar(self, d_b: Quantity, s: Quantity) -> None:
        """Set distributed horizontal (transverse) reinforcement.

        Bars are placed on each face (E.F.):  ρt = n_curtains · Ab / (t × s_h)
        A zero spacing means no rebar, which clears the horizontal reinforcement.
        """
        self._d_b_h = d_b
        self._s_h = s
        if s == 0 * mm:
            self._rho_t = 0 * dimensionless
        else:
            A_b = math.pi / 4 * d_b**2
            self._rho_t = (self._n_curtains * A_b / (self.thickness * s)).to("")

    def set_vertical_rebar(self, d_b: Quantity, s: Quantity) -> None:
        """Set distributed vertical reinforcement.

        Bars are placed on each face (E.F.):  ρl = n_curtains · Ab / (t × s_v)
        A zero spacing means no rebar, which clears the vertical reinforcement.
        """
        self._d_b_v = d_b
        self._s_v = s
        if s == 0 * mm:
            self._rho_l = 0 * dimensionless
        else:
            A_b = math.pi / 4 * d_b**2
            self._rho_l = (self._n_curtains * A_b / (self.thickness * s)).to("")

    # ------------------------------------------------------------------
    # Shear check and design (override RectangularBeam)
    # ------------------------------------------------------------------

    def check_shear(self, forces: list[Forces]) -> DataFrame:
        self._shear_results_list: list = []
        self._shear_results_detailed_list: Dict = {}
        max_dcr: float = 0.0
        self._limiting_case_shear_details = None

        for force in forces:
            code = design_code(self.concrete)
            state = code.requires("check_shear_wall")(self, force)
            # The report tables read the wall, so the state is applied here
            # and not on a values-only path.
            code.requires("apply_wall_shear_state")(self, state)
            result = wall_reports.build_wall_shear_report(self, force)

            self._shear_results_list.append(result)
            self._shear_results_detailed_list[force.id] = {
                "forces": self._forces_shear_wall.copy(),
                "shear_capacity": self._shear_capacity_wall.copy(),
                "min_max": self._data_min_max_wall.copy(),
                "checks_pass": self._all_wall_shear_checks_passed,
            }

            current_dcr = result["DCR"].iloc[0]
            if current_dcr >= max_dcr:
                max_dcr = current_dcr
                self._limiting_case_shear = result
                self._limiting_case_shear_details = self._shear_results_detailed_list[force.id]

        all_data = pd.concat(self._shear_results_list, ignore_index=True)
        units_row = self._get_units_row_shear_wall()
        all_results = pd.concat([units_row, all_data], ignore_index=True)
        self.limiting_case_shear = all_data.loc[all_data["DCR"].idxmax()]

        self._shear_wall_checked = True
        self._shear_checked = True
        return all_results

    def design_shear(self, forces: list[Forces]) -> DataFrame:
        """Design the horizontal (shear) mesh and the minimum vertical mesh.

        Selects a bar diameter + spacing for both directions against the
        worst-case force combination, applies them, and returns the
        re-evaluated check results.
        """
        if not forces:
            raise ValueError("design_shear requires at least one Forces object.")

        design_code(self.concrete).requires("design_shear_wall")(self, forces)

        # Re-run the check so the returned DataFrame / detail dicts reflect the mesh.
        return self.check_shear(forces)

    # ------------------------------------------------------------------
    # Units header row
    # ------------------------------------------------------------------

    def _get_units_row_shear_wall(self) -> pd.DataFrame:
        if self.concrete.unit_system == "metric":
            v_unit = "kN"
        else:
            v_unit = "kip"
        return pd.DataFrame(
            [
                {
                    "Label": "",
                    "Comb.": "",
                    "ρt,min": "",
                    "ρt,req": "",
                    "ρt": "",
                    "ρl,min": "",
                    "ρl": "",
                    "Vu": v_unit,
                    "ØVc": v_unit,
                    "ØVs": v_unit,
                    "ØVn": v_unit,
                    "ØVn,max": v_unit,
                    "Vu≤ØVn,max": "",
                    "Vu≤ØVn": "",
                    "DCR": "",
                }
            ]
        )

    # ------------------------------------------------------------------
    # Top-level check / Node integration
    # ------------------------------------------------------------------

    def check(self, forces: list[Forces]) -> None:
        """Complete check for a shear wall: shear only (no flexure in Phase 0)."""
        self.check_shear(forces)

    def design(self, forces: list[Forces]) -> None:
        """Complete design for a shear wall: shear only (no flexure in Phase 0)."""
        self.design_shear(forces)

    def check_flexure(self, forces: list[Forces]) -> DataFrame:  # type: ignore[override]
        raise NotImplementedError("Flexure check is not implemented for ShearWall (Phase 0).")

    def design_flexure(self, forces: list[Forces]) -> DataFrame:  # type: ignore[override]
        raise NotImplementedError("Flexure design is not implemented for ShearWall (Phase 0).")

    def flexure_results_detailed(self, force: Optional[Forces] = None) -> None:
        raise NotImplementedError("Flexure results are not implemented for ShearWall (Phase 0).")

    # ------------------------------------------------------------------
    # Presentation — every one of these delegates, as on RectangularBeam
    # ------------------------------------------------------------------

    def flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        return wall_reports.wall_flexure_results_detailed_doc(self, force)

    @property
    def data(self) -> None:
        """Show the wall's basic data as Markdown."""
        return wall_reports.wall_data(self)

    @property
    def shear_results(self) -> None:  # type: ignore[override]
        """Show a summary of the shear results as Markdown."""
        return wall_reports.wall_shear_results(self)

    def shear_results_detailed(self, force: Optional[Forces] = None) -> None:  # type: ignore[override]
        """Print the detailed shear tables."""
        return wall_reports.wall_shear_results_detailed(self, force)

    def shear_results_detailed_doc(self, force: Optional[Forces] = None) -> None:  # type: ignore[override]
        """Write the detailed shear results to a Word document."""
        return wall_reports.wall_shear_results_detailed_doc(self, force)

    def plot(self, show: bool = False) -> "Figure":  # type: ignore[override]
        """Draw the wall elevation with its reinforcement."""
        return plot_wall_elevation(self, show=show)

    @property
    def results(self) -> None:  # type: ignore[override]
        """Display wall data + shear results in Markdown (no flexure)."""
        self.data
        if self._shear_wall_checked:
            self.shear_results
        return None

    # ------------------------------------------------------------------
    # Detailed results
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Plot — wall plan view: length lw (horizontal) × thickness t
    # ------------------------------------------------------------------
