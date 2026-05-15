from __future__ import annotations

import math
from typing import Dict, Optional

import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import Markdown, display
from matplotlib.patches import Rectangle
from pandas import DataFrame
from pint import Quantity

from mento._version import __version__ as MENTO_VERSION
from mento.beam import RectangularBeam
from mento.forces import Forces
from mento.material import Concrete, SteelBar
from mento.results import CUSTOM_COLORS, DocumentBuilder, Formatter, TablePrinter
from mento.settings import BeamSettings
from mento.units import cm, dimensionless, kN, mm

from mento.codes.ACI_318_19_wall import _check_shear_ACI_318_19_wall


class ShearWall(RectangularBeam):
    """
    Reinforced concrete structural wall — ACI 318-19 Section 11 shear check/design.

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
        settings: Optional[BeamSettings] = None,
    ) -> None:
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

        # Wall shear result quantities (ACI 318-19)
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
        """
        self._d_b_h = d_b
        self._s_h = s
        A_b = math.pi / 4 * d_b**2
        self._rho_t = (self._n_curtains * A_b / (self.thickness * s)).to("")

    def set_vertical_rebar(self, d_b: Quantity, s: Quantity) -> None:
        """Set distributed vertical reinforcement.

        Bars are placed on each face (E.F.):  ρl = n_curtains · Ab / (t × s_v)
        """
        self._d_b_v = d_b
        self._s_v = s
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
            if self.concrete.design_code == "ACI 318-19" or self.concrete.design_code == "CIRSOC 201-25":
                result = _check_shear_ACI_318_19_wall(self, force)
            else:
                raise NotImplementedError(
                    f"Shear wall check not implemented for design code: {self.concrete.design_code}"
                )

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

        code = self.concrete.design_code
        if code == "ACI 318-19" or code == "CIRSOC 201-25":
            # CIRSOC 201-25 reuses the ACI wall design; the code-specific bar
            # catalogue is selected inside _design_shear_ACI_318_19_wall.
            from mento.codes.ACI_318_19_wall import _design_shear_ACI_318_19_wall

            _design_shear_ACI_318_19_wall(self, forces)
        else:
            raise NotImplementedError(f"Shear wall design not implemented for design code: {code}")

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

    def flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        raise NotImplementedError("Flexure results are not implemented for ShearWall (Phase 0).")

    # ------------------------------------------------------------------
    # Jupyter / Markdown summary properties
    # ------------------------------------------------------------------

    @property
    def data(self) -> None:
        """Wall basic info as Markdown (length, thickness, story height, materials)."""
        markdown_content = (
            f"Shear Wall {self.label}, "
            f"$l_w$={self.length.to('cm')}, "
            f"$t$={self.thickness.to('cm')}, "
            f"$h_w$={self.height.to('cm')}, "
            f"$c_c$={self.c_c.to('cm')}, "
            f"Concrete {self.concrete.name}, Rebar {self.steel_bar.name}."
        )
        self._md_data = markdown_content
        display(Markdown(markdown_content))
        return None

    @property
    def shear_results(self) -> None:  # type: ignore[override]
        if not self._shear_wall_checked:
            self._md_shear_results = "Shear results are not available."
            return None

        details = self._limiting_case_shear_details or {}
        capacity = details.get("shear_capacity", {})
        min_max = details.get("min_max", {})
        forces_dict = details.get("forces", {})
        if not capacity:
            self._md_shear_results = "No shear to check."
            display(Markdown(self._md_shear_results))
            return None

        formatter = Formatter()
        dcr = capacity["Value"][-1]
        phi_Vn = capacity["Value"][2]
        Vu = forces_dict["Value"][0]
        rho_t = min_max["Value"][0]
        rho_l = min_max["Value"][1]
        checks_pass = details.get("checks_pass", False)
        warning = "⚠️ Some checks failed, see detailed results." if not checks_pass else ""

        rebar_h = (
            f"Ø{self._d_b_h.to('mm').magnitude:.0f}/{self._s_h.to('cm').magnitude:.0f} cm E.F."
            if self._s_h.magnitude > 0
            else "not assigned"
        )
        rebar_v = (
            f"Ø{self._d_b_v.to('mm').magnitude:.0f}/{self._s_v.to('cm').magnitude:.0f} cm E.F."
            if self._s_v.magnitude > 0
            else "not assigned"
        )

        markdown_content = (
            f"Horizontal rebar: {rebar_h}, $\\rho_t$={rho_t}, "
            f"Minimum vertical rebar: {rebar_v}, $\\rho_l$={rho_l}, "
            f"$V_u$={Vu} kN, $\\phi V_n$={phi_Vn} kN → "
            f"{formatter.DCR(dcr)} {warning}"
        )
        self._md_shear_results = markdown_content
        display(Markdown(markdown_content))
        return None

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

    def shear_results_detailed(self, force: Optional[Forces] = None) -> None:  # type: ignore[override]
        if not self._shear_wall_checked:
            self._md_shear_results = "Shear results are not available."
            return None
        if force:
            if force.id not in self._shear_results_detailed_list:
                raise ValueError(f"No results found for Forces object with ID {force.id}.")
            result_data = self._shear_results_detailed_list[force.id]
        else:
            result_data = self._limiting_case_shear_details

        print("===== SHEAR WALL DETAILED RESULTS =====")
        TablePrinter("MATERIALS").print_table_data(self._materials_shear_wall, headers="keys")
        TablePrinter("GEOMETRY").print_table_data(self._geometry_shear_wall, headers="keys")
        TablePrinter("FORCES").print_table_data(result_data["forces"], headers="keys")
        TablePrinter("MAX AND MIN LIMIT CHECKS").print_table_data(result_data["min_max"], headers="keys")
        TablePrinter("SHEAR STRENGTH").print_table_data(result_data["shear_capacity"], headers="keys")

    def shear_results_detailed_doc(self, force: Optional[Forces] = None) -> None:  # type: ignore[override]
        if not self._shear_wall_checked:
            self._md_shear_results = "Shear results are not available."
            return None
        if force:
            if force.id not in self._shear_results_detailed_list:
                raise ValueError(f"No results found for Forces object with ID {force.id}.")
            result_data = self._shear_results_detailed_list[force.id]
        else:
            result_data = self._limiting_case_shear_details

        df_materials = pd.DataFrame(self._materials_shear_wall)
        df_geometry = pd.DataFrame(self._geometry_shear_wall)
        df_forces = pd.DataFrame(result_data["forces"])
        df_min_max = pd.DataFrame(result_data["min_max"])
        df_capacity = pd.DataFrame(result_data["shear_capacity"])

        doc_builder = DocumentBuilder(title="Concrete shear wall check")
        doc_builder.add_heading(f"Shear Wall {self.label} shear check", level=1)
        doc_builder.add_text(f"Made with mento {MENTO_VERSION}. Design code: {self.concrete.design_code}")
        doc_builder.add_heading("Materials", level=2)
        doc_builder.add_table_data(df_materials)
        doc_builder.add_table_data(df_geometry)
        doc_builder.add_table_data(df_forces)
        doc_builder.add_heading("Limit checks", level=2)
        doc_builder.add_table_min_max(df_min_max)
        doc_builder.add_heading("Design checks", level=2)
        doc_builder.add_table_dcr(df_capacity)
        doc_builder.save(f"Shear Wall {self.label} shear check {self.concrete.design_code}.docx")

    # ------------------------------------------------------------------
    # Plot — wall plan view: length lw (horizontal) × thickness t
    # ------------------------------------------------------------------

    def plot(self, show: bool = False) -> plt.Figure:  # type: ignore[override]
        """Plan view of the wall (length × thickness) with dimensions and rebar callouts."""
        lw_cm: float = self.length.to("cm").magnitude
        t_cm: float = self.thickness.to("cm").magnitude

        fig, self._ax = plt.subplots(figsize=(10, max(3, t_cm / lw_cm * 8 + 1.5)))

        # Wall outline (plan view: length horizontal, thickness vertical)
        wall = Rectangle(
            (0, 0),
            lw_cm,
            t_cm,
            linewidth=1.3,
            edgecolor=CUSTOM_COLORS["dark_gray"],
            facecolor=CUSTOM_COLORS["light_gray"],
        )
        self._ax.add_patch(wall)

        # Decoupled horizontal/vertical padding so thin walls aren't squashed.
        h_pad = lw_cm * 0.08
        v_pad = max(t_cm * 1.2, lw_cm * 0.04)

        dim_color = CUSTOM_COLORS["dark_blue"]
        text_color = CUSTOM_COLORS["dark_gray"]

        # ---- Length dimension (below the wall) — engineering style ----
        dim_y = -v_pad * 0.6
        tick_v = v_pad * 0.18
        self._ax.annotate(
            "",
            xy=(0, dim_y),
            xytext=(lw_cm, dim_y),
            arrowprops={
                "arrowstyle": "<->",
                "lw": 1,
                "color": dim_color,
                "shrinkA": 0,
                "shrinkB": 0,
            },
        )
        # Extension ticks at both endpoints
        self._ax.plot([0, 0], [dim_y - tick_v, dim_y + tick_v], color=dim_color, lw=1)
        self._ax.plot([lw_cm, lw_cm], [dim_y - tick_v, dim_y + tick_v], color=dim_color, lw=1)
        self._ax.text(
            lw_cm / 2,
            dim_y - tick_v * 1.6,
            "{:.0f~P}".format(self.length.to("cm")),
            ha="center",
            va="top",
            color=text_color,
        )

        # ---- Thickness dimension (right of the wall) — engineering style ----
        x_dim = lw_cm + h_pad * 0.6
        tick_h = h_pad * 0.07
        self._ax.annotate(
            "",
            xy=(x_dim, 0),
            xytext=(x_dim, t_cm),
            arrowprops={
                "arrowstyle": "<->",
                "lw": 1,
                "color": dim_color,
                "shrinkA": 0,
                "shrinkB": 0,
            },
        )
        # Extension ticks at both endpoints
        self._ax.plot([x_dim - tick_h, x_dim + tick_h], [0, 0], color=dim_color, lw=1)
        self._ax.plot([x_dim - tick_h, x_dim + tick_h], [t_cm, t_cm], color=dim_color, lw=1)
        self._ax.text(
            x_dim + tick_h * 4,
            t_cm / 2,
            "{:.0f~P}".format(self.thickness.to("cm")),
            ha="left",
            va="center",
            color=text_color,
        )

        # ---- Rebar callout above the wall (same font as dimensions) ----
        def _fmt_rebar(d_b: Quantity, s: Quantity) -> str:
            if s.magnitude <= 0:
                return "not assigned"
            return f"Ø{d_b.to('mm').magnitude:.0f}/{s.to('cm').magnitude:.0f} cm E.F."

        rebar_text = (
            f"Horizontal rebar: {_fmt_rebar(self._d_b_h, self._s_h)}\n"
            f"Minimum vertical rebar: {_fmt_rebar(self._d_b_v, self._s_v)}"
        )
        self._ax.text(
            lw_cm / 2,
            t_cm + v_pad * 0.35,
            rebar_text,
            ha="center",
            va="bottom",
            color=text_color,
        )

        # Compact limits — no large empty band above or below.
        self._ax.set_xlim(-h_pad * 0.5, lw_cm + h_pad * 2.0)
        self._ax.set_ylim(dim_y - tick_v * 5, t_cm + v_pad * 1.8)
        self._ax.axis("off")
        self._ax.set_title(f"Shear Wall {self.label} — plan view")

        self._fig = fig
        if show:
            plt.show()
        # Prevent Jupyter's inline backend from auto-displaying the figure twice
        # (once when plt.subplots creates it, and again as the cell's return value).
        plt.close(fig)
        return fig
