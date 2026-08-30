"""Notebook views and Word reports for a shear wall.

The wall's counterpart to :mod:`mento.reports.views` and
:mod:`mento.reports.documents`. Both media live in one module here because the
wall's presentation is small enough that splitting it would cost more than it
explains.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import pandas as pd
from IPython.display import Markdown, display

from mento._version import __version__ as MENTO_VERSION
from mento.i18n import get_language, translate
from mento.results import DocumentBuilder, Formatter, TablePrinter

if TYPE_CHECKING:
    from mento.forces import Forces
    from mento.shear_wall import ShearWall


def _show(markdown: str) -> None:
    """Render Markdown in a notebook; IPython ships no type information."""
    display(Markdown(markdown))  # type: ignore[no-untyped-call]


def wall_flexure_results_detailed_doc(self: "ShearWall", force: Optional[Forces] = None) -> None:
    raise NotImplementedError("Flexure results are not implemented for ShearWall (Phase 0).")


def wall_data(self: "ShearWall") -> None:
    """Wall basic info as Markdown (length, thickness, story height, materials)."""
    level_str = f"Level {self.level}, " if self.level else ""
    markdown_content = (
        f"{level_str}Shear Wall {self.label}, "
        f"$l_w$={self.length.to('cm')}, "
        f"$t$={self.thickness.to('cm')}, "
        f"$h_w$={self.height.to('cm')}, "
        f"$c_c$={self.c_c.to('cm')}, "
        f"Concrete {self.concrete.name}, Rebar {self.steel_bar.name}."
    )
    self._md_data = markdown_content
    _show(markdown_content)
    return None


def wall_shear_results(self: "ShearWall") -> None:
    if not self._shear_wall_checked:
        self._md_shear_results = "Shear results are not available."
        return None

    details = self._limiting_case_shear_details or {}
    capacity = details.get("shear_capacity", {})
    min_max = details.get("min_max", {})
    forces_dict = details.get("forces", {})
    if not capacity:
        self._md_shear_results = "No shear to check."
        _show(self._md_shear_results)
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
    _show(markdown_content)
    return None


def wall_shear_results_detailed(self: "ShearWall", force: Optional[Forces] = None) -> None:
    if not self._shear_wall_checked:
        self._md_shear_results = "Shear results are not available."
        return None
    if force:
        if force.id not in self._shear_results_detailed_list:
            raise ValueError(f"No results found for Forces object with ID {force.id}.")
        result_data = self._shear_results_detailed_list[force.id]
    else:
        result_data = self._limiting_case_shear_details

    language = get_language()
    print(translate("===== SHEAR WALL DETAILED RESULTS =====", language))
    TablePrinter("MATERIALS", language).print_table_data(self._materials_shear_wall, headers="keys")
    TablePrinter("GEOMETRY", language).print_table_data(self._geometry_shear_wall, headers="keys")
    TablePrinter("FORCES", language).print_table_data(result_data["forces"], headers="keys")
    TablePrinter("MAX AND MIN LIMIT CHECKS", language).print_table_data(result_data["min_max"], headers="keys")
    TablePrinter("SHEAR STRENGTH", language).print_table_data(result_data["shear_capacity"], headers="keys")


def wall_shear_results_detailed_doc(self: "ShearWall", force: Optional[Forces] = None) -> None:
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

    doc_builder = DocumentBuilder(title="Concrete shear wall check", language=get_language())
    doc_builder.add_heading("Shear Wall {label} shear check", level=1, label=self.label)
    doc_builder.add_text(
        "Made with mento {version}. Design code: {design_code}",
        version=MENTO_VERSION,
        design_code=self.concrete.design_code,
    )
    doc_builder.add_heading("Materials", level=2)
    doc_builder.add_table_data(df_materials)
    doc_builder.add_table_data(df_geometry)
    doc_builder.add_table_data(df_forces)
    doc_builder.add_heading("Limit checks", level=2)
    doc_builder.add_table_min_max(df_min_max)
    doc_builder.add_heading("Design checks", level=2)
    doc_builder.add_table_dcr(df_capacity)
    doc_builder.save(f"Shear Wall {self.label} shear check {self.concrete.design_code}.docx")
