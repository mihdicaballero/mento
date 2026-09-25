"""ACI 318-19, and CIRSOC 201-25 which reprints it under the same numbering.

The registry entries. CIRSOC 201-25 restates the ACI 318-19 provisions clause
by clause, so it reuses every ACI hook except where its printed table gives a
different number. What actually differs, with the clause on each side, is:

* ACI 318-19 Table 9.7.6.2.2 / CIRSOC 201-25 Tabla 9.7.6.2.2 -- the absolute
  caps on stirrup spacing: 600/300 mm (24/12 in.) against 400/200 mm. Hook
  ``stirrup_spacing_caps``.
* ACI 318-19 §9.7.6.4.2 / CIRSOC 201-25 §9.7.6.4.2 with Tabla 9.7.6.4.2 -- the
  smallest stirrup that may laterally support compression reinforcement: a
  No. 10 bar (No. 13 from No. 36 up) against 6, 8, 10 or 12 mm by longitudinal
  bar size. Hook ``min_stirrup_diameter``.
* ACI 318-19 §7.7.2.3 / CIRSOC 201-25 art. 7.7.2.3 -- largest spacing of the
  flexural bars of a one-way slab: the lesser of 3h and 450 mm (18 in.)
  against the lesser of 3h and 300 mm. Hook ``max_bar_spacing_slab``.
* The bar sizes the transverse selection draws from, which are those of
  CIRSOC 201-25 §20.2.1.3, Tabla 20.2.1. Hook ``transverse_rebar``.
* §9.6.1.2 -- the f_y cap for minimum flexural reinforcement: 550 MPa
  (80 ksi) against 500 MPa. Hook ``flexural_min_fy_cap``.
* §9.6.3.1 -- the beam minimum-shear-reinforcement threshold coefficient:
  0.083 (1.0 in psi) against 0.085. Hook ``min_shear_reinforcement_coefficient``.
* Eq. (11.5.4.4) -- the stress divisor for walls under net axial tension:
  3.45 MPa (500 psi) against 3.5 MPa. Hook ``wall_axial_tension_divisor``.

Two further differences are not reflected in the hooks, and are written down
here rather than left implied:

* ACI 318-19 Table 19.2.1.1 / CIRSOC 201-25 Tabla 19.2.1.1 -- smallest f'c the
  code is written for, 17 MPa (2500 psi) against 20 MPa. Neither is checked.
* ACI 318-19 §8.7.2.2 / CIRSOC 201-25 art. 8.7.2.2 -- the same spacing cap for
  a two-way solid slab: the lesser of 2h and 450 mm at the critical sections
  and of 3h and 450 mm elsewhere, against 300 mm in place of both 450 mm.
  Nothing in mento designs a two-way slab, so no hook reads it; it is what
  makes the 300 mm of a footing spanning two ways a clause rather than
  practice under CIRSOC, reached through art. 13.3.3.1.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from mento.codes.ACI_318_19_beam import (
    _check_flexure_ACI_318_19,
    _check_shear_ACI_318_19,
    _design_flexure_ACI_318_19,
    _design_shear_ACI_318_19,
    _flexure_ductile_ACI_318_19,
)
from mento.codes.ACI_318_19_punching import check_punching_ACI_318_19
from mento.codes.ACI_318_19_wall import _check_shear_ACI_318_19_wall, _design_shear_ACI_318_19_wall
from mento.codes.check_state import (
    apply_flexure_state,
    apply_shear_state,
    apply_wall_shear_state,
)
from mento.codes.registry import DesignCode, register
from mento.material import Concrete_ACI_318_19, Concrete_CIRSOC_201_25
from mento.units import cm, dimensionless, inch, kN, kNm, mm, MPa, psi

if TYPE_CHECKING:
    from mento.beam import RectangularBeam
    from mento.rebar import Rebar

#: ACI names the demand M_u and the capacity phi*Mn.
_FLEXURE_SYMBOLS = {
    "demand": "Mu",
    "demand_top": "Mu,top",
    "demand_bot": "Mu,bot",
    "md_demand": "M_u",
    "md_capacity": r"\phi M_n",
}

_UNITS_ROW_SHEAR = {
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

_UNITS_ROW_FLEXURE = {
    "Label": "",
    "Comb.": "",
    "Position": "",
    "As,min": "cm²",
    "As,req top": "cm²",
    "As,req bot": "cm²",
    "As": "cm²",
    "Mu": "kNm",
    "ØMn": "kNm",
    "Mu≤ØMn": "",
    "DCR": "",
}


#: The capacity ticks come out of the all-beams summaries. The DCR column
#: beside them already says whether the section is enough, and the capacities
#: themselves are right there for a reader who wants the margin. They stay in
#: `check_shear` and `check_flexure`, which is where a caller reads them
#: programmatically. EN drops the same three under its own names, so both
#: codes' summaries end up the same shape.
_SUMMARY_DROP_COLUMNS = ("Vu≤ØVn", "Vu≤ØVmax", "Mu≤ØMn")

_SHEAR_SYMBOLS = {
    "reinforcement": "A_v",
    "demand": "V_u",
    "capacity": r"\phi V_n",
    # Row the reduced capacity sits on in this code's shear detail table.
    "capacity_row": 7,
}


#: How far apart the bars of a member on the ground are detailed, whichever
#: code it is designed to. Neither code states either number for a footing: the
#: upper bound is what keeps the pressure from the soil spread across the
#: reinforcement rather than arching between distant bars, and the lower one is
#: EN practice, applied to both codes because nothing about it is particular to
#: EN. Under CIRSOC 201-25 the 300 mm stops being practice and becomes the
#: clause itself, art. 7.7.2.3 (and art. 8.7.2.2 where the footing spans two
#: ways, by art. 13.3.3.1) -- which is why the CIRSOC hook below reaches the
#: same 300 mm for a slab that is not on the ground at all.
_MAX_BAR_SPACING_ON_SOIL = 300 * mm
_MIN_BAR_SPACING_ON_SOIL = 100 * mm

#: The absolute term of the one-way slab spacing cap, which is the whole of the
#: difference between the two codes' §7.7.2.3: ACI 318-19 prints 450 mm in its
#: SI edition and 18 in. in the in-lb one, CIRSOC 201-25 prints 300 mm and has
#: no imperial edition to differ from.
_MAX_BAR_SPACING_SLAB_ACI = 450 * mm
_MAX_BAR_SPACING_SLAB_ACI_IMPERIAL = 18 * inch
_MAX_BAR_SPACING_SLAB_CIRSOC = 300 * mm


def _min_stirrup_diameter(concrete: Any) -> Any:
    """ACI's smallest stirrup: a #3 bar, or 10 mm in metric practice.

    ACI 318-19 §9.7.6.4.2 / CIRSOC 201-25 §9.7.6.4.2 only size the stirrups
    that laterally support compression reinforcement (§9.7.6.4.1 in both);
    for a stirrup placed purely for shear neither code states a minimum
    diameter, so this number is practice -- the smallest bar a cage is tied
    with -- and not the clause. CIRSOC 201-25 §9.7.6.4.2 differs: its Tabla
    9.7.6.4.2 grades the minimum with the longitudinal bar (6, 8, 10 or
    12 mm) where ACI 318-19 §9.7.6.4.2 gives No. 10, or No. 13 from No. 36
    and bundled bars up; see :func:`_min_stirrup_diameter_cirsoc`.
    """
    return 10 * mm if concrete.unit_system == "metric" else 3 / 8 * inch


def _max_bar_spacing_slab_under(section: "RectangularBeam", limit: Any) -> Any:
    """The lesser of 3h and ``limit``, tightened again for a member on the ground.

    The 3h of ACI 318-19 §7.7.2.3 / CIRSOC 201-25 art. 7.7.2.3 is common to the
    two codes; ``limit`` is the absolute term, which is not. The limit is on the
    flexural reinforcement of a one-way slab, so it is the spacing of the layer
    nearest the face that has to meet it.

    A footing is detailed more tightly than §7.7.2.3 alone would ask: a slab on
    the ground is thick, so 3h stops binding long before the bars are close
    enough to spread the bearing pressure into them, and practice caps them at
    300 mm. Its mat is placed as shrinkage and temperature reinforcement
    (ACI 318-19 §13.3.1.2 → §24.4.3.2 / CIRSOC 201-25 art. 13.3.1.2 →
    art. 24.4.3.2), whose own spacing limit -- the lesser of 5h and 450 mm,
    ACI 318-19 §24.4.3.3 / CIRSOC 201-25 art. 24.4.3.3, the same number in both
    -- is looser than what is applied here. Under CIRSOC the 300 mm is no longer
    only practice: art. 7.7.2.3 states it for every slab, so the two coincide.
    """
    if section.support == "soil":
        limit = min(limit, _MAX_BAR_SPACING_ON_SOIL)
    return min(3 * section.height, limit)


def _max_bar_spacing_slab(section: "RectangularBeam") -> Any:
    """ACI 318-19 §7.7.2.3: the lesser of 3h and 450 mm (18 in.).

    Read off the printed page: ACI 318-19 SI p. 95, "maximum spacing s of
    deformed longitudinal reinforcement shall be the lesser of 3h and 450 mm";
    the in-lb edition prints "3h and 18 in." on the same page.

    CIRSOC 201-25 art. 7.7.2.3 puts 300 mm where this has 450 mm, so each code
    registers its own hook; see :func:`_max_bar_spacing_slab_cirsoc`.
    """
    limit = (
        _MAX_BAR_SPACING_SLAB_ACI if section.concrete.unit_system == "metric" else _MAX_BAR_SPACING_SLAB_ACI_IMPERIAL
    )
    return _max_bar_spacing_slab_under(section, limit)


def _min_bar_spacing_slab(section: "RectangularBeam") -> Any:
    """The closest together the bars of a footing are detailed.

    Not a strength limit -- the clear distance of ACI 318-19 §25.2.1 /
    CIRSOC 201-25 §25.2.1 is, and the settings already carry it -- but the
    spacing below which a footing is not detailed in practice, because placing
    and vibrating over the soil stops being worth it. Neither code states it.
    A slab spanning between supports has no such floor.
    """
    return _MIN_BAR_SPACING_ON_SOIL if section.support == "soil" else None


def _min_effective_depth_on_soil(concrete: Any) -> Any:
    """ACI 318-19 §13.3.1.2 / CIRSOC 201-25 art. 13.3.1.2: d >= 150 mm (6 in.).

    Identical in both codes: *overall depth of foundation shall be selected
    such that the effective depth of bottom reinforcement is at least 150 mm*
    (6 in. in the in-lb edition; *al menos 150 mm* in CIRSOC). The limit is on
    the effective depth of the bottom reinforcement, not on the overall
    thickness, which is why it is registered as ``min_effective_depth_on_soil``
    and the section is measured to the bars rather than face to face.

    The distinction is not academic on a footing, because the cover is large:
    ACI 318-19 Table 20.5.1.3.1 asks 75 mm where the concrete is cast against
    and permanently in contact with the ground, and CIRSOC 201-25
    Tabla 20.5.1.3.1, fila (a) asks 55 mm under control of execution or 60 mm
    otherwise, not counting the blinding layer. A 200 mm section with a Ø16
    bottom mat is left with d = 117 mm under ACI and 132 to 137 mm under
    CIRSOC, both short of the 150 mm the clause asks for.
    """
    return 150 * mm if concrete.unit_system == "metric" else 6 * inch


def _min_stirrup_diameter_cirsoc(concrete: Any) -> Any:
    """CIRSOC's catalogue starts one size below ACI's.

    CIRSOC 201-25 §9.7.6.4.2, Tabla 9.7.6.4.2 puts the smallest stirrup that
    may laterally support compression reinforcement at 6 mm, rising to 8, 10
    and 12 mm as the longitudinal bar passes 16, 25 and 32 mm; ACI 318-19
    §9.7.6.4.2 gives a single No. 10 up to No. 32 and No. 13 above. Neither
    code states a minimum diameter for a stirrup placed purely for shear, so
    what this returns is the bottom of the catalogue rather than the table:
    the graded value would need the section's compression bars, which the hook
    does not receive.
    """
    return 6 * mm


def _max_bar_spacing_slab_cirsoc(section: "RectangularBeam") -> Any:
    """CIRSOC 201-25 art. 7.7.2.3: the lesser of 3h and 300 mm.

    Read off the printed page (Reglamento CIRSOC 201-25, Cap. 7-124): "la
    separacion maxima s de la armadura conformada longitudinal debe ser el
    menor entre 3h y 300 mm", against the 450 mm of ACI 318-19 §7.7.2.3. The
    300 mm therefore governs from 100 mm of thickness up, which is every slab
    worth designing: below that, 3h is the lesser of the two anyway.

    CIRSOC 201-25 is published in metric only, so unlike
    :func:`_max_bar_spacing_slab` there is no second number to choose between.
    """
    return _max_bar_spacing_slab_under(section, _MAX_BAR_SPACING_SLAB_CIRSOC)


def _stirrup_spacing_caps(concrete: Any) -> Any:
    """ACI 318-19 Table 9.7.6.2.2: 600 mm (24 in.) under the Vs threshold, 300 mm (12 in.) over it."""
    if concrete.unit_system == "metric":
        return 600 * mm, 300 * mm
    return 24 * inch, 12 * inch


def _stirrup_spacing_caps_cirsoc(concrete: Any) -> Any:
    """CIRSOC 201-25 Tabla 9.7.6.2.2: 400 mm under the Vs threshold, 200 mm over it.

    The ACI 318-19 Table 9.7.6.2.2 rows -- d/2 and d below the threshold, d/4
    and d/2 above it, the threshold itself 0.33*sqrt(f'c)*bw*d in both codes --
    with the caps of CIRSOC 201-2005 kept.
    """
    return 400 * mm, 200 * mm


def _min_shear_reinforcement_coefficient(concrete: Any) -> float:
    """ACI 318-19 §9.6.3.1: 0.083 in SI, 1.0 in US customary."""
    return 1.0 if concrete.is_imperial else 0.083


def _min_shear_reinforcement_coefficient_cirsoc(concrete: Any) -> float:
    """CIRSOC 201-25 §9.6.3.1: 0.085 in SI."""
    return 0.085


def _wall_axial_tension_divisor(concrete: Any) -> Any:
    """ACI 318-19 Eq. (11.5.4.4): 3.45 MPa (500 psi)."""
    return 500 * psi if concrete.is_imperial else 3.45 * MPa


def _wall_axial_tension_divisor_cirsoc(concrete: Any) -> Any:
    """CIRSOC 201-25 Eq. (11.5.4.4): 3.5 MPa."""
    return 3.5 * MPa


def _flexural_min_fy_cap(concrete: Any) -> Any:
    """ACI 318-19 §9.6.1.2: cap f_y in A_s,min at 550 MPa (80,000 psi)."""
    return 550 * MPa if concrete.unit_system == "metric" else 80_000 * psi


def _flexural_min_fy_cap_cirsoc(concrete: Any) -> Any:
    """CIRSOC 201-25 §9.6.1.2: cap f_y in A_s,min at 500 MPa."""
    return 500 * MPa


_SUMMARY_COLUMNS = {
    "moment_demand": "Mu",
    "shear_demand": "Vu",
    "shear_demand_source": "Vu",
    "axial_demand": "Nu",
    "moment_capacity_top": "ØMn,top",
    "moment_capacity_bot": "ØMn,bot",
    "shear_capacity": "ØVn",
}


def _capacity_columns(section: "RectangularBeam") -> dict:
    """ACI reports a reduced capacity phi*Mn on each face."""
    return {
        "ØMn,top": round(section._phi_M_n_top.to("kN*m").magnitude, 1),
        "ØMn,bot": round(section._phi_M_n_bot.to("kN*m").magnitude, 1),
    }


def _initialize_attributes(section: "RectangularBeam") -> None:
    """The zeroed result attributes the ACI report tables read off the beam.

    Moved here from ``RectangularBeam`` verbatim: they are this code's
    compatibility layer (ADR-0001), so the code owns them and a new code brings
    its own set without editing the element.
    """
    section._phi_V_n = 0 * kN
    section._phi_V_s = 0 * kN
    section._phi_V_c = 0 * kN
    section._phi_V_max = 0 * kN
    section._V_u = 0 * kN
    section._M_u = 0 * kNm
    section._M_u_bot = 0 * kNm
    section._M_u_top = 0 * kNm
    section._N_u = 0 * kN
    section._A_cv = 0 * cm**2
    section._k_c_min = 0 * MPa
    section._sigma_Nu = 0 * MPa
    section.V_c = 0 * kN
    section._rho_w = 0 * dimensionless
    section._lambda_s = 0
    section.f_yt = 0 * MPa
    section._max_shear_ok = False
    section._A_s_min_bot = 0 * cm**2
    section._A_s_min_top = 0 * cm**2
    section._A_s_min_eff_bot = 0 * cm**2
    section._A_s_min_eff_top = 0 * cm**2
    section._A_s_max_bot = 0 * cm**2
    section._A_s_max_top = 0 * cm**2
    section._A_s_max_eff_bot = 0 * cm**2
    section._A_s_max_eff_top = 0 * cm**2
    section._phi_M_n_bot = 0 * kNm
    section._phi_M_n_top = 0 * kNm
    section._d_b_max_bot = 0 * mm
    section._d_b_max_top = 0 * mm
    section.flexure_design_results_bot = None
    section.flexure_design_results_top = None
    section._A_s_bool_bot = False
    section._A_s_bool_top = False


def _transverse_rebar_aci(rebar: "Rebar", V_s_req: Any, alpha: float) -> Any:
    return rebar.transverse_rebar_ACI_318_19(V_s_req)


def _transverse_rebar_cirsoc(rebar: "Rebar", V_s_req: Any, alpha: float) -> Any:
    return rebar.transverse_rebar_CIRSOC_201_25(V_s_req)


def _longitudinal_rebar_aci(rebar: "Rebar", A_s_req: Any, A_s_max: Any, mech_cover: Any) -> Any:
    return rebar.longitudinal_rebar_ACI_318_19(A_s_req, A_s_max, mech_cover)


_COMMON = dict(
    check_shear=_check_shear_ACI_318_19,
    check_flexure=_check_flexure_ACI_318_19,
    apply_shear_state=apply_shear_state,
    apply_flexure_state=apply_flexure_state,
    design_shear=_design_shear_ACI_318_19,
    design_flexure=_design_flexure_ACI_318_19,
    longitudinal_rebar=_longitudinal_rebar_aci,
    initialize_attributes=_initialize_attributes,
    check_shear_wall=_check_shear_ACI_318_19_wall,
    # CIRSOC 201-25 reuses the ACI wall design; its bar catalogue is
    # selected inside the design itself.
    design_shear_wall=_design_shear_ACI_318_19_wall,
    apply_wall_shear_state=apply_wall_shear_state,
    # Two-way shear. `design_punching` is Phase 4; requires() names the code.
    check_punching=check_punching_ACI_318_19,
    flexure_symbols=_FLEXURE_SYMBOLS,
    units_row_shear=_UNITS_ROW_SHEAR,
    units_row_flexure=_UNITS_ROW_FLEXURE,
    summary_columns=_SUMMARY_COLUMNS,
    capacity_columns=_capacity_columns,
    shear_symbols=_SHEAR_SYMBOLS,
    summary_drop_columns=_SUMMARY_DROP_COLUMNS,
    # ``max_bar_spacing_slab`` is deliberately absent: the two codes print a
    # different absolute term in their 7.7.2.3, so each registers its own.
    min_bar_spacing_slab=_min_bar_spacing_slab,
    # §13.3.1.2 is written on d, not on h, so there is no overall-thickness
    # hook for these two codes; ``min_thickness_on_soil`` stays unset.
    min_effective_depth_on_soil=_min_effective_depth_on_soil,
    # A_s,max is the tension-controlled limit of §9.3.3.1 (Table 21.2.2), and
    # that is the limit a layout is held to, compression steel included.
    max_steel_is_ductility_limit=True,
    flexure_admissible=_flexure_ductile_ACI_318_19,
)

ACI_318_19 = register(
    DesignCode(
        title="ACI 318-19",
        year=2019,
        materials=(Concrete_ACI_318_19,),
        transverse_rebar=_transverse_rebar_aci,
        min_stirrup_diameter=_min_stirrup_diameter,
        stirrup_spacing_caps=_stirrup_spacing_caps,
        flexural_min_fy_cap=_flexural_min_fy_cap,
        min_shear_reinforcement_coefficient=_min_shear_reinforcement_coefficient,
        wall_axial_tension_divisor=_wall_axial_tension_divisor,
        max_bar_spacing_slab=_max_bar_spacing_slab,
        **_COMMON,  # type: ignore[arg-type]
    )
)

CIRSOC_201_25 = register(
    DesignCode(
        title="CIRSOC 201-25",
        year=2025,
        materials=(Concrete_CIRSOC_201_25,),
        # What CIRSOC does differently: the bar sizes of art. 20.2.1.3,
        # Tabla 20.2.1, the graded stirrup minimum of art. 9.7.6.4.2,
        # Tabla 9.7.6.4.2, the caps of Tabla 9.7.6.2.2, and the slab bar
        # spacing of art. 7.7.2.3. Everything else is the ACI text under the
        # same numbering; see the module docstring for what is not hooked.
        transverse_rebar=_transverse_rebar_cirsoc,
        min_stirrup_diameter=_min_stirrup_diameter_cirsoc,
        stirrup_spacing_caps=_stirrup_spacing_caps_cirsoc,
        flexural_min_fy_cap=_flexural_min_fy_cap_cirsoc,
        min_shear_reinforcement_coefficient=_min_shear_reinforcement_coefficient_cirsoc,
        wall_axial_tension_divisor=_wall_axial_tension_divisor_cirsoc,
        max_bar_spacing_slab=_max_bar_spacing_slab_cirsoc,
        **_COMMON,  # type: ignore[arg-type]
    )
)
