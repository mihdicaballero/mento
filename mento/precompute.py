"""A section's geometry and materials as plain floats, converted once.

ADR-0005 puts pint at the boundary and floats in the equations, and Phase 1 did
that for the equations themselves. It left the *state* in pint, so every helper
unwrapped its inputs and re-wrapped its output, and the boundary was crossed at
every step rather than once.

Measured on an ACI shear check, that was 19 unit-changing conversions costing
~1.17 ms of the check's 1.43 ms. A conversion that actually changes units costs
~58 us; one that does not costs ~3.9 us. Nothing here is about how the unit is
spelled -- ``.to("cm")`` and ``.to(cm)`` measure the same -- it is about not
converting at all.

So a section publishes its values once, in the unit system its design code
writes its coefficients for, and the check runs entirely in floats. The result
is wrapped back into pint where it leaves the check.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from mento.units import cm, ft, inch, kip, kN, kNm, lbf, m, mm, MPa, N, psi

if TYPE_CHECKING:
    from mento.beam import RectangularBeam

#: The units the design-code equations are written in, per unit system. ACI
#: publishes separately rounded coefficients for SI and US customary (ADR-0005),
#: so the system is a property of the section, not a preference.
CANONICAL = {
    False: {"length": mm, "area": mm**2, "stress": MPa, "force": N, "moment": N * mm},
    True: {"length": inch, "area": inch**2, "stress": psi, "force": lbf, "moment": lbf * inch},
}

#: The units the compatibility layer and the report tables expect to read back.
DISPLAY = {
    False: {"length": cm, "area": cm**2, "stress": MPa, "force": kN, "per_length": cm**2 / m, "moment": kNm},
    True: {
        "length": inch,
        "area": inch**2,
        "stress": psi,
        "force": kip,
        "per_length": inch**2 / ft,
        "moment": kip * ft,
    },
}

# A per-length reinforcement ratio (area per unit length) has the dimension of a
# length, and the equations return it as one.
# A dimensionless ratio has no entry here on purpose: it needs no unit, and
# `to_display` returns it without consulting these tables.
for _imperial in (False, True):
    CANONICAL[_imperial]["per_length"] = CANONICAL[_imperial]["length"]


@dataclass(frozen=True)
class SectionFloats:
    """One beam's geometry and materials, in its design code's units.

    Cached on the section and rebuilt whenever the reinforcement or the
    effective depths change, so a loop over load combinations converts nothing.
    """

    is_imperial: bool
    width: float
    height: float
    c_c: float
    d_shear: float
    d_bot: float
    d_top: float
    A_x: float
    A_s_bot: float
    A_s_top: float
    c_mec_bot: float
    c_mec_top: float
    A_v: float
    stirrup_d_b: float
    stirrup_n: int
    f_c: float
    f_y: float
    E_s: float
    # lambda_factor and phi_v are deliberately absent: they are plain floats
    # on the concrete already, and code-specific. Nothing to convert, so
    # nothing to precompute.


def section_floats(section: "RectangularBeam") -> SectionFloats:
    """The section's float view, as of its last change.

    A plain read, with no lazy fallback: ``__post_init__`` builds the view and
    every later change rebuilds it, so a missing one is a bug in that chain and
    should say so rather than be quietly papered over.
    """
    return section.__dict__["_floats"]


def refresh_section_floats(section: "RectangularBeam") -> SectionFloats:
    """Rebuild the float view from the section as it stands now.

    Called eagerly from the one place every geometry and reinforcement change
    funnels through, rather than lazily on first use. Eagerly, because a check
    must not so much as *create* an attribute on the section -- the guards in
    ``test_beam.py`` compare the section's ``__dict__`` before and after, and a
    lazily-filled memo would show up there and be indistinguishable from a
    check writing its results back. It costs nothing: a design triggers seven
    of these, against thousands of conversions saved per check.
    """
    imperial = section.concrete.is_imperial
    length = CANONICAL[imperial]["length"]
    area = CANONICAL[imperial]["area"]
    stress = CANONICAL[imperial]["stress"]

    floats = SectionFloats(
        is_imperial=imperial,
        width=section.width.to(length).magnitude,
        height=section.height.to(length).magnitude,
        c_c=section.c_c.to(length).magnitude,
        d_shear=section._d_shear.to(length).magnitude,
        d_bot=section._d_bot.to(length).magnitude,
        d_top=section._d_top.to(length).magnitude,
        A_x=section._A_x.to(area).magnitude,
        A_s_bot=section._A_s_bot.to(area).magnitude,
        A_s_top=section._A_s_top.to(area).magnitude,
        c_mec_bot=section._c_mec_bot.to(length).magnitude,
        c_mec_top=section._c_mec_top.to(length).magnitude,
        # An area per unit length converts to a plain length.
        A_v=section._A_v.to(length).magnitude,
        stirrup_d_b=section._stirrup_d_b.to(length).magnitude,
        stirrup_n=section._stirrup_n,
        f_c=section.concrete.f_c.to(stress).magnitude,
        f_y=section.steel_bar.f_y.to(stress).magnitude,
        E_s=section.steel_bar._E_s.to(stress).magnitude,
    )
    section.__dict__["_floats"] = floats
    return floats
