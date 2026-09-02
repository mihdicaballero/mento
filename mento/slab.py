from __future__ import annotations
from typing import Any, Callable, ClassVar, TYPE_CHECKING, cast
from dataclasses import dataclass
import math
import warnings

import numpy as np

from mento.beam import RectangularBeam
from mento.codes.registry import design_code
from mento.units import m, mm, cm, inch

#: Area unit the mat search works in, built once.
_MM2 = mm**2

if TYPE_CHECKING:
    from pint import Quantity


def _bars_at_spacing(spacing: Quantity, width: Quantity) -> int:
    """How many bars a strip of ``width`` carries at a centre-to-centre ``spacing``.

    The one place the count-from-spacing rule lives, so that a spacing chosen
    for a given number of bars can be checked against the count it produces.
    """
    if spacing == 0 * mm:
        return 0  # Default to 0 bars if spacing is zero
    # Round up to ensure full bars (e.g., 3.2 bars → 4 bars)
    return int(np.ceil((width.to("cm").magnitude / spacing.to("cm").magnitude)))  # Returns dimensionless count


@dataclass
class OneWaySlab(RectangularBeam):
    """
    One-way slab section. Inherits all flexure and shear checks from RectangularBeam,
    but uses bar diameters + spacing instead of number of bars for longitudinal reinforcement.
    """

    def __post_init__(self) -> None:
        super().__post_init__()

        self.mode = "slab"

        self._stirrup_d_b = 0 * mm
        self._stirrup_s_l = 0 * mm
        self._A_v = 0 * mm**2 / m

        # Allow slabs to place as many bars as minimum spacing permits
        self.settings.max_bars_per_layer = max(self.settings.max_bars_per_layer, 200)
        # ... and no more than the code's own minimum spacing permits, which is
        # the only way that limit can be enforced: the rebar search chooses a
        # bar count, and the spacing is what that count implies. Capping the
        # count is capping the spacing from below.
        minimum = self._min_bar_spacing()
        if minimum is not None:
            self.settings.max_bars_per_layer = max(1, int(self.width / minimum))
        # Force same diameter for all bars within a layer
        self.settings.max_diameter_diff = 0 * self.settings.max_diameter_diff

        self._initialize_longitudinal_rebar_attributes()

    def _initialize_longitudinal_rebar_attributes(self) -> None:
        """Initialize all rebar-related attributes with default values."""

        # Unit system default starter rebar for initial effective height
        # In this case is 0 mm since the slab might not have bottom or top rebar
        if self.concrete.unit_system == "metric":
            self._n1_b, self._d_b1_b = 0, 0 * mm
            self._n1_t, self._d_b1_t = 0, 0 * mm
        else:
            self._n1_b, self._d_b1_b = 0, 0 * inch
            self._n1_t, self._d_b1_t = 0, 0 * inch

        # Bottom rebar defaults (diameter and spacing)
        self._s_b1_b = 0 * mm
        self._d_b2_b, self._s_b2_b = 0 * mm, 0 * mm
        self._d_b3_b, self._s_b3_b = 0 * mm, 0 * mm
        self._d_b4_b, self._s_b4_b = 0 * mm, 0 * mm
        self._n2_b = 0  # No position 2 in slabs
        self._n4_b = 0  # No position 4 in slabs

        # Top rebar defaults (diameter and spacing)
        self._s_b1_t = 0 * mm
        self._d_b2_t, self._s_b2_t = 0 * mm, 0 * mm
        self._d_b3_t, self._s_b3_t = 0 * mm, 0 * mm
        self._d_b4_t, self._s_b4_t = 0 * mm, 0 * mm
        self._n2_t = 0  # No position 2 in slabs
        self._n4_t = 0  # No position 4 in slabs

        # Update dependent attributes
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def set_slab_transverse_rebar(
        self,
        d_b: Quantity = 0 * mm,
        s_long: Quantity = 0 * cm,
        s_trans: Quantity = 0 * cm,
    ) -> None:
        """Sets the transverse rebar in the slab section.

        A zero diameter or a zero spacing in either direction means no stirrups,
        which clears the transverse reinforcement.
        """
        self._stirrup_s_l = s_long
        if d_b == 0 * mm or s_long == 0 * cm or s_trans == 0 * cm:
            # No stirrups: same state as a section without transverse rebar
            self._A_v = 0 * mm**2 / m
            self._stirrup_n = 0
            self._stirrup_d_b = 0 * mm
            self._stirrup_s_w = 0 * mm
        else:
            n_legs_per_unit_width = self.width / s_trans
            A_db = (d_b**2) * math.pi / 4  # Area of one stirrup leg per unit width
            self._A_v = A_db * n_legs_per_unit_width / s_long  # Legs area per unit length
            self._stirrup_d_b = d_b
            self._stirrup_s_w = s_trans
            # Both shear checks read _stirrup_n to decide whether the section
            # carries shear reinforcement at all, and the drawings read it as a
            # number of stirrups. Leaving it at zero was why a slab given
            # transverse rebar was still checked as if it had none. A strip has
            # no closed stirrup to count -- the legs at s_trans need not divide
            # into it a whole number of times -- so this is the equivalent
            # number of rows; _A_v above is what the capacity is computed from.
            self._stirrup_n = max(1, round(n_legs_per_unit_width.to("dimensionless").magnitude / 2))

        # Update effective heights
        self._update_effective_heights()

    def _leg_spacing_across_width(self) -> Quantity:
        """The transverse spacing the strip was detailed with.

        A beam derives this from the geometry of its stirrup cage. A slab does
        not have one: the legs sit on a grid, and how far apart they are across
        the strip is chosen, not implied. Deriving it the beam's way put the
        legs of a designed slab at a spacing it had never been given.
        """
        return self._stirrup_s_w

    def _apply_transverse_design(self, design: Any) -> None:
        """A slab is detailed by a spacing in each direction, not by a stirrup count."""
        self.set_slab_transverse_rebar(
            d_b=design["d_b"],
            s_long=design["s_l"],
            s_trans=design["s_w"],
        )

    def _zero_diameter(self) -> Quantity:
        """A bar diameter of zero, in the length unit of this slab."""
        return 0 * mm if self.concrete.unit_system == "metric" else 0 * inch

    def _max_bar_spacing(self) -> Quantity | None:
        """The largest spacing the design code allows between the flexural bars.

        ACI 318-19 7.7.2.3 and EN 1992-1-1 9.3.1.1(3) both cap it as a multiple
        of the thickness, which is what keeps a slab from being detailed as a
        handful of widely spaced bars that happen to add up to the area. A code
        that states no such limit returns ``None``.
        """
        limit = design_code(self.concrete).max_bar_spacing_slab
        return None if limit is None else cast("Quantity", limit(self))

    def _min_bar_spacing(self) -> Quantity | None:
        """The smallest spacing the design code asks for between the flexural bars.

        Distinct from the clear distance in :class:`~mento.settings.BeamSettings`,
        which every section has and which is about fitting a bar and a vibrator
        between two others. This is the code's own floor, and only a footing has
        one; a slab spanning between supports returns ``None``.
        """
        limit = design_code(self.concrete).min_bar_spacing_slab
        return None if limit is None else cast("Quantity | None", limit(self))

    def _spacing_for_bars(self, n: int) -> Quantity:
        """The spacing that puts ``n`` bars across the strip, as it would be drawn.

        The exact answer is ``width / n``, which is rarely a number anyone
        details to, so it is rounded to the whole centimetre (inch). Rounding up
        is what usually keeps the layout the search chose -- the count comes
        back as ``ceil(width / s)``, so a slightly wider spacing still asks for
        the same ``n`` bars, while a narrower one would silently add one -- but
        it is checked rather than assumed, because it does not always hold. The
        result is then capped at what the code allows between the bars of a
        slab, which only ever asks for more of them.
        """
        if n <= 0:
            return 0 * cm
        unit = cm if self.concrete.unit_system == "metric" else inch
        exact = (self.width / n).to(unit).magnitude
        spacing = math.ceil(exact) * unit
        if _bars_at_spacing(spacing, self.width) < n:
            # Rounding up costs a bar once the bars are close enough that a whole
            # centimetre spans more than one of them -- from about eleven bars in
            # a metre. Round down there instead: that can only ask for more bars
            # than the design chose, never fewer, so the strip is never detailed
            # with less steel than it needs. (The floor is only ever zero at a
            # spacing below one centimetre, where the bars would already overlap.)
            spacing = max(math.floor(exact), 1) * unit
        limit = self._max_bar_spacing()
        if limit is not None:
            # The area the search asked for is not the only thing the layout has
            # to satisfy: bars far enough apart leave the slab between them
            # unreinforced whatever they add up to. Capping the spacing only
            # ever adds bars, so the area is still covered.
            spacing = min(spacing, math.floor(limit.to(unit).magnitude) * unit)
        return spacing

    def _layer_from_design(self, design: Any, first: int, second: int) -> tuple[Quantity, Quantity]:
        """The diameter and spacing of the layer a design gives as two bar groups."""
        n_first = int(design.get(f"n_{first}", 0) or 0)
        n_second = int(design.get(f"n_{second}", 0) or 0)
        # A slab is searched with max_diameter_diff = 0, so a layer only ever
        # comes back with one diameter: whichever group carries the bars.
        d_b = design.get(f"d_b{first}") if n_first > 0 else design.get(f"d_b{second}")
        if d_b is None or d_b.magnitude == 0:
            return self._zero_diameter(), 0 * cm
        return d_b, self._spacing_for_bars(n_first + n_second)

    def _apply_longitudinal_design(self, design: Any, face: str) -> None:
        """Apply a design given in bar counts to a face, as a spacing per layer.

        The rebar search is the beam's -- a slab is not worth a second one --
        and it answers with groups of bars: ``n_1`` and ``n_2`` are the layer
        nearest the face, ``n_3`` and ``n_4`` the one behind it. A slab is not
        detailed that way; it is a diameter repeated at a spacing. So each
        layer is spread over the strip here, and the spacings are what the slab
        stores. Applying the design the beam's way left the bars split across
        two groups that are really one layer, and left the spacings -- the only
        thing a slab is drawn with -- empty.
        """
        d_b1, s_b1 = self._layer_from_design(design, 1, 2)
        d_b3, s_b3 = self._layer_from_design(design, 3, 4)
        setattr(self, f"_d_b1_{face}", d_b1)
        setattr(self, f"_s_b1_{face}", s_b1)
        setattr(self, f"_d_b3_{face}", d_b3)
        setattr(self, f"_s_b3_{face}", s_b3)
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def _apply_longitudinal_design_bot(self, design: Any) -> None:
        """See :meth:`_apply_longitudinal_design`."""
        self._apply_longitudinal_design(design, "b")

    def _apply_longitudinal_design_top(self, design: Any) -> None:
        """See :meth:`_apply_longitudinal_design`."""
        self._apply_longitudinal_design(design, "t")

    def _clear_top_longitudinal(self) -> None:
        """Clear the top face, its spacings included.

        A beam clears a face by zeroing the bar counts. On a slab those counts
        are derived from the spacings, so a spacing left behind would put the
        bars back the next time any layer is recalculated.
        """
        self._d_b1_t, self._s_b1_t = self._zero_diameter(), 0 * cm
        self._d_b3_t, self._s_b3_t = self._zero_diameter(), 0 * cm
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def set_slab_longitudinal_rebar_bot(
        self,
        d_b1: Quantity = 0 * mm,
        s_b1: Quantity = 0 * mm,
        d_b3: Quantity = 0 * mm,
        s_b3: Quantity = 0 * mm,
    ) -> None:
        """Update the bottom rebar configuration and recalculate attributes.

        Args:
            d_b1: Diameter of position 1 bottom rebar (default: 0 mm).
            s_b1: Spacing of position 1 bottom rebar (default: 0 mm).
            d_b3: Diameter of position 3 bottom rebar (default: 0 mm).
            s_b3: Spacing of position 3 bottom rebar (default: 0 mm).
        """
        self._d_b1_b = d_b1 if d_b1 != 0 * mm else self._d_b1_b
        self._s_b1_b = s_b1 if s_b1 != 0 * mm else self._s_b1_b
        self._d_b3_b = d_b3 if d_b3 != 0 * mm else self._d_b3_b
        self._s_b3_b = s_b3 if s_b3 != 0 * mm else self._s_b3_b
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def set_slab_longitudinal_rebar_top(
        self,
        d_b1: Quantity = 0 * mm,
        s_b1: Quantity = 0 * mm,
        d_b3: Quantity = 0 * mm,
        s_b3: Quantity = 0 * mm,
    ) -> None:
        """Update the top rebar configuration and recalculate attributes.

        Args:
            d_b1: Diameter of position 1 top rebar (default: 0 mm).
            s_b1: Spacing of position 1 top rebar (default: 0 mm).
            d_b3: Diameter of position 3 top rebar (default: 0 mm).
            s_b3: Spacing of position 3 top rebar (default: 0 mm).

        """
        self._d_b1_t = d_b1 if d_b1 != 0 * mm else self._d_b1_t
        self._s_b1_t = s_b1 if s_b1 != 0 * mm else self._s_b1_t
        self._d_b3_t = d_b3 if d_b3 != 0 * mm else self._d_b3_t
        self._s_b3_t = s_b3 if s_b3 != 0 * mm else self._s_b3_t
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def _calculate_longitudinal_rebars(self) -> None:
        """Calculate the total rebar area for a slab, given spacing and slab width."""

        # --- BOTTOM REBAR ---
        self._n1_b = _bars_at_spacing(self._s_b1_b, self.width)
        self._n3_b = _bars_at_spacing(self._s_b3_b, self.width)

        # --- TOP REBAR ---
        self._n1_t = _bars_at_spacing(self._s_b1_t, self.width)
        self._n3_t = _bars_at_spacing(self._s_b3_t, self.width)


@dataclass
class Footing(OneWaySlab):
    """A spread footing or raft: a one-way slab bearing directly on the ground.

    Everything about it is a :class:`OneWaySlab` -- the same flexure and shear
    checks, the same reinforcement given as diameter and spacing -- with one
    difference, and it is a difference the design codes make rather than this
    class: the minimum longitudinal reinforcement.

    A member spanning between supports is given a minimum sized to keep it from
    failing the instant it cracks. A member on the ground cannot fail that way,
    because the soil goes on carrying it, so both codes exempt it and put a
    different rule in place:

    * ACI 318-19 §9.6.1.1(b) grants the exemption and §13.3.1.2 replaces the
      flexural minimum with the shrinkage and temperature reinforcement of
      §24.4.3.2 -- 0.0018*b*h at f_y = 420 MPa, on the gross section. CIRSOC
      201-25 shares the clause.
    * EN 1992-1-1 takes the larger of the halved geometric minimum of a
      foundation and the crack-control minimum of §7.3.2(2). The second governs
      the thin footings, not the thick ones: its ratio goes with k/2, and k
      decays from 1.00 to 0.65 between 300 and 800 mm of depth.

    Both are returned already resolved to the largest applicable minimum, so a
    caller designs a footing and takes the answer as it comes.

    Example::

        footing = Footing(
            label="Z1", concrete=concrete, steel_bar=steel,
            width=1 * m, height=60 * cm, c_c=50 * mm,
        )
        Node(section=footing, forces=[Forces(M_y=120 * kNm)]).design()

    It is also detailed differently, which the same ``support`` tells the codes:
    the bars are kept between 100 and 300 mm apart, and a footing thinner than
    its code allows says so when it is built. And because a footing is placed as
    a single mat rather than as two independently detailed faces, a design ends
    with both faces set out at one spacing -- see
    :meth:`_finalize_longitudinal_design`.

    Note:
        A `Footing` designs the section, not the foundation. It says nothing
        about the soil below it -- bearing pressure, punching (see
        :class:`~mento.punching.PunchingSlab`), sliding and overturning are not
        part of a section design.
    """

    #: What makes it a footing, and the only thing the design codes read.
    support: ClassVar[str] = "soil"

    #: The thinnest bar a footing mat is detailed with, whatever the catalogue
    #: holds below it. Practice rather than code: a mesh finer than this is not
    #: what goes down over the ground, and the search below would otherwise
    #: reach for it to shave the lightly loaded face.
    _MIN_MAT_DIAMETER = 10 * mm

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.concrete.unit_system == "metric":
            self.settings.minimum_longitudinal_diameter = max(
                self.settings.minimum_longitudinal_diameter, self._MIN_MAT_DIAMETER
            )
        self._warn_if_thinner_than_the_code_allows()

    #: The layers a mat is placed in, as (diameter attribute, spacing attribute).
    _MAT_LAYERS = (("_d_b1_b", "_s_b1_b"), ("_d_b3_b", "_s_b3_b"), ("_d_b1_t", "_s_b1_t"), ("_d_b3_t", "_s_b3_t"))

    #: How many candidate mats are verified before falling back. They are tried
    #: lightest first, so the first is almost always the answer; the budget only
    #: bounds the cost when the lightest ones do not verify.
    _MAT_CANDIDATES = 12

    def _finalize_longitudinal_design(
        self,
        A_req_bot: Quantity,
        A_req_top: Quantity,
        resists: Callable[[], bool],
    ) -> None:
        """Detail both faces as one mat, at the least steel that still covers them.

        A footing is not built as two independently detailed faces. The bars go
        down as one grid: the bottom is placed, the top is placed over it, and
        the two are tied. The module has to be one number, and on a drawing the
        top sits either at the bottom spacing or at exactly twice it, so that
        one top bar lands on every second bottom bar.

        Reconciling the two spacings the design arrived at is not enough. The
        face governed by the minimum is detailed as many thin bars, which comes
        out closer than the few thick bars the loaded face needs, so simply
        taking the smaller of the two often closes the *loaded* face up and buys
        steel nobody asked for -- about 20% more than the two faces actually
        need, measured across a range of typical footings.

        So the mat is chosen rather than reconciled: over the spacings the code
        allows and the bars the catalogue holds, the lightest combination of
        module and two diameters -- with the top at ``s`` or ``2s`` -- that
        covers what each face requires. That lands within about 7% (EN) to 11%
        (ACI) of the two faces' own requirement, against the 20% of taking the
        smaller.

        Every candidate is verified before it is kept, because the choice feeds
        back into itself: a thicker bar sits deeper, which lowers the effective
        depth and raises the area the face needed in the first place. Anything
        that fails to resist the moments is passed over, and if nothing within
        the budget verifies, the design falls back to reconciling the two
        spacings it had already proved.

        Args:
            A_req_bot: Steel the bottom face has to carry, minimum included.
            A_req_top: The same for the top; zero when it carries none.
            resists: Whether the layout now on the section resists both design
                moments.
        """
        snapshot = [(name, getattr(self, name)) for pair in self._MAT_LAYERS for name in pair]
        for spacing, diameter in self._mat_candidates(A_req_bot, A_req_top):
            self._apply_mat(spacing, diameter)
            if resists():
                return

        for name, value in snapshot:
            setattr(self, name, value)
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()
        self._match_spacings_downward()

    def _mat_candidates(self, A_req_bot: Quantity, A_req_top: Quantity) -> Any:
        """The mats worth trying, lightest first, as (spacings, diameters)."""
        needed_bot = A_req_bot.to(_MM2).magnitude if A_req_bot is not None else 0.0
        needed_top = A_req_top.to(_MM2).magnitude if A_req_top is not None else 0.0
        if needed_bot <= 0 and needed_top <= 0:
            return

        limit_max = self._max_bar_spacing()
        if limit_max is None:
            return  # a code that caps nothing gives the search no space to work in
        limit_min = self._min_bar_spacing()

        width = self.width.to(mm).magnitude
        clear = self.settings.clear_spacing.to(mm).magnitude
        # The unit a spacing is drawn in, so the module is one a detailer writes.
        step = 10.0 if self.concrete.unit_system == "metric" else 25.4
        s_lo = step if limit_min is None else limit_min.to(mm).magnitude
        s_hi = limit_max.to(mm).magnitude

        rebar = self._create_rebar_designer()
        bars = [
            (d.to(mm).magnitude, rebar.rebar_areas[d].to(_MM2).magnitude)
            for d in rebar.rebar_diameters
            if self.settings.minimum_longitudinal_diameter <= d <= self.settings.max_longitudinal_diameter
        ]

        def grid(spacing: float, needed: float) -> Any:
            """The lightest bar that covers ``needed`` at ``spacing``, or None."""
            if needed <= 0:
                return (0.0, 0.0)
            n = math.ceil(width / spacing)
            covering = [(n * area, d) for d, area in bars if spacing - d >= clear and n * area >= needed]
            return min(covering) if covering else None

        found = []
        for k in range(math.ceil(s_lo / step), math.floor(s_hi / step) + 1):
            s_bot = k * step
            bottom = grid(s_bot, needed_bot)
            if bottom is None:
                continue
            for factor in (1, 2):
                s_top = s_bot * factor
                if s_top > s_hi:
                    continue
                top = grid(s_top, needed_top)
                if top is None:
                    continue
                # Lightest first; among equals the wider module, which is the
                # cheaper mat to place for the same steel.
                found.append((bottom[0] + top[0], -s_bot, (s_bot, s_top), (bottom[1], top[1])))

        for _, _, spacing, diameter in sorted(found)[: self._MAT_CANDIDATES]:
            yield spacing, diameter

    def _apply_mat(self, spacing: Any, diameter: Any) -> None:
        """Place a single-layer grid on each face, clearing any second layer."""
        s_bot, s_top = spacing
        d_bot, d_top = diameter
        # A face that needs no steel gets no grid at all -- not a spacing with a
        # zero bar on it, which reads as a grid everywhere it is reported and
        # would put a row of nothing on the drawing.
        s_bot = s_bot if d_bot > 0 else 0.0
        s_top = s_top if d_top > 0 else 0.0
        self._d_b1_b, self._s_b1_b = d_bot * mm, s_bot * mm
        self._d_b1_t, self._s_b1_t = d_top * mm, s_top * mm
        self._d_b3_b, self._s_b3_b = self._zero_diameter(), 0 * cm
        self._d_b3_t, self._s_b3_t = self._zero_diameter(), 0 * cm
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def _match_spacings_downward(self) -> None:
        """Fallback: reconcile the spacings the design arrived at, taking the
        smaller of the two.

        Safe by construction -- it only ever adds bars to the face that had the
        wider spacing -- and needs no verification for that reason, which is
        what makes it the thing to fall back to when no searched mat verifies.
        """
        placed = [
            (diameter, spacing)
            for diameter, spacing in self._MAT_LAYERS
            if getattr(self, diameter).magnitude > 0 and getattr(self, spacing).magnitude > 0
        ]
        if len(placed) < 2:
            return  # one grid, or none: nothing to reconcile

        module = min(getattr(self, spacing) for _, spacing in placed)
        if all(getattr(self, spacing) == module for _, spacing in placed):
            return  # already one mat

        for _, spacing in placed:
            setattr(self, spacing, module)
        self._calculate_longitudinal_rebars()
        self._update_longitudinal_rebar_attributes()

    def _warn_if_thinner_than_the_code_allows(self) -> None:
        """Say so when the section is thinner than a footing is built.

        A warning and not an error: the thickness is the engineer's to choose,
        and the design that follows is still the right design for the section
        it was given. Raising would refuse to answer a question that has an
        answer.
        """
        minimum = design_code(self.concrete).min_thickness_on_soil
        if minimum is None:
            return
        limit = cast("Quantity", minimum(self.concrete))
        if self.height < limit:
            label = f" {self.label}" if self.label else ""
            warnings.warn(
                f"Footing{label} is {self.height:.4g~P} thick, below the "
                f"{limit:.4g~P} {self.concrete.design_code} asks of a footing on soil. "
                "It is designed as given.",
                UserWarning,
                stacklevel=3,
            )
