"""Internal engine shared by the flexural design routines of each design code.

Private to ``mento.codes``: not part of the public API, and deliberately kept out
of the published documentation.
"""

# ---------------------------------------------------------------------------
# Why this module exists
#
# ACI 318-19 and EN 1992-2004 disagree on the equations -- the stress block, the
# safety format, the minimum reinforcement rules -- but not on the strategy that
# surrounds them. CIRSOC 201-25, the Argentine adoption of ACI 318-19, keeps its
# article numbering and reaches this engine through the very same callbacks, so
# a clause written below as ``ACI 318-19 §9.5.1.1(a) / CIRSOC 201-25 §9.5.1.1(a)``
# is one clause carrying two names. The strategy:
#
#   1. Guess the mechanical covers, hence the effective depths.
#   2. Ask the design code for the steel required on each face.
#   3. Let the discrete rebar selector turn those areas into a buildable layout.
#   4. Reconcile the faces: the layout on one face must also cover the
#      compression the opposite face's moment demands from it.
#   5. Re-read the real centroids, which move the effective depths, and iterate
#      (Picard / fixed point) until the covers stop moving.
#   6. Verify the layout actually resists the moment, and if it does not, fall
#      back to the best layout among those visited.
#
# Steps 1 and 3-6 live here. Step 2, and the capacity evaluation of step 6, are
# supplied by the caller as two callbacks -- the only place the design code
# enters. Each code keeps its own equations in
# ``_calculate_flexural_reinforcement_*`` and ``_determine_nominal_moment_*``.
#
# No design-code equation lives in this module, which is why it cites almost no
# clause: A_s,min, A_s,max and the stress block all arrive already evaluated
# through those two callbacks. The one normative rule the engine applies on its
# own account is the design-strength criterion of step 6 -- ACI 318-19
# §9.5.1.1(a) / CIRSOC 201-25 §9.5.1.1(a), phi*Mn >= Mu, and M_Rd >= M_Ed under
# EN 1992-2004 -- which is the comparison in ``_best_visited_pair`` and in the
# final verification of ``_run_flexure_design``.
#
# Kept as comments rather than a module docstring on purpose: Sphinx autodoc
# publishes docstrings, and this is implementation detail, not reference
# material for users.
# ---------------------------------------------------------------------------

from dataclasses import dataclass

import pandas as pd
from typing import TYPE_CHECKING, Any, Callable, Dict, Optional, Tuple

from mento.units import Quantity

from mento.rebar import RebarDesignInfeasibleError
from mento.units import cm, kNm, mm

if TYPE_CHECKING:
    from ..beam import RectangularBeam


_MAX_FLEXURE_ITERATIONS = 30  # safety net for slow divergence without cycling
# Compression layouts verified against the code's capacity per selection; the
# rest are judged by area alone. Enough to reach past the thin two-layer ones.
_MAX_COMPRESSION_CANDIDATES = 40
# Of those short of the area their depth asks for, the ones tried as well when
# none that covers it resists.
_MAX_SHORT_CANDIDATES = 12


@dataclass
class _FaceDemand:
    """Steel areas a design code requires for one tension face.

    Attributes
    ----------
    A_s_min, A_s_max:
        Code limits for the tension face being solved. Under ACI 318-19 and
        CIRSOC 201-25 they are the minimum flexural reinforcement of §9.6.1.2
        — or, for a member the ground supports, the 0.0018*Ag that §13.3.2.1
        sends to ACI 318-19 §7.6.1.1 / CIRSOC 201-25 §7.6.1 — and the
        tension-controlled cap that §9.3.3.1 with Table 21.2.2 puts on a beam
        in both codes. The numbers come from the design code's own hook; this
        dataclass only carries them.
    A_s_tension:
        Steel required on the tension face itself.
    A_s_compression:
        Steel this face's moment requires on the OPPOSITE face, i.e. the
        compression reinforcement of a doubly reinforced section. Zero when the
        section does not need it.
    compression_per_excess:
        Compression steel on the opposite face that each unit of tension steel
        above ``A_s_max`` calls for to keep the section at its ductility limit
        -- ``f_y / f_s'`` under ACI 318-19 / CIRSOC 201-25. ``None`` when
        ``A_s_max`` is not a ductility limit that compression steel extends
        (EN 1992-2004, where it is the 4 % of 9.2.1.1(3)), and then the selector
        never goes past it.
    """

    A_s_min: Quantity
    A_s_max: Quantity
    A_s_tension: Quantity
    A_s_compression: Quantity
    compression_per_excess: Optional[float] = None


# Ask the design code for the areas required on ``face`` ("bot"/"top") by a
# moment ``M`` (always a positive magnitude), given the tension-side effective
# depth ``d`` and the compression-side mechanical cover ``d_prime``.
_RequiredAreas = Callable[..., _FaceDemand]

# Resisting moment of the layout currently applied to the section, on ``face``,
# under a demand of ``M_demand`` (positive magnitude). Includes the safety
# format of the code: phi*Mn for ACI 318-19 and CIRSOC 201-25, with the
# phi = 0.90 that Table 21.2.2 of both gives a tension-controlled section;
# M_Rd for EN.
_Capacity = Callable[..., Quantity]

# Does the layout currently applied, with ``face`` in tension, keep within the
# reinforcement limits of the code? Strength alone does not say it. Under
# ACI 318-19 / CIRSOC 201-25 §9.3.3.1 a beam has to be tension-controlled, and
# past that limit the capacity only drops, through phi, and may still reach the
# moment. Under EN 1992-1-1 §9.2.1.1(3) no face may carry more than its
# maximum, and bars past it add resistance the design must not rely on.
_Admissible = Callable[[str], bool]


def _rebar_design_fingerprint(rebar_design: Any) -> tuple:
    """Canonical, hashable identifier of a discrete rebar layout.

    Two iterations producing the same fingerprint mean the Picard (fixed-point
    iteration) loop has cycled back to a previous configuration. The
    centroid/spacing are intentionally excluded — they are derived from this
    tuple, not part of its identity.
    """

    def _diam_mm(key: str) -> float:
        q = rebar_design.get(key, 0 * mm)
        return float(q.to("mm").magnitude) if q is not None else 0.0

    return (
        int(rebar_design.get("n_1", 0)),
        _diam_mm("d_b1"),
        int(rebar_design.get("n_2", 0)),
        _diam_mm("d_b2"),
        int(rebar_design.get("n_3", 0)),
        _diam_mm("d_b3"),
        int(rebar_design.get("n_4", 0)),
        _diam_mm("d_b4"),
    )


def _best_visited_pair(
    self: "RectangularBeam",
    bot_visited: Dict[tuple, dict],
    top_visited: Dict[tuple, dict],
    assess: Callable[[], tuple[bool, float]],
) -> tuple[Any, Any]:
    """Apply, and return, the best pair of the layouts visited on each face.

    Every layout the loop left on the bottom is tried with every one it left on
    the top -- including none, if the top never carried one -- because a face
    is only right or wrong together with the other: the compression steel of
    one sets what the other can take. ``assess`` reads the pair on the section
    as ``(admissible, worst demand/capacity)``.

    Of the pairs that are admissible and carry both moments, the one with the
    least steel. If none does, the admissible pair that comes closest, and only
    failing that the closest of all; the check then reports what is missing.
    """
    bots: list = list(bot_visited.values()) or [None]
    tops: list = list(top_visited.values()) or [None]
    zero = 0 * (cm**2)
    best_key: Optional[tuple] = None
    best_pair: tuple[Any, Any] = (None, None)
    for bot in bots:
        for top in tops:
            _apply_pair(self, bot, top)
            within, worst = assess()
            steel: Quantity = (bot.get("total_as", zero) if bot is not None else zero) + (
                top.get("total_as", zero) if top is not None else zero
            )
            passes = within and worst <= 1.0
            key = (0, float(steel.to("cm**2").magnitude)) if passes else (1 if within else 2, worst)
            if best_key is None or key < best_key:
                best_key, best_pair = key, (bot, top)
    _apply_pair(self, *best_pair)
    return best_pair


def _apply_pair(self: "RectangularBeam", bot: Any, top: Any) -> None:
    """Put a pair of layouts on the section; ``None`` leaves a face bare."""
    if bot is not None:
        self._apply_longitudinal_design_bot(bot)
    if top is not None:
        self._apply_longitudinal_design_top(top)
    else:
        self._clear_top_longitudinal()


def _run_flexure_design(
    self: "RectangularBeam",
    max_M_y_bot: Quantity,
    max_M_y_top: Quantity,
    required_areas: _RequiredAreas,
    capacity: _Capacity,
    admissible: Optional[_Admissible] = None,
) -> None:
    """Design the longitudinal reinforcement of ``self`` for the two limiting
    moments, using ``required_areas``, ``capacity`` and ``admissible`` as the
    design-code hooks. A code with no limit to hold the layout to passes no
    ``admissible``, and every layout is admissible.

    Implements 'governing-face + reconciliation' so that each face's final
    layout covers tension on that face OR compression from the opposite face,
    whichever is larger.

    The mechanical cover is solved by a Picard (fixed-point) iteration. The
    loop is bounded by ``_MAX_FLEXURE_ITERATIONS`` and stops early on a cycle,
    when a pair of layouts it has already been through comes back. However it
    ends, the layout it leaves is then verified -- strength and admissibility --
    and if it fails, every pair of the layouts visited on each face is tried
    and the best one kept (see :func:`_best_visited_pair`).

    Leaves the chosen layout applied to the section; returns nothing.
    """

    # --- helpers -----------------------------------------------------------------
    # The ranked table behind every row the search handed back, by the row's
    # fingerprint, per face. A row is only ever applied as the head of its own
    # table, so once the face is settled the table of the row on it is the list
    # of alternatives -- run with the mechanical cover the design finished on.
    tables: Dict[str, Dict[tuple, Any]] = {"bot": {}, "top": {}}
    infeasible: Dict[str, bool] = {"bot": False, "top": False}
    # The faces some load puts in tension: a positive moment pulls the bottom,
    # a negative one the top. Only their bars are held to the crack-control
    # cap of §24.3.2 (see :meth:`Rebar.longitudinal_rebar`).
    pulled = {"bot": max_M_y_bot > 0 * kNm, "top": max_M_y_top < 0 * kNm}

    def _design_longitudinal_for_area(A_req: Quantity, A_max: Any, mech_cover: Quantity, face: str) -> Any:
        """Run discrete design for a target area and return best_design dict, or
        None if the rebar designer cannot fit any combination in the section
        geometry (RebarDesignInfeasibleError). Callers must handle the None
        result — preserves the public contract that design_flexure never
        crashes, delegating the "insufficient section" report to check_flexure
        via DCR>1, and to the ``bars_do_not_fit`` warning."""
        rebar = self._create_rebar_designer()
        _ = rebar.longitudinal_rebar(A_req, A_max, mech_cover, face, tension=pulled[face])
        try:
            best = rebar.longitudinal_rebar_design
        except RebarDesignInfeasibleError:
            infeasible[face] = True
            return None
        infeasible[face] = False
        tables[face][_rebar_design_fingerprint(best)] = getattr(rebar, "_long_combos_df", None)
        return best

    def _design_tension_face(A_req: Quantity, demand: _FaceDemand, mech_cover: Quantity, face: str) -> Any:
        """Discrete design of a tension face, capped at ``A_s_max`` while that
        is enough.

        The cap is ``A_s_max`` only while the request is under it. The
        catalogue can leave no layout between the two: in a 15 cm web the bars
        that fit go from 4Ø12 = 4.52 cm² straight to 2Ø16 + 2Ø12 = 6.28 cm², so
        a face asking for 5.06 cm² under a 5.79 cm² cap fell back to the 4.52
        and missed the moment. Where the code lets compression steel extend the
        cap, the smallest layout that covers the request is taken instead, and
        :func:`_compression_for` asks the opposite face for the steel that
        keeps it tension-controlled.

        Where compression steel does not extend it -- the 4 % of EN 1992-1-1
        §9.2.1.1(3), or an ACI 318-19 / CIRSOC 201-25 section whose compression
        bars would sit too close to the neutral axis to help -- A_s_max is a
        hard cap even past the request: the face gets the most it is allowed,
        not a layout the check would reject.
        """
        extendable = demand.compression_per_excess is not None
        A_cap = demand.A_s_max if (A_req <= demand.A_s_max or not extendable) else None
        row = _design_longitudinal_for_area(A_req, A_cap, mech_cover, face)
        short = row is not None and row.get("total_as", 0 * (cm**2)) < A_req
        if short and A_cap is not None and extendable:
            infeasible_before, tables_before = infeasible[face], dict(tables[face])
            uncapped = _design_longitudinal_for_area(A_req, None, mech_cover, face)
            if uncapped is not None and uncapped.get("total_as", 0 * (cm**2)) >= A_req:
                return uncapped
            # Nothing fits past the cap either: keep the capped layout and the
            # table its options are read from.
            infeasible[face], tables[face] = infeasible_before, tables_before
        return row

    def _compression_for(row: Any, demand: _FaceDemand) -> Quantity:
        """Compression steel the tension steel of ``row`` needs on the opposite
        face to stay within the ductility limit -- zero up to ``A_s_max``.

        ACI 318-19 §9.3.3.1 / CIRSOC 201-25 §9.3.3.1 with Table 21.2.2: past
        A_s_max the section stays tension-controlled only while
        A_s <= A_s_max + A_s' * f_s' / f_y.

        Sized for the steel ``row`` places, not for the area the moment asked
        for. The two differ by whatever the bars round up by, and the hook's
        own compression steel only covers the second: a face given more than
        it asked for was left past the limit, and a check that capped it at
        phi = 0.90 credited it with a strength it does not have.
        """
        if row is None or demand.compression_per_excess is None:
            return 0 * (cm**2)
        excess = row.get("total_as", 0 * (cm**2)) - demand.A_s_max
        if excess <= 0 * (cm**2):
            return 0 * (cm**2)
        return (excess * demand.compression_per_excess).to(demand.A_s_max.units)

    M_top_demand: Quantity = abs(max_M_y_top.to("kN*m"))
    settings = self.settings
    assert settings is not None

    def _mech_cover(row: Any) -> Optional[Quantity]:
        """Cover to the centroid of the bars of ``row``, as the section computes it.

        Read off the row without applying it: the bars of each layer sit at
        half their diameter from the stirrup, and the second layer
        ``layers_spacing`` past the thicker bar of the first. ``None`` for a
        row with no bars.
        """
        zero_mm = 0 * mm

        def diameter(key: str) -> Quantity:
            value = row.get(key)
            return zero_mm if value is None else value

        d1, d2, d3, d4 = (diameter(k) for k in ("d_b1", "d_b2", "d_b3", "d_b4"))
        n1, n2, n3, n4 = (int(row.get(k, 0)) for k in ("n_1", "n_2", "n_3", "n_4"))
        second: Quantity = max(d1, d2) + settings.layers_spacing
        groups: Tuple[Tuple[int, Quantity, Quantity], ...] = (
            (n1, d1, d1 / 2),
            (n2, d2, d2 / 2),
            (n3, d3, second + d3 / 2),
            (n4, d4, second + d4 / 2),
        )
        area: Quantity = 0 * mm**2
        moment: Quantity = 0 * mm**3
        for n, db, y in groups:
            area = area + n * db**2
            moment = moment + n * db**2 * y
        if area.magnitude == 0:
            return None
        centroid: Quantity = moment / area
        return self.c_c + self._stirrup_d_b + centroid

    def _compression_need(tension_face: str, tension_row: Any, d_prime: Quantity) -> Quantity:
        """Compression steel the tension steel of ``tension_row`` needs opposite
        it, with the compression bars at depth ``d_prime``.

        The larger of what the moment asks for and what the steel placed asks
        for, both at that depth: the deeper the compression bars, the less
        stress they reach at the ductility limit and the more of them it
        takes.
        """
        cover = _mech_cover(tension_row)
        d_t: Quantity = self.height - (cover if cover is not None else self._c_mec_bot)
        if tension_face == "bot":
            demand = required_areas("bot", max_M_y_bot, d_t, d_prime)
        else:
            demand = required_areas("top", M_top_demand, d_t, d_prime)
        return max(demand.A_s_compression, _compression_for(tension_row, demand))

    def _design_compression_face(face: str, own_req: Quantity, tension_face: str, tension_row: Any) -> Any:
        """Lay out ``face`` as the compression steel of ``tension_row``.

        An area is not enough to choose compression steel by: two layers of
        thin bars sit deeper than one layer of thick ones, reach less stress,
        and need more area -- which, laid out again in two layers, sits deeper
        still. Left to the area alone, a 15x30 section went from 4Ø12 to
        2Ø25 + 2Ø25 on top and never carried its moment, while 2Ø25 in one
        layer does. So every layout the selector ranks is read at its own
        depth, in order, and the first one that makes the tension face carry
        its moment -- as the code verifies it, its limits included -- is taken.
        The ones that cover the area needed at their depth are tried first;
        only if none of them resists are a few of the rest tried too, since a
        code's sizing formula need not agree with its capacity to the last
        square centimetre (EN 1992-2004's asks far more of deep compression
        bars than its resistance needs). Failing all that, the first that
        covers the area; failing that too, the one that comes closest to it.
        """
        shallowest = self.c_c + self._stirrup_d_b + settings.minimum_longitudinal_diameter / 2
        request = max(own_req, _compression_need(tension_face, tension_row, shallowest))
        rebar = self._create_rebar_designer()
        rebar.longitudinal_rebar(request, None, shallowest, face, tension=pulled[face])
        try:
            rebar.longitudinal_rebar_design
        except RebarDesignInfeasibleError:
            infeasible[face] = True
            return None
        infeasible[face] = False
        table = rebar._long_combos_df
        M_t = max_M_y_bot if tension_face == "bot" else M_top_demand
        apply = self._apply_longitudinal_design_top if face == "top" else self._apply_longitudinal_design_bot

        def verdict(row: Any) -> tuple[bool, float]:
            """(admissible, demand/capacity) of the tension face with ``row`` opposite."""
            apply(row)
            M_R = capacity(tension_face, M_t)
            within = admissible is None or admissible(tension_face)
            ratio = float((M_t / M_R).to("dimensionless").magnitude) if M_R.magnitude > 0 else float("inf")
            return within, ratio

        covering: list = []
        short: list = []
        closest: Optional[Any] = None
        closest_margin: Optional[Quantity] = None
        for index, row in table.iterrows():
            cover = _mech_cover(row)
            need = own_req if cover is None else max(own_req, _compression_need(tension_face, tension_row, cover))
            margin: Quantity = row.get("total_as") - need
            if margin >= 0 * (cm**2):
                covering.append(index)
            else:
                short.append(index)
                if closest_margin is None or margin > closest_margin:
                    closest, closest_margin = index, margin
        # Covering ones first, then a few short ones; the first that passes
        # wins. If none does, the one that came closest -- within the limits
        # before past them -- so a section too small for its moment is left
        # with the best the design saw, not the first it tried.
        resisting: Optional[Any] = None
        nearest: Optional[Any] = None
        nearest_key: Optional[tuple] = None
        for index in covering[:_MAX_COMPRESSION_CANDIDATES] + short[:_MAX_SHORT_CANDIDATES]:
            within, ratio = verdict(table.loc[index])
            if within and ratio <= 1.0:
                resisting = index
                break
            key = (0 if within else 1, ratio)
            if nearest_key is None or key < nearest_key:
                nearest, nearest_key = index, key
        candidates = (resisting, nearest, covering[0] if covering else None, closest, table.index[0])
        pick = next(i for i in candidates if i is not None)
        # The alternatives are the table the row heads, so it goes first.
        table = pd.concat([table.loc[[pick]], table.drop(pick)]).reset_index(drop=True)
        best = table.iloc[0]
        tables[face][_rebar_design_fingerprint(best)] = table
        return best

    # --- initial guesses ----------------------------------------------------------
    # Mechanical cover = clear cover + stirrup diameter + the distance from the
    # stirrup to the centroid of the bars. ``c_c`` is the clear cover to the
    # stirrup, which the code in force fixes — ACI 318-19 Table 20.5.1.3.1 /
    # CIRSOC 201-25 Tabla 20.5.1.3.1 — and which Mento takes as given rather
    # than checking. The 1 cm is only a starting guess for the centroid,
    # replaced by the real one at the end of the first iteration; no clause
    # writes it.
    rec_mec = self.c_c + self._stirrup_d_b + 1 * cm  # bottom mechanical cover to centroid (initial)
    d_prima = self.c_c + self._stirrup_d_b + 1 * cm  # top mechanical cover to centroid (initial)

    tol = 0.01 * cm
    Err: Quantity = 2 * tol

    # Cycle detection — store the layout payload for each fingerprint seen.
    # Using a dict preserves insertion order (3.7+) so we can recover the full
    # cycle later if needed for diagnostics.
    bot_visited: Dict[tuple, dict] = {}
    top_visited: Dict[tuple, dict] = {}
    pairs_visited: set = set()
    cycled = False

    for _iteration_count in range(1, _MAX_FLEXURE_ITERATIONS + 1):
        # Effective depths for this iteration
        d = self.height - rec_mec

        # --- bottom tension case (positive moment on bottom face) ----------------
        demand_bot = required_areas("bot", max_M_y_bot, d, d_prima)
        A_s_final_bot_Positive_M = demand_bot.A_s_tension  # tension req. on bottom
        A_s_comp_top = demand_bot.A_s_compression  # compression req. on top from bottom moment

        # init in case no negative moment branch runs
        A_s_comp_bot = 0 * (cm**2)
        A_s_final_top_Negative_M = 0 * (cm**2)
        self._A_s_top = A_s_comp_top
        demand_top: Optional[_FaceDemand] = None

        # --- top tension case (negative moment on top face) ----------------------
        if max_M_y_top < 0:
            demand_top = required_areas(
                "top",
                abs(max_M_y_top.to("kN*m")),
                self.height - d_prima,
                rec_mec,
            )
            A_s_final_top_Negative_M = demand_top.A_s_tension  # tension req. on top
            A_s_comp_bot = demand_top.A_s_compression  # compression req. on bottom from top moment

        # Governing areas on each face (tension on the face vs. opposite-face compression)
        A_req_bot = max(A_s_final_bot_Positive_M, A_s_comp_bot)
        A_req_top = max(A_s_comp_top, A_s_final_top_Negative_M)

        self._A_s_bot = A_req_bot
        self._A_s_top = A_req_top

        # --- Discrete design for each face (independent first pass) ---------------
        # The cap handed to the selector is A_s_max, the tension-controlled
        # limit — ACI 318-19 §9.3.3.1 with Table 21.2.2 / CIRSOC 201-25
        # §9.3.3.1 with Tabla 21.2.2, the same limit in both. Past it the cap
        # is dropped (``None``), because the area being asked for came out of
        # the code hook itself — the tension steel of a doubly reinforced
        # couple, or the compression the opposite face needs — and capping it
        # here would only leave the selector with nothing to fit.
        if A_req_bot >= 0 * (cm**2):
            self.flexure_design_results_bot = _design_tension_face(A_req_bot, demand_bot, self._c_mec_bot, "bot")

        self.flexure_design_results_top = None
        if A_req_top >= 0 * (cm**2):
            if demand_top is not None:
                self.flexure_design_results_top = _design_tension_face(A_req_top, demand_top, self._c_mec_top, "top")
            else:
                # No combination pulls the top: what it is asked for is the
                # compression the bottom needs, which nothing caps. This read
                # the top's A_s,max off the section, and no round of this
                # design writes it: it was whatever the last reporting check
                # left there, so a design redone with its final stirrup
                # searched the top under a cap its first round never had.
                self.flexure_design_results_top = _design_longitudinal_for_area(A_req_top, None, self._c_mec_top, "top")

        # --- Apply both faces (hard overwrite) -----------------------------------
        if self.flexure_design_results_bot is not None:
            self._apply_longitudinal_design_bot(self.flexure_design_results_bot)
        if self.flexure_design_results_top is not None:
            self._apply_longitudinal_design_top(self.flexure_design_results_top)
        else:
            self._clear_top_longitudinal()

        # --- Reconciliation (override only if opposite-face compression governs) -
        A_prov_bot = (
            self.flexure_design_results_bot.get("total_as", 0 * (cm**2))
            if self.flexure_design_results_bot is not None
            else 0 * (cm**2)
        )
        A_prov_top = (
            self.flexure_design_results_top.get("total_as", 0 * (cm**2))
            if self.flexure_design_results_top is not None
            else 0 * (cm**2)
        )

        # Each face owes the other's tension steel the compression it needs:
        # for the steel placed, not the area asked for, and read at the depth
        # the face's own bars sit at. A face that falls short is laid out again
        # as compression steel (see _design_compression_face).
        shallowest = self.c_c + self._stirrup_d_b + settings.minimum_longitudinal_diameter / 2
        if max_M_y_bot > 0 * kNm and self.flexure_design_results_bot is not None:
            cover_top = (
                _mech_cover(self.flexure_design_results_top) if self.flexure_design_results_top is not None else None
            )
            A_s_comp_top = _compression_need("bot", self.flexure_design_results_bot, cover_top or shallowest)
            if A_s_comp_top > A_prov_top:
                row = _design_compression_face("top", A_s_final_top_Negative_M, "bot", self.flexure_design_results_bot)
                if row is not None:
                    self.flexure_design_results_top = row
                    self._apply_longitudinal_design_top(row)
                    A_prov_top = row.get("total_as", A_prov_top)
                    cover = _mech_cover(row)
                    if cover is not None:
                        A_s_comp_top = _compression_need("bot", self.flexure_design_results_bot, cover)

        if demand_top is not None and self.flexure_design_results_top is not None:
            cover_bot = (
                _mech_cover(self.flexure_design_results_bot) if self.flexure_design_results_bot is not None else None
            )
            A_s_comp_bot = _compression_need("top", self.flexure_design_results_top, cover_bot or shallowest)
            if A_s_comp_bot > A_prov_bot:
                row = _design_compression_face("bot", A_s_final_bot_Positive_M, "top", self.flexure_design_results_top)
                if row is not None:
                    self.flexure_design_results_bot = row
                    self._apply_longitudinal_design_bot(row)
                    A_prov_bot = row.get("total_as", A_prov_bot)
                    cover = _mech_cover(row)
                    if cover is not None:
                        A_s_comp_bot = _compression_need("top", self.flexure_design_results_top, cover)

        # --- Update geometry (centroids) for next iteration ----------------------
        c_mec_calc = self.c_c + self._stirrup_d_b + self._bot_rebar_centroid
        # If there is any top steel, use its centroid; otherwise keep previous d_prima
        has_top = (self.flexure_design_results_top is not None) and (
            int(self.flexure_design_results_top.get("n_1", 0))
            + int(self.flexure_design_results_top.get("n_2", 0))
            + int(self.flexure_design_results_top.get("n_3", 0))
            + int(self.flexure_design_results_top.get("n_4", 0))
            > 0
        )
        d_prima_calc = self.c_c + self._stirrup_d_b + self._top_rebar_centroid if has_top else d_prima

        # --- Cycle detection ------------------------------------------------------
        # Bottom and top are coupled: bottom's layout drives `rec_mec`, which is
        # fed back as the compression-side depth of the top design (and
        # vice-versa), so the state of the loop is the PAIR of layouts. Only a
        # pair seen before is a limit cycle. One face repeating its layout is
        # not: it is what a face that has settled does while the other is still
        # moving -- a top face going from 2Ø25 to 2Ø20 + 2Ø16 over a bottom
        # that stays at 2Ø20 -- and stopping there left the top on a layout
        # designed for a depth it no longer had. On a real cycle we exit, and the
        # final verification below picks among the layouts visited.
        fp_bot = fp_top = None
        if self.flexure_design_results_bot is not None:
            fp_bot = _rebar_design_fingerprint(self.flexure_design_results_bot)
            bot_visited.setdefault(fp_bot, dict(self.flexure_design_results_bot))
        if self.flexure_design_results_top is not None:
            fp_top = _rebar_design_fingerprint(self.flexure_design_results_top)
            top_visited.setdefault(fp_top, dict(self.flexure_design_results_top))
        if (fp_bot, fp_top) in pairs_visited:
            cycled = True
        pairs_visited.add((fp_bot, fp_top))

        # --- Convergence update ---------------------------------------------------
        Err = max(abs(c_mec_calc - rec_mec), abs(d_prima_calc - d_prima))
        rec_mec = c_mec_calc
        d_prima = d_prima_calc

        if Err < tol:
            break
        if cycled:
            break

    # --- Final verification ----------------------------------------------------
    # Whether the loop converged, cycled or ran out, the layout it leaves is
    # not guaranteed to resist: the Picard iteration settles the centroids,
    # not the strength. So it is always verified against the design-strength
    # requirement -- ACI 318-19 §9.5.1.1(a) / CIRSOC 201-25 §9.5.1.1(a),
    # phi*Mn >= Mu, and M_Rd >= M_Ed under EN 1992-2004, whichever the
    # ``capacity`` callback speaks -- and against the limits the code puts on
    # the reinforcement (``admissible``). If it fails, every pair of the layouts visited on
    # each face is tried and the best kept: the faces are coupled -- the
    # compression steel of one sets what the other can carry -- so a face is
    # never chosen on its own against whatever the other happens to hold.
    demands = []
    if max_M_y_bot > 0 * kNm:
        demands.append(("bot", max_M_y_bot))
    if max_M_y_top < 0 * kNm:
        demands.append(("top", M_top_demand))

    def _assess() -> tuple[bool, float]:
        """(admissible, worst demand/capacity) of the layout on the section."""
        worst = 0.0
        within = True
        for face, M in demands:
            M_R = capacity(face, M)
            ratio = float((M / M_R).to("dimensionless").magnitude) if M_R.magnitude > 0 else float("inf")
            worst = max(worst, ratio)
            if admissible is not None and not admissible(face):
                within = False
        return within, worst

    def _layout_resists() -> bool:
        """Does the layout on the section carry both moments, as the code allows?"""
        within, worst = _assess()
        return within and worst <= 1.0

    if not _layout_resists() and (bot_visited or top_visited):
        chosen = _best_visited_pair(self, bot_visited, top_visited, _assess)
        self.flexure_design_results_bot, self.flexure_design_results_top = chosen

    # Both faces are settled. An element whose faces are detailed as one -- a
    # footing mat -- gets the last word here, after the verification above and
    # never before it, so that what it starts from is a layout already known to
    # work. It may re-select the bars rather than only the spacing, so it is
    # handed the means to verify its own choice.
    self._finalize_longitudinal_design(A_req_bot, A_req_top, _layout_resists)

    # A layout that still fails is the most that fits, not a design: no bars
    # that fit the width carry the moment, or keep it tension-controlled. The
    # section says so -- the warning ``As_below_required`` -- on the face that
    # fell short, for as long as it carries what the design left on it. A
    # face is short when it holds less than it was asked for; if neither does,
    # the shortfall is the tension face's, whose requirement is then the most
    # it can take and stay tension-controlled.
    self._short_faces = {}
    if not _layout_resists():
        needed = {"bot": max(A_req_bot, A_s_comp_bot), "top": max(A_req_top, A_s_comp_top)}
        placed = {"bot": self._A_s_bot, "top": self._A_s_top}
        for face in ("bot", "top"):
            if placed[face] < needed[face]:
                self._short_faces[face] = (needed[face], placed[face])
        if not self._short_faces and demands:
            face = demands[0][0]
            self._short_faces[face] = (needed[face], placed[face])

    # The alternatives of each face, headed by what it carries now.
    for face, suffix, row in (
        ("bot", "b", self.flexure_design_results_bot),
        ("top", "t", self.flexure_design_results_top),
    ):
        table = tables[face].get(_rebar_design_fingerprint(row)) if row is not None else None
        self._record_longitudinal_options(suffix, row, table)
        if infeasible[face]:
            self._infeasible_faces.add(face)
        else:
            self._infeasible_faces.discard(face)
