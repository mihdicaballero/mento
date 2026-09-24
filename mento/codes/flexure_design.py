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
# EN 1992-2004 -- which is the comparison in ``_select_safe_design`` and in the
# final verification of ``_run_flexure_design``.
#
# Kept as comments rather than a module docstring on purpose: Sphinx autodoc
# publishes docstrings, and this is implementation detail, not reference
# material for users.
# ---------------------------------------------------------------------------

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, Dict, Optional

from mento.units import Quantity

from mento.rebar import RebarDesignInfeasibleError
from mento.units import cm, kNm, mm

if TYPE_CHECKING:
    from ..beam import RectangularBeam


_MAX_FLEXURE_ITERATIONS = 30  # safety net for slow divergence without cycling


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


def _select_safe_design(
    self: "RectangularBeam",
    candidate_designs: list,
    M_demand: Any,
    face: str,
    capacity: _Capacity,
) -> Any:
    """Among a set of candidate rebar designs — those visited during the
    Picard (fixed-point iteration) loop — return the most appropriate one for
    the given face.

    The pass/fail line is the design-strength requirement itself — ACI 318-19
    §9.5.1.1(a) / CIRSOC 201-25 §9.5.1.1(a), phi*Mn >= Mu, and M_Rd >= M_Ed
    under EN 1992-2004. It is read through the ``capacity`` callback, so each
    code applies its own safety format without this function knowing which.

    Selection priority:

    1. Among candidates whose actual resisting moment satisfies
       ``M_Rd >= M_demand``, return the one with the smallest As (most
       economical valid layout).
    2. If no candidate satisfies the check, return the one with the largest
       resisting moment (closest to passing). Downstream ``check_flexure`` will
       surface DCR > 1 so the user is aware that the section is insufficient.

    Never raises: keeping ``design_flexure`` total preserves the public
    contract — summary, plot, shear design and check must keep working even
    when the section is underdesigned.

    Parameters
    ----------
    candidate_designs : list of dict
        Each dict is a rebar_designer payload (keys n_1..n_4, d_b1..d_b4,
        total_as, ...).
    M_demand : Quantity
        Required design moment on the face (positive magnitude).
    face : {"bot", "top"}
        Which face the layout governs.
    capacity : _Capacity
        Code-specific evaluation of the resisting moment.
    """
    M_demand_abs = abs(M_demand.to("kN*m"))

    evaluated: list = []  # tuples of (As_provided, M_Rd, design_dict)
    for design in candidate_designs:
        if face == "bot":
            self._apply_longitudinal_design_bot(design)
        else:
            self._apply_longitudinal_design_top(design)
        # Recompute the capacity with the just-applied layout (centroid included)
        M_Rd = capacity(face, M_demand_abs)
        evaluated.append(
            (
                design.get("total_as", 0 * (cm**2)),
                M_Rd.to("kN*m"),
                design,
            )
        )

    # Primary criterion: candidates that satisfy the check
    passing = [e for e in evaluated if e[1] >= M_demand_abs]
    if passing:
        passing.sort(key=lambda e: e[0])  # smallest As first
        return passing[0][2]

    # Fallback: no candidate passes — pick the one with the largest resisting
    # moment (closest to passing). The downstream check will report DCR > 1.
    evaluated.sort(key=lambda e: -e[1])  # largest capacity first
    return evaluated[0][2]


def _run_flexure_design(
    self: "RectangularBeam",
    max_M_y_bot: Quantity,
    max_M_y_top: Quantity,
    required_areas: _RequiredAreas,
    capacity: _Capacity,
) -> None:
    """Design the longitudinal reinforcement of ``self`` for the two limiting
    moments, using ``required_areas`` and ``capacity`` as the design-code hooks.

    Implements 'governing-face + reconciliation' so that each face's final
    layout covers tension on that face OR compression from the opposite face,
    whichever is larger.

    The mechanical cover is solved by a Picard (fixed-point) iteration. The
    loop is bounded by ``_MAX_FLEXURE_ITERATIONS`` and includes per-face cycle
    detection: if the same layout fingerprint reappears, the loop exits and a
    safe layout is picked among the visited candidates by re-evaluating the
    resisting moment (see :func:`_select_safe_design`). This avoids infinite
    oscillation when two or more discrete layouts alternate without converging
    in the strict tolerance.

    Leaves the chosen layout applied to the section; returns nothing.
    """

    # --- helpers -----------------------------------------------------------------
    # The ranked table behind every row the search handed back, by the row's
    # fingerprint, per face. A row is only ever applied as the head of its own
    # table, so once the face is settled the table of the row on it is the list
    # of alternatives -- run with the mechanical cover the design finished on.
    tables: Dict[str, Dict[tuple, Any]] = {"bot": {}, "top": {}}
    infeasible: Dict[str, bool] = {"bot": False, "top": False}

    def _design_longitudinal_for_area(A_req: Quantity, A_max: Any, mech_cover: Quantity, face: str) -> Any:
        """Run discrete design for a target area and return best_design dict, or
        None if the rebar designer cannot fit any combination in the section
        geometry (RebarDesignInfeasibleError). Callers must handle the None
        result — preserves the public contract that design_flexure never
        crashes, delegating the "insufficient section" report to check_flexure
        via DCR>1, and to the ``bars_do_not_fit`` warning."""
        rebar = self._create_rebar_designer()
        _ = rebar.longitudinal_rebar(A_req, A_max, mech_cover)
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
        """
        A_cap = demand.A_s_max if A_req <= demand.A_s_max else None
        row = _design_longitudinal_for_area(A_req, A_cap, mech_cover, face)
        short = row is not None and row.get("total_as", 0 * (cm**2)) < A_req
        if short and A_cap is not None and demand.compression_per_excess is not None:
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

        Only for a face the moment alone keeps under A_s_max. Once the moment
        itself asks for more, the code hook has already sized the compression
        steel for the area it asked for, and the check caps whatever the bars
        add on top of it without losing capacity; asking for compression for
        that surplus too only fed the next iteration a heavier opposite face,
        a deeper centroid and a larger demand, and the loop ran away.
        """
        if row is None or demand.compression_per_excess is None:
            return 0 * (cm**2)
        if demand.A_s_tension > demand.A_s_max:
            return 0 * (cm**2)
        excess = row.get("total_as", 0 * (cm**2)) - demand.A_s_max
        if excess <= 0 * (cm**2):
            return 0 * (cm**2)
        return (excess * demand.compression_per_excess).to(demand.A_s_max.units)

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
                A_cap_top = self._A_s_max_top if A_req_top <= self._A_s_max_top else None
                self.flexure_design_results_top = _design_longitudinal_for_area(
                    A_req_top, A_cap_top, self._c_mec_top, "top"
                )

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

        # A tension face detailed past A_s_max -- the bars it can take rarely
        # land on the area asked for -- needs compression steel for all of it,
        # not only for the area the moment asked for.
        if max_M_y_bot > 0 * kNm:
            A_s_comp_top = max(A_s_comp_top, _compression_for(self.flexure_design_results_bot, demand_bot))
        if demand_top is not None:
            A_s_comp_bot = max(A_s_comp_bot, _compression_for(self.flexure_design_results_top, demand_top))

        # If compression from top (A_s_comp_bot) exceeds what bottom provides, re-upgrade bottom
        if A_s_comp_bot > A_prov_bot:
            A_cap_bot = self._A_s_max_bot if A_s_comp_bot <= self._A_s_max_bot else None
            self.flexure_design_results_bot = _design_longitudinal_for_area(
                A_s_comp_bot, A_cap_bot, self._c_mec_bot, "bot"
            )
            if self.flexure_design_results_bot is not None:
                self._apply_longitudinal_design_bot(self.flexure_design_results_bot)
                A_prov_bot = self.flexure_design_results_bot.get("total_as", A_s_comp_bot)

        # If compression from bottom (A_s_comp_top) exceeds what top provides, re-upgrade top.
        # The outer guard implies A_s_comp_top > A_prov_top >= 0, hence A_s_comp_top > 0
        # (mirrors the bottom-face reconciliation above).
        if A_s_comp_top > A_prov_top:
            A_cap_top = self._A_s_max_top if A_s_comp_top <= self._A_s_max_top else None
            self.flexure_design_results_top = _design_longitudinal_for_area(
                A_s_comp_top, A_cap_top, self._c_mec_top, "top"
            )
            if self.flexure_design_results_top is not None:
                self._apply_longitudinal_design_top(self.flexure_design_results_top)
                A_prov_top = self.flexure_design_results_top.get("total_as", A_s_comp_top)

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
        # designed for a depth it no longer had. On a real cycle we exit and let
        # `_select_safe_design` pick the best layout among those visited.
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

    # --- Final capacity verification ----------------------------------------
    # Whether the loop converged, cycled, or hit MAX_ITER, the active layout
    # is NOT guaranteed to resist the demand. The Picard (fixed-point
    # iteration) only ensures centroid consistency, not flexural capacity. So
    # we always run a final check against the design-strength requirement —
    # ACI 318-19 §9.5.1.1(a) / CIRSOC 201-25 §9.5.1.1(a), phi*Mn >= Mu, and
    # M_Rd >= M_Ed under EN 1992-2004, whichever the ``capacity`` callback
    # speaks. If the active layout fails, we pick the
    # safest one among the visited layouts (see `_select_safe_design`). If no
    # visited layout passes either, we pick the closest one — the downstream
    # `check_flexure` will surface DCR > 1 to signal that the section is
    # insufficient.
    # `bot_visited`/`top_visited` are empty only when the rebar selector could
    # not fit a single layout on that face; there is then nothing to fall back
    # to and the active (unchanged) layout is left for check_flexure to report.
    if max_M_y_bot > 0 * kNm:
        if capacity("bot", max_M_y_bot) < max_M_y_bot and bot_visited:
            chosen_bot = _select_safe_design(self, list(bot_visited.values()), max_M_y_bot, "bot", capacity)
            self._apply_longitudinal_design_bot(chosen_bot)
            self.flexure_design_results_bot = chosen_bot

    if max_M_y_top < 0 * kNm:
        M_demand_top: Quantity = abs(max_M_y_top.to("kN*m"))
        if capacity("top", M_demand_top) < M_demand_top and top_visited:
            chosen_top = _select_safe_design(self, list(top_visited.values()), M_demand_top, "top", capacity)
            self._apply_longitudinal_design_top(chosen_top)
            self.flexure_design_results_top = chosen_top

    # Both faces are settled. An element whose faces are detailed as one -- a
    # footing mat -- gets the last word here, after the verification above and
    # never before it, so that what it starts from is a layout already known to
    # work. It may re-select the bars rather than only the spacing, so it is
    # handed the means to verify its own choice.
    def _layout_resists() -> bool:
        """Does the layout currently on the section carry both moments?"""
        if max_M_y_bot > 0 * kNm and capacity("bot", max_M_y_bot) < max_M_y_bot:
            return False
        if max_M_y_top < 0 * kNm:
            M_demand: Quantity = abs(max_M_y_top.to("kN*m"))
            if capacity("top", M_demand) < M_demand:
                return False
        return True

    self._finalize_longitudinal_design(A_req_bot, A_req_top, _layout_resists)

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
