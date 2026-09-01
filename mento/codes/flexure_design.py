"""Internal engine shared by the flexural design routines of each design code.

Private to ``mento.codes``: not part of the public API, and deliberately kept out
of the published documentation.
"""

# ---------------------------------------------------------------------------
# Why this module exists
#
# ACI 318-19 and EN 1992-2004 disagree on the equations -- the stress block, the
# safety format, the minimum reinforcement rules -- but not on the strategy that
# surrounds them:
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
# Kept as comments rather than a module docstring on purpose: Sphinx autodoc
# publishes docstrings, and this is implementation detail, not reference
# material for users.
# ---------------------------------------------------------------------------

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, Dict

from pint import Quantity

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
        Code limits for the tension face being solved.
    A_s_tension:
        Steel required on the tension face itself.
    A_s_compression:
        Steel this face's moment requires on the OPPOSITE face, i.e. the
        compression reinforcement of a doubly reinforced section. Zero when the
        section does not need it.
    """

    A_s_min: Quantity
    A_s_max: Quantity
    A_s_tension: Quantity
    A_s_compression: Quantity


# Ask the design code for the areas required on ``face`` ("bot"/"top") by a
# moment ``M`` (always a positive magnitude), given the tension-side effective
# depth ``d`` and the compression-side mechanical cover ``d_prime``.
_RequiredAreas = Callable[..., _FaceDemand]

# Resisting moment of the layout currently applied to the section, on ``face``,
# under a demand of ``M_demand`` (positive magnitude). Includes the safety
# format of the code (phi*Mn for ACI, M_Rd for EN).
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
    def _design_longitudinal_for_area(A_req: Quantity, A_max: Any, mech_cover: Quantity) -> Any:
        """Run discrete design for a target area and return best_design dict, or
        None if the rebar designer cannot fit any combination in the section
        geometry (RebarDesignInfeasibleError). Callers must handle the None
        result — preserves the public contract that design_flexure never
        crashes, delegating the "insufficient section" report to check_flexure
        via DCR>1."""
        rebar = self._create_rebar_designer()
        _ = rebar.longitudinal_rebar(A_req, A_max, mech_cover)
        try:
            return rebar.longitudinal_rebar_design
        except RebarDesignInfeasibleError:
            return None

    # --- initial guesses ----------------------------------------------------------
    rec_mec = self.c_c + self._stirrup_d_b + 1 * cm  # bottom mechanical cover to centroid (initial)
    d_prima = self.c_c + self._stirrup_d_b + 1 * cm  # top mechanical cover to centroid (initial)

    tol = 0.01 * cm
    Err = 2 * tol

    # Cycle detection — store the layout payload for each fingerprint seen.
    # Using a dict preserves insertion order (3.7+) so we can recover the full
    # cycle later if needed for diagnostics.
    bot_visited: Dict[tuple, dict] = {}
    top_visited: Dict[tuple, dict] = {}
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
        if A_req_bot >= 0 * (cm**2):
            A_cap_bot = self._A_s_max_bot if A_req_bot <= self._A_s_max_bot else None
            self.flexure_design_results_bot = _design_longitudinal_for_area(A_req_bot, A_cap_bot, self._c_mec_bot)

        self.flexure_design_results_top = None
        if A_req_top >= 0 * (cm**2):
            A_cap_top = self._A_s_max_top if A_req_top <= self._A_s_max_top else None
            self.flexure_design_results_top = _design_longitudinal_for_area(A_req_top, A_cap_top, self._c_mec_top)

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

        # If compression from top (A_s_comp_bot) exceeds what bottom provides, re-upgrade bottom
        if A_s_comp_bot > A_prov_bot:
            A_cap_bot = self._A_s_max_bot if A_s_comp_bot <= self._A_s_max_bot else None
            self.flexure_design_results_bot = _design_longitudinal_for_area(A_s_comp_bot, A_cap_bot, self._c_mec_bot)
            if self.flexure_design_results_bot is not None:
                self._apply_longitudinal_design_bot(self.flexure_design_results_bot)
                A_prov_bot = self.flexure_design_results_bot.get("total_as", A_s_comp_bot)

        # If compression from bottom (A_s_comp_top) exceeds what top provides, re-upgrade top.
        # The outer guard implies A_s_comp_top > A_prov_top >= 0, hence A_s_comp_top > 0
        # (mirrors the bottom-face reconciliation above).
        if A_s_comp_top > A_prov_top:
            A_cap_top = self._A_s_max_top if A_s_comp_top <= self._A_s_max_top else None
            self.flexure_design_results_top = _design_longitudinal_for_area(A_s_comp_top, A_cap_top, self._c_mec_top)
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
        # A single `cycled` flag covers both faces on purpose. Bottom and top are
        # coupled: bottom's layout drives `rec_mec`, which is fed back as the
        # compression-side depth of the top design (and vice-versa). If one face
        # repeats a layout it already visited, its centroid is oscillating
        # periodically, so `rec_mec`/`d_prima` are too — and that periodic input
        # forces the OTHER face into a limit cycle as well. Detecting recurrence
        # on either face is enough evidence that the whole system is in a limit
        # cycle; continuing would only waste iterations. We exit and let
        # `_select_safe_design` pick the best layout among those visited.
        if self.flexure_design_results_bot is not None:
            fp_bot = _rebar_design_fingerprint(self.flexure_design_results_bot)
            if fp_bot in bot_visited:
                cycled = True
            else:
                bot_visited[fp_bot] = dict(self.flexure_design_results_bot)
        if self.flexure_design_results_top is not None:
            fp_top = _rebar_design_fingerprint(self.flexure_design_results_top)
            if fp_top in top_visited:
                cycled = True
            else:
                top_visited[fp_top] = dict(self.flexure_design_results_top)

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
    # we always run a final check; if the active layout fails, we pick the
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

    if max_M_y_top < 0 * kNm:
        M_demand_top = abs(max_M_y_top.to("kN*m"))
        if capacity("top", M_demand_top) < M_demand_top and top_visited:
            chosen_top = _select_safe_design(self, list(top_visited.values()), M_demand_top, "top", capacity)
            self._apply_longitudinal_design_top(chosen_top)

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
            M_demand = abs(max_M_y_top.to("kN*m"))
            if capacity("top", M_demand) < M_demand:
                return False
        return True

    self._finalize_longitudinal_design(A_req_bot, A_req_top, _layout_resists)
