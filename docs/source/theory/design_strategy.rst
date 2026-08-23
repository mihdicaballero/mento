Design Strategy
===============

Checking a section is a direct calculation: given a layout, the code equations return
a capacity. **Designing** one is not — the required steel area depends on the
effective depth, the effective depth depends on where the bars end up, and the bars
come from a discrete catalogue rather than a continuum.

This page describes how mento resolves that. The strategy lives in
``mento/codes/flexure_design.py`` and is **shared by ACI 318-19 and EN 1992-1-1**:
the codes differ in their equations, not in how a layout is searched for. Each code
supplies two hooks and the driver does the rest.

The two hooks
-------------

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Hook
     - Contract
   * - ``required_areas(face, M, d, d')``
     - Given a moment on a face and the two depths, return a ``FaceDemand``:
       :math:`A_{s,min}`, :math:`A_{s,max}`, the tension steel required on that
       face, and the compression steel this moment requires on the **opposite**
       face.
   * - ``capacity(face, M_demand)``
     - Return the design resisting moment of the layout currently applied — with
       the code's safety format, so :math:`\phi M_n` for ACI and :math:`M_{Rd}`
       for EN.

Everything else — the iteration, the bar selection, the reconciliation, the
verification — is code-agnostic.

The loop
--------

.. code-block:: text

   rec_mec, d'  <-  first guess: c_c + stirrup diameter + 1 cm
   repeat, up to 30 times:

       (1)  d = h - rec_mec
       (2)  ask the code for the required areas on each face
       (3)  take the governing area per face
       (4)  turn each area into a discrete bar layout
       (5)  apply both layouts to the section
       (6)  reconcile: does opposite-face compression govern?
       (7)  read the real bar centroids -> new rec_mec, d'
       (8)  converged (< 0.1 mm) or cycling?  -> exit

   final:  does the active layout actually resist the demand?
           if not, pick the best among the layouts visited

Step 1 — the fixed point
^^^^^^^^^^^^^^^^^^^^^^^^

The mechanical cover is solved by **Picard iteration**: guess it, design with it,
read back what the resulting bars imply, repeat. It converges quickly because each
round changes :math:`d` by only a few millimetres, and the tolerance is 0.1 mm on
both covers.

Step 3 — governing face
^^^^^^^^^^^^^^^^^^^^^^^

Each face has to cover **either** the tension its own moment demands **or** the
compression the opposite face's moment demands from it, whichever is larger:

.. math::

   A_{req,bot} = \max\left(A_{s,tension}^{M^+},\ A_{s,comp}^{M^-}\right)
   \qquad
   A_{req,top} = \max\left(A_{s,comp}^{M^+},\ A_{s,tension}^{M^-}\right)

A beam with a large positive moment at midspan and a large negative moment at the
support gets both faces sized for the worse of the two roles.

Step 4 — discrete bar selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The rebar designer searches bar-count and diameter combinations from the catalogue
that satisfy :math:`A_{prov} \ge A_{req}`, respect the clear-spacing and layer
limits, and are capped at :math:`A_{s,max}`. It returns the most economical fit.

If nothing fits — a narrow section with a large demand — it raises
``RebarDesignInfeasibleError``, which the driver catches and turns into "no layout
for this face". The design **never** raises; the insufficiency surfaces downstream as
``DCR > 1``.

Step 6 — reconciliation
^^^^^^^^^^^^^^^^^^^^^^^

The two faces are designed independently first, then checked against each other. If
the compression the top moment demands on the bottom face exceeds what the bottom
layout actually provides, the bottom is redesigned for that larger area — and
symmetrically for the top.

Step 8 — cycle detection
^^^^^^^^^^^^^^^^^^^^^^^^

Discrete bar selection can oscillate: layout A implies a centroid that calls for
layout B, whose centroid calls back for layout A, forever. The driver fingerprints
every layout it applies — counts and diameters, not the derived centroid — and stops
as soon as a fingerprint repeats.

A single flag covers both faces on purpose. The faces are coupled: the bottom
layout drives ``rec_mec``, which is fed back as the compression-side depth of the top
design, and vice versa. If either face repeats a layout, its centroid is oscillating
periodically, which forces the other face into a limit cycle too. Detecting
recurrence on one face is enough evidence to stop.

Final verification
------------------

This is the step that makes the whole thing trustworthy, and the one that is easy to
leave out.

Converging on the mechanical cover proves the **geometry** is self-consistent. It
does not prove the layout resists the demand. The required-area formula and the
capacity formula are different equations, and any gap between them — a discrete bar
jump, a lever arm evaluated at a slightly different depth, a branch boundary — can
leave a converged layout short.

So after the loop, whatever the exit reason, the driver evaluates the real capacity:

.. math::

   \text{if } \; \text{capacity(face)} < M_{demand}
   \;\Longrightarrow\; \text{pick a safer layout}

The replacement comes from ``select_safe_design``, over the layouts visited during
the loop:

1. Among candidates that satisfy :math:`\text{capacity} \ge M_{demand}`, take the one
   with the **smallest** :math:`A_s` — the most economical layout that works, not the
   strongest one.
2. If none satisfies it, take the one with the **largest** capacity: the closest to
   passing. ``check_flexure`` then reports ``DCR > 1``.

Each candidate is re-applied to the section before being evaluated, so its own
centroid — and therefore its own :math:`d` — is used, not the one from whichever
iteration happened to produce it.

.. note::

   Without this step, ``design_flexure`` can hand back a layout that ``check_flexure``
   immediately rejects. Over a sweep of 288 EN load combinations, 29 designs failed
   their own check before the two codes shared this driver; none do now.

What the driver guarantees
--------------------------

- **Termination.** Bounded at 30 iterations, with early exit on convergence or on a
  detected cycle.
- **Never raises.** An unfittable section degrades to ``DCR > 1``, so the summary,
  plots, shear design and reports keep working on an underdesigned member.
- **Self-consistency.** A layout returned by ``design_flexure`` passes
  ``check_flexure``, unless no layout in the catalogue can.
- **Economy.** Where several layouts work, the smallest area wins.

Where the codes still differ
----------------------------

The driver deliberately does not paper over genuine differences:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Aspect
     - ACI 318-19
     - EN 1992-1-1
   * - Safety format
     - :math:`\phi M_n \ge M_u`
     - :math:`M_{Rd} \ge M_{Ed}`, partial factors on materials
   * - Compression-steel branch
     - :math:`A_s > A_{s,max}`
     - :math:`x_{eff} > x_{eff,lim}`
   * - Minimum steel at :math:`M = 0`
     - Geometric minimum :math:`1.8\text{‰}\,b\,h`
     - :math:`\rho_{min}` applies unconditionally
   * - Extra reported quantities
     - :math:`c/d`, 4/3-rule flag
     - —

Those live entirely inside the two hooks.

Validation
----------

.. list-table::
   :header-rows: 1
   :widths: 34 40 26

   * - Behaviour
     - Test
     - Notes
   * - Smallest passing layout wins
     - ``test_select_safe_design_prefers_smallest_passing_layout``
     - Also covers the fallback to the largest capacity
   * - Cycle exits to a passing visited layout
     - ``test_design_flexure_ACI_318_19_cycle_adopts_passing_visited_layout``
     - Traced against the real loop
   * - Never raises, ACI
     - ``test_design_flexure_rebar_infeasible_does_not_crash``
     - 15×15 cm section, unfittable demand
   * - Never raises, EN
     - ``test_design_flexure_EN_1992_2004_infeasible_does_not_crash``
     - Same contract, gained by sharing the driver
   * - Design output passes its own check
     - ``test_design_flexure_EN_1992_2004_layout_passes_its_own_check``
     - The regression that motivated the final verification
   * - Opposite-face reconciliation
     - ``test_design_flexure_ACI_318_19_compression_bottom_exceeds_provided_bottom``
     - Step 6
   * - Doubly reinforced design, both signs
     - ``test_design_flexure_ACI_318_19_negative_moment_doubly_reinforced``,
       ``test_design_flexure_ACI_318_19_large_negative_moment_doubly_reinforced``
     - Verified through :math:`\phi M_n \ge M_u`
