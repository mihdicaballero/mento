Footing
=======

``Footing`` subclasses ``OneWaySlab``, which subclasses ``RectangularBeam``.
**The strength provisions are identical** — a footing strip is solved as a wide,
shallow beam, and every equation on the :doc:`ACI <beam_aci_318_19>` and
:doc:`EN <beam_en_1992_2004>` pages applies unchanged, as does everything on the
:doc:`one_way_slab` page.

This page covers only what differs: the minimum reinforcement, and the detailing
limits that go with it.

.. warning::

   mento performs **no geotechnical calculation**. Bearing pressure, settlement,
   sliding, overturning, uplift and the plan dimensions needed to spread the load are
   outside its scope, and so is the choice of thickness. The plan size and the depth
   are decided from the soil first; ``Footing`` then designs the section that
   decision produced, for the sectional forces it is given. Punching shear is a
   separate check — see ``PunchingSlab`` in the API reference.

Why the minimum differs
-----------------------

The minimum flexural reinforcement of a beam exists to prevent one failure mode: a
lightly reinforced section that cracks and, at that instant, finds its steel unable to
carry the moment the uncracked concrete had been carrying. The failure is sudden,
which is why the minimum is written as the steel needed to exceed the cracking moment.

A member bearing on the ground cannot fail that way. The soil under it goes on
carrying it after the section cracks, so there is no sudden loss of support to guard
against. Both codes recognise this and exempt the member, then substitute a rule with
a different purpose — restraining shrinkage and temperature movement under ACI,
controlling crack widths under EN.

mento expresses the distinction as a single attribute on the section,
``Section.support``: ``"free"`` for anything spanning between supports, ``"soil"``
for a member bearing on the ground. ``Footing`` sets it to ``"soil"``; the design code
reads it and picks the clause. It is a class attribute, not a constructor argument —
it says what kind of element this is.

ACI 318-19 and CIRSOC 201-25
----------------------------

§9.6.1.1(b) exempts a member supported on the ground from the flexural minimum of
§9.6.1.2, and §13.3.1.2 sends a footing to the shrinkage and temperature
reinforcement of Table 24.4.3.2 instead. That minimum is written on the **gross**
section, so the effective depth does not enter it:

.. math::

   A_{s,min} = \rho_{st}\, b\, h
   \qquad
   \rho_{st} = \max\!\left(0.0018\,\frac{f_{y,ref}}{f_y},\ 0.0014\right)

with :math:`f_{y,ref} = 420` MPa in SI and 60 000 psi in US customary. The ratio is
0.0018 at the reference grade, scales inversely for a stronger steel, and never falls
below the 0.0014 floor of the table.

For comparison, the beam minimum this replaces is

.. math::

   A_{s,min} = \rho_{min}\, b\, d
   \qquad
   \rho_{min} = \max\!\left(\frac{0.25\sqrt{f'_c}}{f_y},\ \frac{1.4}{f_y}\right)

On a 1 m strip 600 mm deep in H25 with :math:`f_y = 420` MPa, the footing rule gives
10.8 cm² against the beam rule's 18.1 cm².

EN 1992-1-1:2004
----------------

Two rules apply and **the larger governs**. Both are asked only for a face the moment
puts in tension: a face that is not in tension is given no minimum at all, not a smaller
one.

Geometric minimum
^^^^^^^^^^^^^^^^^

The non-fragility minimum of §9.2.1.1 is halved per direction for a foundation and
written on the gross section:

.. math::

   A_{s,min} = \rho_{geo}\, b\, h
   \qquad
   \rho_{geo} =
   \begin{cases}
     0.0010 & f_{yk} \le 400\ \text{MPa} \\
     0.0009 & f_{yk} \ge 500\ \text{MPa}
   \end{cases}

Crack control — §7.3.2(2), Eq. (7.1)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Enough steel to carry, without breaking, the tension the concrete releases at the
instant it cracks:

.. math::

   A_{s,min}\,\sigma_s = k_c\, k\, f_{ct,eff}\, A_{ct}
   \qquad\Longrightarrow\qquad
   A_{s,min} = \frac{k_c\, k\, f_{ct,eff}\, A_{ct}}{\sigma_s}

with

.. list-table::
   :header-rows: 1
   :widths: 14 46 40

   * - Term
     - Value mento uses
     - Source
   * - :math:`k_c`
     - 0.40 — rectangular section in pure bending
     - §7.3.2(2)
   * - :math:`k`
     - :math:`1.00` at :math:`h \le 300` mm, :math:`0.65` at :math:`h \ge 800` mm,
       linear between, saturated outside
     - §7.3.2(2)
   * - :math:`f_{ct,eff}`
     - :math:`f_{ctm} = 0.3\,f_{ck}^{2/3}` for :math:`\le` C50/60
     - §7.3.2(2), Table 3.1
   * - :math:`A_{ct}`
     - :math:`b\,h/2` — see below
     - §7.3.2(2)
   * - :math:`\sigma_s`
     - :math:`f_{yd} = f_{yk}/\gamma_s`
     - §7.3.2(2), §2.4.2.4

Worked reference case — a 1 m strip, :math:`h = 600` mm, C25
(:math:`f_{ctm} = 2.565` MPa), B500S:

.. math::

   k = 1.00 - \frac{600 - 300}{800 - 300}\,(1.00 - 0.65) = 0.79

.. math::

   A_{s,min} = \frac{0.40 \times 0.79 \times 2.565 \times (1000 \times 600 / 2)}
                    {500/1.15} = 559\ \text{mm}^2 = 5.59\ \text{cm}^2

against a geometric minimum of :math:`0.0009 \times 1000 \times 600 = 540` mm², so
crack control governs. The beam rule this replaces,
:math:`\max(0.26 f_{ctm}/f_{yk},\ 0.0013)\, b\, d`, would give 7.26 cm².

Which of the two governs depends only on :math:`k`, since both are proportional to
:math:`h`. With :math:`A_{ct} = b\,h/2` the crack-control ratio goes
with :math:`k/2`, and :math:`k` decays from 1.00 to 0.65 between 300 and 800 mm — so
crack control governs the **thin** footings and the geometric minimum takes over in the
thick ones, the crossover falling around :math:`h = 640` mm for C25 with B500S:

.. list-table::
   :header-rows: 1
   :widths: 20 27 27 26

   * - :math:`h`
     - A_s,min applied
     - Geometric alone
     - Governing rule
   * - 0.40 m
     - 1.097 ‰
     - 0.900 ‰
     - Crack control, +22%
   * - 0.60 m
     - 0.932 ‰
     - 0.900 ‰
     - Crack control, +4%
   * - 0.90 m
     - 0.900 ‰
     - 0.900 ‰
     - Geometric

.. note::

   A footing bending both ways has a tension zone on each face, so each is owed the
   crack-control minimum. A face the moment does not put in tension is given no minimum
   at all: whether a footing carries steel there is the consumer's call — it is the only
   one that sees both orthogonal sections of the element — so a design given no negative
   moment leaves the top face empty.

   The distribution reinforcement perpendicular to the span is the same question asked
   of a different section, and is designed as one. mento never reasons about two
   directions at once.

Detailing
---------

Bar spacing
^^^^^^^^^^^

Kept between **100 and 300 mm**, both bounds supplied by the design code through the
registry:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Bound
     - Value
     - Reason
   * - Maximum
     - 300 mm, in place of the slab's
       :math:`\min(3h,\ 450\,/\,400\ \text{mm})`
     - A footing is thick, so :math:`3h` stops binding long before the bars are close
       enough to spread the bearing pressure into them. ACI §7.7.2.3 and
       EN §9.3.1.1(3) still apply; the 300 mm cap is simply always the smaller.
   * - Minimum
     - 100 mm
     - EN detailing practice for foundations, applied under both codes: nothing about
       it is particular to EN. Below it, placing and vibrating over the soil stops
       being worth the closer bars.

The maximum is applied when a design is written back as a spacing, which only ever
adds bars, so the area stays covered. The minimum cannot be applied that way — the
rebar search chooses a **bar count**, and the spacing is what that count implies — so
it is enforced as a cap on ``max_bars_per_layer``:

.. math::

   n_{max} = \left\lfloor \frac{b}{s_{min}} \right\rfloor

A heavier demand is therefore answered with a larger bar rather than with bars closer
together. Both bounds are reported in the flexure check table, so a footing detailed
outside the range fails the check rather than passing quietly.

One mat, chosen rather than reconciled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A footing is not built as two independently detailed faces. The bars go down as one
grid: the bottom is placed, the top is placed over it, and the two are tied. The module
has to be one number, and on a drawing the top sits either at the bottom's spacing or at
exactly twice it, so that one top bar lands on every second bottom bar.

Reconciling the two spacings the per-face design arrived at — taking the smaller of them
— is safe but expensive. The face governed by the minimum is detailed as many thin bars,
which comes out closer than the few thick bars the loaded face needs, so the smaller of
the two is often the *lightly loaded* face's spacing, and adopting it closes the loaded
face up and buys steel nobody asked for. Measured over twenty typical footings, that
costs about 20% more longitudinal steel than the two faces actually require.

So the mat is chosen, not reconciled. Over the modules the code allows and the bars the
catalogue holds, mento takes the lightest combination that covers both faces:

.. math::

   \min_{s,\ d_{bot},\ d_{top}} \;
   n(s)\,A_b(d_{bot}) + n(s_{top})\,A_b(d_{top})

subject to

.. math::

   s_{top} \in \{s,\ 2s\}, \qquad
   s_{min} \le s \le s_{top} \le s_{max}, \qquad
   n(s) = \left\lceil \frac{b}{s} \right\rceil

with each face covering its own required area and each bar leaving the clear distance
the settings ask for. The module is searched in whole centimetres (inches in imperial),
so it is a number a detailer writes. Allowing a triple module as well was tried and
changes nothing: ``s`` and ``2s`` already carry it.

That lands within about **7% (EN) to 11% (ACI)** of what the two faces require, against
the 20% of taking the smaller spacing — and within 1% (EN) to 4% (ACI) of what mento's
own per-face design produces, which is not even buildable as a single mat.

Every candidate is verified before it is kept, because the choice feeds back into
itself: a thicker bar sits deeper, lowering the effective depth the required area was
computed against. Candidates are tried lightest first and applied to the section, and
the first that resists both design moments is the answer. If none within the budget
verifies — a section too small for its demand — the design falls back to reconciling the
two spacings it had already proved, which needs no verification because it only ever
adds bars.

A face carrying no reinforcement is left with none: a footing reinforced on the bottom
only has no second grid to place, and one is not invented for it.

Practical bar size
^^^^^^^^^^^^^^^^^^

The search does not reach below **Ø10** for a footing mat, whatever the catalogue holds
under it. Practice rather than code: a mesh finer than that is not what goes down over
the ground, and without the floor the optimisation reaches for the thinnest bar
available to shave the lightly loaded face.

Minimum thickness
^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Code
     - Minimum
     - Source
   * - ACI 318-19, CIRSOC 201-25
     - 200 mm (8 in.)
     - §13.3.1.2, written as :math:`d \ge 150` mm above the bottom reinforcement
   * - EN 1992-1-1
     - 250 mm
     - Practice, not a clause — the depth below which bar anchorage and the tolerance
       on a surface cast against the ground stop working out

A thinner section **warns and is still designed**. The thickness is the engineer's to
choose, and the design that follows is the right design for the section given;
refusing would decline to answer a question that has an answer.

Implementation decisions
------------------------

.. _footing-decisions:

Where the codes leave room for judgement, mento takes these readings:

**The tension zone is half the section.** :math:`A_{ct}` in Eq. (7.1) is *"the area of
concrete within the tensile zone... just before formation of the first crack"*. Just
before cracking the section is still uncracked, so for a rectangle in pure bending the
neutral axis sits at mid-depth and :math:`A_{ct} = b\,h/2`. Taking the full
:math:`b\,h` would double the minimum and make the geometric rule irrelevant at every
depth.

**The 4/3 relief is not applied to a footing.** ACI §9.6.1.3 lets a section satisfy
:math:`4/3` of the required steel instead of the minimum. It belongs to the clause it
relieves, §9.6.1.2, which §9.6.1.1(b) has already exempted the member from — there is
no sudden cracking failure here for extra steel to buy off, so the shrinkage and
temperature minimum stands as written.

**The 1.8‰ geometric minimum still applies to a free member.** The custom
:math:`1.8\text{‰}\,b\,h` rule mento adds under ACI (see :ref:`aci-decisions`) is
unchanged for beams and suspended slabs. On a footing it is superseded — and at
:math:`f_y = 420` MPa the two coincide, since Table 24.4.3.2 gives exactly 0.0018.

**Steel grades between the EN anchors are interpolated.** The halved geometric
minimum is tabulated at :math:`f_{yk} = 400` and 500 MPa only. mento interpolates
linearly between them and holds the value flat outside, rather than extrapolating a
rule the source does not state.

**ACI's f_y < 420 MPa branch is not special-cased.** Table 24.4.3.2 gives a flat
0.0020 for deformed bars below the reference grade, where the scaling used here gives
:math:`0.0018 \times 420/f_y`, which is larger — 0.0025 at 300 MPa. mento uses the
scaling in both directions, so the answer is conservative for the lower grades rather
than exact.

**The minimum is split by face, as it is for a beam.** A check reports the minimum on
the face in tension for the combination and zero on the other, following the same
convention every element uses. A footing with reinforcement on both faces meets it on
each face through the combination that puts that face in tension.

**An unreinforced face reports DCR ≫ 1.** A footing whose demand cannot be met within
the 100 mm spacing floor comes back with no layout, as any section does when the bar
catalogue cannot satisfy it. The DCR is then formed against a floor of 0.01 kNm rather
than against zero, so it reports a ratio far above 1 instead of failing to divide.

Not implemented
---------------

Everything the :doc:`one_way_slab` page lists as out of scope, plus:

- **Any geotechnical quantity** — see the warning at the top of this page.
- **Punching shear** — a separate element, ``PunchingSlab``.
- **The transverse direction.** A footing spans both ways; mento designs one strip at
  a time, and the second direction is a second ``Footing``. Distribution steel and the
  reduction some codes allow in the short direction of a rectangular footing are the
  engineer's.
- **Anchorage of the bars beyond the critical section**, which for a footing is often
  what sets the plan dimension.

Validation
----------

.. list-table::
   :header-rows: 1
   :widths: 34 40 26

   * - Check
     - Test (``tests/test_footing.py``)
     - Verified against
   * - :math:`\rho_{st}`, both unit systems and the 0.0014 floor
     - ``test_aci_shrinkage_and_temperature_ratio``,
       ``test_aci_shrinkage_and_temperature_ratio_imperial``
     - ACI 318-19 Table 24.4.3.2
   * - :math:`A_{s,min}` on the gross section, ACI
     - ``test_aci_footing_minimum_is_the_gross_section_rule``,
       ``test_aci_footing_minimum_scales_with_the_steel_grade``,
       ``test_aci_footing_minimum_in_imperial_units``
     - ACI 318-19 §9.6.1.1(b), §13.3.1.2, Table 24.4.3.2
   * - The exemption lowers the minimum
     - ``test_aci_footing_minimum_is_lower_than_the_slab_minimum``,
       ``test_en_footing_minimum_is_lower_than_the_slab_minimum``
     - ACI 318-19 §9.6.1.1(b); EN 1992-1-1 §9.2.1.1
   * - CIRSOC shares the ACI clause
     - ``test_cirsoc_shares_the_aci_clause``
     - CIRSOC 201-25
   * - :math:`k` interpolated with the depth
     - ``test_en_crack_control_coefficient_k``
     - EN 1992-1-1 §7.3.2(2)
   * - :math:`\rho_{geo}` at and between the anchors
     - ``test_en_foundation_min_reinforcement_ratio``
     - EN 1992-1-1 §9.2.1.1, §9.8.1
   * - Eq. (7.1) as written
     - ``test_en_crack_control_min_reinforcement_is_equation_7_1``
     - EN 1992-1-1 §7.3.2(2)
   * - The larger of the two rules governs
     - ``test_en_footing_minimum_is_the_crack_control_rule``,
       ``test_en_footing_minimum_is_the_larger_of_the_two_rules``
     - Worked reference case above
   * - Crack control governs the thin sections, not the thick
     - ``test_the_crack_rule_governs_the_thin_footings_not_the_thick``,
       ``test_the_excess_over_the_geometric_minimum_shrinks_with_depth``
     - EN 1992-1-1 §7.3.2(2); the depth table above
   * - Both faces in tension are both owed it
     - ``test_both_faces_in_tension_both_get_the_crack_minimum``
     - EN 1992-1-1 §7.3.2(2)
   * - A face not in tension is given no minimum
     - ``test_a_footing_with_no_top_moment_is_left_without_top_steel``
     - Internal consistency
   * - A design covers the minimum without correction
     - ``test_design_already_covers_the_minimum``,
       ``test_aci_footing_with_no_moment_still_gets_the_minimum``,
       ``test_a_footing_still_carries_its_moment``
     - Internal consistency
   * - Spacing bounds, and a design inside them
     - ``test_footing_bars_are_capped_at_300_mm``,
       ``test_footing_bars_are_floored_at_100_mm``,
       ``test_a_designed_footing_stays_inside_the_spacing_range``,
       ``test_the_floor_is_what_holds_the_footing_apart``
     - Detailing rules above
   * - Both bounds reported in the check table
     - ``test_the_spacing_row_of_a_footing_reports_both_bounds``
     - Internal consistency
   * - One module across both faces, the top at ``s`` or ``2s``
     - ``test_a_designed_footing_is_one_mat``,
       ``test_the_top_may_be_set_out_at_twice_the_bottom``,
       ``test_each_face_keeps_its_own_diameter``
     - Detailing practice; the rule above
   * - The chosen mat is lighter than reconciling the two spacings
     - ``test_the_mat_is_lighter_than_reconciling_the_two_spacings``
     - The measured comparison above
   * - Ø10 floor on the mat
     - ``test_no_mat_is_detailed_below_the_practical_bar``
     - Detailing practice
   * - The mat never undoes the design
     - ``test_matching_the_mat_never_undoes_the_design``,
       ``test_the_mat_keeps_each_face_within_the_spacing_range``,
       ``test_a_footing_reinforced_on_one_face_gets_no_second_grid``,
       ``test_a_slab_still_details_its_faces_independently``
     - Internal consistency
   * - A candidate that does not verify is passed over
     - ``test_a_mat_that_does_not_verify_is_passed_over``,
       ``test_a_section_that_cannot_carry_the_moment_falls_back``,
       ``test_the_fallback_takes_the_smaller_of_the_two_spacings``,
       ``test_the_fallback_leaves_a_single_grid_alone``,
       ``test_reconciling_an_existing_mat_is_a_no_op``
     - Internal consistency
   * - The edges of the search
     - ``test_a_section_needing_nothing_is_offered_no_mat``,
       ``test_a_face_no_bar_can_cover_drops_that_module``
     - Internal consistency
   * - The bounds are a footing's, not every slab's
     - ``test_a_slab_keeps_the_wider_slab_limits``,
       ``test_support_says_which_clause_applies``
     - ACI 318-19 §7.7.2.3; EN 1992-1-1 §9.3.1.1(3)
   * - Thickness warns and still designs
     - ``test_a_thin_footing_warns``,
       ``test_a_footing_of_the_usual_depth_does_not_warn``,
       ``test_a_thin_slab_does_not_warn``,
       ``test_a_thin_footing_is_still_designed``,
       ``test_the_warning_names_the_footing_and_the_code``
     - ACI 318-19 §13.3.1.2; EN practice
   * - An unreinforced face reports a failing DCR
     - ``test_an_unreinforced_face_reports_a_failing_dcr``
     - Internal consistency, both codes
