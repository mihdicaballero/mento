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

Two rules apply and **the larger governs**.

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

Which rule governs depends only on :math:`k`, since both are proportional to
:math:`h`: crack control governs the shallower sections, and the geometric minimum
takes over once :math:`k` has decayed — around :math:`h = 640` mm for C25 with B500S.

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

One spacing for the whole mat
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A footing is not built as two independently detailed faces. The bars go down as one
mat: the bottom grid is placed, the top grid is placed over it on the same module, and
the two are tied at the same spacing. A design that answers the bottom at 120 mm and
the top at 250 mm is not a drawing anyone details, so once both faces are settled they
are reconciled to one number:

.. math::

   s_{mat} = \min\bigl(s_{bot},\ s_{top}\bigr)

The **smaller** of the two, because it is the only direction that is safe: closing the
wider face up adds bars to it, so every face still covers the area its own moment asked
for.

Which face sets the module is not always the one carrying the moment. A face governed
by the minimum is detailed as many thin bars, and that lands at a closer spacing than
the few thick bars the governing face needs — under ACI, whose minimum is the larger of
the two codes', the top face sets the module about half the time and closes the bottom
up with it. Measured over a range of typical footings, the mat costs about **12% more
longitudinal steel** than detailing the two faces independently would, and in the worst
case seen, close to 50%. That is the price of a single grid, and it is the grid the
drawing would have shown anyway.

Only the spacing is reconciled; the diameters stay as the design chose them, so the two
faces may differ there. Adding bars to a layer does not move its centroid, so the
effective depths are unchanged and the capacity of the closed-up face can only rise —
which is why the reconciliation runs **after** the design's final capacity
verification. A face carrying no reinforcement is left with none: a footing reinforced
on the bottom only has no second grid to match, and one is not invented for it.

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
   * - One spacing across both faces, the smaller of the two
     - ``test_a_designed_footing_is_one_mat``,
       ``test_the_mat_takes_the_smaller_of_the_two_spacings``,
       ``test_the_diameters_are_left_as_designed``,
       ``test_matching_is_idempotent``
     - Detailing practice; the rule above
   * - Matching never undoes the design
     - ``test_matching_the_mat_never_undoes_the_design``,
       ``test_the_mat_keeps_each_face_within_the_spacing_range``,
       ``test_a_footing_reinforced_on_one_face_gets_no_second_grid``,
       ``test_a_slab_still_details_its_faces_independently``
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
