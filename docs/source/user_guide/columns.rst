Columns
=======

``RectangularColumn`` and ``CircularColumn`` describe a column section -- its
shape, its longitudinal bars and its ties -- and check its **shear in both
directions**, with the axial load of each combination, under ACI 318-19 and
CIRSOC 201-25.

.. note::

    The column checks shear only. Flexure with axial load (the P–M
    interaction) is not in mento yet, so ``check_flexure`` raises, and the
    ties are checked but not designed: ``design`` raises too. EN 1992-1-1 is
    the next step. The plan is ``docs/architecture/column-shear-roadmap.md``
    in the repository.

Axes and directions
-------------------

The column runs along x (:doc:`local_axes`). ``width`` is measured along y and
``height`` along z. Each combination is checked twice:

- ``V_z`` acts across ``height``: the shear section is ``width`` wide, and its
  tension steel is the bars of a face of length ``width``.
- ``V_y`` acts across ``width``: the same turned, ``height`` wide.

``N_x`` enters both, positive in compression. The moments are not read, and a
combination that carries them is warned about (``column_force_component_not_checked``).

A rectangular column
--------------------

.. code-block:: python

    from mento import Concrete_ACI_318_19, SteelBar, RectangularColumn, Forces, Node
    from mento import MPa, cm, mm, kN

    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)

    column = RectangularColumn(
        label="C1", concrete=concrete, steel_bar=steel,
        width=30 * cm, height=50 * cm, c_c=25 * mm,
    )
    # Ø20 in the corners; between them, 1Ø16 on each face of length width
    # and 2Ø16 on each face of length height.
    column.set_longitudinal_rebar(d_b_corner=20 * mm, n_width=1, d_b_width=16 * mm,
                                  n_height=2, d_b_height=16 * mm)
    # A Ø10 tie every 15 cm, and one Ø8 crosstie holding a bar of each face of length height.
    column.set_ties(d_b=10 * mm, s=15 * cm, crossties_height=1, d_b_crosstie=8 * mm)

    node = Node(column, [Forces(label="1.2D+1.6L", N_x=900 * kN, V_y=60 * kN, V_z=140 * kN)])
    node.check()

    for check in column.shear_checks:
        print(check.direction, check.DCR, check.reinforcement_required, check.A_v_req)

**The legs.** The perimeter tie gives two legs in each direction. A crosstie
is one leg, in the shear it crosses: ``crossties_width`` hold bars of the faces
of length ``width``, span the height and count in ``V_z``; ``crossties_height``
count in ``V_y``. A crosstie may be thinner than the tie (``d_b_crosstie``) and
still counts. Whether it is thick enough to brace its bar is a detailing rule,
which this check does not read.

**d.** From the compressed face to the centroid of the bars of the far face,
as for a beam: a column checked in one direction gives the same numbers as a
``RectangularBeam`` of the same rectangle where the clauses are the same.

A circular column
-----------------

.. code-block:: python

    from mento import CircularColumn

    column = CircularColumn(label="C2", concrete=concrete, steel_bar=steel,
                            diameter=50 * cm, c_c=30 * mm)
    column.set_longitudinal_rebar(n=8, d_b=20 * mm)
    column.set_hoops(d_b=10 * mm, s=15 * cm)

A circle is checked as ACI 318-19 §22.5.2.2 allows: ``b_w = D`` and
``d = 0.8·D``, with the gross area of the circle in the axial term, so a
circle in tension is checked as well as one in compression. A hoop counts as
two legs. The tension steel of ρ_w is the bars past the centre line,
``(n - 1) // 2`` of them.

What a column is checked for
----------------------------

The resistance is the beam's -- V_c of Table 22.5.5.1 with ``N_u/(6·A_g)``,
V_s of §22.5.8.5.3, the section limit of §22.5.1.2. What the column reads
from its own chapter, the same clauses in both codes:

- **When shear reinforcement is required**: §10.6.2.1, ``Vu > 0.5·ØVc``, with
  the axial load in V_c. A beam is held to §9.6.3.1 instead, which ignores the
  axial load; a compressed column needs shear reinforcement later than a beam
  of the same section, a column in tension sooner.
- **How much**: there, the area the shear asks for and at least A_v,min of
  §10.6.2.2. Elsewhere ``A_v_req`` is zero.
- **How far apart**: there, Table 10.7.6.5.2 along the column (``s_max_l``)
  and between legs across it (``s_max_w``) -- the rows of Table 9.7.6.2.2,
  with the caps of each code: 600/300 mm under ACI 318-19, 400/200 mm under
  CIRSOC 201-25.
- The cap of 8.3 MPa on √f'c (§22.5.3.1) stays: §22.5.3.2 lifts it for beams
  and joists only.

Each result, :class:`~mento.column_results.ColumnShearCheck`, carries
``V_concrete`` (ØVc), ``reinforcement_required``, ``A_v``, ``A_v_req``,
``A_v_min``, ``V_capacity``, ``V_max``, the ``DCR`` and the spacings. The ties
of the detailing (ACI 318-19 §25.7.2) are required everywhere and are not
checked here.

Warnings
--------

``column.warnings`` (or ``node.warnings``) lists what the ties miss, per
direction, in the language set with :func:`mento.set_language`:
``column_ties_below_required``, ``column_tie_spacing_exceeds_max``,
``column_leg_spacing_exceeds_max`` and ``column_shear_exceeds_section_limit``,
each with ``values["direction"]``. They are worded when read, so a loop that
only reads ``shear_check_results`` does not pay for them.
