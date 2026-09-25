.. _user_guide/design_results:

Design results
==============

After a check or a design has run, the reinforcement it produced is available as plain
data through two properties of the section: ``flexure_design`` and ``shear_design``.

These return frozen objects, so they are a snapshot of the result rather than a live view
of the section. Every quantity is a pint ``Quantity`` in the unit system of the section's
concrete, and can be converted with ``.to()`` as usual.

.. code-block:: python

    from mento import Concrete_ACI_318_19, SteelBar, RectangularBeam, Node, Forces
    from mento import MPa, cm, mm, kN, kNm

    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="101", concrete=concrete, steel_bar=steel,
        width=20 * cm, height=60 * cm, c_c=25 * mm,
    )

    node = Node(section=beam, forces=[Forces(label="C1", V_z=80 * kN, M_y=100 * kNm)])
    node.design()

Flexure
-------

``flexure_design`` has a ``bottom`` and a ``top`` face:

.. code-block:: python

    flexure = beam.flexure_design

    flexure.bottom.A_s.to("cm**2")      # 5.15 cm², steel provided
    flexure.bottom.A_s_req.to("cm**2")  # 4.96 cm², steel required
    flexure.bottom.A_s_calc.to("cm**2") # 4.96 cm², what the moment alone asks for
    flexure.bottom.A_s_min, flexure.bottom.A_s_max
    flexure.bottom.DCR                  # 0.965
    flexure.bottom.M_capacity           # 103.6 kN·m, ØMn here; MRd under EN 1992
    flexure.bottom.n_bars               # 3

    flexure.DCR                         # worst of the two faces
    str(flexure.bottom)                 # '2Ø16 mm + 1Ø12 mm'

Each face carries its layers, in order, and only the layers that hold bars:

.. code-block:: python

    for layer in flexure.bottom.layers:
        print(layer.n, layer.d_b, layer.A_s.to("cm**2"))

A face with no reinforcement has an empty ``layers`` tuple, which is what
``str(flexure.top)`` reports as ``'no reinforcement'``.

``A_s_req`` and ``A_s_calc`` are the same number on this beam because the moment
governs. They part on a face governed by its minimum: ``A_s_req`` is the steel to detail,
``max(A_s_calc, A_s_min)`` — or the 4/3 rule under ACI — while ``A_s_calc`` is what the
moment alone asked for, and is zero with no moment. Choose bars from the first. Scale an
anchorage by ``A_s,nec / A_s,prov`` from the second: a footing mat governed by its
minimum carries little of the stress that minimum is sized for, and charging it the full
development length anchors a force that is not there.

.. code-block:: python

    face = footing.flexure_design.bottom      # a 1 m × 0.40 m strip under 20 kN·m, ACI

    face.A_s_calc.to("cm**2")                 # 1.5 cm², the moment's own demand
    face.A_s_min.to("cm**2")                  # 7.2 cm², the minimum on the ground
    face.A_s_req == face.A_s_min              # True: the minimum governs

A slab is detailed by a spacing rather than by a bar count, so each of its layers also
carries the spacing it was designed with, and reads as one bar repeated across the strip.
The bars are still there to be counted when what you need is the steel actually placed:

.. code-block:: python

    layer = slab.flexure_design.bottom.layers[0]

    layer.d_b                            # 12 mm
    layer.s.to("cm")                     # 17 cm, None on a beam
    layer.n                              # 5.88, the bars a metre carries at that spacing: 100/17
    str(layer)                           # 'Ø12 mm/17 cm'

    slab.flexure_design.bottom.n_bars    # 5.88, every layer of the face

The count of a slab layer is ``width / s`` and is not a whole number: the strip is a
slice of a slab that goes on past its edges, and its steel is the bar area times the
bars per metre, ``layer.A_s == layer.n * π d_b² / 4`` on a slab as on a beam.

Shear
-----

.. code-block:: python

    shear = beam.shear_design

    shear.n_stirrups        # 1, number of stirrups
    shear.n_legs            # 2, legs crossing the shear plane
    shear.d_b               # 10 mm
    shear.s_l.to("cm")      # 27 cm, longitudinal spacing
    shear.A_v.to("cm**2/m") # 5.82 cm²/m, provided
    shear.A_v_req, shear.A_v_min
    shear.DCR               # 0.462
    shear.V_capacity        # 173 kN, ØVn here; VRd under EN 1992

    str(shear)              # '1eØ10 mm/27 cm'

Several load combinations
-------------------------

A section is normally checked against a list of combinations, and each face is often
governed by a different one. The required areas and the DCRs are the envelope over the
whole list, so they describe the combination that governs — not whichever one happened to
be checked last:

.. code-block:: python

    node = Node(section=beam, forces=[f1, f2, f3])
    node.design()

    beam.flexure_design.bottom.DCR  # worst of the three on the bottom face
    beam.shear_design.A_v_req       # largest stirrup requirement of the three

The provided reinforcement — ``A_s``, the layers and the stirrup layout — describes the
section itself and does not depend on the combination.

The capacity follows the DCR: it is the one of the combination that governs, so the two
remain the ratio they were. Under ACI 318-19 the shear resistance can differ between
combinations, since ``Vc`` depends on the axial load and on which face is in tension; the
per-combination results keep each one's own.

Capacities
----------

Each result also carries the resistance its ``DCR`` was formed from: ``M_capacity`` on a
flexure face and ``V_capacity`` on the shear result. The name is the same under every
design code — it holds ``ØMn`` and ``ØVn`` under ACI 318-19 and CIRSOC 201-25, ``MRd``
and ``VRd`` under EN 1992-1-1 — so a report that prints the symbol names it per code and
the data does not have to.

For any combination with a demand, dividing the demand by its ``DCR`` gives the capacity
back, to within the rounding the design code applies to the ratio. The field is there for
the cases the ratio cannot cover: a face carrying minimum reinforcement against no demand
has a ``DCR`` of zero and a real capacity, and two faces with identical reinforcement
report exactly the same one.

The per-combination results are available too, one per combination of the last check:

.. code-block:: python

    for check in beam.shear_checks:
        check.label, check.DCR, check.V_capacity

    for check in beam.flexure_checks:
        check.label, check.bottom.DCR, check.bottom.M_capacity

    beam.shear_design.V_capacity            # the governing combination's

Design alternatives
-------------------

A design ranks every layout that fits and applies the best one. The runners-up are kept
too, best first, with the applied layout always in first place -- but only the ones the
finished section passes with. Each longitudinal alternative is built on the beam as the
design left it, with the stirrups it ended with and the other face as applied, and kept
if the beam carries both moments with it within the code's limits on its reinforcement
(tension-controlled under ACI 318-19 / CIRSOC 201-25, the 4 % of EN 1992-1-1); its
``DCR`` says at what ratio. A footing offers none, because its mat is chosen as a whole.

.. code-block:: python

    node.design()

    for option in beam.flexure_design.bottom.options:
        str(option), option.A_s, option.functional   # '2Ø16 mm + 1Ø12 mm', ...

    beam.flexure_design.top.options                  # the same for the top face
    beam.shear_design.options                        # StirrupOption: n_stirrups, d_b, s_l, s_w, A_v

A longitudinal option (``RebarOption``) carries its ``layers`` — the same ``RebarLayer``
objects the applied reinforcement is read as — its area and the ``functional`` the search
ranked it by.

The stirrup alternatives are one layout per other bar diameter the code offers, lighter and
heavier alike, in order of diameter: each is the widest spacing with the fewest legs that
covers the demand read at the depth that bar gives the section. Where the spacing limit
governs they share one spacing (``1eØ10/13``, ``1eØ12/13``, ``1eØ16/13``); where the demand
governs, a lighter bar sits closer and a heavier one further apart. Every alternative is
built on the finished section and checked there -- shear and flexure, since a heavier
stirrup lowers the effective depth -- and only the ones the section passes with are kept, so
the list answers "what if I use the bar I have". Each option carries its ``DCR``, the worst
ratio of the section built with it, and its ``functional``, what it adds in steel: the excess
of ``A_v`` over what the section asks for with that bar, plus one per extra closed stirrup.

How many are kept is a setting, three by default:

.. code-block:: python

    beam = RectangularBeam(..., settings=BeamSettings(design_options=5))

The options belong to the design that produced them. A check alone reports none, and
changing the bars by hand afterwards clears the options of what was changed.

A design depends only on its inputs. It starts from the same state every time — the
stirrup diameter the settings assume, the placeholder bars — so running it again, or after
setting reinforcement by hand, gives the same result.

Warnings
--------

A detailing limit can be missed while the strength is fine, or met while it is not, so
the limits are reported beside the ``DCR`` instead of inside it. ``beam.warnings`` (or
``node.warnings``) lists them after a check or a design, one per limit and face:

.. code-block:: python

    node.check()
    for warning in node.warnings:
        warning.code          # 'stirrup_spacing_exceeds_max'
        warning.message       # 'Stirrup spacing along the member: 35 cm exceeds the maximum 13.9 cm.'
        warning.values        # {'s': 35 cm, 's_max': 13.9 cm}
        warning.combinations  # ('1.2D+1.6L', '1.4D')

The ``code`` is stable and is what a program should compare against. The ``message`` is
written in the language set with ``mento.set_language`` when the warnings are read. The
codes and what triggers each are listed in :mod:`mento.design_warnings`.

Reading results too early
-------------------------

Both properties raise ``DesignNotRunError`` if the corresponding check or design has not
been run, rather than returning zeros that could be mistaken for a real result:

.. code-block:: python

    beam = RectangularBeam(...)
    beam.flexure_design
    # DesignNotRunError: No flexure results yet. Run node.design() or
    # node.check_flexure() before reading flexure_design.

Relation to the private attributes
----------------------------------

Sections still keep their results in private attributes such as ``_A_s_bot`` and
``_stirrup_s_l``. Those remain in place and keep working, but they are implementation
details: their names, units and meaning can change between releases. The properties
described here are the supported way to read a result from code.

One difference worth noting: ``_stirrup_n`` counts stirrups, while the area ``A_v`` is
computed from the legs that cross the shear plane. The public object exposes both, as
``n_stirrups`` and ``n_legs``, so there is nothing to infer.
