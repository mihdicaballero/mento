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
    flexure.bottom.A_s_min, flexure.bottom.A_s_max
    flexure.bottom.DCR                  # 0.965
    flexure.bottom.n_bars               # 3

    flexure.DCR                         # worst of the two faces
    str(flexure.bottom)                 # '2Ø16 mm + 1Ø12 mm'

Each face carries its layers, in order, and only the layers that hold bars:

.. code-block:: python

    for layer in flexure.bottom.layers:
        print(layer.n, layer.d_b, layer.A_s.to("cm**2"))

A face with no reinforcement has an empty ``layers`` tuple, which is what
``str(flexure.top)`` reports as ``'no reinforcement'``.

A slab is detailed by a spacing rather than by a bar count, so each of its layers also
carries the spacing it was designed with, and reads as one bar repeated across the strip.
The bars are still there to be counted when what you need is the steel actually placed:

.. code-block:: python

    layer = slab.flexure_design.bottom.layers[0]

    layer.d_b                            # 12 mm
    layer.s.to("cm")                     # 17 cm, None on a beam
    layer.n                              # 6, the bars that spacing puts on the strip
    str(layer)                           # 'Ø12 mm/17 cm'

    slab.flexure_design.bottom.n_bars    # 6, every layer of the face

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
