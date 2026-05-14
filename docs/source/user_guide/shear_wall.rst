Shear Wall
==========

The `ShearWall` class models a reinforced-concrete structural wall for **in-plane
shear analysis and design** per ACI 318-19 Chapter 11. It is a child of
`RectangularBeam` and reuses the underlying materials, geometry, settings, and
units infrastructure. Compared to a beam, the shear formulas, the minimum
reinforcement, and the reinforcement layout differ:

- Concrete shear capacity follows ACI 318-19 §11.5.4.6 with the aspect-ratio
  factor ``α_c``, instead of a longitudinal-reinforcement term.
- Reinforcement is **distributed mesh** in two orthogonal directions
  (``ρt`` horizontal, ``ρl`` vertical), not stirrups.
- Minimum ratios follow §11.6.1 (``ρt,min = ρl,min = 0.0025``), with the
  additional rule ``ρl ≥ ρt`` when ``hw / lw ≤ 2.0``.

.. note::

    Phase 0 of the module covers **shear check and design only**. Flexure design
    for shear walls is not implemented yet. Inherited flexure methods from
    ``RectangularBeam`` are not validated for wall geometry and should not be
    used.

Key Concepts
------------

- **Geometry** (wall-friendly naming, mapped to the parent ``RectangularBeam`` fields):

  - ``thickness`` (*t*) — out-of-plane dimension
  - ``length`` (*lw*) — in-plane length, resists in-plane shear
  - ``height`` (*hw*) — story / overall wall height, used for the ``hw / lw`` aspect ratio

- **Material Properties**: requires a ``Concrete`` object (currently
  ``Concrete_ACI_318_19`` or ``Concrete_CIRSOC_201_25``) and a ``SteelBar``
  object.

- **Reinforcement**: distributed bars defined by **bar diameter + spacing**
  in each direction.

- **Forces**: a list of ``Forces`` objects passed directly to ``check_shear()``
  or ``design_shear()``. ``V_z`` is the in-plane lateral shear demand at the
  section. (Walls do not currently route through a ``Node``.)

Usage
-----

Below is a step-by-step guide on how to use the ``ShearWall`` class.

1. Creating a Shear Wall Object
*******************************

Specify the geometry, materials, and clear cover using the wall-friendly
constructor parameters ``thickness``, ``length``, and ``height``.

.. code-block:: python

    from mento import ShearWall, Concrete_ACI_318_19, SteelBar
    from mento import MPa, cm, mm, m

    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel    = SteelBar(name="ADN 420", f_y=420 * MPa)

    wall = ShearWall(
        label="W1",
        concrete=concrete,
        steel_bar=steel,
        thickness=25 * cm,   # t
        length=4.0 * m,      # lw (in-plane length)
        height=3.5 * m,      # hw (story height; hw/lw = 0.875)
        c_c=20 * mm,
    )

After construction, the dimensions are accessible as ``wall.thickness``,
``wall.length``, and ``wall.height``.

2. Setting the Distributed Reinforcement
****************************************

Use ``set_horizontal_rebar(d_b, s)`` for the horizontal (transverse) mesh that
resists in-plane shear, and ``set_vertical_rebar(d_b, s)`` for the vertical
(longitudinal) mesh. The corresponding ratios are computed internally as
``ρ = Ab / (t × s)``.

.. code-block:: python

    # Ø12 @ 150 mm in both directions
    wall.set_horizontal_rebar(d_b=12 * mm, s=150 * mm)
    wall.set_vertical_rebar(d_b=12 * mm, s=150 * mm)

3. Performing the Shear Check
*****************************

Call ``check_shear()`` directly on the wall with a list of ``Forces``. The
returned DataFrame has one header row (units) followed by one row per
combination.

.. code-block:: python

    from mento import Forces, kN

    f1 = Forces(label="1.2D+1.0E", V_z=800  * kN)
    f2 = Forces(label="0.9D+1.0E", V_z=1200 * kN)

    wall.check_shear([f1, f2])

Returned columns:

+--------------+-----------------------------------------------------+
| Column       | Meaning                                             |
+==============+=====================================================+
| ``ρt,min``   | Minimum horizontal reinforcement ratio (0.0025)     |
+--------------+-----------------------------------------------------+
| ``ρt,req``   | Required horizontal reinforcement ratio             |
+--------------+-----------------------------------------------------+
| ``ρt``       | Provided horizontal reinforcement ratio             |
+--------------+-----------------------------------------------------+
| ``ρl,min``   | Minimum vertical reinforcement ratio                |
+--------------+-----------------------------------------------------+
| ``ρl``       | Provided vertical reinforcement ratio               |
+--------------+-----------------------------------------------------+
| ``Vu``       | Demand shear                                        |
+--------------+-----------------------------------------------------+
| ``ØVc``      | Concrete shear strength (factored)                  |
+--------------+-----------------------------------------------------+
| ``ØVs``      | Steel shear strength (factored)                     |
+--------------+-----------------------------------------------------+
| ``ØVn``      | Total nominal shear strength (factored)             |
+--------------+-----------------------------------------------------+
| ``ØVn,max``  | Section-crushing limit (factored)                   |
+--------------+-----------------------------------------------------+
| ``DCR``      | Demand-to-capacity ratio                            |
+--------------+-----------------------------------------------------+

4. Designing the Shear Reinforcement
************************************

``design_shear()`` runs the check for every combination and stores the worst
case ``ρt,req`` on the wall instance. In Phase 0 the caller picks a bar
diameter and spacing satisfying both ``ρt,req`` and the §11.7.3 spacing
limits, then calls ``set_horizontal_rebar()``.

.. code-block:: python

    wall.design_shear([f1, f2])
    print(wall._rho_t_req.to("").magnitude)   # required horizontal ratio
    print(wall._s_h_max.to("mm").magnitude)    # spacing limit, horizontal bars
    print(wall._s_v_max.to("mm").magnitude)    # spacing limit, vertical bars

5. Inspecting Intermediate Quantities
*************************************

After running a check, the relevant ACI 318-19 quantities are available as
attributes on the wall:

+----------------------+----------------------------------------------+
| Attribute            | Meaning                                      |
+======================+==============================================+
| ``_hw_lw``           | Aspect ratio ``hw / lw``                     |
+----------------------+----------------------------------------------+
| ``_alpha_c``         | ``α_c`` factor (0.25 → 0.17 metric)          |
+----------------------+----------------------------------------------+
| ``_Acv``             | Gross shear area ``lw × t``                  |
+----------------------+----------------------------------------------+
| ``_V_c_wall``        | Nominal concrete shear strength              |
+----------------------+----------------------------------------------+
| ``_V_s_wall``        | Nominal steel shear strength                 |
+----------------------+----------------------------------------------+
| ``_phi_V_n_wall``    | Factored total shear strength                |
+----------------------+----------------------------------------------+
| ``_phi_V_n_max_wall``| Factored crushing-limit shear strength       |
+----------------------+----------------------------------------------+
| ``_DCRv_wall``       | Demand-to-capacity ratio for shear           |
+----------------------+----------------------------------------------+

A worked example with full output is available in the
:doc:`Shear Wall ACI 318-19 example </examples/shear_wall_check_ACI_318-19>`.
