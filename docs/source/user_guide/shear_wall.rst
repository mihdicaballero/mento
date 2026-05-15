Shear Wall
==========

The `ShearWall` class models a reinforced-concrete structural wall for **in-plane
shear analysis and design** per ACI 318-19 Chapter 11.

- Concrete shear capacity follows ACI 318-19 §11.5.4.6 with the aspect-ratio
  factor ``α_c``, instead of a longitudinal-reinforcement term.
- Reinforcement is **distributed mesh** in two orthogonal directions
  (``ρt`` horizontal, ``ρl`` vertical), placed on **both faces** of the wall
  (E.F. — each face), not stirrups.
- The minimum horizontal ratio is ``ρt,min = 0.0025`` (§11.6.1). The minimum
  vertical ratio follows the §11.6.2 interpolation
  ``ρl,min = max(0.0025, 0.0025 + 0.5·(2.5 − hw/lw)·(ρt,req − 0.0025))``,
  with ``hw/lw`` clamped to ``[0.5, 2.5]``.

The same shear provisions serve both **ACI 318-19** and **CIRSOC 201-25**;
CIRSOC differs only in the reinforcing-bar catalogue used for design (it
allows Ø6 mm for the transverse mesh and Ø10 mm minimum for the vertical mesh).

.. note::

    The module covers **shear check and design only**. Flexure design
    for shear walls is not implemented yet. Inherited flexure methods from
    ``RectangularBeam`` are not validated for wall geometry and should not be
    used.

Key Concepts
------------

- **Geometry**:

  - ``thickness`` (*t*) — out-of-plane dimension
  - ``length`` (*lw*) — in-plane length, resists in-plane shear
  - ``height`` (*hw*) — story / overall wall height, used for the ``hw / lw`` aspect ratio

- **Material Properties**: requires a ``Concrete`` object (currently
  ``Concrete_ACI_318_19`` or ``Concrete_CIRSOC_201_25``) and a ``SteelBar``
  object.

- **Reinforcement**: distributed bars defined by **bar diameter + spacing**
  in each direction.

- **Forces**: a list of ``Forces`` objects passed to ``check_shear()`` or
  ``design_shear()`` — either directly on the wall or through a ``Node``.
  ``V_z`` is the in-plane lateral shear demand at the section.

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
(longitudinal) mesh. The mesh is placed on **both faces** of the wall (E.F. —
each face), so the reinforcement ratio counts both curtains:
``ρ = 2 · Ab / (t × s)``.

.. code-block:: python

    # Ø12 @ 150 mm in both directions
    wall.set_horizontal_rebar(d_b=12 * mm, s=150 * mm)
    wall.set_vertical_rebar(d_b=12 * mm, s=150 * mm)

3. Performing the Shear Check
*****************************

Call ``check_shear()`` with a list of ``Forces`` — either directly on the wall
or through a ``Node``. The returned DataFrame has one header row (units)
followed by one row per combination.

.. code-block:: python

    from mento import Forces, Node, kN

    f1 = Forces(label="1.2D+1.0E", V_z=800  * kN)
    f2 = Forces(label="0.9D+1.0E", V_z=1200 * kN)

    # Directly on the wall ...
    wall.check_shear([f1, f2])

    # ... or via a Node (recommended — enables results / detailed reports)
    node = Node(section=wall, forces=[f1, f2])
    node.check()
    node.results

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

``design()`` is fully automatic: it sizes **both** meshes against the
worst-case force combination, applies them to the wall, and returns the
re-evaluated check DataFrame.

.. code-block:: python

    result = wall.design([f1, f2])
    # The wall now carries a designed mesh:
    print(wall._d_b_h, wall._s_h)   # horizontal (shear) bar + spacing
    print(wall._d_b_v, wall._s_v)   # vertical (minimum) bar + spacing

What it does:

1. Runs the check for every force and tracks the worst-case ``ρt,req``.
2. Derives the worst-case ``ρl,min`` from §11.6.2.
3. Selects a bar diameter and spacing for the **horizontal mesh** (against
   ``ρt,req``) and the **vertical mesh** (against ``ρl,min``).
4. Applies both via ``set_horizontal_rebar`` / ``set_vertical_rebar`` and
   re-runs the check.


**Vertical mesh.** Because for now *mento* does not check flexure, the vertical mesh
is always sized to the §11.6.2 *minimum* — it is reported and plotted as
"Minimum vertical rebar".

.. note::

    For **CIRSOC 201-25**, ``design`` automatically uses the CIRSOC bar
    catalogue (Ø6 mm minimum transverse, Ø10 mm minimum vertical). No extra
    configuration is needed — the design code is read from the ``Concrete``
    object.

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
| ``_rho_t_req``       | Required horizontal reinforcement ratio      |
+----------------------+----------------------------------------------+
| ``_rho_l_min``       | Minimum vertical reinforcement ratio         |
+----------------------+----------------------------------------------+
| ``_s_h_max``         | §11.7.3 horizontal spacing limit             |
+----------------------+----------------------------------------------+
| ``_s_v_max``         | §11.7.3 vertical spacing limit               |
+----------------------+----------------------------------------------+
| ``_d_b_h`` / ``_s_h``| Designed horizontal bar diameter / spacing   |
+----------------------+----------------------------------------------+
| ``_d_b_v`` / ``_s_v``| Designed vertical bar diameter / spacing     |
+----------------------+----------------------------------------------+

All reinforcement ratios (``ρt``, ``ρl``) account for the mesh on **both
faces** — ``ρ = 2 · Ab / (t · s)``.

A worked example with full output is available in the
:doc:`Shear Wall ACI 318-19 example </examples/shear_wall_check_ACI_318-19>`.
