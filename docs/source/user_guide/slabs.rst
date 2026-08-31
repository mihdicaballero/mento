Slab
====

The `OneWaySlab` class models a one-way spanning slab for flexure and shear analysis and design.
It inherits all check and design logic from `RectangularBeam`, but uses bar **diameter + spacing**
instead of bar count for longitudinal reinforcement — consistent with standard slab detailing practice.

Key Concepts
------------

- **Geometry**: Defined by a `width` and a slab `height` (thickness). The width is free:
  1 m is the usual choice for a floor slab, but any width works — a footing might use 150 cm.
  Results are the totals across that width, not per-metre values.
- **Material Properties**: Requires a `Concrete` object and a `SteelBar` object.
- **Reinforcement**: Specified by bar diameter and spacing, not bar count.
- **Checks and design**: Performed through a `Node` object, exactly as for beams.

Usage
-----

Below is a step-by-step guide on how to use the `OneWaySlab` class.

1. Creating a Slab Object
*************************

To define a slab, specify its geometry, material properties, and clear cover.
The `width` is the width being designed; `height` is the slab thickness.

.. code-block:: python

    from mento import Concrete_ACI_318_19, SteelBar, mm, cm, MPa
    from mento.slab import OneWaySlab

    # Define materials
    concrete = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)

    # Define slab geometry (1 m wide, 20 cm thick)
    slab = OneWaySlab(label="S101", concrete=concrete, steel_bar=steel,
                      width=100 * cm, height=20 * cm, c_c=25 * mm)

2. Setting Reinforcement
************************

Slab reinforcement is defined by bar diameter and spacing. Each layer accepts a primary
bar group (position 1) and an optional secondary bar group (position 3, second layer).

- **Bottom Longitudinal Reinforcement**: Use `set_slab_longitudinal_rebar_bot`.
- **Top Longitudinal Reinforcement**: Use `set_slab_longitudinal_rebar_top`.
- **Transverse Reinforcement**: Use `set_slab_transverse_rebar` (shear stirrups, rarely needed for slabs).

.. note::

    Unlike beams, slabs do not use positions 2 and 4 (inner bars per layer).
    All bars in a layer must have the same diameter. If no reinforcement is defined,
    *Mento* will assume no rebar and will not perform a check — provide at least a minimum
    rebar before running checks.

.. code-block:: python

    # Bottom reinforcement: Ø12 every 15 cm (layer 1)
    slab.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=15 * cm)

    # Top reinforcement: Ø10 every 20 cm (layer 1)
    slab.set_slab_longitudinal_rebar_top(d_b1=10 * mm, s_b1=20 * cm)

    # Two-layer bottom reinforcement: Ø12 @ 15 cm (layer 1) + Ø10 @ 20 cm (layer 2)
    slab.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=15 * cm,
                                         d_b3=10 * mm, s_b3=20 * cm)

3. Assigning Forces to the Slab
*******************************

Forces are applied through a `Node` object, the same way as for beams.
See the `Node` section for full details.

.. code-block:: python

    from mento import Forces, Node, kN, kNm

    f1 = Forces(label="ELU 1", V_z=30 * kN, M_y=25 * kNm)

    node = Node(section=slab, forces=[f1])

4. Performing Checks
********************

Once forces are assigned, call `check_shear()` and `check_flexure()` on the node.
*Mento* applies the appropriate design code formulas based on the `Concrete` type.

.. code-block:: python

    # Shear check
    node.check_shear()

    # Flexure check
    node.check_flexure()

See the `Node` section for details on interpreting the returned DataFrames.

5. Design the Section
*********************

If no reinforcement is assigned, *Mento* can design flexure and shear automatically.

.. code-block:: python

    node.design_flexure()
    node.design_shear()

The design is written back in the slab's own terms: a diameter and a spacing per layer,
per face. Read it through ``reinforcement`` (what the slab carries, at any time) or
``flexure_design`` (what a check demanded of it), rather than from the private attributes:

.. code-block:: python

    bottom = slab.reinforcement.bottom

    str(bottom)                 # 'Ø12 mm/17 cm'
    bottom.layers[0].d_b        # 12 mm
    bottom.layers[0].s.to("cm") # 17 cm
    bottom.n_bars               # 6 bars across the strip
    bottom.A_s.to("cm**2")      # 6.79 cm² placed

The spacing is rounded to the whole centimetre (inch, in imperial), so it is one that can
be detailed, and never to fewer bars than the design called for. See
:ref:`user_guide/design_results` for the rest of that view.

6. Jupyter Notebook Results
***************************

After performing checks or design, view formatted results with:

.. code-block:: python

    slab.results

The output includes longitudinal and shear reinforcement, applied forces, and DCR values.
See the `Node` section for more information on detailed results and Word report generation.

7. Detailed Results
*******************

See the `Node` section for how to display and save detailed per-load-case results using
`shear_results_detailed()`, `flexure_results_detailed()`, and their `_doc()` variants.
