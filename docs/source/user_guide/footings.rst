Footing
=======

The `Footing` class models the **section** of a spread footing or raft — a one-way
slab bearing directly on the ground.

It is a `OneWaySlab` in every respect but the rules the design codes write differently
for a member on the ground: the minimum longitudinal reinforcement, the spacing the bars
are detailed at, and the thickness the section is expected to have. Everything else —
geometry, materials, reinforcement given as diameter and spacing, flexure and shear
checks through a `Node` — is the slab's, unchanged.

.. warning::

    **Sectional design only. Mento does no geotechnical calculation.**

    `Footing` answers one question: what reinforcement does this section need for
    the forces it is given. It knows nothing about the soil under it and does not
    check bearing pressure, settlement, sliding, overturning, uplift, or the plan
    dimensions the footing needs to spread its load. It does not size the footing.

    Those are the engineer's, and they come first: the plan size and the thickness
    are decided from the soil, and only then is the resulting section handed to
    `Footing`. Punching shear is a separate check — see
    :class:`~mento.punching.PunchingSlab`.

Key Concepts
------------

- **Geometry**: a strip of the footing, given as a `width` and a `height` (the total
  thickness). Results are the totals across that width, so a 1 m strip reads directly
  as per-metre values.
- **Forces**: the sectional forces at the face being designed — the moment and shear
  the soil pressure produces in the strip. Mento does not derive them from the soil;
  they are an input.
- **Reinforcement**: bar diameter and spacing, one direction at a time. A footing spans
  both ways, so it is designed as two strips, one per direction.
- **Design code**: taken from the `Concrete` object, as everywhere else in mento.
  ACI 318-19, CIRSOC 201-25 and EN 1992-2004 are all supported.

Usage
-----

1. Creating a Footing
*********************

.. code-block:: python

    from mento import Concrete_ACI_318_19, Footing, SteelBar, MPa, cm, mm, m

    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)

    # A 1 m strip of a footing 60 cm thick
    footing = Footing(label="Z1", concrete=concrete, steel_bar=steel,
                      width=1 * m, height=60 * cm, c_c=50 * mm)

2. Checking and Designing
*************************

Exactly as for a slab: forces are attached through a `Node`, which drives the check
and the design.

.. code-block:: python

    from mento import Forces, Node, kN, kNm

    node = Node(section=footing, forces=[Forces(label="ELU 1", M_y=120 * kNm, V_z=90 * kN)])

    node.design()          # flexure + shear
    node.check_flexure()   # per-combination table
    node.check_shear()

Reinforcement can also be assigned by hand with `set_slab_longitudinal_rebar_bot()` and
`set_slab_longitudinal_rebar_top()`, and read back afterwards, the same way as for a
`OneWaySlab`. See :doc:`slabs` for the full reinforcement and results API.

What the codes do differently
-----------------------------

Minimum reinforcement
*********************

A member spanning between supports is given a minimum sized so that it cannot fail the
instant it cracks. A member on the ground cannot fail that way — the soil goes on
carrying it — so both codes exempt it and put a different rule in its place:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Code
     - A\ :sub:`s,min`
   * - ACI 318-19,
       CIRSOC 201-25
     - §9.6.1.1(b) grants the exemption and §13.3.1.2 replaces the flexural minimum
       with the shrinkage and temperature reinforcement of Table 24.4.3.2:
       :math:`0.0018\,b\,h` at :math:`f_y = 420` MPa, scaling as
       :math:`0.0018 \cdot 420 / f_y` and never below :math:`0.0014`. Written on the
       **gross** section, so the effective depth does not enter it.
   * - EN 1992-2004
     - The larger of the halved geometric minimum of a foundation
       (:math:`0.0010\,b\,h` at :math:`f_{yk} = 400` MPa, :math:`0.0009\,b\,h` at
       500 MPa) and the crack-control minimum of §7.3.2(2), Eq. (7.1),
       :math:`k_c\,k\,f_{ct,eff}\,A_{ct} / \sigma_s`. The second usually governs a
       thick footing.

The design returns the largest applicable minimum **already applied**, so the
reinforcement it reports needs no correction afterwards.

Bar spacing
***********

Kept between **100 and 300 mm**. The 300 mm cap replaces the 3h of a suspended slab,
which stops binding on a section this thick; the 100 mm floor is the spacing below which
a footing is not detailed in practice. A design answers a heavier demand with a larger
bar rather than with bars closer than 100 mm, and the flexure check table reports the
spacing against both bounds.

Thickness
*********

A section thinner than the code asks of a footing on soil — 200 mm in ACI 318-19
§13.3.1.2, 250 mm in EN practice — warns when it is built:

.. code-block:: python

    Footing(label="Z2", concrete=concrete, steel_bar=steel,
            width=1 * m, height=18 * cm, c_c=50 * mm)
    # UserWarning: Footing Z2 is 18 cm thick, below the 200 mm ACI 318-19 asks
    # of a footing on soil. It is designed as given.

It is advice, not a limit. The section is still designed as given, because the thickness
is the engineer's to choose.
