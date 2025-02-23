.. _examples/index:

Examples
==========

In this folder, you will find detailed examples of how to use mento.


Flexure examples
==========
In this example, we will determine the flexural reinforcement required for a section subjected to various flexural moments in accordance with ACI 318-19.

.. code-block:: python
    from mento import Concrete_ACI_318_19, SteelBar, RectangularConcreteBeam
    from mento import psi, inch, ksi
    from mento.section import Section
    from mento.forces import Forces
    from mento.node import Node

    # Define concrete and steel materials
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steel = SteelBar(name="fy=6000", f_y=60 * ksi)

    # Initialize section using default settings
    section = RectangularConcreteBeam(
        label="B-12x24",
        concrete=concrete,
        steel_bar=steel,
        width=12 * inch,
        height=24 * inch
    )

    # Initialize the flexural loads on the section
    f1 = Forces(label='Combo_01', M_y=400*kip*ft)
    f2 = Forces(label='Combo_02', M_y=-400*kip*ft)
    f3 = Forces(label='Combo_03', M_y=37*kip*ft)
    f4 = Forces(label='Combo_04', M_y=150*kip*ft)
    f5 = Forces(label='Combo_05', M_y=-1*kip*ft)

    # Create the node that associates the section with the loads
    Node(section=beam, forces=[f1,f2,f2,f4,f5])

    # Determine the required steel reinforcement for the beam:
    # Since the concrete was defined using Concrete_ACI_318_19, the design will follow ACI 318-19.
    flexure_results = beam.design_flexure()
    print(flexure_results)
