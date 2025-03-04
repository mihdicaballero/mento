Node
==========

The `Node` class is a fundamental component in the **Mento** library. It represents a node in a structural model, associating a `Section` object with one or more `Forces` objects. This allows you to apply and manage forces on a specific section of your structure.

A `Node` connects a `Section` (e.g., a rectangular concrete beam or column) with the forces acting on it. You can:

- Add one or more forces to the node.
- Retrieve the list of forces applied to the node.
- Reset all forces to zero.
- Check for shear or flexure.
- Design for shear or flexure.

This guide will walk you through how to use the `Node` class effectively.
A `Section` could be a `Beam`, `Slab` or `Column` object since it just connects the section with the forces. 

1. Attributes
********************

The `Node` class has the following attributes:

- **section**: The `Section` object associated with this node. For now, a `RectangularConcreteBeam` is available.
- **forces**: A list of `Forces` objects applied to this node.

To create a `Node` you must define the concrete and steel materials, the concrete section based on those materials and the forces object.

.. code-block:: python

    from mento import Concrete_ACI_318_19, SteelBar, Forces, RectangularBeam, Node
    from mento import mm, cm, kN, MPa, kNm

    # Defines materials and beam section
    conc  = Concrete_ACI_318_19(name="C25",f_c=25*MPa)
    steel = SteelBar(name="ADN 420", f_y=420*MPa)
    beam  = RectangularBeam(label="101",concrete=conc,steel_bar=steel,width=20*cm, height=60*cm)

    # Define forces
    f1 = Forces(label='1.4D', V_z=50*kN, M_y=90*kNm)
    f2 = Forces(label='1.2D+1.6L', V_z=55*kN, M_y=-80*kNm)
    # Assign transverse and longitudinal reinforcement
    beam.set_transverse_rebar(n_stirrups=1, d_b=10*mm, s_l=20*cm)
    beam.set_longitudinal_rebar_bot(n1=3,d_b1=16*mm)
    beam.set_longitudinal_rebar_top(n1=3,d_b1=12*mm)

    # Create node and assign beam section and list of forces
    node_1 = Node(section=beam, forces=[f1, f2])

2. Performing Checks
********************

The `Node` class provides the following methods for checking, designing and printing results. They all apply for beams and columns.
Once the `Section` is defined and forces are assigned in a `Node` object, you can perform checks for shear and flexure.

*Mento* will apply corresponding design code formulas depending on the type of `Concrete` object 
created for all the forces assigned and store the limiting case for shear and top and bottom bending moment. 

- **Shear Check**: Use `check_shear()`.
- **Flexure Check**: Use `check_flexure()`.

.. code-block:: python

    # Perform shear and flexure checks in the Node section
    node_1.check_shear()
    node_1.check_flexure()

3. Design the section
********************

If you don't assign transverse or longitudinal rebar, you can ask *Mento* to design for shear and flexure.
*Mento* will apply corresponding design code formulas depending on the type of `Concrete` object 
created for all the forces assigned and store the limiting case for shear and top and bottom bending moment. 

After performing the numerical design, *Mento* will assign the optimal rebar combination of stirrups and longitudinal reinforcement,
considering all the longitudinal positions in the beam. This optimizations takes into account the default settings 
for longitudinal rebar limitation (vibrator size, maximum rebar diameter difference) and engineering criteria to suggest
the best rebar configuration balancing the amount of rebars, layers and different diamaters.

- **Shear Design**: Use `design_shear()`.
- **Flexure Design**: Use `design_flexure()`.

.. code-block:: python

    # Perform shear and flexure checks
    beam.design_shear()
    beam.design_flexure()

*Mento* will also create a Pandas DataFrame with all the check results for each Load Case in the Force object assigned to the Node, both for shear and flexure analysis.
You can print those results from the previous method.

.. code-block:: python

    # Perform shear and flexure checks
    shear_results = beam.design_shear()
    print(shear_results)
    flexure_results = beam.design_flexure()
    print(flexure_results)

4. Jupyter Notebook Results
******************

After performing the checks, you can view the results in a formatted way in a Notebook.

When you run `node.results`, the output includes for a `Beam` object:

- **Top and bottom longitudinal reinforcement**.
- **Shear reinforcement**.
- **Applied moments and shear forces**.
- **Design capacity ratios (DCR)**.
- **Warnings** (if any).

The output is formatted using LaTeX math notation for clarity and precision.
See the `Beam` or `Column` section for more information on how to display results.
The results are presented in a user-friendly format, with color-coded Demand-Capacity Ratios (DCR) for quick assessment.

5. Detailed Results
*******************

For more detailed results, you can use the following methods:

- **Shear Results**: Use `shear_results_detailed()`.
- **Flexure Results**: Use `flexure_results_detailed()`.

These methods provide a comprehensive breakdown of the calculations, which can be useful for reporting or further analysis. 
This reuslts will print in the Terminal or in a Jupyter Notebook the same way.

.. code-block:: python

    # View detailed shear results
    node_1.shear_results_detailed()
    # View detailed flexure results
    node_1.flexure_results_detailed()

You can also print the detailed results for the analysis of a specific corce if you pass it as an input:
.. code-block:: python

    # Print detailed results for shear, for specific force
    node_1.shear_results_detailed(f1)
    # Print detailed results for flexure, for specific force
    node_1.flexure_results_detailed(f1)

If you want to save the detailed results as a report in Microsoft Word, just run the following methods instead:

.. code-block:: python

    # View detailed shear results
    node_1.shear_results_detailed_doc()
    # View detailed flexure results
    node_1.flexure_results_detailed_doc()

6. Force methods
********************

The `Node` class provides the following methods for changing the Forces assigned to it:

- **Add a force**: Use `add_forces(forces)`.
- **Get list of forces**: Use `get_forces_list()`.
- **Reset forces to zero**: Use `reset_forces()`.

.. code-block:: python

    # Add a single force
    node_1.add_forces(force1) 
    # Add multiple forces
    node_1.add_forces([force2, force3])  
    # Get forces list.
    forces_list = node_1.get_forces_list()
    print(forces_list)
    # Reset forces
    node_1.reset_forces()