Forces Class
============

The `Forces` class in the `mento` package allows users to define and manipulate the primary forces acting on a structural element, such as axial force, shear force, and bending moment. This class is designed for flexibility in defining forces, with the ability to adjust and retrieve them in different stages of a structural analysis or design workflow.

Key Concepts
------------

- **Axial Force (`N_x`)**: Force applied along the axis of the element, along the x-x axis.
- **Shear Force (`V_z`)**: Force acting perpendicular to the axis of the element, along the z-z local axis.
- **Bending Moment (`M_y`)**: The moment caused by forces that induce bending about the y-y axis.

These forces are defined using compatible units from the `Pint` library, specifically `kN` (kilonewtons) for forces and `kNm` (kilonewton meters) for moments.

Usage
-----

Below is a step-by-step guide on how to use the `Forces` class in your structural analysis workflows.

1. Creating a Forces Object
**********

To define the forces acting on a structural element, instantiate a `Forces` object. You can provide initial values for axial force (`N_x`), shear force (`V_z`), and bending moment (`M_y`). If no values are provided, the forces default to zero.

.. code-block:: python

    from mento import Forces, kN, kNm

    # Define a forces object with specific values
    forces = Forces(N_x=2*kN, V_z=10*kN, M_y=5*kNm)

2. Accessing Forces
**********

Once a `Forces` object is created, you can access individual forces using the respective properties `N_x`, `V_z`, and `M_y`. These properties will return the forces in kilonewtons (kN) or kilonewton meters (kN*m), making it easy to inspect the current state of forces in your structure.

.. code-block:: python

    print(forces.N_x)  # Output: 2.00 kN
    print(forces.V_z)  # Output: 10.00 kN
    print(forces.M_y)  # Output: 5.00 kN*m

3. Modifying Forces
**********

Forces can be modified at any point by calling the `set_forces()` method. This method allows you to update the values of `N_x`, `V_z`, and `M_y`.

.. code-block:: python

    # Update the axial and moment forces
    forces.set_forces(N_x=3*kN, M_y=7*kNm)

4. Retrieving Forces as a Dictionary
**********

You can retrieve the forces in the form of a dictionary for easy manipulation, storage, or reporting. The `get_forces()` method returns a dictionary where the keys are `N_x`, `V_z`, and `M_y`, with values corresponding to the respective forces.

.. code-block:: python

    forces_dict = forces.get_forces()
    print(forces_dict)
    # Output: {'N_x': 3.00 kN, 'V_z': 10.00 kN, 'M_y': 7.00 kN*m}

5. Assigning a Label to a Force
**********

Optionally, you can assign a label to a force object to describe the specific load condition or scenario (e.g., "Crane load", "Wind load"). This is useful in complex models where multiple forces are acting on different elements.

.. code-block:: python

    forces.label = "Crane load"
    print(forces.label)  # Output: Crane load

6. Force Object ID
**********

Each `Forces` object is automatically assigned a unique ID, which can be accessed through the `id` property. This is helpful when tracking multiple force objects in more complex analyses.

.. code-block:: python

    print(forces.id)  # Output: Unique ID (e.g., 1, 2, etc.)

Example Workflow
----------------

Here's a full example of how the `Forces` class could be used in a typical workflow:

.. code-block:: python

    from mento import Forces, kN, kNm

    # Create a new Forces object
    forces = Forces(N_x=2*kN, V_z=10*kN, M_y=5*kNm)

    # Check current values of forces
    print(forces.N_x)  # Output: 2.00 kN
    print(forces.V_z)  # Output: 10.00 kN
    print(forces.M_y)  # Output: 5.00 kN*m

    # Modify the forces
    forces.set_forces(N_x=3*kN, M_y=7*kNm)

    # Retrieve forces as a dictionary
    forces_dict = forces.get_forces()
    print(forces_dict)

    # Assign a label to the forces object
    forces.label = "Crane load"
    print(forces.label)

    # Check the unique ID assigned to this object
    print(forces.id)

This flexible interface ensures that you can easily manage forces during the design and analysis of structural elements, while maintaining clear and consistent units.
