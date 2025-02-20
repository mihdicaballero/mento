Node
==========

The `Node` class is a fundamental component in the **Mento** library. It represents a node in a structural model, associating a `Section` object with one or more `Forces` objects. This allows you to apply and manage forces on a specific section of your structure.

A `Node` connects a `Section` (e.g., a beam or column) with the forces acting on it. You can:
- Add one or more forces to the node.
- Retrieve the list of forces applied to the node.
- Reset all forces to zero.

This guide will walk you through how to use the `Node` class effectively.
A `Section` could be a `Beam`, `Column`, `Slab` or `Footing` object since ti just connects the section with the forces. 

Creating a Node
---------------

To create a `Node`, you need a `Section` object and one or more `Forces` objects. Here's how you can do it:

.. code-block:: python

    from mento.section import Section
    from mento.forces import Forces
    from mento.node import Node

    # Create a Section and Forces object
    section = Section()
    force1 = Forces()
    force2 = Forces()

    # Create a Node with a single force
    node = Node(section, force1)

    # Add another force to the node
    node.add_forces(force2)


Attributes
----------

The `Node` class has the following attributes:

- **section**: The `Section` object associated with this node.
- **forces**: A list of `Forces` objects applied to this node.


Methods
-------

The `Node` class provides the following methods:

### `add_forces(forces)`
Adds one or more `Forces` objects to the node.

- **Parameters**:
  - `forces`: A single `Forces` object or a list of `Forces` objects.
- **Example**:

  .. code-block:: python

      node.add_forces(force1)  # Add a single force
      node.add_forces([force2, force3])  # Add multiple forces

### `get_forces_list()`
Returns the list of `Forces` objects applied to this node.

- **Returns**: A list of `Forces` objects.
- **Example**:

  .. code-block:: python

      forces_list = node.get_forces_list()
      print(forces_list)

### `reset_forces()`
Resets all forces applied to this node to zero.

- **Example**:

  .. code-block:: python

      node.reset_forces()


Example: Working with Nodes
---------------------------

Here’s a complete example that demonstrates how to use the `Node` class:

.. code-block:: python

    from mento.section import Section
    from mento.forces import Forces
    from mento.node import Node

    # Create a Section and Forces objects
    section = Section()
    force1 = Forces()
    force2 = Forces()

    # Create a Node with a single force
    node = Node(section, force1)

    # Add another force to the node
    node.add_forces(force2)

    # Get the list of forces applied to the node
    forces_list = node.get_forces_list()
    print("Forces applied to the node:", forces_list)

    # Reset all forces in the node
    node.reset_forces()
    print("Forces after reset:", node.get_forces_list())


Tips and Best Practices
----------------------

- **Reuse Nodes**: If you have multiple sections with similar forces, consider reusing `Node` objects to simplify your code.
- **Check Input Types**: Always ensure that the `section` and `forces` passed to the `Node` constructor are of the correct types (`Section` and `Forces`, respectively).
- **Reset Forces When Needed**: Use the `reset_forces()` method to clear forces before applying new ones.