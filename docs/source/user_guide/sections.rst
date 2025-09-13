Sections
==============================

The `mento` package provides a structured way to model structural sections. This documentation focuses on
the foundational `Section` class and its specialized subclass `RectangularSection`.

These classes represent the geometric and material properties of structural elements, forming the basis
for specific implementations like beams and columns.

Section
-------

The `Section` class serves as the base class for all structural sections. It includes core attributes
and methods for managing materials, settings, and unique identification for each section.

**Key Attributes:**

- **id** (*int*): A unique identifier automatically assigned to each section instance.
- **label** (*String*): The label of the section.
- **concrete** (*Concrete*): The concrete material used in the section.
- **steel_bar** (*SteelBar*): The steel reinforcement used in the section.
- **c_c** (*Quantity*): The clear cover of the section.

RectangularSection
------------------

The `RectangularSection` class extends `Section` to represent rectangular-shaped structural elements, such
as beams or columns. It adds specific geometric attributes and properties relevant to rectangular sections.

**Attributes:**

- **width** (*Quantity*): The width of the section.
- **height** (*Quantity*): The height of the section.

**Properties:**

- **A_x** (*Quantity*): Cross-sectional area.
- **I_y** (*Quantity*): Moment of inertia about the Y-axis.
- **I_z** (*Quantity*): Moment of inertia about the Z-axis.

**Hierarchy:**

`RectangularSection` inherits all attributes and methods from `Section` while introducing geometry-specific properties.

Usage Notes
-----------

- **Section and RectangularSection Usage**:
  Since these are base classes, you won't be using them to create objects per se but to access its atrributes
  which other classes inherit, like RectangularBeam.

You can access the section properties directly in a `RectangularBeam` object.

.. code-block:: python
  section = RectangularBeam(
      label="101",
      concrete=concrete,
      steel_bar=steel,
      width=20 * cm,
      height=50 * cm,
      c_c = 2.5 * cm
  )
  print('b =', section.width)
  print('h =', section.height)
  print('A_x =', section.A_x)
  print('I_y =', section.I_y)
  print('I_z =', section.I_z)


The `Section` class is a foundational building block in the `mento` package. It can serve as the base for specialized structural elements such as:

- **Rectangular Beams**: See the `RectangularBeam` documentation for its extended features like shear and flexural checks.
