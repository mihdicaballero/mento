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
- **concrete** (*Concrete*): The concrete material used in the section.
- **steel_bar** (*SteelBar*): The steel reinforcement used in the section.
- **settings** (*Settings*): A dictionary-like object containing default and user-defined settings.
  - Includes settings like:
    - `clear_cover`: Clear cover for reinforcement.
    - `stirrup_diameter_ini`: Initial diameter of stirrups.
    - `longitudinal_diameter_ini`: Initial diameter of longitudinal reinforcement.

**Key Methods:**

- **update_settings(new_settings)**: Updates the section settings with a new dictionary of values.
- **get_settings**: Returns the current settings as a dictionary.

RectangularSection
------------------

The `RectangularSection` class extends `Section` to represent rectangular-shaped structural elements, such
as beams or columns. It adds specific geometric properties and methods relevant to rectangular sections.

**Key Attributes:**

- **width** (*PlainQuantity*): The width of the section.
- **height** (*PlainQuantity*): The height of the section.
- **A_x** (*PlainQuantity*): Cross-sectional area.
- **I_y** (*PlainQuantity*): Moment of inertia about the Y-axis.
- **I_z** (*PlainQuantity*): Moment of inertia about the Z-axis.

**Hierarchy:**

`RectangularSection` inherits all attributes and methods from `Section` while introducing geometry-specific properties.

Usage Notes
-----------

- **Settings Initialization**: 
  The `Section` class uses default settings unless explicitly updated. Ensure that any custom settings are 
  provided as a dictionary during initialization or by calling `update_settings`.
- **Section and RectangularSection Usage**:
  Since these are base classes, you won't be using them to create objects per se but to access its atrributes
  which other classes inherit, like RectangularBeam, RectangularColumn or CircularColumn.


Future Extensions
-----------------

The `Section` class is a foundational building block in the `mento` package. It can serve as the base for specialized structural elements such as:

- **Rectangular Beams**: See the `RectangularBeam` documentation for its extended features like shear and flexural checks.
- **Rectangular Columns**: Planned extension to handle axial and biaxial bending interactions for rectangular sections.
- **Circular Columns**: Planned extension to handle axial and biaxial bending interactions for circular sections.
