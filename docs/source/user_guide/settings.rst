Settings
==============

The `BeamSettings` class is used to define and manage configuration options
for structural sections such as beams.
It holds default values for various parameters and allows users to
modify these settings to suit specific design needs.
Each `Section` can have its own unique settings, either using
the default settings or customized values provided by the user.

Overview
--------

When designing concrete elements, certain design parameters (e.g., maximum clear spacing) are commonly required. The `BeamSettings` class
provides a structured way to define these parameters, ensuring consistency
across sections while allowing flexibility for customization.

Key Features:

- **Default settings**: Pre-defined values for typical design cases.
- **Customization**: Users can override default settings with specific
  values for individual `Sections`.
- **Integration with section classes**: Each structural `Section` object
  has its own instance of `BeamSettings`, allowing for unique configurations
  per `Section`.

Usage Scenarios
---------------

1. **Using Default Settings**:

   If no custom settings are provided, a section will use default settings.
   These settings are suitable for common design cases and ensure design
   calculations can proceed without requiring additional input from the user.
   All these settings are used when designing the reinforcement of beams.

   Defaut settings change based on the metric system selected.

   **Default settings for metric system include**:

   * Clear spacing: 20 mm
   * Stirrup diameter: 8 mm
   * Minimum longitudinal bar diameter: 8 mm
   * Vibrator size: 30 mm
   * Layers spacing: 25 mm
   * Maximum diameter difference: 5 mm
   * Maximum bars per layer: 5

   **Default settings for imperial system include**:

   * Clear spacing: 1 inch
   * Stirrup diameter: 3/8 inch
   * Minimum longitudinal bar diameter: 3/8 inch
   * Vibrator size: 1.25 inch
   * Layers spacing: 1 inch
   * Maximum diameter difference: 2/8 inch
   * Maximum bars per layer: 5

.. note::
    The default values are used to calculate effective heights when designing a concrete section,
    before knowing which diameters bars will be final. Also for designing user-tailored rebar layouts in beams.
    The metric or imperial initial settings are selected based on the unit of the concrete resistance `f_c` provided when defining the concrete material.

2. **Customizing Settings**:

   Users can provide custom settings when initializing a section or later by
   updating specific values. For example, if the vibrator size
   differs from the defaults, these can be adjusted accordingly.

Example 1: Default Settings
---------------------------

In this example, we create a section using the default settings:

.. code-block:: python

  from mento import Concrete_ACI_318_19, SteelBar, RectangularBeam
  from mento import psi, inch, ksi, mm

  # Define concrete and steel materials
  concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
  steel = SteelBar(name="ADN 420", f_y=60 * ksi)

  # Initialize section using default settings
  section = RectangularBeam(
      label="101",
      concrete=concrete,
      steel_bar=steel,
      width=10 * inch,
      height=16 * inch,
      c_c = 25 * mm
  )

  # Check default settings
  print(section.settings)

This section will automatically use the default settings for
spacing, stirrup diameters, etc. for the corresponding unit system depending on the unit of the concrete resistance.
If `f_c` is indicated in psi or ksi, the default unit system will be `imperial` and if it is indicated in
Pa or MPa, it will be `metric`.

Example 2: Custom Settings
--------------------------

You can customize specific settings by defining specific attributes during
BeamSettings class initialization. In this example, we increase the clear spacing and
modify the minimum longitudinal bar diameter.

.. code-block:: python

  from mento import BeamSettings
  settings = BeamSettings(clear_spacing= 40 * mm, minimum_longitudinal_diameter = 12 * mm)

  # Create section with custom settings
  section = RectangularBeam(
      label="101",
      concrete=concrete,
      steel_bar=steel,
      width = 20 * cm,
      height = 60 * cm,
      c_c = 25 * mm,
      settings=settings
  )

  # Print the updated settings
  print(section.settings)

Attributes
----------

The attributes of the settings class are as follows:
  - clear_spacing
  - stirrup_diameter_ini
  - vibrator_size
  - layers_spacing
  - max_diameter_diff
  - minimum_longitudinal_diameter
  - max_bars_per_layer
