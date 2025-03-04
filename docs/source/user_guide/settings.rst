Settings
==============

The `Settings` class is used to define and manage configuration options
for structural sections such as beams.
It holds default values for various parameters and allows users to
modify these settings to suit specific design needs.
Each `Section` can have its own unique settings, either using
the default settings or customized values provided by the user.

Overview
--------

When designing concrete elements, certain design parameters (e.g., clear
cover, stirrup diameter) are commonly required. The `Settings` class
provides a structured way to define these parameters, ensuring consistency
across sections while allowing flexibility for customization.

Key Features:

- **Default settings**: Pre-defined values for typical design cases.
- **Customization**: Users can override default settings with specific
  values for individual `Sections`.
- **Integration with section classes**: Each structural `Section` object
  has its own instance of `Settings`, allowing for unique configurations
  per `Section`.

Usage Scenarios
---------------

1. **Using Default Settings**:

   If no custom settings are provided, a section will use default settings.
   These settings are suitable for common design cases and ensure
   calculations can proceed without requiring additional input from the user.

   Defaut settings change based on the metric system selected.

   **Default settings for metric system include**:

   * Clear cover: 25 mm
   * Clear spacing: 20 mm
   * Stirrup diameter: 8 mm
   * Minimum longitudinal bar diameter: 8 mm
   * Vibrator size: 30 mm
   * Layers spacing: 25 mm
   * Maximum diameter difference: 5 mm
   * Maximum bars per layer: 5

   **Default settings for imperial system include**:

   * Clear cover: 1 inch
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


2. **Customizing Settings**:

   Users can provide custom settings when initializing a section or later by
   updating specific values. For example, if the clear cover or longitudinal
   bar diameter differs from the defaults, these can be adjusted accordingly.

3. **Loading Design Code Specific Settings**:

   The `Settings` class includes predefined configuration values that adhere
   to a specific design code. These settings, for example, can be loaded
   by calling `load_aci_318_19_settings()` on the settings object.
   This is done internally when a check or design method is called for
   a `Section`.

Example 1: Default Settings
---------------------------

In this example, we create a section using the default settings:

.. code-block:: python

  from mento import Concrete_ACI_318_19, SteelBar, RectangularBeam
  from mento import psi, inch, ksi

  # Define concrete and steel materials
  concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
  steel = SteelBar(name="ADN 420", f_y=60 * ksi)

  # Initialize section using default settings
  section = RectangularBeam(
      label="101",
      concrete=concrete,
      steel_bar=steel,
      width=10 * inch,
      height=16 * inch
  )

  # Check default settings
  print(section.settings.default_settings_imperial)
  print(section.settings.default_settings_metric)

This section will automatically use the default settings for clear cover,
spacing, stirrup diameters, etc. for the corresponding unit system depending on the unit of the concrete resistance.
If `f_c` is indicated in psi or ksi, the default unit system will be `imperial` and if it is indicated in
Pa or MPa, it will be `metric`.

Example 2: Custom Settings
--------------------------

You can customize specific settings by passing a dictionary of values during
section initialization. In this example, we increase the clear cover and
modify the minimum longitudinal bar diameter. You must call the properties by their
name correctly for this to work.

.. code-block:: python

  custom_settings = {'clear_cover': 50 * mm, 'minimum_longitudinal_diameter': 12 * mm}

  # Create section with custom settings
  section = RectangularBeam(
      label="101",
      concrete=concrete,
      steel_bar=steel,
      width=12 * inch,
      height=18 * inch,
      settings=custom_settings
  )

  # Print the updated settings
  print(section.settings)

Attributes
----------

- **default_settings_imperial**: Contains default design parameters such
  as clear cover, spacing, stirrup diameter, etc. for imperial units.
- **default_settings_metric**: Contains default design parameters such
  as clear cover, spacing, stirrup diameter, etc. for metric units.
- **ACI_318_19_settings**: Contains ACI 318-19 specific settings
  (e.g., reduction factors, minimum reinforcement considerations).
- **EN_1992_2004_settings**: Contains EN 1992-1-1:2004 specific settings.
- **settings**: Current instance settings, which can be a mix of
  defaults and user-defined values.

Methods
-------

- **load_aci_318_19_settings()**: Loads the ACI 318-19 design settings.
- **load_en_1992_2004_settings()**: Loads the EN 1992-1-1:2004 design settings.
- **get_setting(key)**: Retrieves the value of a specific setting.
- **set_setting(key, value)**: Sets the value of a specific setting.
