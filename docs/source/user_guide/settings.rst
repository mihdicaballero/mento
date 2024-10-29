Settings Class
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

   **Default settings include**:

   * Clear cover: 25 mm
   * Clear spacing: 20 mm
   * Stirrup diameter: 8 mm
   * Longitudinal bar diameter: 16 mm
   * Vibrator size: 30 mm
   * Layers spacing: 25 mm

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

    from mento import Concrete_ACI_318_19, SteelBar, RectangularConcreteBeam
    from mento import psi, inch, ksi, mm

    # Define concrete and steel materials
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steel = SteelBar(name="ADN 420", f_y=60 * ksi)

    # Initialize section using default settings
    section = RectangularConcreteBeam(
        label="V-10x16",
        concrete=concrete,
        steel_bar=steel,
        width=10 * inch,
        height=16 * inch
    )

    # Check default settings
    print(section.settings.default_settings)

This section will automatically use the default settings for clear cover,
spacing, stirrup diameters, etc.

Example 2: Custom Settings
--------------------------

You can customize specific settings by passing a dictionary of values during
section initialization. In this example, we increase the clear cover and
modify the longitudinal bar diameter. You must call the properties by their
name correctly for this to work.

.. code-block:: python

    custom_settings = {'clear_cover': 50 * mm, 'longitudinal_diameter_ini': 25 * mm}

    # Create section with custom settings
    section = RectangularConcreteBeam(
        label="V-12x18",
        concrete=concrete,
        steel_bar=steel,
        width=12 * inch,
        height=18 * inch,
        settings=custom_settings
    )

    # Print the updated settings
    print(section.settings.settings)

Attributes
----------

- **default_settings**: Contains default design parameters such
  as clear cover, spacing, stirrup diameter, etc.
- **aci_318_19_settings**: Contains ACI 318-19 specific settings
  (e.g., reduction factors, minimum reinforcement considerations).
- **settings**: Current instance settings, which can be a mix of
  defaults and user-defined values.

Methods
-------

- **load_aci_318_19_settings()**: Loads the ACI 318-19 design settings.
- **get_setting(key)**: Retrieves the value of a specific setting.
- **set_setting(key, value)**: Sets the value of a specific setting.
- **add_settings(new_settings)**: Adds or updates multiple settings at once.
- **update(new_settings)**: Updates the current settings with a
  dictionary of new values.

Conclusion
----------

The `Settings` class provides a robust way to manage design parameters,
allowing users to work with defaults or customize their section properties.
Whether you're following standard design codes like ACI 318-19 or using
unique configurations, the `Settings` class ensures flexibility and
clarity in structural calculations.
