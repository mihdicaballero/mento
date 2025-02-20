Materials
===================

The `mento` package offers a wide range of material models,
allowing structural engineers to easily define and work with
common materials such as concrete and steel. Each material is
equipped with key mechanical properties and is compliant with
various international design codes.

Concrete Models
---------------

Concrete is one of the primary materials used in structural
engineering. In `mento`, you can work with various concrete
types, each compliant with different design codes:

* **Concrete_ACI_318_19**: Complies with the American Concrete
  Institute (ACI) 318-19 design code.
* **Concrete_EN_1992_2004**: Complies with the Eurocode 2 (EN 1992) design code.
* **Concrete_CIRSOC_201_25**: Complies with the CIRSOC 201-25 design code.

.. note::
   Create a material considering the design code that you want to use. The method for checking or designing elements won't change it's name. So just change the code and you are done.


Concrete properties depend in the design code selected. For example:

* **Compressive strength (f_c)**: The characteristic strength of
  concrete in compression for ACI, in MPa or psi units.
* **Density**: Concrete density (typically around 2500 kg/m³) for
  every code.
* **Elastic Modulus (E_c or E_cm)**: Secant modulus of elasticity
  derived from concrete strength for ACI or EC2.
* **Tensile strength (f_r or f_ctm)**: Concrete’s tensile strength
  calculated as per the design code ACI or EC2.
* **Beta_1**: A factor for the equivalent rectangular stress block,
  specific to ACI.

.. note::
   The metric system of units or imperial system of units will be automatically set based on the unit for the concrete strength. 

Example: Creating ACI Concrete
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To create concrete compliant with the ACI 318-19 standard with a
25 MPa strength:

.. code-block:: python

    from mento import Concrete_ACI_318_19, MPa

    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    print(concrete)

Steel Models
------------

Steel materials are commonly used for reinforcement bars and
prestressed tendons. `mento` provides models for these steel types:

- **SteelBar**: Reinforcing steel bars.
- **SteelStrand**: Prestressed steel strands.

Key properties of steel include:

* **Yield strength (f_y)**: The stress at which the material
  begins to deform plastically.
* **Modulus of elasticity (E_s)**: Steel’s elastic modulus,
  typically around 200 GPa.
* **Density**: Steel’s density, generally taken as 7850 kg/m³.

Example: Creating Reinforcing Steel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To create reinforcing steel with a yield strength of 420 MPa:

.. code-block:: python

    from mento import SteelBar, MPa

    steel_bar = SteelBar(name="ADN 420", f_y=420 * MPa)
    print(steel_bar)

Accessing Material Properties
-----------------------------

Each material class provides a `get_properties` method that
returns a dictionary of key properties such as:

- **Concrete**: Compressive strength, tensile strength, modulus of elasticity.
- **Steel**: Yield strength, elastic modulus, and density.

Simply call this method to access material attributes. If you `print()` 
a material you will get a string output with all it's properties.

.. note::
   All material properties in `mento` are automatically converted to appropriate units (MPa, kg/m³, etc.) based on the design code selected.

Getting Started
---------------

To get started with the materials in `mento`, simply import
the material class, create an instance, and access its properties.
Make sure to select the correct design code for your project.

For more details, consult the relevant structural design codes.
