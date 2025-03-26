What is Mento?
--------------

Mento is a Python package designed to simplify the structural
analysis and design of concrete elements.
It provides tools that are intuitive for structural engineers,
enabling efficient design checks and generation of detailed reports.
The package supports various design codes, allowing for flexible
application in different regions and standards.

Mento can handle the design and analysis of:

- **Rectangular concrete beams** for flexure and shear.

In the future Mento will support:
- **Circular and rectangular concrete columns** for different structural needs.

Some key features of Mento include:

- **Unit-sensitive design**: Variables can be input with their respective units for accurate calculations.
- **Interactive usage**: Mento integrates seamlessly with Jupyter Notebooks, allowing engineers to build custom workflows using its modules.
- **Results in Markdown and DataFrames**: Results are provided in markdown format and as Pandas DataFrames, facilitating the handling and presentation of multiple design checks.
- **Report generation**: Mento can generate detailed reports in LaTeX format, making it easy to document the results of the analysis.

Mento is thoroughly tested for compliance with major design codes such as **ACI 318-19**, **EN 1992**, and **CIRSOC 201-2005**, ensuring reliable results that meet industry standards.

Using Mento is easy and intuitive:

.. code-block:: python

    from mento import Concrete_ACI_318_19, SteelBar, RectangularBeam
    from mento import cm, MPa

    # Define concrete and steel properties
    conc= Concrete_ACI_318_19(name="C25",f_c=25*MPa)
    steel= SteelBar(name="ADN 420", f_y=420*MPa)
    # Define beam section
    section = RectangularBeam(
            label="101",concrete=conc,steel_bar=steel,width=20*cm, height=40*cm)
    # Display data with LaTeX formatting in a Jupyter Notebook
    section.results

Expected Output
---------------

In a Jupyter notebook, this will display the beam data with LaTeX-style formatting.
Below is the expected output:

.. math::

   \textsf{Beam 101}, \, b = 20.00 \, \textsf{cm}, \, h = 40.00 \, \textsf{cm}, \, c_{\text{c}} = 2.50 \, \textsf{cm}, \, \textsf{Concrete C25}, \, \textsf{Rebar ADN 420}.

This is an ideal way to present structural data, making the results clear and easy to read.
The use of Jupyter Notebooks and LaTeX ensures that all units and parameters are well-formatted for structural engineering reports.

Design Principles
-----------------

Mento was developed to meet the needs of structural engineers for a flexible, code-compliant, and user-friendly design tool. The package is built with the following principles in mind:

- **Code-compliant checks**: Mento supports multiple design codes and is easily extensible to add more in the future.
- **Modular design**: The package’s design ensures that engineers can mix and match different sections for custom analyses.
- **Integration with Pandas and Word**: Mento generates reports and tables that can be exported into different formats, supporting data analysis and documentation needs.

For more detailed help getting started, see the :ref:`User Guide <user_guide/index>` and explore :ref:`Examples <examples/index>`.
