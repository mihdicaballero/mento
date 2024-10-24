Units in `mento` Package
=========================

The `mento` package is designed to work with both Metric (SI) and
Imperial unit systems. The units are managed using the `Pint` library,
which allows for seamless unit conversions and arithmetic.

Importing units
--------------------

Units are imported from mento directly, you just type what you
want to use. For example:

.. code-block:: python

    from mento import mm, cm, kN, MPa, inch, ft, ksi

Usage within `mento`
--------------------

The units in `mento` are fully compatible with the `Pint` unit
registry, allowing for easy conversions between different systems.
You can add values with different units of the same type, and change
to other units using the `.to()` method.

The example below adds attributes in different unit imported
from mento and then changes it to feet.

.. code-block:: python

    from mento import cm, m
    a = 2*m
    b = 15*cm
    c = a + b
    print(a, b, c)
    d = c.to('ft')
    print(d)

Available Units
---------------

The following units are available in `mento`:

* **Metric (SI)**:

  * Length: `m`, `cm`, `mm`
  * Force: `kN`
  * Moment: `kNm`
  * Stress: `MPa`, `GPa`
  * Mass: `kg`
  * Time: `sec`

* **Imperial**:

  * Length: `inch`, `ft`
  * Force: `lb`, `kip`
  * Stress: `psi`, `ksi`

Unit Formatting
---------------

The units in `mento` follow a standardized formatting scheme to
ensure clarity in outputs, showing usually 2 decimals for any attribute.
You can format and display units using the `~P` specifier in
the `Pint` library. Here’s an example:

.. code-block:: python

    from mento import kN
    F = 15.2354*kN
    print(F)
    F2 = f"Force is {F:.4f~P}" # Output: "Force is 15.2354 kN"
    print(F2)

For more advanced usage, refer to the `Pint` library documentation
to explore additional formatting options.
