.. _user_guide/table_style:

Table style
===========

The tables in a Word report are banded grey with a bold header, ruled above, below and
under the header row and nowhere else. ``set_table_style`` replaces that look for the rest
of the session, so it is written once, before the documents are produced:

.. code-block:: python

    import mento
    from mento import TableStyle

    mento.set_table_style(TableStyle(band_fill="F4F7FB", header_fill="2C6DC4",
                                     header_color="FFFFFF"))

Every Word document produced afterwards uses it:

.. code-block:: python

    from mento import Concrete_ACI_318_19, SteelBar, RectangularBeam, Node, Forces
    from mento import MPa, cm, mm, kN, kNm

    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="101", concrete=concrete, steel_bar=steel,
        width=20 * cm, height=60 * cm, c_c=25 * mm,
    )

    node = Node(section=beam, forces=[Forces(label="C1", V_z=80 * kN, M_y=100 * kNm)])
    node.design()
    node.shear_results_detailed_doc()   # a Word document in the new style

``mento.get_table_style()`` returns the style in force. One document can depart from it
without changing the setting, by building its own ``DocumentBuilder``:

.. code-block:: python

    from mento import DocumentBuilder, TableStyle

    builder = DocumentBuilder(title="Report", table_style=TableStyle(band_fill=""))

What can be set
---------------

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Field
     - Default
     - Meaning
   * - ``name``
     - ``"Mento Table"``
     - The style's name in Word's gallery.
   * - ``band_fill``
     - ``"F2F2F2"``
     - Fill of the banded rows. Empty for no banding.
   * - ``band_size``
     - ``1``
     - How many rows a band spans.
   * - ``header_fill``
     - ``""``
     - Fill of the header row. Empty for no fill.
   * - ``header_color``
     - ``""``
     - Text colour of the header row. Empty to use the body's.
   * - ``header_bold``
     - ``True``
     - Whether the header row is bold.
   * - ``text_color``
     - ``"1A1A1A"``
     - Text colour of the table.
   * - ``border_color``
     - ``"404040"``
     - Colour of every rule the style draws.
   * - ``top_pt``
     - ``1.0``
     - Rule above the header.
   * - ``header_rule_pt``
     - ``1.0``
     - Rule between the header and the first data row.
   * - ``bottom_pt``
     - ``1.0``
     - Rule closing the table.
   * - ``inside_pt``
     - ``0.0``
     - Rule between data rows.
   * - ``vertical_pt``
     - ``0.0``
     - Rule between columns.

Colours and thicknesses
-----------------------

A **colour** is six hexadecimal digits with no ``#``. Word does not object to ``"#F2F2F2"``;
it discards the fill and renders the table without it, so ``TableStyle`` raises
``ValueError`` instead:

.. code-block:: python

    TableStyle(band_fill="#F2F2F2")   # ValueError
    TableStyle(band_fill="F2F2F2")    # a light grey

An empty colour means the field is left unset, and whatever is behind it shows through:
``header_fill=""`` is a header with no fill, ``text_color=""`` takes the document's text
colour.

A **thickness** is given in points. Word draws no rule thinner than 0.25 pt and none
thicker than 12 pt, and anything outside that range is drawn at the nearest end of it. A
thickness of exactly ``0`` is not a hairline but the absence of a rule, which is why the
default table has no lines between its rows: the banding already separates them.

Why a style and not formatting
------------------------------

The look is written into the document once, as a table style, and every table points at
it. Two things follow from that. A table keeps its look when it is edited — a row added in
Word is banded, because the banding is a rule in the document rather than a fill painted
on the rows that existed when it was written. And a report of forty tables carries one
definition, so the choice can be changed in one place, including by hand in Word.

The pass/fail column of a summary is the exception: green and red are applied to those
cells directly, which outranks the style, so the verdict keeps its colour whatever the
table style says.
