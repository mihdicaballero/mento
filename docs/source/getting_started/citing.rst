.. _getting_started/citing:

Citing mento
============

If mento supported a paper, a thesis, a technical note or a design report, please cite it.
Citations are what make an open source engineering tool visible to the people deciding
whether to trust it, and they are the main way this project gets credit for the work.

Which version to cite
---------------------

Always cite the version you actually ran. mento reports it:

.. code-block:: python

   import mento

   print(mento.__version__)

Calculations change between versions — a corrected clause or a new code edition can move a
result — so the version is part of the reproducibility of whatever you report.

Citation metadata
-----------------

The authoritative metadata lives in `CITATION.cff
<https://github.com/mihdicaballero/mento/blob/main/CITATION.cff>`_ at the root of the
repository. On GitHub, the **Cite this repository** button in the sidebar renders it as APA
or BibTeX for you, so you rarely need to write the entry by hand.

BibTeX
------

.. code-block:: bibtex

   @software{mento,
     author  = {Caballero, Mehdí and Romaris, Juan Pablo},
     title   = {mento: an intuitive tool for structural engineers to design concrete
                elements efficiently},
     version = {0.5.0},
     year    = {2026},
     url     = {https://github.com/mihdicaballero/mento}
   }

Replace ``version`` and ``year`` with the release you used.

.. note::

   Releases are being archived on `Zenodo <https://zenodo.org/>`_, which issues a DOI for
   each one. Once the first archived release is out, the DOI will be shown here and in
   ``CITATION.cff``, and it is preferable to the plain repository URL in a formal citation.

Citing a design code
--------------------

mento implements published standards; it does not replace them. When your work depends on a
particular clause, cite the standard itself — ACI 318-19, EN 1992-2004 or CIRSOC 201-2025 —
alongside mento, which is the implementation you used to apply it.
