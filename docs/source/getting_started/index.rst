.. _getting_started/index:

Getting Started
===============

.. toctree::
   :hidden:
   :maxdepth: 2

   what_is_mento

The getting started guide aims to get you using mento productively as quickly as possible.


Installation
------------

Install with pip
^^^^^^^^^^^^^^^

You can install Mento via ``pip``, the Python package manager. Simply run the following command in your terminal:

.. code-block:: bash

   pip install mento

This will install the latest stable version of the Mento package, along with any necessary dependencies.

Try Mento without installation (Google Colab)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you don't have Python installed or want to quickly try Mento, you can use Google Colab:

1. **Open Google Colab**: Go to https://colab.research.google.com/
2. **Create a new notebook**: Click "File" > "New notebook"
3. **Install Mento**: In the first cell, paste and run this command:

   .. code-block:: bash

      !pip install mento

4. **Use Mento**: In the next cell, import the package and run an example from the documentation:

   You can find detailed examples in the `Examples` section of the documentation:

   - :ref:`Beam Summary <examples/beam_summary>`
   - :ref:`Rectangular Beam Check <examples/rectangular_beam_check>`

5. **Run your code**: Press :kbd:`Shift` + :kbd:`Enter` to execute each cell

.. note::
   Google Colab provides free access to Python runtime with all dependencies,
   so you can test Mento without any local setup.
