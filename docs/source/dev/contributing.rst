.. _dev/contributing:

Contributing to Mento
=====================

Contributions of every size are welcome, from fixing a typo in these docs to implementing a
whole new design code.

The full contributing guide lives in the repository, at `CONTRIBUTING.md
<https://github.com/mihdicaballero/mento/blob/main/CONTRIBUTING.md>`__. This page summarises
what you need to get started locally. By participating you agree to abide by our
:ref:`Code of Conduct <dev/codeofconduct>`.

Mento uses (and thanks):

- GitHub to host the code
- GitHub Actions to test all commits and PRs
- codecov to monitor test coverage
- ReadTheDocs to host the documentation
- ruff for linting and formatting, mypy for type checking, and pre-commit to enforce them
- pytest to write tests
- Sphinx to write docs
- pint for units
- `CalcPad`_ for calculation validation
- GitHub Discussions for community support.

.. _CalcPad: https://github.com/Proektsoftbg/Calcpad

Ways to contribute
------------------

**Report an issue.** Bugs, wrong results, documentation problems and feature requests all go
to the `issue tracker <https://github.com/mihdicaballero/mento/issues>`__, which has
templates for each. If a design result looks wrong, include the inputs, the value Mento
produced, and the expected value with its source.

**Ask a question.** Use `Discussions
<https://github.com/mihdicaballero/mento/discussions>`__.

**Contribute code.** Issues labelled ``good first issue`` and ``help wanted`` are the easiest
place to start. For anything substantial, open a Discussion first so we can agree on the
approach before you write code.

Setting up your environment
---------------------------

Mento requires Python 3.10 or newer.

.. code-block:: bash

    $ git clone https://github.com/mihdicaballero/mento.git
    $ cd mento
    $ python -m venv venv
    $ source venv/bin/activate      # Windows: venv\Scripts\activate
    $ pip install -e .
    $ pip install -r requirements_dev.txt
    $ pre-commit install

Development workflow
--------------------

1. Create a branch off ``main``.
2. Make your changes, with tests and documentation alongside the code.
3. Run the checks below until they are clean.
4. Open a pull request against ``main``, writing ``Closes #<issue number>`` in the
   description.

We will not merge a pull request if tests fail, if calculations are not validated against a
recognised source, if linting or type checking does not pass, or if the documentation builds
with errors.

Running the checks
------------------

.. code-block:: bash

    $ pytest                        # test suite, with coverage
    $ ruff check . --fix            # lint
    $ ruff format .                 # format
    $ mypy mento/                   # strict type checking
    $ pre-commit run --all-files    # everything the hooks enforce

To build the documentation:

.. code-block:: bash

    $ pip install -r docs/requirements_docs.txt
    $ cd docs
    $ make html

Writing tests
-------------

- For bug fixes, add a regression test that fails before your change and passes after it.
- For new features, add tests in the test module matching the source module — code in
  ``mento/beam.py`` is tested in ``tests/test_beam.py``.
- Prefer plain functions to classes for tests.
- Use ``pytest.mark.parametrize`` for families of similar cases, and fixtures instead of
  constructing the same beam or material repeatedly.
- Compare physical quantities with an explicit tolerance rather than exact equality.

Calculation validation
----------------------

Every contribution that touches a design calculation must be validated before it is merged.
This is what makes Mento trustworthy.

1. Find a worked reference example: the design code itself, an official design guide, or a
   recognised structural engineering textbook. Name it in your issue or pull request.
2. Reproduce it in `CalcPad`_, with units, including the edge cases the implementation needs
   to handle.
3. Include the CalcPad file with your pull request.
4. Reference the validation source in your test module, so the next reader can trace where
   the expected numbers came from.

If no published example exists for what you are implementing, say so in the pull request and
we will work out an acceptable validation path together.
