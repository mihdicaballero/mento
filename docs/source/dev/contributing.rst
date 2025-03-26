.. _dev/contributing:

Contributing to Mento
====================

Mento uses (and thanks):

- GitHub to host the code
- GitHub Actions to test all commits and PRs
- codecov to monitor test coverage
- ReadTheDocs to host the documentation
- black, ruff and mypy as code linter and pre-commit to enforce them.
- pytest to write tests
- Sphinx to write docs
- `CalcPad`_ for calculation validation
- GitHub Discussions for community support.

.. _CalcPad: https://github.com/Proektsoftbg/Calcpad

You can contribute in different ways:

Report Issues
-------------
You can report any issues with the package or documentation to the Mento issue tracker.
Also feel free to submit feature requests, comments, bugs or questions.

Contribute Code
---------------
To contribute fixes, code or documentation to Mento:

1. Open a Discussions thread to discuss the changes you want to make.
2. Fork Mento on GitHub
3. Submit changes using a pull request against the master branch
4. If submitting new code:

   - Add tests (see below)
   - Add documentation
   - Validate all calculations in CalcPad
   - Include the CalcPad validation files with your PR
5. Write "Closes #<bug number>" in the PR description or comment
6. Execute ``pre-commit run --all-files`` and resolve any issues

We won't merge a PR if:
- Tests are failing
- Calculations aren't validated in CalcPad against design code guides
- Code doesn't adhere to linting standards
- Documentation builds with errors

Setting Up Your Environment
---------------------------
For first-time contributors on Linux or OSX:

.. code-block:: bash

    $ git clone git@github.com:mihdicaballero/mento.git
    $ cd mento
    $ python -m virtualenv venv
    $ source venv/bin/activate
    $ pip install -e .
    $ pip install -r requirements_dev.txt
    $ pip install pre-commit
    $ pre-commit install

Writing Tests
------------
We use pytest for testing. When contributing code:

- For bug fixes: add a test to test_issues.py or amend existing tests
- For new features: add tests in the appropriate test file
- Prefer functions to classes for tests
- Use parametrize as much as possible
- Use fixtures instead of instantiating components directly
- Validate all calculations in CalcPad first
- Include the CalcPad validation files with your tests

Calculation Validation
---------------------
All new calculations must be validated in CalcPad before implementation:

1. Present and confirm a Validation example from design code guies or recognized books for structural engineering design
2. Create a CalcPad file demonstrating the calculation, with units
3. Include all edge cases and validation scenarios
4. Include both the CalcPad file with your PR
5. Reference these validations in your pytest modules

Running Tests and Building Docs
------------------------------

To review code:

.. code-block:: bash

    $ cd mento
    $ pre-commit run --all-files

To run tests:

.. code-block:: bash

    $ cd mento
    $ pytest

To build documentation:

.. code-block:: bash

    $ cd docs
    $ make html
