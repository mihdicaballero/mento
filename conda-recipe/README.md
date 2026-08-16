# conda-forge recipe

This directory holds the conda recipe for mento. It is **not** used to build anything in
this repository: conda-forge builds packages from its own feedstock, and this is the copy we
keep in step with `pyproject.toml` so that submitting or updating the recipe is a copy, not a
rewrite.

- `meta.yaml` — the recipe itself, pinned to a released version on PyPI and its sha256.
- `conda_build_config.yaml` — raises `python_min` to 3.10, mento's minimum.

The submission and update procedure is in
[CONTRIBUTING.md](../CONTRIBUTING.md#publishing-on-conda-forge).
