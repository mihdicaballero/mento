# `_wip/` — temporary, delete before merging

Working material shared on a branch so it can be read in review. **Nothing in this
directory is part of the package**, and the directory is expected to be deleted in the
last commit of the branch, before the PR is merged. `pyproject.toml` excludes it from
ruff; that exclusion goes with it.

## Contents

- **`punching_phase_1.ipynb`** — the punching module after Phase 1 of
  [`docs/architecture/punching-roadmap.md`](../docs/architecture/punching-roadmap.md):
  declaring the slab's top reinforcement and reading back the derived `d` and ρ, which
  force components a punching node reads, and the four column cases. `check()` and
  `design()` are still `NotImplementedError`.

  Outputs are committed on purpose — the point is that it reads in a browser or a diff
  without anyone having to run a kernel.

- **`build_punching_notebook.py`** — regenerates the notebook from source. The notebook is
  written by this script rather than by hand so the prose and the code stay in one
  reviewable file:

  ```
  python _wip/build_punching_notebook.py
  jupyter nbconvert --to notebook --execute --inplace _wip/punching_phase_1.ipynb
  ```

The permanent home for a worked example is `docs/source/examples/`, once the check exists
and the notebook can show a result rather than a placeholder.
