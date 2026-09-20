---
name: prepare-release
description: Prepare a mento release branch — bump the version in pyproject.toml and CITATION.cff, move the CHANGELOG Unreleased entries under the new version, and verify the suite. Use when asked to cut, prepare or bump a release. Publishing itself happens on GitHub, not here.
---

# Preparing a mento release

Merging a pull request publishes nothing. PyPI is triggered only by publishing a GitHub
Release, which starts `.github/workflows/publish.yml`. This skill covers the branch that
precedes it. The full procedure, including Zenodo and conda-forge, is in CONTRIBUTING.md
under *Making a release*.

## Steps

1. **Decide the number.** mento is past 1.0, so strict SemVer applies
   (`docs/architecture/adr/0003-pre-1.0-versioning-policy.md`): a breaking change is a
   major, a backwards-compatible addition a minor, a fix a patch. Every break gets a
   migration note in the CHANGELOG. Ask the user if the entries do not make the level
   obvious.

2. **Bump three files, and only these three:**
   - `pyproject.toml` → `version`
   - `CITATION.cff` → `version` and `date-released` (today, `YYYY-MM-DD`)
   - `CHANGELOG.md` → move every entry under `Unreleased` into a new `## [X.Y.Z] - YYYY-MM-DD`
     section, leaving `Unreleased` in place but empty

   `mento/_version.py` reads the installed distribution metadata, so it needs no edit.

3. **Verify before committing.** Run the suite, ruff and mypy with the platform interpreter
   from CLAUDE.md (`ruff check .`, `ruff format --check .`, `mypy mento/`, `pytest tests/`). A release branch that fails CI wastes a tag.

4. **Commit and push the branch.** Do not create the tag or the GitHub Release — the tag
   must point at the merge commit on `main`, and `publish.yml` refuses to continue if the
   tag and `pyproject.toml` disagree.

5. **Hand back** the version, the CHANGELOG summary, and a reminder of what the user still
   does by hand: merge the branch, then draft the release with tag `vX.Y.Z` (the leading
   `v` matters) against `main`.

## Do not

- Do not publish to PyPI from a local machine. Trusted publishing expects the workflow.
- Do not invent CHANGELOG entries. If `Unreleased` is empty, say so and stop.
