<img width="800" alt="mento_github" src="https://github.com/user-attachments/assets/32128ec0-a8c2-4782-b210-1afc633f0da6" />

*An intuitive tool for structural engineers to design concrete elements efficiently.*

[![Tests](https://github.com/mihdicaballero/mento/actions/workflows/tests.yml/badge.svg)][tests]
[![Docs](https://readthedocs.org/projects/mento-docs/badge/?version=latest)](https://mento-docs.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/github/mihdicaballero/mento/graph/badge.svg?token=9X81ZRKMCX)](https://codecov.io/github/mihdicaballero/mento)
[![PyPI](https://img.shields.io/pypi/v/mento.svg)](https://pypi.org/project/mento/)
[![Python versions](https://img.shields.io/pypi/pyversions/mento.svg)](https://pypi.org/project/mento/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)][ruff]

[tests]: https://github.com/mihdicaballero/mento/actions/workflows/tests.yml
[ruff]: https://github.com/charliermarsh/ruff

This repository provides a comprehensive package for the design and structural analysis of concrete sections, beams and columns. The package adheres to multiple design codes, ensuring broad applicability in structural engineering projects.

#### Installation

```bash
pip install mento
```

Requires Python 3.10 or newer.

#### Features
- Check and design for flexure and shear of:
    - Rectangular Concrete Beam
    - Rectangular One way Slab
- Unit-sensitive design, allowing users to input variables with their respective units for accurate calculations.
- Interactive usage in Jupyter Notebooks, allowing users to create custom calculations using package modules.
- Results are presented in markdown format within the notebook and as pandas DataFrames for easy handling of multiple checks.
- Ability to generate detailed calculation reports in Word.
- Comprehensive testing for design code compliance, including ACI 318-19, EN 1992-2004, and CIRSOC 201-2025.

#### Where mento fits

mento designs and verifies **structural members against a design code**, and documents the
result: given a section, its materials and a set of load combinations, it returns the
required reinforcement, the demand-to-capacity ratios, and a calculation report you can hand
to a reviewer.

mento is also, as far as we know, the only open source package that implements
**CIRSOC 201-2025**, the Argentinian concrete design standard.

#### Using mento at your company?

mento is free and open-source under MIT. If your team wants help going further, we offer a few things on top:

- **Team onboarding** — short remote sessions to get your engineers productive with Python, Jupyter, VS Code and Git, using mento as the starting point.
- **Custom tools** — building specific workflows on top of mento (custom reports, integrations with Revit / ETABS / Robot, batch design tools, internal company libraries).
- **Sponsored features** — if your firm needs a specific element type, design code, or check that's on the roadmap (or not), we can scope it as paid work and contribute it back.

If any of that is useful for your team, fill out [this form](https://forms.gle/QoDzczQToLa78jMo7) and we'll follow up.

Not looking for anything? A ⭐ on the repo or feedback in [Discussions](https://github.com/mihdicaballero/mento/discussions) is also genuinely appreciated.

#### Roadmap
The development is structured around key milestones, with ongoing tasks that aim to enhance functionality and compliance with design standards:
- [x] Rectangular concrete beam section check and design for ACI 318-19 and CIRSOC 201-25.
- [x] Rectangular concrete beam section check and design for EN 1992-2004.
- [x] One way concrete slab check and design for ACI 318-19 and CIRSOC 201-25.
- [x] One way concrete slab check and design for EN 1992-2004.
- [x] Shear wall shear check and design for ACI 318-19 and CIRSOC 201-25.
- [ ] Slab shear punching check and design for ACI 318-19 and CIRSOC 201-25. (in progress)
- [ ] Slab shear punching check and design for EN 1992-2004.
- [ ] Shear wall shear check and design for EN 1992-2004.

Each milestone incorporates rigorous testing and continuous integration to ensure code quality and reliability.

#### Documentation
You can read the official documentation in this link: [Mento Docs](https://mento-docs.readthedocs.io/)

#### Contributing
We welcome contributions from the community to expand and enhance the package. Start with the [contributing guide](CONTRIBUTING.md), then look through the [open issues](https://github.com/mihdicaballero/mento/issues) — those labelled `good first issue` are a good entry point. Every calculation contributed to mento is validated against a published worked example before it is merged; the guide explains how.

Participation in this project is governed by our [Code of Conduct](CODE_OF_CONDUCT.md).

#### Disclaimer
mento is a tool to assist structural engineers, not a replacement for engineering judgement. Results must be reviewed and accepted by a qualified engineer who takes responsibility for the design. The software is provided "as is", without warranty of any kind, and the authors accept no liability for its use. Verify the output against the applicable design code before relying on it for any real structure.

#### Citing mento
If mento supports your research or professional work, citation metadata is in [CITATION.cff](CITATION.cff), and GitHub's "Cite this repository" button will format it for you. The [citing guide](https://mento-docs.readthedocs.io/en/latest/getting_started/citing.html) explains which version to cite and gives a BibTeX entry.

#### License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for the full text.
