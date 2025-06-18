# Mento
*An intuitive tool for structural engineers to design concrete elements efficiently.*

[![Tests](https://github.com/mihdicaballero/mento/actions/workflows/tests.yml/badge.svg)][tests]
[![Docs](https://readthedocs.org/projects/mento-docs/badge/?version=latest)](https://mento-docs.readthedocs.io/en/latest/?badge=latest)
[![Codecov](https://codecov.io/gh/robbievanleeuwen/concrete-properties/branch/master/graph/badge.svg)][codecov]
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)][ruff]

[tests]: https://github.com/mihdicaballero/mento/actions/workflows/tests.yml
[ruff]: https://github.com/charliermarsh/ruff
[codecov]: https://app.codecov.io/github/mihdicaballero/mento

This repository provides a comprehensive package for the design and structural analysis of concrete sections, beams and columns. The package adheres to multiple design codes, ensuring broad applicability in structural engineering projects.

#### Features
- Check and design for flexure and shear of:
    - Rectangular Concrete Beam
    - Circular Concrete Column
    - Rectangular Concrete Column
- Unit-sensitive design, allowing users to input variables with their respective units for accurate calculations.
- Interactive usage in Jupyter Notebooks, allowing users to create custom calculations using package modules.
- Results are presented in markdown format within the notebook and as pandas DataFrames for easy handling of multiple checks.
- Ability to generate detailed calculation reports in Word.
- Comprehensive testing for design code compliance, including ACI 318-19, EN 1992-2004, and CIRSOC 201-2025.

#### Roadmap
The development is structured around key milestones, with ongoing tasks that aim to enhance functionality and compliance with design standards:
- [x] Rectangular concrete beam section check and design for ACI 318-19 and CIRSOC 201-25.
- [ ] Rectangular concrete beam section check and design for EN 1992-2004.
- [ ] Circular concrete column section check and design for ACI 318-19 and CIRSOC 201-25.
- [ ] Circular concrete column section check and design for EN 1992-2004.
- [ ] Rectangular concrete column section check and design for ACI 318-19 and CIRSOC 201-25.
- [ ] Rectangular concrete column section check and design for EN 1992-2004.

Each milestone incorporates rigorous testing and continuous integration to ensure code quality and reliability.

#### Documentation
You can read the official documentation in this link: [Mento Docs](https://mento-docs.readthedocs.io/)

#### Contributing
We welcome contributions from the community to expand and enhance the package. Please check the roadmap for current milestones and open issues for collaboration opportunities.

#### License
This project is licensed under the MIT License.
