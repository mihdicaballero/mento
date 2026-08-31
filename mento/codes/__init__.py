"""Design code implementations.

Deliberately empty. This module used to re-export the private check and design
functions of every code, which meant that touching any submodule — including a
leaf like ``codes.aci_318_19.equations.shear``, which imports nothing from
mento — pulled the entire code layer in with it. Nothing imported those
re-exports; what it did buy was an import cycle, since the code modules import
``mento.rebar`` and ``mento.rebar`` now reaches the equation modules.

Import the modules directly (``from mento.codes.ACI_318_19_beam import ...``).
"""
