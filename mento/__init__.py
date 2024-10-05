from pint import UnitRegistry

# Initialize the UnitRegistry
ureg = UnitRegistry(system='mks')
Quantity = ureg.Quantity
ureg.formatter.default_format = '.2f~P'