from structurelab.common.materials import *

import math
import forallpeople
forallpeople.environment('structural',top_level=True)
#Useful extra units
cm = 1e-2*m


concrete_aci = Concrete("H-25","ACI",25*MPa)
properties_aci = concrete_aci.get_concrete_properties()
print(properties_aci)