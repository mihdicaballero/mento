from dataclasses import dataclass, field
from mento.units import kN, kNm
from pint import Quantity
from devtools import debug
from typing import Dict, Optional

@dataclass
class Forces:
    """
    A class to represent the forces acting on a structural element.

    Attributes
    ----------
    N_x : float
        Axial force along the x-axis (default is 0 kN).
    V_z : float
        Shear force along the z-axis (default is 0 kN).
    M_y : float
        Bending moment about the y-axis (default is 0 kN*m).

    Methods
    -------
    get_forces() -> dict
        Returns the forces as a dictionary with keys 'N_x', 'V_z', and 'M_y'.
    set_forces() -> None
        Sets the forces of the object with the provided values.
    """
    _id: int = field(init=False, repr=False)  # Instance ID, assigned internally
    _last_id: int = field(default=0, init=False, repr=False)  # Class variable to keep track of last assigned ID
    label: Optional[str] = None

    def __init__(self, label: Optional[str] = None, N_x: Quantity = 0 * kN, V_z: Quantity = 0 * kN,
                  M_y: Quantity = 0 * kNm) -> None:
 
        # Increment the class variable for the next unique ID
        Forces._last_id += 1
        self._id = Forces._last_id  # Private ID assigned internally, unique per instance

        # Initialize label (user-defined or None by default)
        self.label = label

        # Set the forces
        self._N_x = N_x
        self._V_z = V_z
        self._M_y = M_y

    @property
    def id(self) -> int:
        """Read-only property for accessing the unique ID of the instance."""
        return self._id

    @property
    def N_x(self) -> Quantity:
        return self._N_x

    @property
    def V_z(self) -> Quantity:
        return self._V_z

    @property
    def M_y(self) -> Quantity:
        return self._M_y

    def get_forces(self) -> Dict[str, Quantity]:
        """Returns the forces as a dictionary."""
        return  {
            'N_x': self._N_x,
            'V_z': self._V_z,
            'M_y': self._M_y
        }

    def set_forces(self, N_x: Quantity = 0*kN, V_z: Quantity = 0*kN, M_y: Quantity = 0*kNm) -> None:
        """Sets the forces in the object."""
        self._N_x = N_x
        self._V_z = V_z
        self._M_y = M_y


def main() -> None:
    f = Forces(M_y=10 * kNm, N_x=2 * kN)
    debug(f) 
    print(f.M_y, f.id)
    debug(f.N_x)
    debug(f.get_forces())
    f.set_forces(N_x=10*kN,M_y=15*kNm)
    f.label = "Crane load"
    debug(f)
    debug(f.get_forces())

if __name__ == "__main__":
    main()
