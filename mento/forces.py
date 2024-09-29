from dataclasses import dataclass, field
from mento.units import kN, m
from devtools import debug
from typing import Dict, Any, Optional

@dataclass
class Forces:
    """
    A class to represent the forces acting on a structural element.

    Attributes
    ----------
    Nx : float
        Axial force along the x-axis (default is 0 kN).
    Vz : float
        Shear force along the z-axis (default is 0 kN).
    My : float
        Bending moment about the y-axis (default is 0 kN*m).

    Methods
    -------
    get_forces() -> dict
        Returns the forces as a dictionary with keys 'Nx', 'Vz', and 'My'.
    set_forces() -> None
        Sets the forces of the object with the provided values.
    """
    _last_id: int = field(default=0, init=False, repr=False)  # Class variable to keep track of last assigned ID
    id: int = field(init=False)  # Instance variable for the ID
    Nx: float = field(default=0*kN)
    Vz: float = field(default=0*kN)
    My: float = field(default=0*kN*m)
    _forces: Dict[str, float] = field(init=False, repr=False)

    def __init__(self, id: Optional[int] = None, Nx: float = 0*kN, Vz: float = 0*kN, My: float = 0*kN*m, **kwargs: Any) -> None:
        # Initialize the default fields
        if id is None:
            Forces._last_id += 1  # Increment the class variable for the next ID
            self.id = Forces._last_id  # Assign the next available ID
        else:
            self.id = id  # Use the specified ID
        self.Nx = Nx
        self.Vz = Vz
        self.My = My
        # Initialize the internal dictionary to hold the forces
        super().__setattr__('_forces', {'Nx': Nx, 'Vz': Vz, 'My': My})
        # Update the dictionary with additional keyword arguments
        for key, value in kwargs.items():
            self._forces[key] = value
    
    def __getattr__(self, name: str) -> float:
        if name in self._forces:
            return self._forces[name]
        else:
            raise AttributeError(f"'Forces' object has no attribute '{name}'")
    
    def __setattr__(self, name: str, value: float) -> None:
        """Override setattr to only allow setting known force attributes."""
        if name in ('id', 'Nx', 'Vz', 'My'):
            super().__setattr__(name, value)
            # Update _forces if it has been initialized
            if '_forces' in self.__dict__:
                self._forces[name] = value
        else:
            raise AttributeError(f"'Forces' object has no attribute '{name}'")


    def get_forces(self) -> Dict[str, float]:
        """
        Get the current force values as a dictionary.

        Returns
        -------
        Dict[str, float]
            A dictionary containing the current values of 'Nx', 'Vz', and 'My'.
        """
        return self._forces

    def set_forces(self, **kwargs: Any) -> None:
        """
        Set the forces to new values.

        Parameters
        ----------
        Nx : float
            New axial force along the x-axis (default is 0*kN).
        Vz : float
            New shear force along the z-axis (default is 0*kN).
        My : float
            New bending moment about the y-axis (default is 0*kN*m).

        Returns
        -------
        None
        """
        for key, value in kwargs.items():
            if key in self._forces:
                setattr(self, key, value)
            else:
                raise AttributeError(f"'Forces' object has no attribute '{key}'")


def main() -> None:
    f1 = Forces(My=10 * kN * m, Nx=2 * kN)
    debug(f1) 
    debug(f1.My, f1.id)
    f = Forces(id=10)
    f.Nx = 20 * kN
    debug(f.Nx)
    debug(f.get_forces())
    f.set_forces(Nx=10*kN,My=15*kN*m)
    debug(f)

if __name__ == "__main__":
    main()
