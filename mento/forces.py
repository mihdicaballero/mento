from dataclasses import dataclass, field
from typing import Dict, Optional
from pint import Quantity
from pint.facets.plain import PlainQuantity

from mento.units import kN, kNm


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
    unit_system : str
        The unit system to use for displaying forces ('metric' or 'imperial').

    Methods
    -------
    get_forces() -> dict
        Returns the forces as a dictionary with keys 'N_x', 'V_z', and 'M_y'.
    set_forces() -> None
        Sets the forces of the object with the provided values.
    compare_to(other: 'Forces', by: str = 'V_z') -> bool
        Compares this force with another force based on a given attribute.
    """

    _id: int = field(init=False, repr=False)  # Instance ID, assigned internally
    _last_id: int = field(
        default=0, init=False, repr=False
    )  # Class variable to keep track of last assigned ID
    label: Optional[str] = None
    _N_x: Quantity = field(default=0 * kN)  # Ensure you use the correct unit type here
    _V_z: Quantity = field(default=0 * kN)
    _M_y: Quantity = field(default=0 * kNm)
    unit_system: str = field(default="metric")  # Add unit_system as a field

    def __init__(
        self,
        label: Optional[str] = None,
        N_x: Quantity = 0 * kN,
        V_z: Quantity = 0 * kN,
        M_y: Quantity = 0 * kNm,
        unit_system: str = "metric",
    ) -> None:
        # Increment the class variable for the next unique ID
        Forces._last_id += 1
        self._id = (
            Forces._last_id
        )  # Private ID assigned internally, unique per instance

        # Initialize the label
        self.label = label
        self.unit_system = unit_system  # Set the unit system

        # Set the forces upon initialization
        self.set_forces(N_x, V_z, M_y)

    @property
    def id(self) -> int:
        """Read-only property for accessing the unique ID of the instance."""
        return self._id
    
    @id.setter
    def id(self, value) -> None:  # type: ignore[no-untyped-def]
        # Normalize error message across Python versions (3.10 vs 3.11+)
        raise AttributeError("property 'id' of 'Forces' object has no setter")

    @property
    def N_x(self) -> Quantity:
        """Axial force along the x-axis (default is 0 kN)."""
        if self.unit_system == "metric":
            return self._N_x.to("kN")
        else:
            return self._N_x.to("kip")

    @property
    def V_z(self) -> Quantity:
        """Shear force along the z-axis (default is 0 kN)."""
        if self.unit_system == "metric":
            return self._V_z.to("kN")
        else:
            return self._V_z.to("kip")

    @property
    def M_y(self) -> Quantity:
        """Bending moment about the y-axis (default is 0 kN*m)."""
        if self.unit_system == "metric":
            return self._M_y.to("kN*m")
        else:
            return self._M_y.to("ft*kip")

    def get_forces(self) -> Dict[str, Quantity]:
        """Returns the forces as a dictionary with keys 'N_x', 'V_z', and 'M_y'."""
        if self.unit_system == "metric":
            return {
                "N_x": self._N_x.to("kN"),
                "V_z": self._V_z.to("kN"),
                "M_y": self._M_y.to("kN*m"),
            }
        else:
            return {
                "N_x": self._N_x.to("kip"),
                "V_z": self._V_z.to("kip"),
                "M_y": self._M_y.to("ft*kip"),
            }

    def set_forces(
        self, N_x: Quantity = 0 * kN, V_z: Quantity = 0 * kN, M_y: Quantity = 0 * kNm
    ) -> None:
        """Sets the forces in the object."""
        self._N_x = N_x
        self._V_z = V_z
        self._M_y = M_y

    def compare_to(self, other: "Forces", by: str = "V_z") -> bool:
        """Compares this force with another force based on a selected attribute.

        Parameters
        ----------
        other : Forces
            Another Forces instance to compare with.
        by : str
            The attribute to compare by ('N_x', 'V_z', or 'M_y').

        Returns
        -------
        bool
            True if this force is greater than the other force by the selected attribute.
        """
        if by not in ["N_x", "V_z", "M_y"]:
            raise ValueError(
                "Comparison attribute must be one of 'N_x', 'V_z', or 'M_y'"
            )
        return getattr(self, by).magnitude > getattr(other, by).magnitude

    def __str__(self) -> str:
        return (
            f"Force ID: {self.id}, Label: {self.label}, "
            f"N_x: {self.N_x}, V_z: {self.V_z}, M_y: {self.M_y}"
        )
