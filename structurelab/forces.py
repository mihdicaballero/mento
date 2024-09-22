from dataclasses import dataclass, field
from structurelab.units import kN, m

@dataclass
class Forces:
    Nx: float = field(default=0*kN) #type: ignore
    Vz: float = field(default=0*kN) #type: ignore
    My: float = field(default=0*kN*m) #type: ignore

    def __init__(self, Nx: float = 0*kN, Vz: float = 0*kN, My: float = 0*kN*m, **kwargs): #type: ignore
        # Initialize the default fields
        super().__setattr__('Nx', Nx)
        super().__setattr__('Vz', Vz)
        super().__setattr__('My', My)
        # Initialize the internal dictionary to hold the forces
        super().__setattr__('_forces', {'Nx': Nx, 'Vz': Vz, 'My': My})
        # Update the dictionary with additional keyword arguments
        for key, value in kwargs.items():
            self._forces[key] = value

    def __getattr__(self, name):
        if name in self._forces:
            return self._forces[name]
        else:
            raise AttributeError(f"'Forces' object has no attribute '{name}'")

    def __setattr__(self, name, value):
        if '_forces' in self.__dict__:
            if name in self._forces:
                self._forces[name] = value
                super().__setattr__(name, value)
            else:
                raise AttributeError(f"'Forces' object has no attribute '{name}'")
        else:
            super().__setattr__(name, value)

def main():
    f1 = Forces(My=10*kN*m,Nx=2*kN) # type: ignore
    print(f1.My,f1.Nx) 
    f = Forces()
    f.Nx = 20*kN #type: ignore
    print(f)

if __name__ == "__main__":
    main()