from devtools import debug
from typing import List
from mento.beam import RectangularBeam
from mento.material import Concrete, SteelBar, Concrete_ACI_318_19, Concrete_EN_1992_2004, Concrete_CIRSOC_201_25
from mento.section import Section
from mento.forces import Forces
from typing import Union
import pandas as pd
from mento.units import MPa, ksi, psi, kip, mm, inch, kN, m, cm, kNm, ft, dimensionless

class Node:
    def __init__(self, section: Section, forces: Union[Forces, List[Forces]]) -> None:
        if not isinstance(section, Section):
            raise TypeError("section must be an instance of Section")
        
        # Ensure forces is either a single Forces object or a list of Forces
        if isinstance(forces, Forces):
            self.forces = [forces]  # Wrap single Forces object in a list
        elif isinstance(forces, list) and all(isinstance(force, Forces) for force in forces):
            self.forces = forces  # Assign the list directly if valid
        else:
            raise TypeError("forces must be an instance of Forces or a list of Forces")
        
        self.section = section
        # Set the node reference in the section (beam)
        self.section.node = self  # Associate the Node with the RectangularBeam

    def add_forces(self, forces: Union[Forces, List[Forces]]) -> None:
        if isinstance(forces, Forces):
            self.forces.append(forces)  # Append single Forces object
        elif isinstance(forces, list) and all(isinstance(force, Forces) for force in forces):
            self.forces.extend(forces)  # Extend list with multiple Forces objects
        else:
            raise TypeError("forces must be an instance of Forces or a list of Forces")

    def get_forces_list(self) -> List[Forces]:
        """Returns the list of forces applied to this node."""
        return self.forces

    def reset_forces(self) -> None:
        """Reset all forces in the node to zero."""   
        if isinstance(self.forces, Forces):
            # Replace the single Forces object with a zero-initialized Forces object
            self.forces = [Forces()]
        elif isinstance(self.forces, list):
            # Replace each force in the list with a zero-initialized Forces object
            self.forces = [Forces() for _ in self.forces]


    def check_flexure(self) -> pd.DataFrame:
        """
        """        
        try:
            return self.section._check_flexure(self.forces)
        except AttributeError:
            print("La sección actual no tiene el método 'check_flexure' definido")


    def design_flexure(self) -> pd.DataFrame:
        """
        """        
        try:
            return self.section._design_flexure(self.forces)
        except AttributeError:
            print("La sección actual no tiene el método 'design_flexure' definido")



    def check_shear(self) -> pd.DataFrame:
        """
        """        
        try:
            return self.section._check_shear(self.forces)
        except AttributeError:
            print("La sección actual no tiene el método 'design_flexure' definido")


    def design_shear(self) -> pd.DataFrame:
        """
        """        
        try:
            return self.section._design_shear(self.forces)
        except AttributeError:
            print("La sección actual no tiene el método 'design_flexure' definido")




##########################################################################
# DEBUG ZONE
##########################################################################

def flexure_design_test() -> None:
    concrete = Concrete_ACI_318_19(name="fc 4000", f_c=4000*psi)  
    steelBar = SteelBar(name="fy 60000", f_y=60*ksi)  
    custom_settings = {'clear_cover': 1.5*inch} 
    beam = RectangularBeam(
        label="B-12x24",
        concrete=concrete,
        steel_bar=steelBar,
        width=12*inch,  
        height=24*inch,
        settings=custom_settings  
    )
    
    f = Forces(label='Test_01', V_z = 40*kip, M_y=400*kip*ft)
    f2 = Forces(label='Test_01', V_z = 100*kip, M_y=-400*kip*ft)
    N1=Node(section=beam, forces=[f,f2])

    flexure_results = N1.design_flexure()
    print(flexure_results)

def shear_ACI_imperial() -> None:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000*psi)  
    steelBar = SteelBar(name="ADN 420", f_y=60*ksi)  
    custom_settings = {'clear_cover': 1.5*inch} 
    beam = RectangularBeam(
        label="102",
        concrete=concrete,
        steel_bar=steelBar,
        width=10*inch,  
        height=16*inch,
        settings=custom_settings  
    )

    # f1 = Forces(label='D', V_z=37.727*kip, N_x=20*kip)
    f1 = Forces(label='D', V_z=37.727*kip)
    # f1 = Forces(label='D', V_z=8*kip)
    # f2 = Forces(label='L', V_z=6*kip) # No shear reinforcing
    forces=[f1]
    beam.set_transverse_rebar(n_stirrups=1, d_b=0.5*inch, s_l=6*inch)
    N1=Node(beam,forces)
    # beam.set_longitudinal_rebar_bot(n1=2, d_b1=0.625*inch)
    print(beam._A_v)
    results = N1.check_shear()
    # results = beam.design_shear()
    print(results)
    # section.design_shear(f, A_s=0.847*inch**2)
    N1.section.shear_results_detailed()  
    # section.shear_results_detailed_doc()

if __name__ == "__main__":
    flexure_design_test()
    shear_ACI_imperial()