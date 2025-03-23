from typing import Union, Optional, List
import pandas as pd
from dataclasses import field

from mento.beam import RectangularBeam
from mento.material import SteelBar, Concrete_ACI_318_19, Concrete_EN_1992_2004
from mento.section import Section
from mento.forces import Forces
from mento.units import ksi, psi, kip, inch, ft, MPa, cm, mm, kN

class Node:
    _id: int
    _last_id: int = 0  # Class variable to keep track of last assigned ID
    
    def __init__(self, section: Section, forces: Union[Forces, List[Forces]]) -> None:
        # Initialize the private ID automatically
        Node._last_id += 1  # Increment the class variable for the next ID
        self._id = Node._last_id  # Assign the next available ID 

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

    def clear_forces(self) -> None:
        """Remove all forces from the node."""
        self.forces.clear()

    def check_flexure(self) -> pd.DataFrame:
        return self.section.check_flexure(self.forces)

    def design_flexure(self) -> pd.DataFrame:
        return self.section.design_flexure(self.forces)

    def check_shear(self) -> pd.DataFrame:
        return self.section.check_shear(self.forces)

    def design_shear(self) -> pd.DataFrame:
        return self.section.design_shear(self.forces)
    
    def shear_results_detailed(self, force: Optional[Forces] = None) -> None:
        return self.section.shear_results_detailed(force)
    
    def shear_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        return self.section.shear_results_detailed_doc(force)

    def flexure_results_detailed(self, force: Optional[Forces] = None) -> None:
        return self.section.flexure_results_detailed(force)

    def flexure_results_detailed_doc(self, force: Optional[Forces] = None) -> None:
        return self.section.flexure_results_detailed_doc(force)
    
    def __repr__(self) -> str:
        # Get the section label (assuming Section has a `_id` attribute)
        section_label = f"Section label: {self.section.label}"
        
        # Get the list of forces applied
        forces_list = [str(force) for force in self.forces]  # Assuming Forces has a `__repr__` or `__str__` method
        forces_str = "\n  - " + "\n  - ".join(forces_list) if forces_list else "No forces applied"
        
        # Combine section label and forces into a single string
        return f"Node ID: {self.id} - {section_label}\nForces Applied:{forces_str}"

    # Beam results for Jupyter Notebook
    @property
    def results(self) -> None:
        return self.section.results
    @property
    def id(self) -> int:
        """Read-only property to access the private _id."""
        return self._id

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
    node_1=Node(section=beam, forces=[f,f2])

    flexure_results = node_1.design_flexure()
    node_1.flexure_results_detailed()
    # node_1.flexure_results_detailed_doc()
    print(node_1)

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
    node_1=Node(beam,forces)
    # beam.set_longitudinal_rebar_bot(n1=2, d_b1=0.625*inch)
    print(beam._A_v)
    results = node_1.check_shear()
    # results = beam.design_shear()
    print(results)
    # section.design_shear(f, A_s=0.847*inch**2)
    node_1.shear_results_detailed()  
    node_1.shear_results_detailed_doc()

def shear_ACI_metric() -> None:
    concrete= Concrete_ACI_318_19(name="C30",f_c=30*MPa) 
    steelBar= SteelBar(name="ADN 420", f_y=420*MPa)
    custom_settings = {'clear_cover': 30*mm}
    beam = RectangularBeam(label="101",
                                      concrete=concrete,steel_bar=steelBar,width=20*cm, height=50*cm,
                                       settings=custom_settings)
    f1 = Forces(label='1.4D', V_z=100*kN)
    f2 = Forces(label='1.2D+1.6L', V_z=1.55*kN)
    f3 = Forces(label='W', V_z=2.20*kN)
    f4 = Forces(label='S', V_z=8.0*kN)
    f5 = Forces(label='E', V_z=1.0*kN)
    forces=[f1, f2, f3, f4, f5]
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)
    beam.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=20*cm) 
    node = Node(section=beam, forces=forces)
    results = node.check_shear()
    # results = beam.design_shear(forces)
    print(node)
    # print(beam.shear_design_results)
    # beam.shear_results_detailed()
    # print(beam.shear_design_results)
    # print(beam.results)
    # beam.shear_results_detailed_doc()

def shear_EN_1992() -> None:
    concrete= Concrete_EN_1992_2004(name="C25",f_ck=25*MPa) 
    steelBar= SteelBar(name="B500S", f_y=500*MPa)
    custom_settings = {'clear_cover': 2.6*cm}
    beam = RectangularBeam(label="101",
                                      concrete=concrete,steel_bar=steelBar,width=20*cm, height=60*cm,
                                       settings=custom_settings)
    # f = Forces(V_z=100*kN, M_y=100*kNm)
    f = Forces(V_z=30*kN, N_x=0*kN)
    forces=[f]
    beam.set_transverse_rebar(n_stirrups=1, d_b=6*mm, s_l=25*cm)
    beam.set_longitudinal_rebar_bot(n1=4,d_b1=16*mm)
    # results = beam.check_shear(forces)
    results = beam.design_shear(forces)
    print(results)
    beam.shear_results_detailed()
    # beam.shear_results_detailed_doc()
if __name__ == "__main__":
    # flexure_design_test()
    # shear_ACI_imperial()
    shear_ACI_metric()
    # shear_EN_1992()