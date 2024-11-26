from typing import List, Dict, Optional
from pandas import DataFrame
import pandas as pd

from mento.material import Concrete, SteelBar, Concrete_ACI_318_19
from mento.forces import Forces
from mento.beam import RectangularBeam
from mento import mm, cm, kN, MPa, m, inch, ft, kNm
from mento.node import Node

pd.set_option('future.no_silent_downcasting', True)

class BeamSummary:
    def __init__(self, concrete: Concrete, steel_bar: SteelBar, beam_list: DataFrame) -> None:
        self.concrete: Concrete = concrete
        self.steel_bar: SteelBar = steel_bar
        self.beam_list: DataFrame = beam_list
        self.units_row: List[str] = []
        self.data: DataFrame = None
        self.beams: List[RectangularBeam] = []
        self._beam_summary: List = []
        self.check_and_process_input()
        self.convert_to_beams()

    def check_and_process_input(self) -> None:
        # Separate the header, units, and data
        self.units_row = self.beam_list.iloc[0].tolist()  # Second row (units)
        data = self.beam_list.iloc[1:].copy()  # Data rows (after removing the units row)

        # Convert NaN in units to "dimensionless"
        self.units_row = ['' if pd.isna(unit) else unit for unit in self.units_row]

        # Validate the units row
        self.validate_units(self.units_row)
        
        # Convert NaN to 0 in the data rows
        data.iloc[:, 1:] = data.iloc[:, 1:].fillna(0).astype(float)
        # Convert specific columns to int and others to float
        columns_to_int = ['ns', 'n1', 'n2', 'n3', 'n4']
        for col in columns_to_int:
            if col in data.columns:
                data[col] = data[col].astype(int)

        # Apply units to the corresponding columns, skipping the first column
        for i in range(1, len(self.units_row)):  # Start from the second column (index 1)
            unit_str = self.units_row[i]
            if unit_str != '':
                unit = self.get_unit_variable(unit_str)
                if isinstance(data.iloc[:, i], pd.Series):
                    # data.iloc[:, i] = data.iloc[:, i].apply(lambda x: x * unit)
                    data.iloc[:, i] = data.iloc[:, i].apply(lambda x: x * unit)

        # Store the processed data
        self.data = data
        # print("Processed Data: Ok")

    def validate_units(self, units_row: List) -> None:
        valid_units = {"m", "mm", "cm", "inch", "ft", "kN", "kNm", ''}
        for unit_str in units_row:
            if unit_str and unit_str not in valid_units:
                raise ValueError(f"Invalid unit '{unit_str}' detected. Allowed units: {valid_units}")
        # print("Processed Units: Ok")
    
    def get_unit_variable(self, unit_str: str) -> Dict:
        # Map strings to actual unit variables (predefined in the script)
        unit_map = {
            'mm': mm,
            'cm': cm,
            'm': m,
            'in': inch,
            'ft': ft,
            'kN': kN,
            'kNm': kNm,
            'MPa': MPa,
        }
        if unit_str in unit_map:
            return unit_map[unit_str]
        else:
            raise ValueError(f"Unit '{unit_str}' is not recognized.")
        
    def convert_to_beams(self) -> None:
        self.beams = []

        for index, row in self.data.iterrows():
            # Extract forces for each row
            M_y = row['My']  # Example: My in kNm
            N_x = row['Nx']  # Example: Nx in kN
            V_z = row['Vz']  # Example: Vz in kN

            # Ensure these are pint.Quantity objects with correct units
            forces = Forces(M_y=M_y, N_x=N_x, V_z=V_z)

            # Extract geometric properties of the beam (width and height)
            width = row['b']  # Example: width in cm
            height = row['h']  # Example: height in cm

            # Create a rectangular concrete beam using the extracted values
            beam = RectangularBeam(label=row['Label'],
                                              concrete=self.concrete,
                                              steel_bar=self.steel_bar,
                                              width=width,
                                              height=height)        
            # Create a Node for each pair of beam and forces
            Node(section=beam, forces=forces)
            # Set transverse rebar (stirrups) for the beam
            n_stirrups = row['ns']  # Number of stirrups
            d_b = row['dbs']  # Diameter of rebar (mm)
            s_l = row['sl']  # Spacing of stirrups (cm)

            beam.set_transverse_rebar(n_stirrups=n_stirrups, d_b=d_b, s_l=s_l)

            # Store the section and its corresponding forces
            self.beams.append(beam)
    
    def check(self, capacity_check: bool = False) -> DataFrame:
        """
        Perform a check on all beams in the summary.
        
        Parameters
        ----------
        capacity_check : bool, optional
            If True, resets all forces in the node to zero to perform a capacity check. 
            Otherwise, uses the forces currently assigned to the node.
        
        Returns
        -------
        DataFrame
            A DataFrame with the results of the check.
        """
        results_list = []

        for beam in self.beams:
            # Save the original forces
            if beam.node is None:
                raise ValueError(f"Node is not set for the section {beam.label}")
            assert beam.node is not None  # Inform MyPy that node is not None
            node = beam.node # get the node associated with the beam
            original_forces = node.get_forces_list()

            if capacity_check:
                # Reset forces to zero for capacity check
                node.reset_forces()
                # Perform the shear check
                shear_results = beam.check_shear()
                
                if beam._stirrup_n ==0:
                    rebar_v = '-'
                else:
                    rebar_v = f"{int(beam._stirrup_n)}eØ{int(beam._stirrup_d_b.to('mm').magnitude)}/{beam._stirrup_s_l.to('cm').magnitude}"  # noqa: E501
                # Design results
                results_dict = {
                    'Beam': beam.label,  # Minimum shear reinforcement area
                    'b': beam.width.magnitude, # Required shear reinforcing area
                    'h': beam.height.magnitude,  # Provided stirrup reinforcement per unit length
                    'Av': rebar_v,  # Reinforcement contribution to shear capacity
                    'Av,real': round(shear_results['Av'][0].magnitude,2),  #Reinforcement contribution to shear capacity
                    'ØVn': round(shear_results['ØVn'][0].magnitude,2),  #Check if applied shear is within total capacity
                }
                # Create a units row as a DataFrame
                units_row = pd.DataFrame([{
                    'Beam': '',
                    'b': 'cm',
                    'h': 'cm',
                    'Av': '',
                    'Av,real': 'cm²/m',
                    'ØVn': 'kN',
                }])
            else:
                # Perform the shear check
                shear_results = beam.check_shear()
                # Design results
                results_dict = {
                    'Beam': beam.label,  # Minimum shear reinforcement area
                    'b': beam.width.magnitude, # Required shear reinforcing area
                    'h': beam.height.magnitude,  # Provided stirrup reinforcement per unit length
                    'Vu': round(shear_results['Vu'][0].magnitude, 2),  # Max Vu for the design
                    'Nu': round(shear_results['Nu'][0].magnitude, 2),  # Max Nu for the design
                    'Av,req': round(shear_results['Av,req'][0].magnitude, 2), # Reinforcement contribution to shear capacity  # noqa: E501
                    'Av,real': round(shear_results['Av'][0].magnitude, 2),  # Reinforcement contribution to shear capacity # noqa: E501
                    'ØVn': round(shear_results['ØVn'][0].magnitude, 2),  # Check if applied shear is within total capacity # noqa: E501
                    'DCRv': round(shear_results['DCR'][0].magnitude,2),  # Check if applied shear is within total capacity # noqa: E501
                }
                # Create a units row as a DataFrame
                units_row = pd.DataFrame([{
                    'Beam': '',
                    'b': 'cm',
                    'h': 'cm',
                    'Vu': 'kN',
                    'Nu': 'kN',
                    'Av,req': 'cm²/m',
                    'Av,real': 'cm²/m',
                    'ØVn': 'kN',
                    'DCRv': '',
                }])
            # Restore the original forces after capacity check
            node.forces = original_forces
            # Add the results to the list
            results_list.append(results_dict)

        # Convert results list into a DataFrame
        results_df = pd.DataFrame(results_list)
        
        # Combine the units row with the results DataFrame
        final_df = pd.concat([units_row, results_df], ignore_index=True)
        return final_df
    
    def shear_results(self, index: Optional[int] = None, capacity_check: bool = False) -> DataFrame:
        """
        Access the beam by its index from the beam_summary list and retrieve 
        detailed results and beam data. If no index is provided, calculate
        shear results for all beams and return a complete DataFrame.

        :param index: Optional index of the beam in beam_summary.beams list.
        :param capacity_check: If True, performs a capacity check (resets forces to zero).
        :return: A DataFrame of detailed results for the specific beam, 
                or a DataFrame of all beams if no index is provided.
        """
        def process_beam(beam: RectangularBeam) -> DataFrame:
            """
            Process a single beam to calculate shear results.

            :param beam: A RectangularBeam object.
            :return: A DataFrame with shear results for the beam.
            """
            if beam.node is None:
                raise ValueError(f"Node is not set for the section {beam.label}")
            assert beam.node is not None  # Inform MyPy that node is not None
            node = beam.node # get the node associated with the beam
            original_forces = node.get_forces_list()

            # Perform capacity check by resetting forces to zero
            if capacity_check:
                beam.node.reset_forces()

            # Run shear check
            results = beam.check_shear()

            # Restore original forces after capacity check
            if capacity_check:
                beam.node.forces = original_forces

            return results

        # If index is provided, return results for the specific beam
        if index is not None:
            if index - 1 >= len(self.beams):
                raise IndexError(f"Index {index} is out of range for the beam list.")

            # Access the specific beam item
            item = self.beams[max(index - 1, 0)]
            return process_beam(item)  # Return results for the specific beam

        # If no index is provided, calculate results for all beams
        all_shear_results = []

        for item in self.beams:
            results = process_beam(item)
            all_shear_results.append(results)

        # Combine all shear results into a single DataFrame
        return pd.concat(all_shear_results, ignore_index=True)

def main() -> None:
    conc = Concrete_ACI_318_19(name="C25", f_c=25*MPa)
    steel = SteelBar(name="ADN 420", f_y=420*MPa)
    input_df = pd.read_excel(r'.\tests\examples\Mento-Input.xlsx', sheet_name='Beams', usecols='B:R', skiprows=4)
    # # data = {'Label': ['', 'V101', 'V102', 'V103', 'V104'],
    #         'b': ['cm', 20, 20, 20, 20],
    #         'h': ['cm', 50, 50, 50, 50],
    #         'Nx': ['kN', 0, 0, 0, 50],
    #         'Vz': ['kN', 20, 50, 100, 100],
    #         'My': ['kNm', 0, 35, 40, 45],
    #         'ns': ['', 0, 1.0, 1.0, 1.0],
    #         'dbs': ['mm', 0, 6, 6, 6],
    #         'sl': ['cm', 0, 20, 20, 20],
    #         'n1': ['', 2.0, 2.0, 2.0, 2.0],
    #         'db1': ['mm', 12, 12, 12, 12],
    #         'n2': ['', 1.0, 0.0, 1.0, 0.0],
    #         'db2': ['mm', 10, 0, 10, 0],
    #         'n3': ['', 2.0, 0.0, 2.0, 0.0],
    #         'db3': ['mm', 12, 0, 0, 0],
    #         'n4': ['', 1.0, 0.0, 1.0, 0.0],
    #         'db4': ['mm', 10, 0, 0, 0]}
    # input_df = pd.DataFrame(data)
    # print(input_df)
    beam_summary = BeamSummary(concrete=conc, steel_bar=steel, beam_list=input_df)
    # print(beam_summary.data)
    capacity = beam_summary.check(capacity_check=True)
    print(capacity)
    check = beam_summary.check()
    print(check)
    # beam_summary.check().to_excel('hola.xlsx', index=False)
    results = beam_summary.shear_results(capacity_check=False)
    print(results)
    # beam_summary.beams[0].shear_results_detailed()

if __name__ == "__main__":
    main()
