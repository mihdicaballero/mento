from typing import List, Dict, Optional
from pandas import DataFrame
import pandas as pd
import copy
from collections import OrderedDict

from mento.material import (
    Concrete,
    SteelBar,
    Concrete_EN_1992_2004,
)
from mento.forces import Forces
from mento.beam import RectangularBeam
from mento import mm, cm, kN, MPa, m, inch, ft, kNm
from mento.node import Node

pd.set_option("future.no_silent_downcasting", True)


class BeamSummary:
    def __init__(
        self, concrete: Concrete, steel_bar: SteelBar, beam_list: DataFrame
    ) -> None:
        self.concrete: Concrete = concrete
        self.steel_bar: SteelBar = steel_bar
        self.beam_list: DataFrame = beam_list
        self.units_row: List[str] = []
        self.data: DataFrame = None
        self.nodes: List[Node] = []
        self._beam_summary: List = []
        self.check_and_process_input()
        self.convert_to_nodes()

    def check_and_process_input(self) -> None:
        # Separate the header, units, and data
        self.units_row = self.beam_list.iloc[0].tolist()  # Second row (units)
        data = self.beam_list.iloc[
            1:
        ].copy()  # Data rows (after removing the units row)

        # Convert NaN in units to "dimensionless"
        self.units_row = ["" if pd.isna(unit) else unit for unit in self.units_row]

        # Validate the units row
        self.validate_units(self.units_row)

        # Convert NaN to 0 in the data rows
        data.iloc[:, 2:] = data.iloc[:, 2:].fillna(0).astype(float)
        # Convert specific columns to int and others to float
        columns_to_int = ["ns", "n1", "n2", "n3", "n4"]
        for col in columns_to_int:
            if col in data.columns:
                data[col] = data[col].astype(int)

        # Apply units to the corresponding columns, skipping the first column
        for i in range(
            1, len(self.units_row)
        ):  # Start from the second column (index 1)
            unit_str = self.units_row[i]
            if unit_str != "":
                unit = self.get_unit_variable(unit_str)
                if isinstance(data.iloc[:, i], pd.Series):
                    # data.iloc[:, i] = data.iloc[:, i].apply(lambda x: x * unit)
                    data.iloc[:, i] = data.iloc[:, i].apply(lambda x: x * unit)

        # Store the processed data
        self.data = data
        # print("Processed Data: Ok")

    def validate_units(self, units_row: List) -> None:
        valid_units = {"m", "mm", "cm", "inch", "ft", "kN", "kNm", ""}
        for unit_str in units_row:
            if unit_str and unit_str not in valid_units:
                raise ValueError(
                    f"Invalid unit '{unit_str}' detected. Allowed units: {valid_units}"
                )
        # print("Processed Units: Ok")

    def get_unit_variable(self, unit_str: str) -> Dict:
        # Map strings to actual unit variables (predefined in the script)
        unit_map = {
            "mm": mm,
            "cm": cm,
            "m": m,
            "in": inch,
            "ft": ft,
            "kN": kN,
            "kNm": kNm,
            "MPa": MPa,
        }
        if unit_str in unit_map:
            return unit_map[unit_str]
        else:
            raise ValueError(f"Unit '{unit_str}' is not recognized.")

    def convert_to_nodes(self) -> None:
        self.nodes = []

        for index, row in self.data.iterrows():
            # Extract forces for each row
            M_y = row["My"]
            N_x = row["Nx"]  # Positive for compression
            V_z = row["Vz"]
            comb = row["Comb."]

            # Ensure these are pint.Quantity objects with correct units
            forces = Forces(label=comb, M_y=M_y, N_x=N_x, V_z=V_z)

            # Extract geometric properties of the beam (width and height)
            width = row["b"]
            height = row["h"]
            c_c = row["cc"]

            # Create a rectangular concrete beam using the extracted values
            beam = RectangularBeam(
                label=row["Label"],
                concrete=self.concrete,
                steel_bar=self.steel_bar,
                width=width,
                height=height,
                c_c=c_c,
            )
            # Set transverse rebar (stirrups) for the beam
            n_stirrups = row["ns"]  # Number of stirrups
            d_b = row["dbs"]  # Diameter of rebar (mm)
            s_l = row["sl"]  # Spacing of stirrups (cm)

            if n_stirrups != 0:
                beam.set_transverse_rebar(n_stirrups=n_stirrups, d_b=d_b, s_l=s_l)

            # Set longitudinal rebar at the bottom if n1 is not 0
            n1 = row["n1"]
            d_b1 = row["db1"]  # Diameter in mm
            n2 = row["n2"]
            d_b2 = row["db2"]  # Diameter in mm
            n3 = row["n3"]
            d_b3 = row["db3"]  # Diameter in mm
            n4 = row["n4"]
            d_b4 = row["db4"]  # Diameter in mm
            if n1 != 0:
                if M_y >= 0 * kNm:
                    beam.set_longitudinal_rebar_bot(
                        n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4
                    )
                else:
                    beam.set_longitudinal_rebar_top(
                        n1, d_b1, n2, d_b2, n3, d_b3, n4, d_b4
                    )
            # Create a Node for each pair of beam and forces
            node = Node(section=beam, forces=forces)

            # Store the section and its corresponding forces
            self.nodes.append(node)

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
        for node in self.nodes:
            beam: RectangularBeam = node.section  # type: ignore
            original_forces = [copy.deepcopy(force) for force in node.get_forces_list()]

            rebar_v = (
                "-"
                if beam._stirrup_n == 0
                else f"{int(beam._stirrup_n)}eØ{int(beam._stirrup_d_b.to('mm').magnitude)}/{int(beam._stirrup_s_l.to('cm').magnitude)}"
            )  # noqa: E501
            rebar_f_top = (
                "-"
                if beam._n1_t == 0
                else (
                    f"{beam._format_longitudinal_rebar_string(beam._n1_t, beam._d_b1_t, beam._n2_t, beam._d_b2_t)}"
                    + (
                        f" ++ {beam._format_longitudinal_rebar_string(beam._n3_t, beam._d_b3_t, beam._n4_t, beam._d_b4_t)}"
                        if beam._n3_t != 0
                        else ""
                    )
                )
            )
            rebar_f_bot = (
                "-"
                if beam._n1_b == 0
                else (
                    f"{beam._format_longitudinal_rebar_string(beam._n1_b, beam._d_b1_b, beam._n2_b, beam._d_b2_b)}"
                    + (
                        f" ++ {beam._format_longitudinal_rebar_string(beam._n3_b, beam._d_b3_b, beam._n4_b, beam._d_b4_b)}"
                        if beam._n3_b != 0
                        else ""
                    )
                )
            )

            if capacity_check:
                # Remove all forces assignments
                node.clear_forces()
                # Create empty force
                empty_force = Forces()
                node.add_forces(empty_force)
                # Perform the shear check
                shear_results = node.check_shear()
                flexure_results = node.check_flexure()
                # Common data
                common_data = {
                    "Beam": beam.label,
                    "b": beam.width.magnitude,
                    "h": beam.height.magnitude,
                    "As,top": rebar_f_top,
                    "As,bot": rebar_f_bot,
                    "Av": rebar_v,
                    "As,top,real": round(beam._A_s_top.to("cm**2").magnitude, 1),
                    "As,bot,real": round(beam._A_s_bot.to("cm**2").magnitude, 1),
                    "Av,real": round(shear_results["Av"][0], 1),
                }

                common_units = {
                    "Beam": "",
                    "b": "cm",
                    "h": "cm",
                    "As,top": "",
                    "As,bot": "",
                    "Av": "",
                    "As,top,real": "cm²",
                    "As,bot,real": "cm²",
                    "Av,real": "cm²/m",
                }

                # Code-specific data #TODO
                if isinstance(self.concrete, Concrete_EN_1992_2004):
                    code_specific_data = {
                        "MRd,top": round(beam._M_Rd_top.to("kN*m").magnitude, 1),
                        "MRd,bot": round(beam._M_Rd_bot.to("kN*m").magnitude, 1),
                        "VRd": shear_results["V_Rd"][0],
                    }
                    code_specific_units = {
                        "MRd,top": "kNm",
                        "MRd,bot": "kNm",
                        "VRd": "kN",
                    }
                else:  # ACI 318-19 or CIRSOC 201-25 design code
                    code_specific_data = {
                        "ØMn,top": round(beam._phi_M_n_top.to("kN*m").magnitude, 1),
                        "ØMn,bot": round(beam._phi_M_n_bot.to("kN*m").magnitude, 1),
                        "ØVn": shear_results["ØVn"][0],
                    }
                    code_specific_units = {
                        "ØMn,top": "kNm",
                        "ØMn,bot": "kNm",
                        "ØVn": "kN",
                    }

                # Merge data dictionaries
                merged_data = {**common_data, **code_specific_data}
                # Extract DCR values and remove them from the merged dictionary
                dcr_keys = ["DCRb,top", "DCRb,bot", "DCRv"]
                dcr_values = {key: merged_data.pop(key) for key in dcr_keys}
                # Rebuild results_dict with DCR values placed at the end
                results_dict = OrderedDict({**merged_data, **dcr_values})

                # Merge units dictionaries
                merged_units = {**common_units, **code_specific_units}

                # Extract and reorder DCR unit keys
                dcr_units = {key: merged_units.pop(key) for key in dcr_keys}

                # Build the units row with DCRs at the end
                units_row = pd.DataFrame([OrderedDict({**merged_units, **dcr_units})])
                # Restore the original forces after capacity check
                node.clear_forces()
                node.add_forces(original_forces)
            else:
                # Perform the shear check
                shear_results = (
                    node.check_shear().iloc[1:].reset_index(drop=True)
                )  # Skip the first row (units)
                flexure_results = (
                    node.check_flexure().iloc[1:].reset_index(drop=True)
                )  # Skip the first row (units)
                # Common data that doesn't change between codes
                common_data = {
                    "Beam": beam.label,
                    "b": int(beam.width.magnitude),
                    "h": int(beam.height.magnitude),
                    "cc": int(beam.c_c.magnitude),
                    "As,top": rebar_f_top,
                    "As,bot": rebar_f_bot,
                    "Av": rebar_v,
                    "As,req,top": round(beam._A_s_req_top.to("cm**2").magnitude, 1),
                    "As,req,bot": round(beam._A_s_req_bot.to("cm**2").magnitude, 1),
                    "Av,req": round(shear_results["Av,req"][0], 2),
                    "Av,real": round(shear_results["Av"][0], 2),
                    "DCRb,top": round(beam._DCRb_top, 3),
                    "DCRb,bot": round(beam._DCRb_bot, 3),
                    "DCRv": shear_results["DCR"][0],
                }

                common_units = {
                    "Beam": "",
                    "b": "cm",
                    "h": "cm",
                    "cc": "mm",
                    "As,top": "",
                    "As,bot": "",
                    "Av": "",
                    "As,req,top": "cm²/m",
                    "As,req,bot": "cm²/m",
                    "Av,req": "cm²/m",
                    "Av,real": "cm²/m",
                    "DCRb,top": "",
                    "DCRb,bot": "",
                    "DCRv": "",
                }

                # Code-specific data
                if isinstance(self.concrete, Concrete_EN_1992_2004):  # TODO
                    code_specific_data = {
                        "MEd": round(flexure_results["Mu"][0], 1),
                        "VEd": round(shear_results["Vu"][0], 1),
                        "NEd": round(shear_results["Nu"][0], 1),
                        "MRd,top": round(beam._M_Rd_top.to("kN*m").magnitude, 1),
                        "MRd,bot": round(beam._M_Rd_bot.to("kN*m").magnitude, 1),
                        "VRd": round(shear_results["V_Rd"][0], 1),
                    }
                    code_specific_units = {
                        "MEd": "kNm",
                        "VEd": "kN",
                        "NEd": "kN",
                        "MRd,top": "kNm",
                        "MRd,bot": "kNm",
                        "VRd": "kN",
                    }
                else:  # ACI 318-19 or CIRSOC 201-25 design code
                    code_specific_data = {
                        "Mu": round(flexure_results["Mu"][0], 1),
                        "Vu": round(shear_results["Vu"][0], 1),
                        "Nu": round(shear_results["Nu"][0], 1),
                        "ØMn,top": round(beam._phi_M_n_top.to("kN*m").magnitude, 1),
                        "ØMn,bot": round(beam._phi_M_n_bot.to("kN*m").magnitude, 1),
                        "ØVn": round(shear_results["ØVn"][0], 1),
                    }
                    code_specific_units = {
                        "Mu": "kNm",
                        "Vu": "kN",
                        "Nu": "kN",
                        "ØMn,top": "kNm",
                        "ØMn,bot": "kNm",
                        "ØVn": "kN",
                    }

                # Merge data dictionaries
                merged_data = {**common_data, **code_specific_data}
                # Extract DCR values and remove them from the merged dictionary
                dcr_keys = ["DCRb,top", "DCRb,bot", "DCRv"]
                dcr_values = {key: merged_data.pop(key) for key in dcr_keys}
                # Rebuild results_dict with DCR values placed at the end
                results_dict = OrderedDict({**merged_data, **dcr_values})

                # Merge units dictionaries
                merged_units = {**common_units, **code_specific_units}

                # Extract and reorder DCR unit keys
                dcr_units = {key: merged_units.pop(key) for key in dcr_keys}

                # Build the units row with DCRs at the end
                units_row = pd.DataFrame([OrderedDict({**merged_units, **dcr_units})])

            # Add the results to the list
            results_list.append(results_dict)

        # Convert results list into a DataFrame
        results_df = pd.DataFrame(results_list)

        # Combine the units row with the results DataFrame
        final_df = pd.concat([units_row, results_df], ignore_index=True)
        return final_df

    def shear_results(
        self, index: Optional[int] = None, capacity_check: bool = False
    ) -> DataFrame:
        """
        Get shear results for one or all beams.
        Includes a units row only once at the top.
        """
        if index is not None:
            if index - 1 >= len(self.nodes):
                raise IndexError(f"Index {index} is out of range for the beam list.")
            node = self.nodes[max(index - 1, 0)]
            df = self._process_beam_for_check(node, "shear", capacity_check)
            units_row = node.section._get_units_row_shear()  # type: ignore
            df_all = pd.concat([units_row, df], ignore_index=True)
        else:
            # For all nodes, only include the units row once
            results = []
            units_row_added = False
            for item in self.nodes:
                df = self._process_beam_for_check(item, "shear", capacity_check)
                if not units_row_added:
                    units_row = item.section._get_units_row_shear()  # type: ignore
                    results.append(units_row)
                    units_row_added = True
                results.append(df)

            df_all = pd.concat(results, ignore_index=True)

        # Separate the units row (first row) from the data
        units_row = df_all.iloc[[0]]  # DataFrame with 1 row
        data_rows = df_all.iloc[1:].copy()

        # Recombine units + data
        df_final = pd.concat([units_row, data_rows], ignore_index=True)
        return df_final

    def flexure_results(
        self, index: Optional[int] = None, capacity_check: bool = False
    ) -> DataFrame:
        """
        Get flexure results for one or all beams.
        Includes a units row only once at the top.
        """
        if index is not None:
            if index - 1 >= len(self.nodes):
                raise IndexError(f"Index {index} is out of range for the beam list.")
            node = self.nodes[max(index - 1, 0)]
            df = self._process_beam_for_check(node, "flexure", capacity_check)
            units_row = node.section._get_units_row_flexure()  # type: ignore
            df_all = pd.concat([units_row, df], ignore_index=True)
        else:
            # For all nodes, only include the units row once
            results = []
            units_row_added = False
            for item in self.nodes:
                df = self._process_beam_for_check(item, "flexure", capacity_check)
                if not units_row_added:
                    units_row = item.section._get_units_row_flexure()  # type: ignore
                    results.append(units_row)
                    units_row_added = True
                results.append(df)

            df_all = pd.concat(results, ignore_index=True)

        # Separate the units row (first row) from the data
        units_row = df_all.iloc[[0]]  # DataFrame with 1 row
        data_rows = df_all.iloc[1:].copy()

        # Recombine units + data
        df_final = pd.concat([units_row, data_rows], ignore_index=True)
        return df_final

    def _process_beam_for_check(
        self, node: Node, check_type: str, capacity_check: bool
    ) -> DataFrame:
        """
        Shared method to process beam for either shear or flexure checks.

        :param node: Node object to check
        :param check_type: Either 'shear' or 'flexure'
        :param capacity_check: If True, performs capacity check (resets forces)
        :return: Results DataFrame
        """
        original_forces = [copy.deepcopy(force) for force in node.get_forces_list()]

        if capacity_check:
            node.clear_forces()
            node.add_forces(Forces())  # Add empty force

        # Run the appropriate check
        if check_type == "shear":
            results = node.check_shear().iloc[1:].reset_index(drop=True)
        elif check_type == "flexure":
            results = node.check_flexure().iloc[1:].reset_index(drop=True)
        else:
            raise ValueError("check_type must be either 'shear' or 'flexure'")

        # Restore original forces if we did a capacity check
        if capacity_check:
            node.clear_forces()
            node.add_forces(original_forces)
        return results
