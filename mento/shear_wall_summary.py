from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
from collections import OrderedDict

import pandas as pd
from pandas import DataFrame

from mento.material import Concrete, SteelBar
from mento.forces import Forces
from mento.shear_wall import ShearWall
from mento import mm, cm, kN, m, kNm, MPa, inch, ft
from mento.i18n import translate_dataframe
from mento.node import Node
from mento.reports.summaries import wall_summary_doc


class ShearWallSummary:
    def __init__(self, concrete: Concrete, steel_bar: SteelBar, wall_list: DataFrame) -> None:
        self.concrete: Concrete = concrete
        self.steel_bar: SteelBar = steel_bar
        self.wall_list: DataFrame = wall_list
        self.units_row: List[str] = []
        self.data: DataFrame = DataFrame()
        self.nodes: List[Node] = []
        self.wall_keys: List[Tuple[str, str]] = []
        self.check_and_process_input()
        self.convert_to_walls()

    # ------------------------------------------------------------------
    # Input processing
    # ------------------------------------------------------------------

    def check_and_process_input(self) -> None:
        self.units_row = self.wall_list.iloc[0].tolist()
        data = self.wall_list.iloc[1:].copy()

        self.units_row = ["" if pd.isna(unit) else unit for unit in self.units_row]
        self.validate_units(self.units_row)

        # Columns 3 onward are numeric (after Level, Label, Comb.)
        data.iloc[:, 3:] = data.iloc[:, 3:].astype(float).fillna(0)

        for i in range(1, len(self.units_row)):
            unit_str = self.units_row[i]
            if unit_str != "":
                unit = self.get_unit_variable(unit_str)
                if isinstance(data.iloc[:, i], pd.Series):
                    data.iloc[:, i] = data.iloc[:, i].apply(lambda x, u=unit: x * u)

        self.data = data

    def validate_units(self, units_row: List[str]) -> None:
        valid_units = {"m", "mm", "cm", "inch", "ft", "kN", "kNm", ""}
        for unit_str in units_row:
            if unit_str and unit_str not in valid_units:
                raise ValueError(f"Invalid unit '{unit_str}' detected. Allowed units: {valid_units}")

    def get_unit_variable(self, unit_str: str) -> Any:
        unit_map: Dict[str, Any] = {
            "mm": mm,
            "cm": cm,
            "m": m,
            "in": inch,
            "inch": inch,
            "ft": ft,
            "kN": kN,
            "kNm": kNm,
            "MPa": MPa,
        }
        if unit_str in unit_map:
            return unit_map[unit_str]
        raise ValueError(f"Unit '{unit_str}' is not recognized.")

    # ------------------------------------------------------------------
    # Grouping rows → walls + nodes
    # ------------------------------------------------------------------

    def convert_to_walls(self) -> None:
        self.nodes = []
        self.wall_keys = []
        groups: OrderedDict[Tuple[str, str], Dict[str, Any]] = OrderedDict()

        for _, row in self.data.iterrows():
            level = str(row["Level"])
            label = str(row["Label"])
            key = (level, label)

            forces = Forces(
                label=row["Comb."],
                N_x=row["Nx"],
                V_z=row["Vz"],
                M_y=row["My"],
            )

            if key not in groups:
                groups[key] = {
                    "level": level,
                    "label": label,
                    "t": row["t"],
                    "lw": row["lw"],
                    "hw": row["hw"],
                    "cc": row["cc"],
                    "dbh": row["dbh"],
                    "sh": row["sh"],
                    "dbv": row["dbv"],
                    "sv": row["sv"],
                    "forces": [forces],
                }
            else:
                self._validate_geometry_consistency(groups[key], row, key)
                groups[key]["forces"].append(forces)

        for key, g in groups.items():
            wall = ShearWall(
                level=g["level"],
                label=g["label"],
                concrete=self.concrete,
                steel_bar=self.steel_bar,
                thickness=g["t"],
                length=g["lw"],
                height=g["hw"],
                c_c=g["cc"],
            )

            if g["dbh"].magnitude != 0 and g["sh"].magnitude != 0:
                wall.set_horizontal_rebar(d_b=g["dbh"], s=g["sh"])
            if g["dbv"].magnitude != 0 and g["sv"].magnitude != 0:
                wall.set_vertical_rebar(d_b=g["dbv"], s=g["sv"])

            node = Node(section=wall, forces=g["forces"])
            self.nodes.append(node)
            self.wall_keys.append(key)

    def _validate_geometry_consistency(
        self, first: Dict[str, Any], row: "pd.Series[Any]", key: Tuple[str, str]
    ) -> None:
        for col, field in [("t", "t"), ("lw", "lw"), ("hw", "hw"), ("cc", "cc")]:
            val_first = first[field]
            val_row = row[col]
            if abs(val_first.magnitude - val_row.magnitude) > 1e-6:
                raise ValueError(
                    f"Geometry mismatch for wall {key}: '{col}' differs between rows "
                    f"({val_first} vs {val_row}). All rows in the same (Level, Label) "
                    f"group must have identical geometry."
                )

    # ------------------------------------------------------------------
    # Check
    # ------------------------------------------------------------------

    def check(self) -> DataFrame:
        results_list = []

        for node in self.nodes:
            wall: ShearWall = node.section  # type: ignore

            if wall._d_b_h.magnitude == 0 or wall._s_h.magnitude == 0:
                raise ValueError(
                    f"Wall '{wall.level} - {wall.label}' has no horizontal rebar assigned. "
                    f"All walls must have rebar for check(). "
                    f"Either run design() first or provide rebar in the input."
                )
            if wall._d_b_v.magnitude == 0 or wall._s_v.magnitude == 0:
                raise ValueError(
                    f"Wall '{wall.level} - {wall.label}' has no vertical rebar assigned. "
                    f"All walls must have rebar for check(). "
                    f"Either run design() first or provide rebar in the input."
                )

            node.check_shear()

            limiting = wall.limiting_case_shear
            dcr = limiting["DCR"]

            rebar_h = f"Ø{wall._d_b_h.to('mm').magnitude:.0f}/{wall._s_h.to('cm').magnitude:.0f}"
            rebar_v = f"Ø{wall._d_b_v.to('mm').magnitude:.0f}/{wall._s_v.to('cm').magnitude:.0f}"

            all_checks_pass = wall._all_wall_shear_checks_passed
            status = "✅" if dcr <= 1 and all_checks_pass else "❌"

            results_dict = OrderedDict(
                {
                    "Level": wall.level,
                    "Label": wall.label,
                    "t": int(wall.thickness.to("cm").magnitude),
                    "lw": round(wall.length.to("m").magnitude, 2),
                    "hw": round(wall.height.to("m").magnitude, 2),
                    "Horiz.": rebar_h,
                    "Vert.": rebar_v,
                    "ρt": round(float(wall._rho_t.magnitude), 5),
                    "ρl": round(float(wall._rho_l.magnitude), 5),
                    "Vu,max": round(limiting["Vu"], 1),
                    "ØVn": round(limiting["ØVn"], 1),
                    "DCR": round(dcr, 3),
                    "Status": status,
                }
            )
            results_list.append(results_dict)

        if self.concrete.unit_system == "metric":
            v_unit = "kN"
        else:
            v_unit = "kip"

        units_row = pd.DataFrame(
            [
                OrderedDict(
                    {
                        "Level": "",
                        "Label": "",
                        "t": "cm",
                        "lw": "m",
                        "hw": "m",
                        "Horiz.": "",
                        "Vert.": "",
                        "ρt": "",
                        "ρl": "",
                        "Vu,max": v_unit,
                        "ØVn": v_unit,
                        "DCR": "",
                        "Status": "",
                    }
                )
            ]
        )

        results_df = pd.DataFrame(results_list)
        return translate_dataframe(pd.concat([units_row, results_df], ignore_index=True))

    # ------------------------------------------------------------------
    # Design
    # ------------------------------------------------------------------

    def design(self) -> DataFrame:
        design_df: DataFrame = self.data.reset_index(drop=True).copy()

        for i, node in enumerate(self.nodes):
            wall: ShearWall = node.section  # type: ignore
            node.design_shear()

            key = self.wall_keys[i]
            mask = design_df.apply(
                lambda r: str(r["Level"]) == key[0] and str(r["Label"]) == key[1],
                axis=1,
            )
            design_df.loc[mask, "dbh"] = wall._d_b_h  # type: ignore
            design_df.loc[mask, "sh"] = wall._s_h  # type: ignore
            design_df.loc[mask, "dbv"] = wall._d_b_v  # type: ignore
            design_df.loc[mask, "sv"] = wall._s_v  # type: ignore

        self.design_data = design_df
        print("✅ Shear wall design completed for all walls in Summary.")
        return design_df

    # ------------------------------------------------------------------
    # Shear results
    # ------------------------------------------------------------------

    def shear_results(self, index: Optional[int] = None) -> DataFrame:
        if index is not None:
            if index < 1 or index > len(self.nodes):
                raise IndexError(f"Index {index} is out of range. Valid: 1 to {len(self.nodes)}")
            node = self.nodes[index - 1]
            return translate_dataframe(node.check_shear())

        results = []
        units_row_added = False
        for node in self.nodes:
            df = node.check_shear()
            if not units_row_added:
                results.append(df.iloc[[0]])
                units_row_added = True
            results.append(df.iloc[1:])
        out: DataFrame = pd.concat(results, ignore_index=True)
        return translate_dataframe(out)

    # ------------------------------------------------------------------
    # Excel I/O
    # ------------------------------------------------------------------

    def export_design(self, path: str) -> None:
        if not hasattr(self, "design_data"):
            raise AttributeError("No design data found. Run .design() before exporting.")

        df_numeric = self.design_data.copy()
        for col in df_numeric.columns:
            df_numeric[col] = df_numeric[col].apply(lambda x: x.magnitude if hasattr(x, "magnitude") else x)

        df_export = pd.concat(
            [
                pd.DataFrame([self.units_row], columns=self.wall_list.columns),
                df_numeric,
            ],
            ignore_index=True,
        )
        df_export.to_excel(path, index=False)
        print(f"✅ Shear wall design exported to {path}")

    def import_design(self, path: str) -> None:
        wall_df = pd.read_excel(path)
        self.wall_list = wall_df
        self.check_and_process_input()
        self.convert_to_walls()
        print("✅ Shear wall design imported and summary data updated.")

    # ------------------------------------------------------------------
    # Word export
    # ------------------------------------------------------------------

    def results_detailed_doc(self, index: int = 1) -> None:
        """Export detailed results for one wall, plus summary tables for all, to Word.

        The assembly lives in :mod:`mento.reports.summaries`.
        """
        wall_summary_doc(self, index)
