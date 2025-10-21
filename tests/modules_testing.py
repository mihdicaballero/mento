from devtools import debug
import numpy as np
import pandas as pd

from mento import cm, MPa, m, kg, kN, kNm, inch, ft, kip, lb, psi, mm, ksi
from mento.settings import BeamSettings
from mento.material import (
    Concrete,
    SteelBar,
    Concrete_EN_1992_2004,
    Concrete_ACI_318_19,
    Concrete_CIRSOC_201_25,
)
from mento.section import Section
from mento.rectangular import RectangularSection
from mento import Forces, RectangularBeam, Node
from mento.rebar import Rebar
from mento.summary import BeamSummary


def units() -> None:
    debug(2 * cm, 3 * MPa, 4 * kg, 1 * mm, 1 * m, 3 * kN, 2 * kNm)
    debug(1 * psi, 1 * lb, 1 * kip, 1 * psi, 1 * ksi, 1 * inch, 1 * ft)
    N = 3 * kN
    A = 3 * m**2
    s = N / A
    a = 3.25 * m
    print(N.to_compact(), N.to("kip"), A, s)
    # Get only the value
    debug(N.magnitude, N.to("kN").magnitude)
    debug(N, A, s)

    # How to format specific output
    debug("Length is {:~P}".format(a))


def settings() -> None:
    # Metric defaults
    settings = BeamSettings()
    print(settings, "\n")

    # Metric defaults with override
    settings = BeamSettings(clear_spacing=35 * mm)
    print(settings, "\n")

    # Imperial defaults with override
    settings = BeamSettings(unit_system="imperial", clear_spacing=2 * inch)
    print(settings, "\n")

    # Create with imperial defaults
    # settings = BeamSettings(unit_system="imperial")

    # Update settings
    # settings.update(
    #     max_bars_per_layer=4
    # )


def section() -> None:
    concrete = Concrete("C25")
    steel_bar = SteelBar(name="B500S", f_y=500 * MPa)
    section = Section(
        c_c=25 * mm, concrete=concrete, steel_bar=steel_bar, label="Test Section"
    )
    print(section)


def rectangular() -> None:
    concrete = Concrete("C25", f_c=25 * MPa)
    steel_bar = SteelBar(name="ADN 420", f_y=420 * MPa)
    section = RectangularSection(
        concrete=concrete,
        steel_bar=steel_bar,
        c_c=25 * mm,
        width=20 * cm,
        height=50 * cm,
        label="Test Rectangular Section",
    )
    debug(section)
    debug(section.A_x, section.I_y, section.I_z)


def material() -> None:
    # Test cases
    concrete = Concrete_ACI_318_19(name="H25", f_c=4 * ksi)
    # concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    print(concrete)
    debug(concrete.get_properties())
    concrete_eurocode = Concrete_EN_1992_2004(name="C25/30", f_c=25 * MPa)
    print(concrete_eurocode)
    debug(concrete_eurocode.get_properties())
    steelbar = SteelBar(name="ADN 500", f_y=500 * MPa)
    debug(steelbar.get_properties())
    print(steelbar)
    # steelstrand = SteelStrand(name='Y1860',f_y=1700*MPa)
    # debug(steelstrand.get_properties())
    # print(concrete.f_c.to('MPa'), concrete.f_c.to('MPa').magnitude)


def shear_EN_1992() -> None:
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    steelBar = SteelBar(name="B500S", f_y=500 * MPa)
    beam = RectangularBeam(
        label="101",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=60 * cm,
        c_c=2.6 * cm,
    )
    # f = Forces(V_z=100*kN, M_y=100*kNm)
    # f = Forces(V_z=100 * kN, N_x=0 * kN)
    f = Forces(V_z=130 * kN, N_x=0 * kN)
    beam.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=25 * cm)
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam, forces=f)
    # results = node.check_shear()
    results = node.design_shear()
    print(results)
    print(node.shear_results_detailed())
    # beam.plot()


def shear_ACI_318_19() -> None:
    concrete = Concrete_ACI_318_19(name="C4", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam = RectangularBeam(
        concrete=concrete,
        steel_bar=steelBar,
        c_c=2.5 * cm,
        width=20 * cm,
        height=40 * cm,
        label="Test",
    )
    # section.set_transverse_rebar(
    #     n_stirrups=0, d_b=0.375 * inch, s_l=6 * inch
    # )
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)

    # f = Forces(label="ELU_01", V_z=30 * kip)
    f = Forces(label="ELU_01", V_z=37.727 * kip)
    node = Node(section=beam, forces=f)
    # section.plot()
    # results = node.check_shear()
    results = node.design_shear()
    print(results)
    print(beam.shear_design_results)
    # print(node.shear_results_detailed())
    # node.shear_results_detailed_doc()


def shear_CIRSOC_201_2025() -> None:
    concrete = Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam = RectangularBeam(
        concrete=concrete,
        steel_bar=steelBar,
        c_c=2.5 * cm,
        width=20 * cm,
        height=40 * cm,
        label="Test",
    )
    # section.set_transverse_rebar(
    #     n_stirrups=0, d_b=0.375 * inch, s_l=6 * inch
    # )

    f1 = Forces(V_z=73 * kN, M_y=109 * kNm)
    node = Node(section=beam, forces=f1)
    # section.plot()
    # results = node.check_shear()
    # results = node.design_shear()
    results = node.design_flexure()
    print(results)
    print(beam.flexure_design_results_bot)
    # print(node.shear_results_detailed())


def flexure_ACI_318_19_imperial() -> None:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steelBar = SteelBar(name="ADN 420", f_y=60 * ksi)
    section = RectangularBeam(
        concrete=concrete,
        steel_bar=steelBar,
        c_c=1.5 * inch,
        width=12 * inch,
        height=24 * inch,
        label="Test Rectangular Beam",
    )
    section.set_longitudinal_rebar_bot(
        n1=2,
        d_b1=1.41 * inch,
        n2=0,
        d_b2=1.128 * inch,
        n3=2,
        d_b3=1.27 * inch,
        n4=0,
        d_b4=1 * inch,
    )
    section.set_longitudinal_rebar_top(
        n1=2,
        d_b1=0.75 * inch,
        n2=0,
        d_b2=1.128 * inch,
        n3=0,
        d_b3=1 * inch,
        n4=0,
        d_b4=1 * inch,
    )

    f = Forces(label="Test_01", M_y=400 * kip * ft, V_z=10 * kip)
    node = Node(section=section, forces=f)
    # section.plot()
    results = node.check_flexure()
    print(results)


def flexure_ACI_318_19_metric() -> None:
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    section = RectangularBeam(
        concrete=concrete,
        steel_bar=steelBar,
        c_c=2.5 * cm,
        width=20 * cm,
        height=60 * cm,
        label="Test",
    )

    f = Forces(label="Test_01", M_y=50 * kNm)
    node = Node(section=section, forces=f)
    # section.plot()
    # results = node.check_flexure()
    results = node.design_flexure()
    # print(results)
    print(section.flexure_design_results_bot)
    # node.flexure_results_detailed_doc()


def rebar() -> None:
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam = RectangularBeam(
        label="V 20x50",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=40 * cm,
        c_c=25 * mm,
    )
    as_req = 9.26 * cm**2

    beam_rebar = Rebar(beam)
    long_rebar_df = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=as_req)
    best_design = beam_rebar.longitudinal_rebar_design
    print(long_rebar_df, best_design)


def rebar_df() -> None:
    concrete = Concrete_ACI_318_19(name="H30", f_c=30 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    beam_settings = BeamSettings(stirrup_diameter_ini=8 * mm)
    beam = RectangularBeam(
        label="V 20x50",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=50 * cm,
        c_c=30 * mm,
    )
    beam_rebar = Rebar(beam)
    # Create a list of required steel areas from 0.5 to 10 with a step of 0.5 cm²
    as_req_list = np.arange(0.5, 10.5, 0.5) * cm**2
    # Initialize an empty DataFrame to store the results
    results_df = pd.DataFrame()

    # Loop through each required steel area
    for as_req in as_req_list:
        # Run the longitudinal_rebar_ACI_19 method
        long_rebar_df = beam_rebar.longitudinal_rebar_ACI_318_19(A_s_req=as_req)

        # Extract the first row of the resulting DataFrame
        first_row = long_rebar_df.iloc[0:1].copy()

        # Add a column to store the required steel area
        first_row[
            "as_req"
        ] = as_req.magnitude  # Store the magnitude (value without units)

        # Append the first row to the results DataFrame
        results_df = pd.concat([results_df, first_row], ignore_index=True)

    # Display the results DataFrame
    print(results_df)


def summary() -> None:
    # --- Step 1: Define materiales and load input file ---
    # conc = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    conc = Concrete_CIRSOC_201_25(name="C25", f_c=25 * MPa)
    # conc = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 500", f_y=500 * MPa)
    from pathlib import Path

    # path = Path(__file__).parents[1] / "_notebooks" / "Mento-Input.xlsx"
    path = "G:\Mi unidad\mento\Tools\Mento-Input.xlsx"
    input_df = pd.read_excel(path, sheet_name="Design", usecols="B:T", skiprows=4)
    beam_summary = BeamSummary(concrete=conc, steel_bar=steel, beam_list=input_df)
    # print(beam_summary.data)

    # # --- Step 2: Run design and export suggested rebars ---
    # print(beam_summary.design())
    # beam_summary.export_design("Beam-Design.xlsx")

    # --- Step 3: Import edited design and run checks ---
    beam_summary.import_design("Beam-Design.xlsx")
    # print(beam_summary.data)
    # capacity = beam_summary.check(capacity_check=True)
    # print(capacity)
    check = beam_summary.check(capacity_check=False)
    print(check)
    # # beam_summary.check().to_excel('hola.xlsx', index=False)
    # results_shear = beam_summary.shear_results(capacity_check=False)

    # results = beam_summary.shear_results(capacity_check=True)
    # results_flexure = beam_summary.flexure_results(capacity_check=False)
    # print(results_shear, "\n", results_flexure)
    # beam_summary.nodes[2].shear_results_detailed()


if __name__ == "__main__":
    # units()
    # settings()
    # material()
    # section()
    # rectangular()
    # shear_EN_1992()
    # shear_ACI_318_19()
    # shear_CIRSOC_201_2025()
    # flexure_ACI_318_19_imperial()
    # flexure_ACI_318_19_metric()
    # rebar()
    summary()
