from devtools import debug
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

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
from mento import Forces, RectangularBeam, Node, OneWaySlab
from mento.rebar import Rebar
from mento.summary import BeamSummary


#######################################################################################
def clear_console() -> None:
    """
    Auxiliar function for clear console
    Clears the console based on the operating system.
    """
    if os.name == "nt":  # For Windows
        os.system("cls")
    else:  # For macOS and Linux
        os.system("clear")

#######################################################################################

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



#######################################################################################


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
    f = Forces(V_z=100 * kN, N_x=0 * kN)
    # f = Forces(V_z=130 * kN, N_x=0 * kN)
    beam.set_transverse_rebar(n_stirrups=1, d_b=6 * mm, s_l=25 * cm)
    beam.set_longitudinal_rebar_bot(n1=4, d_b1=16 * mm)
    node = Node(section=beam, forces=f)
    results = node.check_shear()
    # results = node.design_shear()
    print(results)
    print(node.shear_results_detailed())
    # beam.plot()


def shear_ACI_318_19_imperial() -> None:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4000 * psi)
    steelBar = SteelBar(name="ADN 420", f_y=60 * ksi)
    beam = RectangularBeam(
        label="V101",
        concrete=concrete,
        steel_bar=steelBar,
        width=10 * inch,
        height=16 * inch,
        c_c=1.5 * inch,
    )
    # section.set_transverse_rebar(
    #     n_stirrups=0, d_b=0.375 * inch, s_l=6 * inch
    # )
    beam.set_longitudinal_rebar_bot(n1=2, d_b1=0.625 * inch)

    # f = Forces(label="ELU_01", V_z = 6 * kip)
    f = Forces(label="ELU_01", V_z=37.727 * kip)
    node = Node(section=beam, forces=f)
    # section.plot()
    # results = node.check_shear()
    results = node.design_shear()
    # print(results)
    # print(beam.shear_design_results)
    node.shear_results_detailed()
    # node.shear_results_detailed_doc()


def shear_ACI_318_19_slab() -> None:
    concrete = Concrete_ACI_318_19(name="C4", f_c=4 * ksi)
    steelBar = SteelBar(name="ADN 420", f_y=60 * ksi)
    slab = OneWaySlab(
        concrete=concrete,
        steel_bar=steelBar,
        c_c=0.75 * inch,
        width=12 * inch,
        height=7 * inch,
        label="Test",
    )

    # f = Forces(label="ELU_01", V_z=30 * kip)
    f = Forces(label="ELU_01", V_z=1.52 * kip)
    node = Node(section=slab, forces=f)
    # section.plot()
    slab.set_slab_longitudinal_rebar_bot(d_b1=0.5 * inch, s_b1=10 * inch)
    results = node.check_shear()
    # results = node.design_shear()
    print(results)
    # print(slab.shear_design_results)
    print(node.shear_results_detailed())
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



#######################################################################################
# TEST ASSOCIATED WITH FLEXURAL FUNCTIONS #############################################
#######################################################################################


def check_flexure_ACI_318_19_imperial() -> None:
    clear_console()
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
    debug(section.settings)
    section.set_longitudinal_rebar_bot(
        n1=2,
        d_b1=1.41 * inch,
        n2=0,
        d_b2=0 * inch,
        n3=2,
        d_b3=1.27 * inch,
        n4=0,
        d_b4=0 * inch,
    )
    section.set_transverse_rebar(n_stirrups=1, d_b=(3/8) *inch, s_l=6 *inch)
    debug(section._d_b1_b, section._d_b3_b)
    section.set_longitudinal_rebar_top(
        n1=2,
        d_b1=0.75 * inch,
        n2=0,
        d_b2=0 * inch,
        n3=0,
        d_b3=1 * inch,
        n4=0,
        d_b4=0 * inch,
    )

    f = Forces(label="Test_01", M_y=400 * kip * ft, V_z=10 * kip)
    node = Node(section=section, forces=f)
    # section.plot()
    results = node.check_flexure()
    print(results)


def flexure_ACI_318_19_metric() -> None:
    # concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    concrete = Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    section = RectangularBeam(
        concrete=concrete,
        steel_bar=steelBar,
        c_c=2.5 * cm,
        width=20 * cm,
        height=60 * cm,
        label="Test",
    )

    f = Forces(label="Test_01", M_y=100 * kNm, V_z=50 * kN)
    node = Node(section=section, forces=f)
    # results = node.check_flexure()
    section.set_longitudinal_rebar_bot(n1=2, d_b1=16 * mm)
    section.set_transverse_rebar(n_stirrups=1, d_b=8 * mm, s_l=10 * cm)
    # node.design_shear()
    section.plot()
    # print(section.shear_results_detailed())
    # print(section.flexure_design_results_bot)
    # node.flexure_results_detailed_doc()


def flexure_ACI_318_19_metric_slab() -> None:
    # concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    concrete = Concrete_ACI_318_19(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)
    slab = OneWaySlab(
        concrete=concrete,
        steel_bar=steelBar,
        c_c=2.5 * cm,
        width=100 * cm,
        height=20 * cm,
        label="Test",
    )

    f = Forces(label="Test_01", M_y=20 * kNm)
    node = Node(section=slab, forces=f)
    slab.set_slab_longitudinal_rebar_bot(d_b1=12 * mm, s_b1=25 * cm)
    print(node.check_flexure())
    slab.flexure_results_detailed()


def check_flexure_EN_1992_2004_TEST_01() -> None:
    clear_console()
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    steelBar = SteelBar(name="B500S", f_y=500 * MPa)
    BEAM_TEST_01 = RectangularBeam(
        label="BEAM_TEST_01",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=60 * cm,
        c_c=2.6 * cm,
    )
    BEAM_TEST_01._stirrup_d_b=6*mm
    BEAM_TEST_01.set_longitudinal_rebar_bot(
        n1=4,
        d_b1=16 * mm,
        n2=0,
        d_b2=10 *mm,
        n3=0,
        d_b3=10 * mm,
        n4=0,
        d_b4=10 * mm,
    )
    BEAM_TEST_01.set_longitudinal_rebar_top(
        n1=0,
        d_b1=10 * mm,
        n2=0,
        d_b2=10 * mm,
        n3=0,
        d_b3=10 * mm,
        n4=0,
        d_b4=10 * mm,
    )
    f = Forces(M_y=150 * kNm)
    node = Node(section=BEAM_TEST_01, forces=f)
    results = node.check_flexure()
    print(results)
    debug(results.iloc[1]["M_Rd"].to(kNm).magnitude )

def check_flexure_EN_1992_2004_TEST_02() -> None:
    clear_console()
    concrete = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    steelBar = SteelBar(name="B500S", f_y=500 * MPa)
    BEAM_TEST_01 = RectangularBeam(
        label="BEAM_TEST_01",
        concrete=concrete,
        steel_bar=steelBar,
        width=20 * cm,
        height=60 * cm,
        c_c=2.6 * cm,
    )
    BEAM_TEST_01._stirrup_d_b=6*mm
    BEAM_TEST_01.set_longitudinal_rebar_bot(
        n1=5,
        d_b1=25 * mm,
        n2=0,
        d_b2=10 *mm,
        n3=0,
        d_b3=10 * mm,
        n4=0,
        d_b4=10 * mm,
    )
    BEAM_TEST_01.set_longitudinal_rebar_top(
        n1=5,
        d_b1=12 * mm,
        n2=0,
        d_b2=10 * mm,
        n3=0,
        d_b3=10 * mm,
        n4=0,
        d_b4=10 * mm,
    )
    f = Forces(M_y=450 * kNm)
    node = Node(section=BEAM_TEST_01, forces=f)
    results = node.check_flexure()
    print(results)
    debug(results.iloc[1]["M_Rd"].to(kNm).magnitude )

def check_flexure_EN_1992_2004_TEST_04() -> None:
    clear_console()
    concrete = Concrete_EN_1992_2004(name="C60", f_c=60 * MPa)
    steelBar = SteelBar(name="B400S", f_y=400 * MPa)
    BEAM_TEST_03 = RectangularBeam(
        label="BEAM_TEST_03",
        concrete=concrete,
        steel_bar=steelBar,
        width=30 * cm,
        height=50 * cm,
        c_c=3 * cm,
    )
    BEAM_TEST_03._stirrup_d_b=6*mm
    BEAM_TEST_03.set_longitudinal_rebar_top(
        n1=6,
        d_b1=25 * mm,
        n2=0,
        d_b2=10 *mm,
        n3=0,
        d_b3=10 * mm,
        n4=0,
        d_b4=10 * mm,
    )
    f = Forces(M_y=-370 * kNm)
    node = Node(section=BEAM_TEST_03, forces=f)
    results = node.check_flexure()
    debug(BEAM_TEST_03._d_bot)
    debug(BEAM_TEST_03._d_top)
    print(results)
    debug(results.iloc[1]["M_Rd"].to(kNm).magnitude )

def design_flexure_EN_1992_2004_test_01() -> None:
    clear_console()
    concrete = Concrete_EN_1992_2004(name="C60", f_c=60 * MPa)
    steelBar = SteelBar(name="B400S", f_y=400 * MPa)
    f1 = Forces(label="C1", M_y=20 * kNm)
    f2 = Forces(label="C2", M_y=-20 * kNm)
    f3 = Forces(label="C3", M_y=-25 * kNm)
    f4 = Forces(label="C4", M_y=150 * kNm)
    f5 = Forces(label="C5", M_y=120 * kNm)
    f6 = Forces(label="C6", M_y=-28 * kNm)
    BEAM_TEST_DESIGN_FLEXURE_01 = RectangularBeam(
        label="test 01",
        concrete=concrete,
        steel_bar=steelBar,
        width=30 * cm,
        height=50 * cm,
        c_c=3 * cm,
    )
    forces = [f1, f2, f3, f4, f5, f6]
    results = BEAM_TEST_DESIGN_FLEXURE_01.design_flexure(forces)
    # print(beam.flexure_design_results_bot,'\n', beam.flexure_design_results_top)
    print(results)

def design_flexure_ACI_318_test_01() -> None:
    concrete = Concrete_ACI_318_19(name="fc 4000", f_c=4000 * psi)
    steelBar = SteelBar(name="fy 60000", f_y=60 * ksi)
    beam = RectangularBeam(
        label="B-12x24",
        concrete=concrete,
        steel_bar=steelBar,
        c_c=1.5*inch,
        width=12 * inch,
        height=24 * inch
    )

    f1 = Forces(label="Test_01", V_z=40 * kip, M_y=400 * kip * ft)
    f2 = Forces(label="Test_02", V_z=100 * kip, M_y=-400 * kip * ft)
    forces = [f1, f2]

    flexure_results = beam.design_flexure(forces)
    print(flexure_results)
    
def flexure_design_test_EN() -> None:
    #TODO Sale mal el DCR
    concrete = Concrete_EN_1992_2004(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="fy 60000", f_y=500*MPa)
    beam = RectangularBeam(
        label="V-20x60",
        concrete=concrete,
        steel_bar=steelBar,
        c_c=2.6*cm,
        width=20 * cm,
        height=60 * cm
    )
    beam._stirrup_d_b = 6*mm

    f = Forces(label="Test_01", M_y=100 * kNm)
    forces = [f]

    flexure_results = beam.design_flexure(forces)
    print(flexure_results)

def flexure_check_test() -> None:
    clear_console()
    concrete = Concrete_ACI_318_19(name="H-25", f_c=25 * MPa)
    steelBar = SteelBar(name="420", f_y=420 * MPa)

    beam = RectangularBeam(
        label="101",
        concrete=concrete,
        steel_bar=steelBar,
        c_c=2.5*cm,
        width=20 * cm,
        height=60 * cm
    )


    beam.set_longitudinal_rebar_bot(n1=2, d_b1=12 * mm, n2=1, d_b2=12 * mm)
    beam.set_longitudinal_rebar_top(n1=2, d_b1=16 * mm)
    f1 = Forces(label="D", M_y=0 * kNm, V_z=50 * kN)
    f2 = Forces(label="L", M_y=-100 * kNm)
    f3 = Forces(label="W", M_y=-50 * kNm)
    f4 = Forces(label="S", M_y=110 * kNm)
    forces = [f1, f2, f3, f4]
    beam.check_flexure(forces)

    # beam.check_shear()
    beam.flexure_results_detailed()
    # beam.flexure_results_detailed_doc()
    # beam.shear_results_detailed_doc()


def flexure_ACI_318_19_metric_single(
    d_long_mm: float, d_stir_mm: float, out_path: str
) -> None:
    """
    Build one RectangularBeam with given bar diameters and save plot as PNG.
    d_long_mm  = diameter (mm) of bottom longitudinal corner bars (n1)
    d_stir_mm  = diameter (mm) of stirrup
    out_path   = full path to png
    """

    # Materials
    concrete = Concrete_CIRSOC_201_25(name="H25", f_c=25 * MPa)
    steelBar = SteelBar(name="ADN 420", f_y=420 * MPa)

    # Section
    section = RectangularBeam(
        concrete=concrete,
        steel_bar=steelBar,
        c_c=2.5 * cm,
        width=20 * cm,
        height=60 * cm,
        label=f"b20_h60_long{int(d_long_mm)}_stir{int(d_stir_mm)}",
    )

    # dummy load etc. (not strictly needed for plotting but keeping consistent)
    f = Forces(label="Test_01", M_y=100 * kNm, V_z=50 * kN)
    node = Node(
        section=section, forces=f
    )  # noqa: F841 (we don't use `node` directly here)

    # Assign bottom longitudinal reinforcement
    # n1=2 bars of diameter d_long_mm
    section.set_longitudinal_rebar_bot(
        n1=2,
        d_b1=d_long_mm * mm,
        # if you also define n2 layer you can add it here,
        # but I’ll leave only n1 since corners are what we’re testing
        n2=0,
        d_b2=None,
    )

    # Assign stirrup
    section.set_transverse_rebar(
        n_stirrups=1,
        d_b=d_stir_mm * mm,
        s_l=10 * cm,
    )

    # Now plot and save instead of just showing
    fig = plt.figure()  # create a figure handle so we can close later
    plt.close(fig)  # we'll let section.plot() create its own fig

    section.plot(show=False)  # this calls plt.show() in your code currently

    # Instead of relying on show(), grab current figure and save:
    current_fig = plt.gcf()
    current_fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(current_fig)

#######################################################################################

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



#######################################################################################



def summary() -> None:
    # --- Step 1: Define materiales and load input file ---
    # conc = Concrete_ACI_318_19(name="C25", f_c=25 * MPa)
    conc = Concrete_CIRSOC_201_25(name="C25", f_c=25 * MPa)
    # conc = Concrete_EN_1992_2004(name="C25", f_c=25 * MPa)
    steel = SteelBar(name="ADN 500", f_y=500 * MPa)
    from pathlib import Path

    # path = Path(__file__).parents[1] / "_notebooks" / "Mento-Input.xlsx"
    # path = "G:\Mi unidad\mento\Tools\Mento-Input.xlsx"
    # input_df = pd.read_excel(path, sheet_name="Design", usecols="B:T", skiprows=4)
    # beam_summary = BeamSummary(concrete=conc, steel_bar=steel, beam_list=input_df)
    # print(beam_summary.data)

    # # --- Step 2: Run design and export suggested rebars ---
    # print(beam_summary.design())
    # beam_summary.export_design("Beam-Design.xlsx")

    # --- Step 3: Import edited design and run checks ---
    # beam_summary.import_design("Beam-Design.xlsx")
    # print(beam_summary.data)
    # capacity = beam_summary.check(capacity_check=True)
    # print(capacity)
    # check = beam_summary.check(capacity_check=False)
    # print(check)
    # # beam_summary.check().to_excel('hola.xlsx', index=False)
    # results_shear = beam_summary.shear_results(capacity_check=False)

    # results = beam_summary.shear_results(capacity_check=True)
    # results_flexure = beam_summary.flexure_results(capacity_check=False)
    # print(results_shear, "\n", results_flexure)
    # beam_summary.nodes[2].shear_results_detailed()




def batch_test_beam_plots() -> None:
    """
    Iterate over combinations of longitudinal bar diameters and stirrup diameters,
    generate plots, and save them to ~/Desktop/BeamPlotTests.
    """

    # diameters to test
    long_diams_mm = [8, 10, 12, 16, 20, 25]  # longitudinal bottom corner bars
    stirrup_diams_mm = [6, 8, 10, 12]  # stirrups

    # output folder on Desktop
    desktop_dir = os.path.join(os.path.expanduser("~"), "Desktop")
    out_dir = os.path.join(desktop_dir, "BeamPlotTests")
    os.makedirs(out_dir, exist_ok=True)

    for d_long in long_diams_mm:
        for d_stir in stirrup_diams_mm:
            filename = f"beam_b20_h60_long{d_long}mm_stir{d_stir}mm.png"
            out_path = os.path.join(out_dir, filename)
            print(f"Generating {filename} ...")
            flexure_ACI_318_19_metric_single(
                d_long_mm=d_long,
                d_stir_mm=d_stir,
                out_path=out_path,
            )

    print(f"Finished. Images saved in: {out_dir}")



############################################################################################




if __name__ == "__main__":
    # units()
    # settings()
    # material()
    # section()
    # rectangular()
    # shear_EN_1992()
    # shear_ACI_318_19_imperial()
    # shear_ACI_318_19_slab()
    flexure_ACI_318_19_metric_slab()
    # shear_CIRSOC_201_2025()
    # check_flexure_ACI_318_19_imperial()
    # check_flexure_EN_1992_2004_TEST_01()
    # check_flexure_EN_1992_2004_TEST_02()
    # check_flexure_EN_1992_2004_TEST_04()
    # design_flexure_EN_1992_2004_test_01()
    # design_flexure_ACI_318_test_01()
    # rebar()
    # summary()
    # batch_test_beam_plots()
