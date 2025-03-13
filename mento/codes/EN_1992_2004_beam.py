import math
from pint.facets.plain import PlainQuantity
from typing import TYPE_CHECKING
import pandas as pd
from pandas import DataFrame

from mento.material import Concrete_EN_1992_2004
from mento.rebar import Rebar
from mento.units import MPa, mm, kNm, dimensionless, kN, inch
from mento.forces import Forces  


if TYPE_CHECKING:
    from ..beam import RectangularBeam  # Import Beam for type checking only

def _initialize_variables_EN_1992_2004(self: 'RectangularBeam', force: Forces) -> None:
    """
    Initialize variables for EN 1992-2004 design code.
    """
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Set the initial variables
        self._N_Ed = force.N_x
        self._V_Ed_1 = force.V_z  # Consider the same shear at the edge of support and in d
        self._V_Ed_2 = force.V_z  # Consider the same shear at the edge of support and in d

        # Load settings for gamma factors
        self.settings.load_EN_1992_2004_settings()
        self._alpha_cc = self.settings.get_setting('alpha_cc')
        self._gamma_c = self.settings.get_setting('gamma_c')
        self._gamma_s = self.settings.get_setting('gamma_s')
        self._f_ywk = self._f_yk
        self._f_ywd = self._f_ywk/self._gamma_s
        self._f_cd = self._alpha_cc*self._f_ck/self._gamma_c

        # Consider bottom or top tension reinforcement
        self._A_s_tension = self._A_s_bot if force._M_y >= 0*kNm else self._A_s_top

        # Minimum shear reinforcement calculation
        self._alpha = math.radians(90)
        rho_min = 0.08*math.sqrt(self._f_ck.to('MPa').magnitude) / (self._f_ywk)*MPa
        self._A_v_min = rho_min * self.width * math.sin(self._alpha)
        # Compression stress, positive
        if force._M_y >= 0*kNm:
            self._rho_l_bot = min((self._A_s_tension + self._A_p) / (self.width * self._d_shear), 0.02*dimensionless)
        else:
            self._rho_l_top = min((self._A_s_tension + self._A_p) / (self.width * self._d_shear), 0.02*dimensionless)
        # Shear calculation for sections without rebar
        self._k_value = min(1 + math.sqrt(200 * mm / self._d_shear), 2)
        # Positive of compression
        self._sigma_cp = min(self._N_Ed / self.A_x, 0.2*self._f_cd)

def _shear_without_rebar_EN_1992_2004(self: 'RectangularBeam') -> PlainQuantity:
    self._stirrup_d_b = 0*mm
    self._theta = 0
    # Total shear capacity without rebar
    C_rdc = 0.18/self._gamma_c
    v_min = 0.035*self._k_value**(3/2)*math.sqrt(self._f_ck.to('MPa').magnitude)
    k_1 = 0.15
    V_Rd_c_min = ((v_min+k_1*self._sigma_cp.to('MPa').magnitude)* self.width * self._d_shear * MPa).to('kN')
    rho_l  = self._rho_l_bot if self._M_Ed >= 0*kNm else self._rho_l_top
    V_Rd_c = ((C_rdc*self._k_value*(100*rho_l*self._f_ck.to('MPa').magnitude)**(1/3)*MPa\
                +k_1*self._sigma_cp.to('MPa'))* self.width * self._d_shear).to('kN')        
    return max(V_Rd_c_min, V_Rd_c)
    
def _check_shear_EN_1992_2004(self: 'RectangularBeam', Force:Forces) -> DataFrame:
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        # Initialize all the code related variables
        _initialize_variables_EN_1992_2004(self, Force)            
        
        if self._stirrup_n == 0:
            # Calculate V_Rd_c
            self._V_Rd_c = _shear_without_rebar_EN_1992_2004(self)
            # According to EN1992-1-1 §6.2.1(4) minimum shear reinforcement should nevertheless be provided
            # according to EN1992-1-1 §9.2.2. The minimum shear reinforcement may be omitted in members where
            # transverse redistribution of loads is possible (such as slabs) and members of minor importance
            # which do not contribute significantly to the overall resistance and stability of the structure.
            self._A_v_req = self._A_v_min
            #Maximum shear capacity is the same as the concrete capacity
            self._V_Rd = self._V_Rd_c
            self._V_Rd_max = self._V_Rd
            self._max_shear_ok = self._V_Ed_1 <= self._V_Rd_max

        else:
            # Shear reinforcement calculations
            d_bs = self._stirrup_d_b
            s_l = self._stirrup_s_l
            n_legs = self._stirrup_n*2
            A_db = (d_bs ** 2) * math.pi / 4  # Area of one stirrup leg
            A_vs = n_legs * A_db  # Total area of stirrups
            self._A_v = A_vs / s_l  # Stirrup area per unit length
            # Total shear strength with rebar
            alpha_cw = 1  # Non-prestressed members or members subject to tensile stress due to axial force
            v_1 = 0.6 * (1 - self._f_ck.to('MPa').magnitude / 250)  # Strength reduction factor for concrete struts
            z = 0.9 * self._d_shear  # Lever arm

            # The θ angle is lmited between 21,8° ≤ θ ≤ 45°(1 ≤ cot(θ) ≤ 2.5)
            # Check the minimum strut angle θ = 21.8° (cot(θ) = 2.5)
            theta_min: float = math.radians(21.8)
            cot_theta_min: float = 1 / math.tan(theta_min)

            V_Rd_max_min_angle = (alpha_cw * self.width * z * v_1 * self._f_cd / (cot_theta_min +
                                                                                    math.tan(theta_min))).to('kN')

            if self._V_Ed_1 <= V_Rd_max_min_angle:
                # If within the minimum angle
                self._theta = theta_min
                self._cot_theta = cot_theta_min
                self._V_Rd_max = V_Rd_max_min_angle
                self._max_shear_ok = True
            else:
                # Check the maximum strut angle θ = 45° (cot(θ) = 1.0)
                theta_max: float = math.radians(45)
                cot_theta_max = 1 / math.tan(theta_max)
                V_Rd_max_max_angle: PlainQuantity = (alpha_cw * self.width * z * v_1 * self._f_cd / (cot_theta_max +
                                                                                        math.tan(theta_max))).to('kN')

                if self._V_Ed_1 > V_Rd_max_max_angle:
                    self._theta = theta_max
                    self._cot_theta = 1 / math.tan(self._theta)
                    self._V_Rd_max = V_Rd_max_max_angle
                    self._max_shear_ok = False
                else:
                    self._max_shear_ok = True
                    # Determine the angle θ of the strut based on the shear force
                    self._theta = 0.5 * math.asin((self._V_Ed_1 / V_Rd_max_max_angle))
                    self._cot_theta = 1 / math.tan(self._theta)
                    self._V_Rd_max = (alpha_cw * self.width * z * v_1 * self._f_cd / (self._cot_theta +
                                                                                        math.tan(self._theta))).to('kN')
            # Required shear reinforcing area
            self._A_v_req = max((self._V_Ed_2 / (z * self._f_ywd * self._cot_theta)), self._A_v_min)
            self._V_Rd_s = (self._A_v * z * self._f_ywd * self._cot_theta)
            #Maximum shear capacity is the same as the steel capacity
            self._V_Rd = self._V_Rd_s

            # Rebar spacing checks
            section_rebar = Rebar(self)
            n_legs_actual = self._stirrup_n * 2      # Ensure legs are even
            self._stirrup_s_l = self._stirrup_s_l
            self._stirrup_s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs_actual - 1) 
            self._stirrup_s_max_l, self._stirrup_s_max_w =\
                    section_rebar.calculate_max_spacing_EN_1992_2004(self._alpha)

        self._DCRv = round(abs((self._V_Ed_2.to('kN').magnitude / self._V_Rd.to('kN').magnitude)),3)
        # Design results
        results = {
            'Label': self.label, #Beam label
            'Av,min': self._A_v_min.to('cm ** 2 / m'),  # Minimum shear reinforcement area
            'Av,req': self._A_v_req.to('cm ** 2 / m'), # Required shear reinforcing area
            'Av': self._A_v.to('cm ** 2 / m'),  # Provided stirrup reinforcement per unit length
            'VEd,1': self._V_Ed_1.to('kN'), # Max Vu for the design at the support
            'VEd,2': self._V_Ed_2.to('kN'), # Max Vu for the design at d from the support
            'VRd,c': self._V_Rd_c.to('kN'),  # Concrete contribution to shear capacity
            'VRd,s': self._V_Rd_s.to('kN'),  # Reinforcement contribution to shear capacity
            'VRd': self._V_Rd.to('kN'),  # Total shear capacity
            'VRd,max': self._V_Rd_max.to('kN'),  # Maximum shear capacity
            'VEd,1<VRd,max': self._max_shear_ok,  # Check if applied shear is within max shear capacity
            'VEd,2<VRd': self._V_Ed_2 <= self._V_Rd,  # Check if applied shear is within total capacity
            "DCR" :  self._DCRv
        }
        _initialize_dicts_EN_1992_2004_shear(self)
        return pd.DataFrame([results], index=[0])
    else:
        raise ValueError("Concrete type is not compatible with EN 1992 shear check.")

def _design_shear_EN_1992_2004(self: 'RectangularBeam', Force:Forces) -> None:
    return None

def _design_flexure_EN_1992_2004(self: 'RectangularBeam', force: Forces) -> None:
    pass


##########################################################
# RESULTS
##########################################################

def _initialize_dicts_EN_1992_2004_shear(self: 'RectangularBeam') -> None:
    if isinstance(self.concrete, Concrete_EN_1992_2004):
        """Initialize the dictionaries used in check and design methods."""
        self._materials_shear = {
            "Materials": [
            "Section Label",
            "Concrete strength",
            "Steel reinforcement yield strength",
            "Safety factor for concrete",
            "Safety factor for steel",
            "Coefficient for long term effects and loading effects"
            ],
        "Variable": ["","fck", "fywk", "γc", "γs", "αcc"],
            "Value": [self.label, round(self.concrete.f_ck.to('MPa').magnitude,2),
                    round(self.steel_bar.f_y.to('MPa').magnitude,2),
                    self._gamma_c, self._gamma_s, self._alpha_cc
                    ],
            "Unit": ["", "MPa", "MPa", "","",""]
        }
        self._geometry_shear = {
                "Geometry": [
                "Section height",
                "Section width",
                "Clear cover",
                "Longitudinal tension rebar",
            ],
            "Variable": ["h", "b", "rgeom", "As"],
            "Value": [self.height.to('cm').magnitude, self.width.to('cm').magnitude, self.c_c.to('cm').magnitude,
                    round(self._A_s_tension.to('cm**2').magnitude,2)],
            "Unit": ["cm", "cm", "cm", "cm²"]
        }
        self._forces_shear = {
            "Design forces": [
                "Axial, positive for compression",
                "Shear",
            ],
            "Variable": ["NEd", "VEd,2"],
            "Value": [round(self._N_Ed.to('kN').magnitude,2), round(self._V_Ed_2.to('kN').magnitude,2)],
            "Unit": ["kN", "kN"]
        }
        # Min max lists
        if self._V_Rd_s == 0*kN:
            self._stirrup_d_b = 0*mm if self.concrete.unit_system == "metric" else 0 * inch  
        # Min max lists
        min_values = [None, None, self._A_v_min]   # Use None for items without a minimum constraint
        max_values = [self._stirrup_s_max_l, self._stirrup_s_max_w, None]  # Use None for items without a maximum constraint
        current_values = [self._stirrup_s_l, self._stirrup_s_w, self._A_v]  # Current values to check

        # Generate check marks based on the range conditions
        checks = [
            '✔️' if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else '❌'
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_shear_checks_passed = all(check == '✔️' for check in checks)
        self._data_min_max_shear = {
            'Check': ['Stirrup spacing along length', 'Stirrup spacing along width', 'Minimum shear reinforcement'],
            'Unit': ['cm', 'cm', 'cm²/m'],
            'Valor': [round(self._stirrup_s_l.to('cm').magnitude,2), round(self._stirrup_s_w.to('cm').magnitude,2),
            round(self._A_v.to('cm**2/m').magnitude,2)],
            'Min.': ["", "", round(self._A_v_min.to('cm**2/m').magnitude,2)],
            'Max.': [round(self._stirrup_s_max_l.to('cm').magnitude,2), round(self._stirrup_s_max_w.to('cm').magnitude,2), ""],
            'Ok?': checks
        }
        self._shear_reinforcement = {
                "Shear reinforcement strength": [
                "Number of stirrups",
                "Stirrup diameter",
                "Stirrup spacing",
                "Effective height",
                "Minimum shear reinforcing",
                "Required shear reinforcing",
                "Defined shear reinforcing",
                "Shear rebar strength"
            ],
            "Variable": ["ns", "db", "s", "d", "Asw,min","Asw,req","Asw", "VRd,s"],
            "Value": [self._stirrup_n, self._stirrup_d_b.to('mm').magnitude, self._stirrup_s_l.to('cm').magnitude,
                    self._d_shear.to('cm').magnitude, round(self._A_v_min.to('cm**2/m').magnitude,2),
                    round(self._A_v_req.to('cm**2/m').magnitude,2),
                    round(self._A_v.to('cm**2/m').magnitude,2),
                    round(self._V_Rd_s.to('kN').magnitude,2)],
            "Unit": ["", "mm", "cm", "cm", "cm²/m","cm²/m", "cm²/m","kN"]
        }
        check_max = '✔️' if self._max_shear_ok else '❌'
        check_DCR = '✔️' if self._DCRv < 1 else '❌'
        rho_l  = self._rho_l_bot if self._M_Ed >= 0*kNm else self._rho_l_top
        self._shear_concrete = {
                "Shear strength": [
                "Longitudinal reinforcement ratio",
                'k value',
                "Axial stress",
                "Concrete strut angle",
                "Concrete strength",
                "Maximum shear strength",
                "Total shear strength", 
                "Max shear check",
                "Demand Capacity Ratio"
            ],
            "Variable": ["ρl", "k", "σcd","Θ","VRd,c", "VRd,max", "VRd","" ,"DCR"],
            "Value": [round(rho_l.magnitude,4),
                    round(self._k_value,2),round(self._sigma_cp.to('MPa').magnitude,2),
                    round(math.degrees(self._theta),1),
                    round(self._V_Rd_c.to('kN').magnitude,2), round(self._V_Rd_max.to('kN').magnitude,2), 
                    round(self._V_Rd.to('kN').magnitude,2), check_max, round(self._DCRv,2)],
            "Unit": ["", "", "MPa", "deg", "kN", "kN", "kN", "", check_DCR]
        }


# # Helper function to convert degrees to radians
# def degrees_to_radians(degrees: Union[int, float]) -> float:
#     return degrees * (math.pi / 180.0)

# *** (A) ***
# class Concrete:
#     def __init__(self, fck_value):
#         self._f_ck = fck_value  # Characteristic compressive cylinder strength of concrete at 28 days [Pa]
#         self._gamma_c = 1.5      # Concrete partial material safety factor
#         self._f_cd = self._f_ck / self._gamma_c  # Design value of concrete compressive strength
#         self._f_ck_MPa = self._f_ck * 1e-6  # Convert to MPa
#         self._f_cm_MPa = self._f_ck_MPa + 8.0  # Mean value of concrete cylinder compressive strength [MPa]

#         if self._f_ck <= 50 * 1e6:
#             self._f_ctm_MPa = 0.3 * (self._f_ck_MPa ** (2.0 / 3.0))  # Mean value of axial tensile strength [MPa]
#             self._f_ctm = self._f_ctm_MPa * 1e6  # Convert back to Pa
#             self._epsilon_cu3 = 0.0035  # Ultimate compressive strain
#             self._lambda = 0.8  # EN 1992-1-1 (3.19)
#             self._eta = 1.0     # EN 1992-1-1 (3.21)
#         else:
#             self._f_ctm_MPa = 2.11 * math.log(1 + 0.1 * self._f_cm_MPa)  # Mean value of axial tensile strength [MPa]
#             self._f_ctm = self._f_ctm_MPa * 1e6  # Convert back to Pa
#             self._epsilon_cu3 = 0.001 * (2.6 + 35 * ((90 - self._f_ck_MPa) / 100) ** 4.0)  # Ultimate compressive strain
#             self._lambda = 0.8 - ((self._f_ck_MPa - 50) / 400)  # EN 1992-1-1 (3.20)
#             self._eta = 1.0 - ((self._f_ck_MPa - 50) / 200)  # EN 1992-1-1 (3.22)

#         self._eta_x_fcd = self._eta * self._f_cd  # eta * f_cd

#     # Getters
#     def f_ck(self):
#         return self._f_ck

#     def f_ck_MPa(self):
#         return self._f_ck_MPa

#     def f_ctm(self):
#         return self._f_ctm

#     def f_ctm_MPa(self):
#         return self._f_ctm_MPa

#     def gamma_c(self):
#         return self._gamma_c

#     def f_cd(self):
#         return self._f_cd

#     def epsilon_cu3(self):
#         return self._epsilon_cu3

#     def lambda_(self):
#         return self._lambda

#     def eta(self):
#         return self._eta

#     def eta_x_fcd(self):
#         return self._eta_x_fcd


# # *** part of (A) ***
# class Steel:
#     def __init__(self, fyk_value, fyk_value_trans):
#         self._f_yk = fyk_value  # Characteristic yield strength of reinforcement [Pa]
#         self._gamma_s = 1.15    # Partial factor for reinforcing steel
#         self._f_yd = self._f_yk / self._gamma_s  # Design yield strength of reinforcement
#         self._E_s = 200 * 1e9   # Modulus of elasticity of reinforcing steel [Pa]
#         self._epsilon_yd = self._f_yd / self._E_s  # Yield strain
#         self._f_ywk = fyk_value_trans  # Characteristic yield of shear reinforcement [Pa]
#         self._f_ywd = self._f_ywk / self._gamma_s  # Design yield of shear reinforcement
#         self._f_ywk_MPa = self._f_ywk * 1e-6  # Convert to MPa

#     # Getters
#     def f_yk(self):
#         return self._f_yk

#     def f_yd(self):
#         return self._f_yd

#     def gamma_s(self):
#         return self._gamma_s

#     def E_s(self):
#         return self._E_s

#     def epsilon_yd(self):
#         return self._epsilon_yd

#     def f_ywk(self):
#         return self._f_ywk

#     def f_ywd(self):
#         return self._f_ywd

#     def f_ywk_MPa(self):
#         return self._f_ywk_MPa


# # *** part of (A) ***
# class Geometry:
#     def __init__(self, c_nom, phi, phi_s, height, width):
#         self._c_nom = c_nom  # Nominal cover [m]
#         self._phi = phi      # Diameter of longitudinal bars [m]
#         self._phi_s = phi_s  # Diameter of transverse reinforcement (stirrups) [m]
#         self._height = height  # Height of the beam [m]
#         self._width = width    # Width of the beam [m]
#         self._a = self._c_nom + self._phi_s + 0.5 * self._phi  # Effective depth [m]
#         self._d = self._height - self._a  # Effective depth [m]
#         self._A_c = self._width * self._height  # Concrete cross-sectional area [m²]
#         self._b_w = self._width  # Smallest width of the cross-section in the tensile area [m]

#     # Getters
#     def width(self):
#         return self._width

#     def height(self):
#         return self._height

#     def a(self):
#         return self._a

#     def d(self):
#         return self._d

#     def A_c(self):
#         return self._A_c

#     def b_w(self):
#         return self._b_w


# # *** part of (A) ***
# class Reinforcement:
#     def __init__(self, concrete, steel, geometry, phi, phi_s, n_st):
#         f_ctm = concrete.f_ctm()  # Mean value of axial tensile strength [Pa]
#         f_yk = steel.f_yk()       # Characteristic yield strength of reinforcement [Pa]
#         width = geometry.width()  # Width of the beam [m]
#         d = geometry.d()          # Effective depth [m]
#         A_c = geometry.A_c()      # Concrete cross-sectional area [m²]

#         A_s_min1 = 0.26 * (f_ctm / f_yk) * width * d
#         A_s_min2 = 0.0013 * width * d
#         self._A_s_min = max(A_s_min1, A_s_min2)  # Minimum reinforcement area [m²]
#         self._A_s_max = 0.04 * A_c               # Maximum reinforcement area [m²]
#         self._A_phi = (math.pi * (phi ** 2)) / 4  # Cross-sectional area of one bar [m²]
#         self._A_phi_s = (math.pi * (phi_s ** 2)) / 4  # Cross-sectional area of one stirrup [m²]
#         self._rho_l_min = self._A_s_min / (width * d)  # Minimum reinforcement ratio
#         self._rho_l_max = self._A_s_max / (width * d)  # Maximum reinforcement ratio
#         self._n_st = n_st  # Number of stirrups' legs
#         self._A_sw = self._n_st * self._A_phi_s  # Area of shear reinforcement [m²]

#     # Getters
#     def A_s_min(self):
#         return self._A_s_min

#     def A_s_max(self):
#         return self._A_s_max

#     def A_phi(self):
#         return self._A_phi

#     def A_phi_s(self):
#         return self._A_phi_s

#     def rho_l_min(self):
#         return self._rho_l_min

#     def rho_l_max(self):
#         return self._rho_l_max

#     def A_sw(self):
#         return self._A_sw

#     def n_st(self):
#         return self._n_st


# # *** (K) ***
# class Forces:
#     def __init__(self, M_Ed, V_Ed, N_Ed):
#         self.M_Ed = M_Ed  # Design moment [Nm]
#         self.V_Ed = V_Ed  # Design shear force [N]
#         self.N_Ed = N_Ed  # Design axial force [N]


# class Angle:
#     def __init__(self, theta=45, alpha=90):
#         self.theta = degrees_to_radians(theta)  # Angle between concrete compression strut and beam axis [rad]
#         self.alpha = degrees_to_radians(alpha)  # Angle between shear reinforcement and beam axis [rad]
#         self.cot_theta = 1 / math.tan(self.theta)
#         self.tan_theta = math.tan(self.theta)
#         self.sin_alpha = math.sin(self.alpha)
#         self.cot_alpha = 1 / math.tan(self.alpha)


# class ReinforcementAreas:
#     def __init__(self, bottom=0, top=0, is_correct=True):
#         self.bottom = bottom  # Bottom reinforcement area [m²]
#         self.top = top        # Top reinforcement area [m²]
#         self.is_correct = is_correct  # Flag to indicate if the reinforcement is within limits


# class TransverseReinforcement:
#     def __init__(self, s_req=0, density=0, is_correct=True):
#         self.s_req = s_req      # Required spacing of stirrups [m]
#         self.density = density  # Transverse reinforcement density
#         self.is_correct = is_correct  # Flag to indicate if the reinforcement is within limits


# # *** (B), (C), (D), (E) ***
# def calculate_required_reinforcement_area(M_Ed, concrete, steel, geometry, reinforcement):
#     if M_Ed > 0:
#         # *** (B) ***
#         xi_eff_lim = concrete.lambda_() * (concrete.epsilon_cu3() / (concrete.epsilon_cu3() + steel.epsilon_yd()))
#         x_eff_lim = xi_eff_lim * geometry.d()
#         M_lim = concrete.eta_x_fcd() * geometry.width() * x_eff_lim * (geometry.d() - 0.5 * x_eff_lim)

#         if M_Ed <= M_lim:
#             # *** (D) ***
#             S_c_eff = M_Ed / (geometry.width() * (geometry.d() ** 2) * concrete.eta_x_fcd())
#             xi_eff = 1 - math.sqrt(1 - 2 * S_c_eff)
#             x_eff = xi_eff * geometry.d()
#             A_s = (concrete.eta_x_fcd() * x_eff * geometry.width()) / steel.f_yd()
#             A_s_req = max(A_s, reinforcement.A_s_min())
#             return ReinforcementAreas(A_s_req, 0)
#         else:
#             # *** (E) ***
#             a_2 = geometry.a()
#             A_s1 = M_lim / ((geometry.d() - 0.5 * x_eff_lim) * steel.f_yd())
#             delta_M = M_Ed - M_lim
#             A_s2 = delta_M / ((geometry.d() - a_2) * steel.f_yd())
#             A_s_req1 = max(A_s1 + A_s2, reinforcement.A_s_min())
#             return ReinforcementAreas(A_s_req1, A_s2)
#     else:
#         return ReinforcementAreas(0, 0)


# # *** (F), (G), (H), (I), (J) ***
# def calculate_provided_reinforcement_area(required_reinforcement_areas, reinforcement):
#     A_s_req1 = required_reinforcement_areas.bottom
#     A_s_req2 = required_reinforcement_areas.top
#     A_s_req = A_s_req1 + A_s_req2

#     n1 = math.ceil(A_s_req1 / reinforcement.A_phi())
#     A_s_prov1 = n1 * reinforcement.A_phi()

#     n2 = math.ceil(A_s_req2 / reinforcement.A_phi())
#     A_s_prov2 = n2 * reinforcement.A_phi()

#     A_s_prov = A_s_prov1 + A_s_prov2

#     if reinforcement.A_s_min() <= A_s_req <= reinforcement.A_s_max():
#         if A_s_prov <= reinforcement.A_s_max():
#             return ReinforcementAreas(A_s_prov1, A_s_prov2)
#         else:
#             return ReinforcementAreas(A_s_prov1, A_s_prov2, False)
#     else:
#         return ReinforcementAreas(A_s_prov1, A_s_prov2, False)