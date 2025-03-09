import math
from typing import Union

# Helper function to convert degrees to radians
def degrees_to_radians(degrees: Union[int, float]) -> float:
    return degrees * (math.pi / 180.0)

# *** (A) ***
class Concrete:
    def __init__(self, fck_value):
        self._f_ck = fck_value  # Characteristic compressive cylinder strength of concrete at 28 days [Pa]
        self._gamma_c = 1.5      # Concrete partial material safety factor
        self._f_cd = self._f_ck / self._gamma_c  # Design value of concrete compressive strength
        self._f_ck_MPa = self._f_ck * 1e-6  # Convert to MPa
        self._f_cm_MPa = self._f_ck_MPa + 8.0  # Mean value of concrete cylinder compressive strength [MPa]

        if self._f_ck <= 50 * 1e6:
            self._f_ctm_MPa = 0.3 * (self._f_ck_MPa ** (2.0 / 3.0))  # Mean value of axial tensile strength [MPa]
            self._f_ctm = self._f_ctm_MPa * 1e6  # Convert back to Pa
            self._epsilon_cu3 = 0.0035  # Ultimate compressive strain
            self._lambda = 0.8  # EN 1992-1-1 (3.19)
            self._eta = 1.0     # EN 1992-1-1 (3.21)
        else:
            self._f_ctm_MPa = 2.11 * math.log(1 + 0.1 * self._f_cm_MPa)  # Mean value of axial tensile strength [MPa]
            self._f_ctm = self._f_ctm_MPa * 1e6  # Convert back to Pa
            self._epsilon_cu3 = 0.001 * (2.6 + 35 * ((90 - self._f_ck_MPa) / 100) ** 4.0)  # Ultimate compressive strain
            self._lambda = 0.8 - ((self._f_ck_MPa - 50) / 400)  # EN 1992-1-1 (3.20)
            self._eta = 1.0 - ((self._f_ck_MPa - 50) / 200)  # EN 1992-1-1 (3.22)

        self._eta_x_fcd = self._eta * self._f_cd  # eta * f_cd

    # Getters
    def f_ck(self):
        return self._f_ck

    def f_ck_MPa(self):
        return self._f_ck_MPa

    def f_ctm(self):
        return self._f_ctm

    def f_ctm_MPa(self):
        return self._f_ctm_MPa

    def gamma_c(self):
        return self._gamma_c

    def f_cd(self):
        return self._f_cd

    def epsilon_cu3(self):
        return self._epsilon_cu3

    def lambda_(self):
        return self._lambda

    def eta(self):
        return self._eta

    def eta_x_fcd(self):
        return self._eta_x_fcd


# *** part of (A) ***
class Steel:
    def __init__(self, fyk_value, fyk_value_trans):
        self._f_yk = fyk_value  # Characteristic yield strength of reinforcement [Pa]
        self._gamma_s = 1.15    # Partial factor for reinforcing steel
        self._f_yd = self._f_yk / self._gamma_s  # Design yield strength of reinforcement
        self._E_s = 200 * 1e9   # Modulus of elasticity of reinforcing steel [Pa]
        self._epsilon_yd = self._f_yd / self._E_s  # Yield strain
        self._f_ywk = fyk_value_trans  # Characteristic yield of shear reinforcement [Pa]
        self._f_ywd = self._f_ywk / self._gamma_s  # Design yield of shear reinforcement
        self._f_ywk_MPa = self._f_ywk * 1e-6  # Convert to MPa

    # Getters
    def f_yk(self):
        return self._f_yk

    def f_yd(self):
        return self._f_yd

    def gamma_s(self):
        return self._gamma_s

    def E_s(self):
        return self._E_s

    def epsilon_yd(self):
        return self._epsilon_yd

    def f_ywk(self):
        return self._f_ywk

    def f_ywd(self):
        return self._f_ywd

    def f_ywk_MPa(self):
        return self._f_ywk_MPa


# *** part of (A) ***
class Geometry:
    def __init__(self, c_nom, phi, phi_s, height, width):
        self._c_nom = c_nom  # Nominal cover [m]
        self._phi = phi      # Diameter of longitudinal bars [m]
        self._phi_s = phi_s  # Diameter of transverse reinforcement (stirrups) [m]
        self._height = height  # Height of the beam [m]
        self._width = width    # Width of the beam [m]
        self._a = self._c_nom + self._phi_s + 0.5 * self._phi  # Effective depth [m]
        self._d = self._height - self._a  # Effective depth [m]
        self._A_c = self._width * self._height  # Concrete cross-sectional area [m²]
        self._b_w = self._width  # Smallest width of the cross-section in the tensile area [m]

    # Getters
    def width(self):
        return self._width

    def height(self):
        return self._height

    def a(self):
        return self._a

    def d(self):
        return self._d

    def A_c(self):
        return self._A_c

    def b_w(self):
        return self._b_w


# *** part of (A) ***
class Reinforcement:
    def __init__(self, concrete, steel, geometry, phi, phi_s, n_st):
        f_ctm = concrete.f_ctm()  # Mean value of axial tensile strength [Pa]
        f_yk = steel.f_yk()       # Characteristic yield strength of reinforcement [Pa]
        width = geometry.width()  # Width of the beam [m]
        d = geometry.d()          # Effective depth [m]
        A_c = geometry.A_c()      # Concrete cross-sectional area [m²]

        A_s_min1 = 0.26 * (f_ctm / f_yk) * width * d
        A_s_min2 = 0.0013 * width * d
        self._A_s_min = max(A_s_min1, A_s_min2)  # Minimum reinforcement area [m²]
        self._A_s_max = 0.04 * A_c               # Maximum reinforcement area [m²]
        self._A_phi = (math.pi * (phi ** 2)) / 4  # Cross-sectional area of one bar [m²]
        self._A_phi_s = (math.pi * (phi_s ** 2)) / 4  # Cross-sectional area of one stirrup [m²]
        self._rho_l_min = self._A_s_min / (width * d)  # Minimum reinforcement ratio
        self._rho_l_max = self._A_s_max / (width * d)  # Maximum reinforcement ratio
        self._n_st = n_st  # Number of stirrups' legs
        self._A_sw = self._n_st * self._A_phi_s  # Area of shear reinforcement [m²]

    # Getters
    def A_s_min(self):
        return self._A_s_min

    def A_s_max(self):
        return self._A_s_max

    def A_phi(self):
        return self._A_phi

    def A_phi_s(self):
        return self._A_phi_s

    def rho_l_min(self):
        return self._rho_l_min

    def rho_l_max(self):
        return self._rho_l_max

    def A_sw(self):
        return self._A_sw

    def n_st(self):
        return self._n_st


# *** (K) ***
class Forces:
    def __init__(self, M_Ed, V_Ed, N_Ed):
        self.M_Ed = M_Ed  # Design moment [Nm]
        self.V_Ed = V_Ed  # Design shear force [N]
        self.N_Ed = N_Ed  # Design axial force [N]


class Angle:
    def __init__(self, theta=45, alpha=90):
        self.theta = degrees_to_radians(theta)  # Angle between concrete compression strut and beam axis [rad]
        self.alpha = degrees_to_radians(alpha)  # Angle between shear reinforcement and beam axis [rad]
        self.cot_theta = 1 / math.tan(self.theta)
        self.tan_theta = math.tan(self.theta)
        self.sin_alpha = math.sin(self.alpha)
        self.cot_alpha = 1 / math.tan(self.alpha)


class ReinforcementAreas:
    def __init__(self, bottom=0, top=0, is_correct=True):
        self.bottom = bottom  # Bottom reinforcement area [m²]
        self.top = top        # Top reinforcement area [m²]
        self.is_correct = is_correct  # Flag to indicate if the reinforcement is within limits


class TransverseReinforcement:
    def __init__(self, s_req=0, density=0, is_correct=True):
        self.s_req = s_req      # Required spacing of stirrups [m]
        self.density = density  # Transverse reinforcement density
        self.is_correct = is_correct  # Flag to indicate if the reinforcement is within limits


# *** (B), (C), (D), (E) ***
def calculate_required_reinforcement_area(M_Ed, concrete, steel, geometry, reinforcement):
    if M_Ed > 0:
        # *** (B) ***
        xi_eff_lim = concrete.lambda_() * (concrete.epsilon_cu3() / (concrete.epsilon_cu3() + steel.epsilon_yd()))
        x_eff_lim = xi_eff_lim * geometry.d()
        M_lim = concrete.eta_x_fcd() * geometry.width() * x_eff_lim * (geometry.d() - 0.5 * x_eff_lim)

        if M_Ed <= M_lim:
            # *** (D) ***
            S_c_eff = M_Ed / (geometry.width() * (geometry.d() ** 2) * concrete.eta_x_fcd())
            xi_eff = 1 - math.sqrt(1 - 2 * S_c_eff)
            x_eff = xi_eff * geometry.d()
            A_s = (concrete.eta_x_fcd() * x_eff * geometry.width()) / steel.f_yd()
            A_s_req = max(A_s, reinforcement.A_s_min())
            return ReinforcementAreas(A_s_req, 0)
        else:
            # *** (E) ***
            a_2 = geometry.a()
            A_s1 = M_lim / ((geometry.d() - 0.5 * x_eff_lim) * steel.f_yd())
            delta_M = M_Ed - M_lim
            A_s2 = delta_M / ((geometry.d() - a_2) * steel.f_yd())
            A_s_req1 = max(A_s1 + A_s2, reinforcement.A_s_min())
            return ReinforcementAreas(A_s_req1, A_s2)
    else:
        return ReinforcementAreas(0, 0)


# *** (F), (G), (H), (I), (J) ***
def calculate_provided_reinforcement_area(required_reinforcement_areas, reinforcement):
    A_s_req1 = required_reinforcement_areas.bottom
    A_s_req2 = required_reinforcement_areas.top
    A_s_req = A_s_req1 + A_s_req2

    n1 = math.ceil(A_s_req1 / reinforcement.A_phi())
    A_s_prov1 = n1 * reinforcement.A_phi()

    n2 = math.ceil(A_s_req2 / reinforcement.A_phi())
    A_s_prov2 = n2 * reinforcement.A_phi()

    A_s_prov = A_s_prov1 + A_s_prov2

    if reinforcement.A_s_min() <= A_s_req <= reinforcement.A_s_max():
        if A_s_prov <= reinforcement.A_s_max():
            return ReinforcementAreas(A_s_prov1, A_s_prov2)
        else:
            return ReinforcementAreas(A_s_prov1, A_s_prov2, False)
    else:
        return ReinforcementAreas(A_s_prov1, A_s_prov2, False)