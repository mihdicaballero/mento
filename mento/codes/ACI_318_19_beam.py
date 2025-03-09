import math
from pint.facets.plain import PlainQuantity
import pandas as pd
from typing import TYPE_CHECKING, Dict, Any
import warnings

from mento.material import Concrete_ACI_318_19
from mento.rebar import Rebar
from mento.units import MPa, mm, kN, inch, ksi, psi, cm, m, kNm, ft, kip
from mento.forces import Forces  


if TYPE_CHECKING:
    from ..beam import RectangularBeam  # Import Beam for type checking only

def _initialize_variables_ACI_318_19(self: 'RectangularBeam', force: Forces) -> None:
        self._N_u = force._N_x
        self._V_u = force._V_z
        self._M_u = force._M_y
        if self._M_u > 0*kNm: 
            self._M_u_bot = self._M_u
            self._M_u_top = 0*kNm
        else:
            self._M_u_bot = 0*kNm
            self._M_u_top = self._M_u
        self.settings.load_ACI_318_19_settings()
        self.phi_v = self.settings.get_setting('phi_v')
        self.phi_t = self.settings.get_setting('phi_t')
        self.lambda_factor = self.settings.get_setting('lambda')
        self.f_yt = _calculate_f_yt_aci(self)
        # Consider bottom or top tension reinforcement
        self._A_s_tension = self._A_s_bot if force._M_y >= 0*kNm else self._A_s_top

def _calculate_shear_reinforcement_aci(self: 'RectangularBeam') -> None:
    V_s = self._A_v * self.f_yt * self._d_shear  # Shear contribution of reinforcement
    self._phi_V_s = self.phi_v * V_s  # Reduced shear contribution of reinforcement

def _calculate_effective_shear_area_aci(self: 'RectangularBeam') -> None:
    self._A_cv = self.width * self._d_shear  # Effective shear area
    self._rho_w = self._A_s_tension.to('cm**2') / self._A_cv.to('cm**2')  # Longitudinal reinforcement ratio
    if self.concrete.unit_system == "metric":
        self._lambda_s = math.sqrt(2 / (1 + 0.004*self._d_shear/mm))
    else:
        self._lambda_s = math.sqrt(2 / (1 + self._d_shear / (10 * inch)))

def _calculate_concrete_shear_strength_aci(self: 'RectangularBeam') -> None:
    f_c = self.concrete.f_c
    self._sigma_Nu = min(self._N_u / (6 * self.A_x), 0.05 * f_c)  # Axial stress influence
    if self.concrete.unit_system == "metric":
        V_cmin = 0 * kN
        if self._A_v < self._A_v_min:
            if self._A_s_tension == 0*cm**2:
                warnings.warn("Longitudinal rebar As cannot be zero if A_v is less than A_v_min.", UserWarning)
            self._k_c_min = 0.66 * self._lambda_s * self.lambda_factor * self._rho_w ** (1/3)\
                * math.sqrt(f_c / MPa)* MPa + self._sigma_Nu
        else:
            self._k_c_min = max(0.17 * self.lambda_factor * math.sqrt(f_c / MPa) * MPa + self._sigma_Nu,
                            0.66 * self.lambda_factor * self._rho_w ** (1/3) * math.sqrt(f_c / MPa) * MPa\
                                    + self._sigma_Nu)
    else:
        V_cmin = 0 * kip
        if self._A_v < self._A_v_min:
            self._k_c_min = 8 * self._lambda_s * self.lambda_factor * self._rho_w ** (1/3) * math.sqrt(f_c / psi)\
                    * psi + self._sigma_Nu
        else:
            self._k_c_min = max(2 * self.lambda_factor * math.sqrt(f_c / psi) * psi + self._sigma_Nu,
                            8 * self.lambda_factor * self._rho_w ** (1/3) * math.sqrt(f_c / psi) * psi +\
                                    self._sigma_Nu)
    # Maximum concrete shear strength
    if self.concrete.unit_system == "metric":
        V_cmax = (0.42 * self.lambda_factor * math.sqrt(self.concrete.f_c / MPa) * MPa) * self._A_cv
    else:
        V_cmax = (5 * self.lambda_factor * math.sqrt(self.concrete.f_c / psi) * psi) * self._A_cv
    self.V_c = min(V_cmax, max(V_cmin, self._k_c_min * self._A_cv))
    self._phi_V_c = self.phi_v * self.V_c

def _calculate_max_shear_capacity_aci(self: 'RectangularBeam') -> None:
    "Formula for maximum total shear capacity (V_max)"
    if self.concrete.unit_system == "metric":
        V_max = self.V_c + (0.66 * self.lambda_factor * math.sqrt(self.concrete.f_c / MPa) * MPa) * self._A_cv
    else:
        V_max = self.V_c + (8 * self.lambda_factor * math.sqrt(self.concrete.f_c / psi) * psi) * self._A_cv
    self._phi_V_max = self.phi_v * V_max
    self._max_shear_ok = self._V_u < self._phi_V_max  

def _calculate_A_v_min_ACI(self: 'RectangularBeam', f_c: PlainQuantity) -> None:
    """Calculate the minimum shear reinforcement based on unit system."""
    # 'Minimum reinforcement should be placed if the factored shear Vu 
    # is greater than half the shear capacity of the concrete,
    # reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
    # Rebar needed, V_u > φ_v*V_c/2 for Imperial system
    f_yt = _calculate_f_yt_aci(self)
    
    if self.concrete.unit_system == "metric":
        self._A_v_min = max(
            (0.062 * math.sqrt(f_c / MPa) * MPa / f_yt) * self.width,
            (0.35 * MPa / f_yt) * self.width
        )
    else:
        self._A_v_min = max(
            (0.75 * math.sqrt(f_c / psi) * psi / f_yt) * self.width,
            (50 * psi / f_yt) * self.width
        )   
    
def _calculate_f_yt_aci(self: 'RectangularBeam') -> PlainQuantity:
    """Determine the yield strength of steel based on unit system."""
    if self.concrete.unit_system == "metric":
        return min(self.steel_bar.f_y, 420 * MPa)
    else:
        return min(self.steel_bar.f_y, 60 * ksi)

def _check_minimum_reinforcement_requirement_aci(self: 'RectangularBeam') -> None:
    if self.concrete.unit_system == "metric":
        if self._V_u < 0.083*self.phi_v*self.lambda_factor*math.sqrt(self.concrete.f_c / MPa) * MPa*self._A_cv:
            self._A_v_req = 0 * cm**2/m
            self._A_v_min = 0 * cm**2/m
            self._max_shear_ok = True
        elif 0.083*self.phi_v*self.lambda_factor*math.sqrt(self.concrete.f_c / MPa) * MPa*self._A_cv\
            < self._V_u < self._phi_V_max:
            _calculate_A_v_min_ACI(self, self.concrete.f_c)
            self._max_shear_ok = True
        else:
            _calculate_A_v_min_ACI(self, self.concrete.f_c)
            self._max_shear_ok = False
    else:
        if self._V_u < self.phi_v*self.lambda_factor*math.sqrt(self.concrete.f_c / psi) * psi*self._A_cv:
            self._A_v_req = 0 * inch**2/ft
            self._A_v_min = 0 * inch**2/ft
            self._max_shear_ok = True
            # if self._V_u != 0*kN:
            #     self._stirrup_d_b = 0*inch
        elif self.phi_v*self.lambda_factor*math.sqrt(self.concrete.f_c / psi) * psi*self._A_cv\
            < self._V_u < self._phi_V_max:
            _calculate_A_v_min_ACI(self, self.concrete.f_c)
            self._max_shear_ok = True
        else:
            _calculate_A_v_min_ACI(self, self.concrete.f_c)
            self._max_shear_ok = False

def _calculate_V_s_req(self: 'RectangularBeam') -> None:
    self._V_s_req = self._V_u - self._phi_V_c
    self._A_v_req = max(self._V_s_req / (self.phi_v * self.f_yt * self._d_shear), self._A_v_min).to('cm ** 2 / m')

def _calculate_total_shear_strength_aci(self: 'RectangularBeam') -> None:
    self._phi_V_n = self.phi_v * (self.V_c + self._A_v * self.f_yt * self._d_shear)
    V_d_max = min(self._phi_V_n, self._phi_V_max)
    self._DCRv = round(abs((self._V_u.to('kN').magnitude / V_d_max.to('kN').magnitude)),3)

def _calculate_rebar_spacing_aci(self: 'RectangularBeam') -> None:
    section_rebar = Rebar(self)
    n_legs_actual = self._stirrup_n * 2  # Ensure legs are even
    self._stirrup_s_w = (self.width - 2 * self.c_c - self._stirrup_d_b) / (n_legs_actual - 1)
    self._stirrup_s_max_l, self._stirrup_s_max_w = section_rebar.calculate_max_spacing_ACI_318_19(self._V_u - self._phi_V_c, self._A_cv)
    self._stirrup_s_l = max(self._stirrup_s_l, 0 * inch)
    self._stirrup_s_w = max(self._stirrup_s_w, 0 * inch)

def _check_shear_ACI_318_19(self: 'RectangularBeam', force: Forces) -> pd.DataFrame:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        # Set the initial variables
        _initialize_variables_ACI_318_19(self, force)

        # Minimum shear reinforcement calculation
        _calculate_A_v_min_ACI(self, self.concrete.f_c)
        if self._stirrup_n > 0:
            # Shear reinforcement calculations
            _calculate_shear_reinforcement_aci(self)

        # Effective shear area and longitudinal reinforcement ratio
        _calculate_effective_shear_area_aci(self)

        # Check if minimum reinforcement is required
        _check_minimum_reinforcement_requirement_aci(self)

        # Concrete shear strength calculation
        _calculate_concrete_shear_strength_aci(self)

        # Maximum total shear capacity
        _calculate_max_shear_capacity_aci(self)

        # Calculate required shear reinforcement
        _calculate_V_s_req(self)

        # Total shear strength
        _calculate_total_shear_strength_aci(self)

        # Rebar spacing checks
        _calculate_rebar_spacing_aci(self)

        # Check results and return DataFrame
        results = _compile_results_ACI_shear(self, force)
        _initialize_dicts_ACI_318_19_shear(self)
        return pd.DataFrame([results], index=[0])
    else:
        raise ValueError("Concrete type is not compatible with EN 1992 shear check.")

def _design_shear_ACI_318_19(self: 'RectangularBeam', force: Forces) -> None:
        # Set the initial variables
        _initialize_variables_ACI_318_19(self, force)
        # Minimum shear reinforcement calculation
        _calculate_A_v_min_ACI(self, self.concrete.f_c)
        # Consider that the beam has minimum reinforcement
        self._A_v = self._A_v_min
        # Effective shear area and longitudinal reinforcement ratio
        _calculate_effective_shear_area_aci(self)
        # Concrete shear strength calculation
        _calculate_concrete_shear_strength_aci(self)
        # Maximum total shear capacity
        _calculate_max_shear_capacity_aci(self)
        # Check if minimum reinforcement is required
        _check_minimum_reinforcement_requirement_aci(self)
        # Calculate required shear reinforcement
        _calculate_V_s_req(self)

        return None

def _compile_results_ACI_shear(self: 'RectangularBeam', force: Forces) -> Dict[str, Any]:
    return {
        'Section Label': self.label,
        'Load Combo': force.label,
        'Av,min': self._A_v_min.to('cm ** 2 / m'),
        'Av,req': self._A_v_req.to('cm ** 2 / m'),
        'Av': self._A_v.to('cm ** 2 / m'),
        'Vu': self._V_u.to('kN'),
        'Nu': self._N_u.to('kN'),
        'ØVc': self._phi_V_c.to('kN'),
        'ØVs': self._phi_V_s.to('kN'),
        'ØVn': self._phi_V_n.to('kN'),
        'ØVmax': self._phi_V_max.to('kN'),
        'Vu<ØVmax': self._max_shear_ok,
        'Vu<ØVn': self._V_u <= self._phi_V_n,
        'DCR': self._DCRv
    }

def _initialize_dicts_ACI_318_19_shear(self: 'RectangularBeam') -> None:
    if isinstance(self.concrete, Concrete_ACI_318_19):
        """Initialize the dictionaries used in check and design methods."""
        self._materials_shear = {
            "Materials": [
                "Section Label",
                "Concrete strength",
                "Steel reinforcement yield strength",
                "Concrete density",
                "Normalweight concrete",
                "Safety factor for shear"
            ],
            "Variable": ["","fc", "fy", "wc", "λ", "Øv"],
            "Value": [self.label, round(self.concrete.f_c.to('MPa').magnitude,2), 
                        round(self.steel_bar.f_y.to('MPa').magnitude,2),round(self.concrete.density.to('kg/m**3').magnitude,1),
                        self.settings.get_setting('lambda'), self.settings.get_setting('phi_v')],
            "Unit": ["", "MPa", "MPa", "kg/m³", "", ""]
        }
        self._geometry_shear = {
            "Geometry": [
                "Section height",
                "Section width",
                "Clear cover",
                "Longitudinal tension rebar",
            ],
            "Variable": ["h", "b", "cc", "As"],
            "Value": [self.height.to('cm').magnitude, self.width.to('cm').magnitude, self.c_c.to('cm').magnitude,
                        round(self._A_s_tension.to('cm**2').magnitude,2)],
            "Unit": ["cm", "cm", "cm", "cm²"]
        }
        self._forces_shear = {
            "Design forces": [
                "Axial, positive for compression",
                "Shear",
            ],
            "Variable": ["Nu", "Vu"],
            "Value": [round(self._N_u.to('kN').magnitude,2), round(self._V_u.to('kN').magnitude,2) ],
            "Unit": ["kN", "kN"]
        }
        # Min max lists
        if self._phi_V_s == 0*kN:
            db_min = 0*mm if self.concrete.unit_system == "metric" else 0 * inch
            self._stirrup_d_b = 0*mm if self.concrete.unit_system == "metric" else 0 * inch
        else:
            db_min = 10 * mm if self.concrete.unit_system == "metric" else 3 / 8 * inch
        min_values = [None, None, self._A_v_min, db_min]   # Use None for items without a minimum constraint
        max_values = [self._stirrup_s_max_l, self._stirrup_s_max_w, None, None]  # Use None for items without a maximum constraint
        current_values = [self._stirrup_s_l, self._stirrup_s_w, self._A_v, self._stirrup_d_b]  # Current values to check

        # Generate check marks based on the range conditions
        checks = [
            '✔️' if (min_val is None or curr >= min_val) and (max_val is None or curr <= max_val) else '❌'
            for curr, min_val, max_val in zip(current_values, min_values, max_values)
        ]
        self._all_shear_checks_passed = all(check == '✔️' for check in checks)
        self._data_min_max_shear = {
            'Check': ['Stirrup spacing along length', 'Stirrup spacing along width', 'Minimum shear reinforcement',
                        'Minimum rebar diameter'],
            'Unit': ['cm', 'cm', 'cm²/m', 'mm'],
            'Value': [round(self._stirrup_s_l.to('cm').magnitude,2), round(self._stirrup_s_w.to('cm').magnitude,2),
            round(self._A_v.to('cm**2/m').magnitude,2), round(self._stirrup_d_b.magnitude,0)],
            'Min.': ["", "", round(self._A_v_min.to('cm**2/m').magnitude,2), round(db_min.magnitude,0)],
            'Max.': [round(self._stirrup_s_max_l.to('cm').magnitude,2), round(self._stirrup_s_max_w.to('cm').magnitude,2), "", ""],
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
            "Variable": ["ns", "db", "s", "d", "Av,min","Av,req","Av", "ØVs"],
            "Value": [self._stirrup_n, self._stirrup_d_b.to('mm').magnitude, self._stirrup_s_l.to('cm').magnitude,
                    self._d_shear.to('cm').magnitude, round(self._A_v_min.to('cm**2/m').magnitude,2),
                    round(self._A_v_req.to('cm**2/m').magnitude,2),
                    round(self._A_v.to('cm**2/m').magnitude,2),
                    round(self._phi_V_s.to('kN').magnitude,2)],
            "Unit": ["", "mm", "cm", "cm", "cm²/m","cm²/m", "cm²/m","kN"]
        }
        check_max = '✔️' if self._max_shear_ok else '❌'
        check_FU = '✔️' if self._DCRv < 1 else '❌'
        self._shear_concrete = {
            "Shear strength": [
                "Effective shear area",
                "Longitudinal reinforcement ratio",
                "Size modification factor",
                "Axial stress",
                "Concrete effective shear stress",
                "Concrete strength",
                "Maximum shear strength",
                "Total shear strength", 
                "Max shear check",
                "Demand Capacity Ratio"
            ],
            "Variable": ["Acv", "ρw", "λs", "σNu", "kc", "ØVc", "ØVmax", "ØVn", "" ,"DCR"],
            "Value": [round(self._A_cv.to('cm**2').magnitude,2),
                        round(self._rho_w.magnitude,5),round(self._lambda_s,3),
                        round(self._sigma_Nu.to('MPa').magnitude,2), round(self._k_c_min.to('MPa').magnitude,2),
                        round(self._phi_V_c.to('kN').magnitude,2), round(self._phi_V_max.to('kN').magnitude,2), 
                        round(self._phi_V_n.to('kN').magnitude,2), check_max, round(self._DCRv,2)],
            "Unit": ["cm²", "", "", "MPa", "MPa", "kN", "kN", "kN", "", check_FU]
        }
        self._shear_all_checks = self._all_shear_checks_passed and (check_max == '✔️') and (check_FU == '✔️')