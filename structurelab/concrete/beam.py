from dataclasses import dataclass
from structurelab.concrete.rectangular import RectangularConcreteSection
import numpy as np
from structurelab import material
import math
import forallpeople
forallpeople.environment('structural', top_level=True)
# Definir algunas unidades adicionales útiles
cm = 1e-2 * m  # type: ignore


@dataclass
class Beam(RectangularConcreteSection):
    def __init__(self, name: str, concrete: material.Concrete, steelBar: material.SteelBar, width: float, depth: float):  # type: ignore
        super().__init__(name, concrete, steelBar, width, depth)

    def __determine_maximum_flexural_reinforcement_ratio_ACI_318_19(self):
        # Determination of maximum reinforcement ratio
        concrete_properties=self.concrete.get_properties()
        beta_1=concrete_properties["beta_1"]
        f_c=concrete_properties["f_c"]
        epsilon_c=concrete_properties["epsilon_c"]
        rebar_properties=self.steelBar.get_properties()
        f_y=rebar_properties["f_y"]
        epsilon_y=rebar_properties["epsilon_y"]

        epsilon_min_rebar_ACI_318_19=epsilon_y+epsilon_c # ESTO CHEQUEARLO BIEN, CREO QUE ES ASI, PERO REVISAR

        rho_max=0.85*beta_1*f_c/f_y*(epsilon_c/(epsilon_c+epsilon_min_rebar_ACI_318_19))
        return rho_max

    def __calculate_phi_ACI_318_19(self, epsilon_mas_deformado:float):
        # CREO QUE NO LA USO PARA NADA, PERO OJO, NO SE. 
        rebar_properties=self._steelBar.get_properties()
        concrete_properties=self._concrete.get_properties()
        epsilon_c=concrete_properties["epsilon_c"]
        rebar_properties=self._steelBar.get_properties()
        epsilon_y=rebar_properties["epsilon_y"]


        if epsilon_mas_deformado<=epsilon_y:
            return 0.65
        elif epsilon_mas_deformado<=epsilon_y+epsilon_c:
            return (0.9-0.65)*(epsilon_mas_deformado-epsilon_y)/epsilon_c+0.65
        else:
            return 0.9

    def design_flexure_ACI_318_19(self, M_u:float):
        a_assumed=0.18*self._depth # Partimos de un valor propuesto para la altura del bloque de compresion de Whitney.

        concrete_properties=self.concrete.get_properties()
        f_c=concrete_properties["f_c"]
        rebar_properties=self.steelBar.get_properties()
        f_y=rebar_properties["f_y"]

        tol=0.01
        error=tol*2

        phi_assumed=0.9 # Asumimos que trabaja dominada a traccion (en realidad mas que asumirlo, lo forzamso)
        #while phi_assumed==0.9:
        while(error>=tol):
            A_s_calc=M_u/phi_assumed/(f_y*(0.9*self._depth-a_assumed/2)) # DONDE PUSE 0.9*self._depth en realidad va d, pero no lo tenemos bien el cover aun :(
            a=A_s_calc*f_y/(0.85*f_c*self._width)
            error=abs(a_assumed-a)/a_assumed
            a_assumed=a

        A_s_min=max((0.25*np.sqrt(f_c / MPa)*MPa/f_y*self._depth*self._width/cm**2) , (1.4*MPa/f_y*self._depth*self._width/cm**2))*cm**2# type: ignore

        #ACA HAY QUE HACER ALGO CON LOS SETTINGS PARA QUE POR DEFECTO NO HAGA LO DE REDUCIR CUANTIA MINIMA
        # PORQUE ES MEDIO PELIGROSO, PORQUE AFECTA A LA CAPACIDAD POR CORTE :( A DONINI NO LE GUSTA :(
        # PODEMOS PONER ESO EN LOS SETTINGS :)
        # VERLOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

        settingMISTERIOSOS= False # ACA VER ESTO DESPUES, ES PARA LO DE LA REDUCCIOn

        if A_s_calc>A_s_min:
            self.__A_s_calculated=A_s_calc
        elif  4*A_s_calc/3 > A_s_min:
            self.__A_s_calculated=A_s_min
        else: 
            if settingMISTERIOSOS:
                self.__A_s_calculated=4*A_s_calc/3
            else:
                self.__A_s_calculated=A_s_min
            
        rho_max=self.__determine_maximum_flexural_reinforcement_ratio_ACI_318_19()
        A_s_max=rho_max*self._depth*self._width
        if self.__A_s_calculated > A_s_max:
            print(f"El armado requerido es mayor al maximo")
        else:
            print(f"PARA EL MOMENTO {M_u}, EL ARMADO DE CALCULO ES {self.__A_s_calculated}")


    '''
    def desgin_flexure():
        if isinstance(self.concrete, concreteACI):
            self.design_flexure_ACI()
        else:
            raise ValueError("Concrete type not supported")
    '''

    def check_flexure_ACI_318_19(self):
        max_rho=self.__determine_maximum_flexural_reinforcement_ratio_ACI_318_19()
        print(max_rho)
        pass   

    def check_shear_ACI_318_19(self, V_u:float, N_u:float, A_s, d_b=0.5, s=8, n_legs=2):
        concrete_properties=self.concrete.get_properties()
        f_c=concrete_properties["f_c"]
        rebar_properties=self.steelBar.get_properties()
        f_y=rebar_properties["f_y"]
        self._settings.load_aci_318_19_settings()
        phi_v = self._settings.get_setting('phi_v')
        f_yt = min(f_y,400*MPa) #type: ignore
        lambda_factor = self._settings.get_setting('lambda')

        # Minimum shear reinforcement calculation
        # 'Minimum reinforcement should be placed if the factored shear Vu is greater than half the shear capacity of the concrete, reduced by 0.5ϕVc. It is assumed that minimum reinforcement is required.
        # Rebar needed, V_u > φ_v*V_c/2
        A_vmin = max((0.75 * math.sqrt(f_c / psi) * psi) / (f_yt * self._width) , (50 * psi) / (f_yt * self._width))  # type: ignore

        # Shear reinforcement calculations
        A_db = (d_b ** 2) * math.pi / 4  # Area of one stirrup leg
        A_v = n_legs * A_db  # Total area of stirrups
        A_vs = A_v / s  # Stirrup area per unit length

        d_bs_ini = self._settings.get_setting('longitudinal_diameter') # Longitudinal diameter if none is defined
        d = self._depth-self.cc-d_bs_ini
        V_s = A_v * f_yt * d / s  # Shear contribution of reinforcement
        phi_V_s = phi_v * V_s  # Reduced shear contribution of reinforcement

        # Effective shear area and longitudinal reinforcement ratio
        A_cv = self._width * d  # Effective shear area
        A_g = self._width * self._depth  # Gross area
        rho_w = A_s / A_cv  # Longitudinal reinforcement ratio
        
        # Size modification factor
        lambda_s = math.sqrt(2 / (1 + d / (10*inch))) #type: ignore

        # Concrete shear strength calculation
        sigma_Nu = min(N_u / (6 * A_g), 0.05 * f_c)  # Axial stress influence

        # Concrete shear capacity depending on whether min rebar is present
        k_c_min_rebar = max(2 * lambda_factor * math.sqrt(f_c / psi) * psi + sigma_Nu, #type: ignore
                            8 * lambda_factor * rho_w ** (1/3) * math.sqrt(f_c / psi) * psi + sigma_Nu) #type: ignore

        V_cmin = 0*kN  # type: ignore
        V_cmax = (5 * lambda_factor * math.sqrt(f_c / psi) * psi) * A_cv  # Maximum concrete shear strength #type: ignore
        
        # Calculate actual concrete shear strength
        V_c = min(V_cmax, max(V_cmin, k_c_min_rebar * A_cv))
        phi_V_c = phi_v * V_c  # Reduced concrete shear strength
        
        # Maximum total shear capacity
        V_max = V_c + (8 * lambda_factor * math.sqrt(f_c / psi) * psi) * A_cv #type: ignore
        phi_V_max = phi_v * V_max  # Reduced maximum shear capacity

        # Total shear strength
        phi_V_n = phi_v * (V_c + V_s)  # Total reduced shear strength (concrete + rebar)

        # Check results
        result = {
            'A_vmin': A_vmin,  # Minimum shear reinforcement area
            'A_vs': A_vs,  # Provided stirrup reinforcement per length
            'phi_V_c': phi_V_c,  # Concrete contribution to shear capacity
            'phi_V_s': phi_V_s,  # Reinforcement contribution to shear capacity
            'phi_V_n': phi_V_n,  # Total shear capacity
            'phi_V_max': phi_V_max,  # Maximum shear capacity
            'shear_ok': V_u <= phi_V_n,  # Check if applied shear is within total capacity
            'max_shear_ok': V_u <= phi_V_max  # Check if applied shear is within max shear capacity
        }

        return result



def main():
    # Ejemplo de uso
    concrete=material.create_concrete(name="H30",f_c=30*MPa, design_code="ACI 318-19") # type: ignore
    steelBar=material.SteelBar(name="ADN 420", f_y=420*MPa) # type: ignore
    section = Beam(
        name="V-40x50",
        concrete=concrete,
        steelBar=steelBar,
        width=400 * mm,  # type: ignore
        depth=500 * mm,  # type: ignore
    )

    print(f"Nombre de la sección: {section.get_name()}")
    section.design_flexure_ACI_318_19(500*kN*m)  # type: ignore

if __name__ == "__main__":
    main()