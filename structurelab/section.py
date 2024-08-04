from settings import Settings
from material import create_concrete, SteelBar, Concrete

from dataclasses import dataclass
import forallpeople
forallpeople.environment('structural', top_level=True)
import numpy as np


# Definir algunas unidades adicionales útiles
cm = 1e-2 * m  # type: ignore

@dataclass
class Section:
    def __init__(self, name:str, settings=None):
        self._name=name
        self._settings = settings if settings is not None else Settings()
    
    def get_name(self):
        return self._name


@dataclass
class ConcreteSection(Section):
    def __init__(self, name: str, concrete, steelBar: str):  # type: ignore
        super().__init__(name)
        self.concrete = concrete
        self.steelBar = steelBar
        self.cc = self._settings.get_setting('clear_cover')

@dataclass
class RectangularConcreteSection(ConcreteSection):
    def __init__(self, name: str, concrete, steelBar: str, width: float, depth: float):
        super().__init__(name, concrete, steelBar)
        self.width = width
        self.depth = depth

    def get_area(self):
        return self._width * self._depth



@dataclass
class Beam(RectangularConcreteSection):
    def __init__(self, name: str, concrete: Concrete, steelBar: SteelBar, width: float, depth: float):  # type: ignore
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
            raise ValueError("Tipo de hormigón no soportado")
    '''

    def check_flexure_ACI_318_19(self):
        max_rho=self.__determine_maximum_flexural_reinforcement_ratio_ACI_318_19()
        print(max_rho)
        pass        



def main():
    # Ejemplo de uso
    concrete=create_concrete(name="H30",f_c=30*MPa, design_code="ACI 318-19") # type: ignore
    steelBar=SteelBar(name="ADN 420", f_y=420*MPa) # type: ignore
    section = Beam(
        name="V-40x50",
        concrete=concrete,
        steelBar=steelBar,
        width=400 * mm,  # type: ignore
        depth=500 * mm,  # type: ignore
    )

    print(f"Nombre de la sección: {section.get_name()}")
    section.design_flexure_ACI_318_19(500*kN*m)  # type: ignore
    section.design_flexure_ACI_318_19(600*kN*m)  # type: ignore

if __name__ == "__main__":
    main()