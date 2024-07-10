from dataclasses import dataclass
import forallpeople

forallpeople.environment('structural', top_level=True)

@dataclass
class Forces:
    My: float = 0 * kN * m  # type: ignore

    def get_moment(self):
        return self.My

    def set_moment(self, My):
        self.My = My


# Ejemplo de uso de la clase Forces
f1 = Forces() #Creando con valor por defecto
print(f1.get_moment())  # Imprimirá 0 kN·m

f2 = Forces(10 * kN * m) # Creando con argumento
print(f2.get_moment())  # Imprimirá 10 kN·m