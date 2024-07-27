import forallpeople
forallpeople.environment('structural',top_level=True)
#Useful extra units
cm = 1e-2*m # type: ignore

class Settings: 
    def __init__(self):
        # Concrete clear cover:
        self.cc=2.5*cm # type: ignore

    