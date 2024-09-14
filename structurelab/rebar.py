from structurelab.material import create_concrete
from structurelab.settings import Settings
from structurelab.concrete.beam import Beam
import forallpeople
forallpeople.environment('structural',top_level=True)
cm = 1e-2*m # type: ignore

class Rebar:
    rebar_diameters = [6*mm, 8*mm, 10*mm, 12*mm, 16*mm, 20*mm, 25*mm, 32*mm] #type: ignore
    rebar_areas = {
                6*mm: 0.283*cm**2,   #type: ignore
                8*mm: 0.503*cm**2,  #type: ignore
                10*mm: 0.785*cm**2, #type: ignore
                12*mm: 1.131*cm**2, #type: ignore
                16*mm: 2.011*cm**2, #type: ignore
                20*mm: 3.142*cm**2, #type: ignore
                25*mm: 4.909*cm**2,#type: ignore
                32*mm: 8.043*cm**2,#type: ignore
            }
# Add beam class of rebar
    def __init__(self, required_as, beam):
            self.required_as = required_as
            self.beam = beam
            self.settings = beam.get_settings()
            self.cc = self.settings.get('clear_cover')
            self.clear_spacing = self.settings.get('clear_spacing')
            self.def_stirrup_db = self.settings.get('stirrup_diameter')
            self.bars = self.calculate_rebars()

    def calculate_rebars(self):
        effective_width = self.beam._width - 2 * (self.cc + self.def_stirrup_db)
        
        best_combination = None
        min_total_area = float('inf')
        layer = 1 # Assuming a single layer for simplicity
        
        # Try all combinations of rebar sizes and numbers
        try:
            best_combination = None  # Initialize best_combination
            
            for diameter in self.rebar_diameters:
                rebar_area = self.rebar_areas[diameter]  # Convert diameter to rebar area
                num_bars = 2  # Start with 2 bars minimum
                total_as = 0 * cm**2
                
                while num_bars * diameter + (num_bars - 1) * self.clear_spacing <= effective_width:
                    total_as = num_bars * rebar_area
                    
                    if total_as >= self.required_as:
                        if best_combination is None or total_as < min_total_area:
                            min_total_area = total_as
                            available_spacing = (effective_width - (num_bars * diameter)) / (num_bars - 1)
                            best_combination = {
                                'layer_1': layer,
                                'num_bars_1': num_bars,
                                'diameter_1': diameter,
                                'total_as': total_as,
                                'available_spacing_1': available_spacing
                            }
                        break  # Exit the loop when the required area is satisfied
                    
                    num_bars += 1

            # After the loop, check if no valid combination was found
            if best_combination is None:
                raise ValueError("Cannot fit the required reinforcement within the beam width considering clear cover and spacing.")

        except ValueError as e:
            # Handle the ValueError and raise it again
            raise ValueError(str(e)) from None
                
        return best_combination

def main():
    concrete=create_concrete(name="H30",f_c=30*MPa, design_code="ACI 318-19") # type: ignore
    section = Beam(
        name="V 20x50",
        concrete=concrete,
        steelBar="Barras Longitudinales",
        width=20*cm,
        depth=50*cm,
    )
    print(f"Nombre de la sección: {section.get_name()}")
    # Custom settings for a specific beam
    custom_settings = {
        'clear_cover': 30 * mm, #type: ignore
        'stirrup_diameter': 8 * mm,#type: ignore
    }
    # Creating a Settings instance with custom settings
    section.update_settings(custom_settings)
    as_nec=150*cm**2
    rebar = Rebar(required_as=as_nec, beam=section)
    # Call the calculate_rebars method
    best_combination = rebar.calculate_rebars()

    print(f"Best combination: {best_combination}")

if __name__ == "__main__":
    main()