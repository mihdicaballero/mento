
import math
import forallpeople
forallpeople.environment('structural',top_level=True)
cm = 1e-2*m # type: ignore

class Rebar:
    rebar_diameters = [6*mm, 8*mm, 10*mm, 12*mm, 16*mm, 20*mm, 25*mm, 32*mm] #type: ignore
    rebar_areas = {
                6*mm: 0.283*cm**2,   
                8*mm: 0.503*cm**2,  
                10*mm: 0.785*cm**2, 
                12*mm: 1.131*cm**2, 
                16*mm: 2.011*cm**2, 
                20*mm: 3.142*cm**2, 
                25*mm: 4.909*cm**2,
                32*mm: 8.043*cm**2,
            }
    def __init__(self, required_as, beam, rebar_settings):
            self.required_as = required_as
            self.beam = beam
            self.clear_cover = rebar_settings['clear_cover']
            self.clear_spacing = rebar_settings['clear_spacing']
            self.stirrup_diameter = rebar_settings['stirrup_diameter']
            self.bars = self.calculate_rebars()

    def calculate_rebars(self):
        effective_width = self.beam.width - 2 * (self.clear_cover + self.stirrup_diameter)
        available_width = effective_width
        
        best_combination = None
        min_total_area = float('inf')
        layer = 1 # Assuming a single layer for simplicity
        
        # Try all combinations of rebar sizes and numbers
        for diameter in self.rebar_diameters:
            rebar_area = self.rebar_areas[diameter]  # Convert diameter to rebar area
            num_bars = 2  # Start with 2 bars minimum
            total_as = 0 * cm**2
            
            while num_bars*diameter + (num_bars - 1) * self.clear_spacing <= effective_width:
                total_as = num_bars * rebar_area
                if total_as >= self.required_as:
                    if total_as < min_total_area:
                        min_total_area = total_as
                        available_spacing = (effective_width-(num_bars*diameter))/(num_bars - 1)
                        best_combination = {
                            'layer_1': layer,
                            'num_bars_1': num_bars,
                            'diameter_1': diameter,
                            'total_as': total_as,
                            'DCR': total_as / self.required_as,
                            'available_spacing_1': available_spacing
                        }
                    break
                num_bars += 1
            
        if best_combination is None:
            raise ValueError("Cannot fit the required reinforcement within the beam width considering clear cover and spacing.")
        
        return best_combination


# Test the class and methods
class BeamSection:
    def __init__(self, name, concrete, steel, width, height):
        self.name = name
        self.concrete = concrete
        self.steel = steel
        self.width = width
        self.height = height

# Example usage:
beam = BeamSection(
    name="V20x50",
    concrete="aci_concrete",
    steel="steel_rebar",
    width=18 * cm,
    height=60 * cm
)
rebar_settings = {
    'clear_cover': 25*mm,       #type: ignore
    'clear_spacing': 20*mm,     #type: ignore
    'stirrup_diameter': 6*mm,  #type: ignore
}

as_nec=5.6*cm**2
rebar = Rebar(required_as=as_nec, beam=beam, rebar_settings=rebar_settings)

# Call the calculate_rebars method
best_combination = rebar.calculate_rebars()

print(f"Best combination: {best_combination}")