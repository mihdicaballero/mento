from typing import Optional, Dict, Any

from mento.units import mm

class Settings:
    default_settings: Dict[str, Any] = {
            # Beam design settings
            'clear_cover': 25 * mm, 
            'clear_spacing': 25 * mm, 
            'stirrup_diameter_ini': 10 * mm,
            # 'longitudinal_diameter_ini': 16*mm, Esto ya no se usa, porque puse 2Ø8 interno
            'vibrator_size': 30 * mm, 
            'layers_spacing': 25 * mm,
            'max_diameter_diff': 5*mm,
            'max_bars_per_layer': 5,
            'minimum_longitudinal_diameter': 12*mm,
        }
    
    def __init__(self, settings: Optional[Dict[str, Any]] = None):
        # Initialize instance settings with default settings
        self.settings: Dict[str, Any] = self.default_settings.copy()
        
        # Update defaults with provided settings if any
        if settings:
            self.update(settings)

        self.aci_318_19_settings : Dict[str, Any] = {
                'lambda': 1, # Normalweight concrete
                'phi_v': 0.75, # Shear strength reduction factor
                'phi_c': 0.65, # Compression controlled strength reduction factor
                'phi_t': 0.90,  # Tension controlled strength reduction factor
                'flexural_min_reduction': "True"  # True selects 4/3 of calculated steel if it's less than minimum
            }
        self.en_1992_2004_settings : Dict[str, Any] = {
                'gamma_c': 1.5,
                'gamma_s': 1.15, 
                'alpha_cc': 1.00, 
            }
    def load_aci_318_19_settings(self) -> None:
        """
        Load settings specific to ACI 318-19.
        This will override only the settings that are different in ACI 318-19.
        """
        # Update current settings with ACI 318-19 specific settings
        self.add_settings(self.aci_318_19_settings)

    def load_en_1992_2004_settings(self) -> None:
        """
        Load settings specific to ACI 318-19.
        This will override only the settings that are different in ACI 318-19.
        """
        # Update current settings with ACI 318-19 specific settings
        self.add_settings(self.en_1992_2004_settings)
    
    def get_setting(self, key: str) -> Any:
        if key in self.settings:
            return self.settings.get(key)
        else:
            raise KeyError(f"Setting '{key}' does not exist.")

    def set_setting(self, key: str, value: Any) -> None:
        if key in self.settings:
            self.settings[key] = value
        else:
            raise KeyError(f"Setting '{key}' does not exist.")
        
    def add_settings(self, new_settings: Dict[str, Any]) -> None:
            """Adds the settings dictionary to the current settings."""
            for key, value in new_settings.items():
                self.settings[key] = value
    
    def update(self, new_settings: Dict[str, Any]) -> None:
        """Updates the settings dictionary with new values."""
        self.add_settings(new_settings)

    def __str__(self) -> str:
        """Customize the string representation for user-friendly display."""
        output = "Settings:\n"
        for key, value in self.settings.items():
            output += f"  {key}: {value}\n"
        return output.strip()

def settings() -> None:
    settings_test = Settings()
    custom_settings = {'clear_cover': 50*mm}
    settings_test.update(custom_settings)

    print(settings_test)

if __name__ == "__main__":
    settings()