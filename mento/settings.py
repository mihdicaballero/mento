from mento.units import mm
from typing import Optional, Dict, Any

class Settings:
    def __init__(self, settings: Optional[Dict[str, Any]] = None):
        # Default settings
        self.settings: Dict[str, Any] = {
            # Beam design settings
            'clear_cover': 25 * mm, 
            'clear_spacing': 20 * mm, 
            'stirrup_diameter': 6 * mm, 
            'vibrator_size': 30 * mm, 
            'layers_spacing': 25 * mm,
        }
        
        # Update defaults with provided settings if any
        if settings:
            self.settings.update(settings)

    def load_aci_318_19_settings(self) -> None:
        """
        Load settings specific to ACI 318-19.
        This will override only the settings that are different in ACI 318-19.
        """
        aci_318_19_settings: Dict[str, Any] = {
            'lambda': 1, # Normalweight concrete
            'phi_v': 0.75, # Shear strength reduction factor
            'phi_c': 0.65, # Compression controlled strength reduction factor
            'phi_t': 0.90,  # Tension controlled strength reduction factor
            'flexural_min_reduction': "True"  # True selects 4/3 of calculated steel if it's less than minimum
        }

        # Update current settings with ACI 318-19 specific settings
        self.settings.update(aci_318_19_settings)
    
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
    
    def update(self, new_settings: Dict[str, Any]) -> None:
        """Updates the settings dictionary with new values."""
        self.settings.update(new_settings)