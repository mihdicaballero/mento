from mento.units import mm

class Settings: 
    def __init__(self, settings_dict=None):
        # Default settings
        self.settings = {
            # Beam design settings
            'clear_cover': 25 * mm, 
            'clear_spacing': 20 * mm, 
            'stirrup_diameter': 6 * mm, 
            'vibrator_size': 30 * mm, 
            'layers_spacing': 25 * mm, 
            'longitudinal_diameter': 12*mm 
            # Add other settings with their default values
            # Safety factors
        }
        
        # Update defaults with provided settings if any
        if settings_dict:
            self.settings.update(settings_dict)

    def load_aci_318_19_settings(self):
        """
        Load settings specific to ACI 318-19.
        This will override only the settings that are different in ACI 318-19.
        """
        aci_318_19_settings = {
            'lambda': 1, # Normalweight concrete
            'phi_v': 0.75, # Shear strength reduction factor
            'phi_c': 0.65, # Compression controlled strength reduction factor
            'phi_t': 0.90,  # Tension controlled strength reduction factor
            'flexural_min_reduction': "True" # If "True", selects 4/3 of the calculated steel if it is less than the minimum. If "False", always uses the code's minimum values.
            # Add other ACI-specific settings
        }

        # Update current settings with ACI 318-19 specific settings
        self.settings.update(aci_318_19_settings)
    
    def get_setting(self, key):
        if key in self.settings:
            return self.settings.get(key)
        else:
            raise KeyError(f"Setting '{key}' does not exist.")

    def set_setting(self, key, value):
        if key in self.settings:
            self.settings[key] = value
        else:
            raise KeyError(f"Setting '{key}' does not exist.")
    
    def update(self, new_settings: dict):
        """Updates the settings dictionary with new values."""
        self.settings.update(new_settings)