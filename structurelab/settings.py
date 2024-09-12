import forallpeople
forallpeople.environment('structural',top_level=True)
#Useful extra units
cm = 1e-2*m # type: ignore

class Settings: 
    def __init__(self, settings_dict=None):
        # Default settings
        self.settings = {
            # Beam design settings
            'clear_cover': 25 * mm, # type: ignore
            'clear_spacing': 20 * mm, # type: ignore
            'stirrup_diameter': 6 * mm, # type: ignore
            'vibrator_size': 30 * mm, # type: ignore
            'layers_spacing': 25 * mm, # type: ignore
            'longitudinal_diameter': 12*mm # type: ignore
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
            'phi_t': 0.90  # Tension controlled strength reduction factor
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

def main():
    # Define custom settings
    custom_settings = {
        'clear_cover': 30 * mm, # type: ignore
        'stirrup_diameter': 8 * mm # type: ignore
    }

    # Initialize the Settings class with custom settings
    settings = Settings(settings_dict=custom_settings)

    # Test getting existing settings
    print(f"Clear Cover: {settings.get_setting('clear_cover')}")  # Output: 30.0 mm
    print(f"Stirrup Diameter: {settings.get_setting('stirrup_diameter')}")  # Output: 8.0 mm

    # Test setting existing settings
    settings.set_setting('clear_cover', 35 * mm) # type: ignore
    print(f"Updated Clear Cover: {settings.get_setting('clear_cover')} ")  # Output: 35.0 mm

    settings = Settings()
    print("Default Settings:", settings.settings)
    settings.load_aci_318_19_settings()
    print("ACI 318-19 Settings:", settings.settings)

if __name__ == "__main__":
    main()