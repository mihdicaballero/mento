import forallpeople
forallpeople.environment('structural',top_level=True)
#Useful extra units
cm = 1e-2*m # type: ignore

class Settings: 
    def __init__(self, settings_dict=None):
        # Default settings
        self.settings = {
            'clear_cover': 25 * mm, # type: ignore
            'clear_spacing': 20 * mm, # type: ignore
            'stirrup_diameter': 6 * mm, # type: ignore
            'vibrator_size': 30 * mm, # type: ignore
            'layers_spacing': 25 * mm, # type: ignore
            # Add other settings with their default values
        }
        
        # Update defaults with provided settings if any
        if settings_dict:
            self.settings.update(settings_dict)

    def get_setting(self, key):
        return self.settings.get(key)

    def set_setting(self, key, value):
        self.settings[key] = value

    