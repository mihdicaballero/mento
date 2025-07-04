from typing import Optional, Dict, Any

from mento.units import mm, inch
from mento.material import Concrete


class Settings:
    default_settings_metric: Dict[str, Any] = {
        # Beam design settings
        "clear_cover": 25 * mm,
        "clear_spacing": 25 * mm,
        "stirrup_diameter_ini": 8 * mm,
        "vibrator_size": 30 * mm,
        "layers_spacing": 25 * mm,
        "max_diameter_diff": 5 * mm,
        "max_bars_per_layer": 5,
        "minimum_longitudinal_diameter": 8 * mm,
    }
    default_settings_imperial: Dict[str, Any] = {
        # Beam design settings
        "clear_cover": 1 * inch,
        "clear_spacing": 1 * inch,
        "stirrup_diameter_ini": 3 / 8 * inch,
        "vibrator_size": 1.25 * inch,
        "layers_spacing": 1 * inch,
        "max_diameter_diff": 2 / 8 * inch,
        "max_bars_per_layer": 5,
        "minimum_longitudinal_diameter": 3 / 8 * inch,
    }

    def __init__(
        self,
        concrete: Optional[Concrete] = None,
        settings: Optional[Dict[str, Any]] = None,
    ):
        """
        Initialize the Settings class with defaults based on unit system.

        Parameters
        ----------
        concrete : Concrete, optional
            A Concrete object to determine the unit system.
        settings : dict, optional
            Custom settings to override defaults.
        """
        # Select default settings based on the concrete's unit system
        if concrete and concrete.unit_system == "imperial":
            self.settings = self.default_settings_imperial.copy()
        else:  # Default to metric
            self.settings = self.default_settings_metric.copy()

        # Update defaults with provided settings if any
        if settings:
            self.update(settings)

        self.ACI_318_19_settings: Dict[str, Any] = {
            "lambda": 1,  # Normalweight concrete
            "phi_v": 0.75,  # Shear strength reduction factor
            "phi_c": 0.65,  # Compression controlled strength reduction factor
            "phi_t": 0.90,  # Tension controlled strength reduction factor
            "flexural_min_reduction": "False",  # True selects 4/3 of calculated steel if it's less than minimum
        }
        self.EN_1992_2004_settings: Dict[str, Any] = {
            "gamma_c": 1.5,
            "gamma_s": 1.15,
            "alpha_cc": 0.85,
            "delta": 0.85
        }

    def load_ACI_318_19_settings(self) -> None:
        """
        Load settings specific to ACI 318-19.
        This will override only the settings that are different in ACI 318-19.
        """
        # Update current settings with ACI 318-19 specific settings
        self.add_settings(self.ACI_318_19_settings)

    def load_EN_1992_2004_settings(self) -> None:
        """
        Load settings specific to ACI 318-19.
        This will override only the settings that are different in ACI 318-19.
        """
        # Update current settings with ACI 318-19 specific settings
        self.add_settings(self.EN_1992_2004_settings)

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
