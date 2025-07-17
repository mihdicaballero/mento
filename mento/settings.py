from dataclasses import dataclass, field, fields
from typing import Any, ClassVar, Dict

from mento.units import mm, inch

# Sentinel to detect if user passed a value
_NOT_SET = object()


@dataclass
class BeamSettings:
    """Settings for beam design with separate metric and imperial defaults.
    Must input values with units (in metric or imperial)
    
    Available Parameters with Default Values:
    --------------------------------------
    Unit system: "metric" or "imperial"

    Metric Defaults:
      - clear_spacing: 25 mm
      - stirrup_diameter_ini: 8 mm
      - vibrator_size: 30 mm
      - layers_spacing: 25 mm
      - max_diameter_diff: 5 mm
      - minimum_longitudinal_diameter: 8 mm
      - max_bars_per_layer: 5
    
    Imperial Defaults:
      - clear_spacing: 1 inch
      - stirrup_diameter_ini: 3/8 inch
      - vibrator_size: 1.25 inch
      - layers_spacing: 1 inch
      - max_diameter_diff: 0.25 inch
      - minimum_longitudinal_diameter: 3/8 inch
      - max_bars_per_layer: 5
    """
    
    # Class-level default values (documented but not shown in hover)
    _metric_defaults: ClassVar[Dict[str, Any]] = {
        "clear_spacing": 25 * mm,
        "stirrup_diameter_ini": 8 * mm,
        "vibrator_size": 30 * mm,
        "layers_spacing": 25 * mm,
        "max_diameter_diff": 5 * mm,
        "minimum_longitudinal_diameter": 8 * mm,
        "max_bars_per_layer": 5
    }
    
    _imperial_defaults: ClassVar[Dict[str, Any]] = {
        "clear_spacing": 1 * inch,
        "stirrup_diameter_ini": 3/8 * inch,
        "vibrator_size": 1.25 * inch,
        "layers_spacing": 1 * inch,
        "max_diameter_diff": 0.25 * inch,
        "minimum_longitudinal_diameter": 3/8 * inch,
        "max_bars_per_layer": 5
    }
    
    unit_system: str = "metric"

    clear_spacing: Any = field(default=_NOT_SET)
    stirrup_diameter_ini: Any = field(default=_NOT_SET)
    vibrator_size: Any = field(default=_NOT_SET)
    layers_spacing: Any = field(default=_NOT_SET)
    max_diameter_diff: Any = field(default=_NOT_SET)
    minimum_longitudinal_diameter: Any = field(default=_NOT_SET)
    max_bars_per_layer: Any = field(default=_NOT_SET)

    def __post_init__(self) -> None:
        defaults = (
            self._imperial_defaults if self.unit_system == "imperial"
            else self._metric_defaults
        )

        for f in fields(self):
            if not f.init or f.name == "unit_system":
                continue

            current_value = getattr(self, f.name)

            if current_value is _NOT_SET:
                setattr(self, f.name, defaults[f.name])

        if self.max_bars_per_layer < 1:
            raise ValueError("max_bars_per_layer must be at least 1")
   
    def __str__(self) -> str:
        """Returns only the current settings, excluding class defaults."""
        settings_list = []
        for name, value in vars(self).items():
            # Skip private and special attributes, and the unit_system field
            if not name.startswith('_') and name != 'unit_system':
                if hasattr(value, 'magnitude') and hasattr(value, 'units'):  # It's a Quantity
                    settings_list.append(f"{name}: {value.magnitude:.2f} {value.units:~}")
                else:
                    settings_list.append(f"{name}: {value}")
        return "\n".join(settings_list)

GLOBAL_BEAM_SETTINGS = BeamSettings()