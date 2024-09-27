from dataclasses import dataclass
from mento.section import Section


@dataclass
class ConcreteSection(Section):
    def __init__(self, name: str, concrete, steelBar: str, settings=None):  # type: ignore
        super().__init__(name, settings)
        self.concrete = concrete
        self.steelBar = steelBar
        