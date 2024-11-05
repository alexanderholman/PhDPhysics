from .IonPosition import IonPosition

class Species:
    NAMES = [
        # todo: atomic symbols
    ]
    name: str = ""
    ion_positions: list[IonPosition] = []
    def __init__(
            self,
            name: str,
            ion_positions: list[IonPosition]
        ) -> None:
        self.name = name
        self.ion_positions = ion_positions