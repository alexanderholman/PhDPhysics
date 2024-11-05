class Position:
    value: float = 0.0
    selective_dynamics: bool = True
    def __init__(
            self,
            value: float,
            selective_dynamics: bool = True
        ) -> None:
        self.value = value
        self.selective_dynamics = selective_dynamics

    def __str__(self) -> str:
        return str(self.value)