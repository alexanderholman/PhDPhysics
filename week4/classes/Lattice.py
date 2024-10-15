class Lattice:
    a1: list[float] = [0.0, 0.0, 0.0]
    a2: list[float] = [0.0, 0.0, 0.0]
    a3: list[float] = [0.0, 0.0, 0.0]
    def __init__(
            self,
            a1: list[float],
            a2: list[float],
            a3: list[float]
        ) -> None:
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3