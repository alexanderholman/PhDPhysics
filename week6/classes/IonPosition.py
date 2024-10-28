from .Position import Position

class IonPosition:
    x: Position
    y: Position
    z: Position
    def __init__(
            self,
            x: Position,
            y: Position,
            z: Position
        ) -> None:
        self.x = x
        self.y = y
        self.z = z

    def __str__(self) -> str:
        return f"({self.x}, {self.y}, {self.z})"