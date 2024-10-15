from .Lattice import Lattice
from .Species import Species

class POSCAR:
    COORDINATES = ["Direct", "Cartesian"]
    comment: str = ""
    scaling_factor: float = 1.0
    lattice: Lattice
    species: list[Species] = []
    selective_dynamics: bool = True
    coordinate_mode: str = "Direct"

    def __init__(
            self,
            comment: str,
            species: list[Species],
            scaling_factor: float = 1.0,
            lattice: Lattice = Lattice(
                a1=[5.43, 0.0, 0.0],
                a2=[0.0, 5.43, 0.0],
                a3=[0.0, 0.0, 5, 43]
            ),
            selective_dynamics: bool = True,
            coordinate_mode: str = "Direct"
    ) -> None:
        self.comment = comment
        self.scaling_factor = scaling_factor
        self.lattice = lattice
        self.species = species
        self.selective_dynamics = selective_dynamics
        self.coordinate_mode = coordinate_mode

    def __str__(self) -> str:
        # generate a VASP POSCAR string

        # add comment to the POSCAR string
        poscar = self.comment + "\n"

        # add scaling factor to the POSCAR string
        poscar += str(self.scaling_factor) + "\n"

        # add lattice vectors to the POSCAR string
        poscar += " ".join(map(str, self.lattice.a1)) + "\n"
        poscar += " ".join(map(str, self.lattice.a2)) + "\n"
        poscar += " ".join(map(str, self.lattice.a3)) + "\n"

        # add species names and quantities to the POSCAR string
        for s in self.species:
            poscar += s.name + " "
        poscar += "\n"
        for s in self.species:
            poscar += str(len(s.ion_positions)) + " "
        poscar += "\n"

        # add selective dynamics and coordinate mode to the POSCAR string
        if self.selective_dynamics:
            poscar += "Selective Dynamics\n"
        poscar += self.coordinate_mode + "\n"

        # add ion positions to the POSCAR string
        for s in self.species:
            for ion_position in s.ion_positions:
                poscar += str(ion_position.x) + " "
                poscar += str(ion_position.y) + " "
                poscar += str(ion_position.z) + " "
                if self.selective_dynamics:
                    poscar += ("T" if ion_position.x.selective_dynamics else "F") + " "
                    poscar += ("T" if ion_position.y.selective_dynamics else "F") + " "
                    poscar += ("T" if ion_position.z.selective_dynamics else "F") + " "
                poscar += "\n"
        return poscar.rstrip()