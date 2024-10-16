import os
from ase.io import read, write
from ase.optimize import FIRE
from ase.visualize import view
from mace.calculators import mace_mp
from copy import deepcopy
from random import randint

from numpy.random import random

from .Helper import Helper
from .Lattice import Lattice
from .Species import Species
from .IonPosition import IonPosition
from .Position import Position

class POSCAR:
    COORDINATES = ["Direct", "Cartesian"]
    comment: str = ""
    scaling_factor: float = 1.0
    lattice: Lattice
    species: list[Species] = []
    selective_dynamics: bool = True
    coordinate_mode: str = "Direct"
    atoms: None
    atom_iterations = []
    written_to: str = None

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

    def append_to_x(self, poscar: "POSCAR") -> None:
        pass

    def append_to_x_negative(self, poscar: "POSCAR") -> None:
        poscar.append_to_x(self)

    def append_to_y(self, poscar: "POSCAR") -> None:
        pass

    def append_to_y_negative(self, poscar: "POSCAR") -> None:
        poscar.append_to_y(self)

    def append_to_z(self, poscar: "POSCAR") -> None:
        pass

    def append_to_z_negative(self, poscar: "POSCAR") -> None:
        poscar.append_to_z(self)

    def expand_to_super_cell(self, x: int = 1, y:int = 1, z: int = 1) -> None:
        if x == 1 and y == 1 and z == 1:
            return
        if x < 1 or y < 1 or z < 1:
            raise ValueError("The super cell dimensions must be greater than 0.")

        # self.comment += f"-expanded_by_{x}x{y}x{z}"

        # expand lattice
        self.lattice.a1 = [x * a for a in self.lattice.a1]
        self.lattice.a2 = [y * a for a in self.lattice.a2]
        self.lattice.a3 = [z * a for a in self.lattice.a3]

        # expand species
        result_species = []
        i = 0
        while i < x:
            j = 0
            while j < y:
                k = 0
                while k < z:
                    for l, species in enumerate(self.species):
                        if (len(result_species) < l + 1):
                            result_species.append(
                                Species(
                                    name=species.name,
                                    ion_positions=[]
                                )
                            )
                        for ion_position in species.ion_positions:
                            old_x, old_y, old_z = ion_position.x.value, ion_position.y.value, ion_position.z.value
                            new_x, new_y, new_z = old_x + i, old_y + j, old_z + k
                            result_species[l].ion_positions.append(
                                IonPosition(
                                    x=Position(
                                        value=new_x
                                    ),
                                    y=Position(
                                        value=new_y
                                    ),
                                    z=Position(
                                        value=new_z
                                    )
                                )
                            )
                    k += 1
                j += 1
            i += 1
        for species in result_species:
            for ion_position in species.ion_positions:
                ion_position.x.value = ion_position.x.value / x
                ion_position.y.value = ion_position.y.value / y
                ion_position.z.value = ion_position.z.value / z
        self.species = result_species

        # move layers of species to fit the new lattice

    def randomise_species(self, species: list[str] = None, weights: list[int] = None) -> None:
        result_species = []
        # self.comment += "-randomised"
        if species is None:
            species = [s.name for s in self.species]
        if weights is None:
            weights = [1 for _ in species]
        if len(weights) != len(species):
            raise ValueError("The number of species and weights must be the same.")
        weighted_species = []
        for i, s in enumerate(species):
            for _ in range(weights[i]):
                weighted_species.append(s)

        self.comment += "_to_" + "_".join(species)
        for s in species:
            result_species.append(
                Species(
                    name=s,
                    ion_positions=[]
                )
            )
        for s in self.species:
            for i, ion in enumerate(s.ion_positions):
                # randomly assign a new species
                result_species[species.index(weighted_species[randint(0, len(weighted_species) - 1)])].ion_positions.append(
                    IonPosition(
                        x=ion.x,
                        y=ion.y,
                        z=ion.z
                    )
                )
        self.species = result_species

    def write(self, filename: str, dirname: str = None) -> str:
        if dirname is None:
            dirname = "./structures/" + Helper.slugify(self.comment) + "/"
        filename = dirname + filename

        # Create the directory if it does not exist.
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        # Write the generated structure to a VASP POSCAR file.
        file = open(filename, "w")
        file.write(str(self))
        file.close()

        self.written_to = filename

        return filename

    def load_into_ace(self):
        temp = False
        if self.written_to is None:
            temp = True
            dirname = "./temp/"
            filename = dirname + "POSCAR"
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            file = open(filename, "w")
            file.write(str(self))
            file.close()
            self.written_to = filename
        self.atoms = read(self.written_to)
        self.atom_iterations.append(deepcopy(self.atoms))
        if temp:
            os.remove(self.written_to)
            self.written_to = None
        return self.atoms


    def relax(self, write: bool = False) -> None:
        # Set the periodic boundary conditions to True.
        self.atoms.set_pbc(True)

        # Initialize the MACE machine learning potential.
        macemp = mace_mp(model="large", dispersion=True, default_dtype="float64")

        # Assign the MACE potential to the atoms object for subsequent calculations.
        self.atoms.calc = macemp

        # Set up and run the optimization.
        optimizer = FIRE(self.atoms)

        # Run the optimization until the maximum force on each atom is less than 0.005 eV/Å.
        optimizer.run(fmax=0.005)

        if write:
            dirname = "./structures/" + Helper.slugify(self.comment) + "/"
            filename = dirname + "relaxed.vasp"
            # Write the relaxed structure to a VASP POSCAR file.
            self.atoms.write(filename, format='vasp', direct=True)

        self.atom_iterations.append(deepcopy(self.atoms))

    def write_energy(self, filename: str) -> None:
        dirname = "./structures/" + Helper.slugify(self.comment) + "/"
        filename = dirname + "energy.txt"

        # Write potential energy to file
        file = open(filename, "w")
        file.write("total: " + str(self.atoms.get_potential_energy()) + "\nper atom: " + str(self.atoms.get_potential_energy() / len(self.atoms)))
        file.close()

    def image(self, filename: str = None) -> None:
        dirname = "./structures/" + Helper.slugify(self.comment) + "/"
        if filename is None:
            filename = "image.png"
        filename = dirname + filename
        write(filename, self.atoms, rotation="45x,45y,45z", scale=150)

    def view(self) -> None:
        view(self.atoms)

    def from_file(self, filename: str) -> None:
        pass

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
            if (len(s.ion_positions) > 0):
                poscar += s.name + " "
        poscar += "\n"
        for s in self.species:
            if (len(s.ion_positions) > 0):
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