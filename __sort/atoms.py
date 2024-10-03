import os
import unicodedata
import re
from ase.io import read, write
from ase.optimize import FIRE
from ase.visualize import view
from mace.calculators import mace_mp

# define a class structures for VASP POSCAR file format
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
                a3=[0.0, 0.0, 5,43]
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

# Define a function to convert a string to a slug, found here: https://stackoverflow.com/questions/295135/turn-a-string-into-a-valid-filename
def slugify(value, allow_unicode=False):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')

from structures_list import *

# Process structures, generate files, and relax structures.
for st in structures:
    poscar = POSCAR(
        comment=st["comment"],
        species=list(map(
            lambda sp: Species(
                name=sp["name"], 
                ion_positions=list(map(
                    lambda p: IonPosition(Position(p[0]), Position(p[1]), Position(p[2])), 
                    sp["positions"]
                ))
            ), 
            st["species"])
        )

    )

    #if lattice is provided, use it
    if "lattice" in st:
        poscar.lattice = st["lattice"]

    # generate filenames
    dirname = "./structures/" + slugify(st["comment"]) + "/"
    generated_file = dirname + "generated.vasp"
    relaxed_file = dirname + "relaxed.vasp"
    image_file = dirname + "image.png"
    energy_file = dirname + "potential_energy.txt"

    # Create the directory if it does not exist.
    os.makedirs(os.path.dirname(generated_file), exist_ok=True)

    # Write the generated structure to a VASP POSCAR file.
    file = open(generated_file, "w")
    file.write(str(poscar))
    file.close()
    
    # Read the generated structure from a VASP POSCAR file.
    atoms = read(generated_file)

    # Set the periodic boundary conditions to True.
    atoms.set_pbc(True)

    # Initialize the MACE machine learning potential.
    macemp = mace_mp(model="large", dispersion=True, default_dtype="float64")

    # Assign the MACE potential to the atoms object for subsequent calculations.
    atoms.calc = macemp

    # Set up and run the optimization.
    optimizer = FIRE(atoms)

    # Run the optimization until the maximum force on each atom is less than 0.005 eV/Å.
    optimizer.run(fmax=0.005)

    # Write the relaxed structure to a VASP POSCAR file.
    atoms.write(relaxed_file, format='vasp', direct=True)

    # Write potential energy to file
    file = open(energy_file, "w")
    file.write("total: " + str(atoms.get_potential_energy()) + "\nper atom: " + str(atoms.get_potential_energy() / len(atoms)))
    file.close()

    # Write the image of the relaxed structure.
    write(image_file, atoms, rotation="45x,45y,45z", scale=150)

    # View the relaxed structure.
    view(atoms)