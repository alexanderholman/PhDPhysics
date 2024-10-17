# Import the classes.
from copy import deepcopy

from classes.POSCAR import POSCAR
from classes.Species import Species
from classes.IonPosition import IonPosition
from classes.Position import Position
from classes.Lattice import Lattice

st = {
    "comment": "si-2",
    "lattice": Lattice(
        a1=[2.715, 0.0, 0.0],
        a2=[0.0, 2.715, 0.0],
        a3=[0.0, 0.0, 2.715]
    ),
    "species": [
        {
            "name": "Si",
            "positions": [
                [0.0, 0.0, 0.0],
                [1/4, 1/4, 1/4],
            ]
        }
    ]
}

# Process structures, generate files, and relax structures.
comment = st["comment"] + "-half-mp-149-lattice"

poscar = POSCAR(
    comment=comment,
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

# if lattice is provided, use it
if "lattice" in st:
    poscar.lattice = st["lattice"]

original_poscar = deepcopy(poscar)

# Write the generated structure to a VASP POSCAR file.
poscar.write(f"generated.vasp")

# Load randomised structure into ASE.
poscar.load_into_ace()

# Load original structure
poscar.species = deepcopy(original_poscar.species)
poscar.lattice = deepcopy(original_poscar.lattice)
poscar.atoms = poscar.atom_iterations[0]

# define number of times to expand the structure to a super-cell
n = 2

# Expand the structure to a super-cell.
poscar.expand_to_super_cell(
    x=n,
    y=n,
    z=n
)

# Write the expanded structure to a VASP POSCAR file.
poscar.write(f"expanded-{n}.vasp")

# Read the generated structure from a VASP POSCAR file.
poscar.load_into_ace()

# Take snapshot of the structure, pre-relaxation.
poscar.image(f"pre-relaxation-expanded-{n}.png")

# Relax the structure.
poscar.relax(
    write=True,
    filename=f"relaxed-expanded-{n}.vasp"
)

# Take snapshot of the structure, post-relaxation.
poscar.image(f"post-relaxation-expanded-{n}.png")

# Write the energy of the relaxed structure.
poscar.write_energy(filename=f"energy-expanded-{n}.txt")

# view all iterations of the structure
for i, atoms in enumerate(poscar.atom_iterations):
    poscar.atoms = atoms
    poscar.view()