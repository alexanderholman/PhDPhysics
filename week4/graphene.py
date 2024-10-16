# Import the classes.
from copy import deepcopy

from classes.POSCAR import POSCAR
from classes.Species import Species
from classes.IonPosition import IonPosition
from classes.Position import Position
from classes.Lattice import Lattice

# Number of times to expand the structure to a super-cell.
N = 5

st = {
    "comment": "graphene",
    "lattice": Lattice(
        a1=[2.4559495449 , 0.0, 0.0],
        a2=[-1.2279240765, 2.1269439648, 0.0],
        a3=[0.0, 0.0, 10.0]
    ),
    "species": [
        {
            "name": "C",
            "positions": [
                [0.0, 0.0, 0.0],
                [2/3, 1/3, 0.0],
            ]
        }
    ]
}

# Process structures, generate files, and relax structures.
comment = st["comment"] # + f"-super-cell-{N}"

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

# Expand the structure to a super-cell N times
for n in range(1, N + 1):
    # Load original structure
    poscar.species = deepcopy(original_poscar.species)
    poscar.lattice = deepcopy(original_poscar.lattice)
    poscar.atoms = deepcopy(poscar.atom_iterations[0])

    # Expand the structure to a super-cell.
    poscar.expand_to_super_cell(
        x=n,
        y=n,
        z=n
    )

    # Write the expanded structure to a VASP POSCAR file.
    poscar.write(f"generated-expanded-{n}.vasp")

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