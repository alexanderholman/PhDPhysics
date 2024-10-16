# Import the classes.
from classes.POSCAR import POSCAR
from classes.Species import Species
from classes.IonPosition import IonPosition
from classes.Position import Position

# Import the list of structures to process.
from structures_list import *

# Number of times to expand the structure to a super-cell.
N = 3

randomise = [
    ["Si", 50],
    ["Ge", 50]
]

# Process structures, generate files, and relax structures.
for st in structures:
    comment = st["comment"] + f"-super-cell-{N}"

    if len(randomise):
        comment += "-randomised"
        # split randomise into species and weights
        random_species = []
        random_weights = []
        for r in randomise:
            comment += f"-{r[0]}-{r[1]}"
            random_species.append(r[0])
            random_weights.append(r[1])

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

    # Write the generated structure to a VASP POSCAR file.
    poscar.write(f"generated.vasp")

    if len(randomise):
        # Randomise the species of the structure.
        poscar.randomise_species(
            species=random_species,
            weights=random_weights
        )

        # Write the randomised structure to a VASP POSCAR file.
        poscar.write(f"generated-random.vasp")

    # Load randomised structure into ASE.
    poscar.load_into_ace()

    # Expand the structure to a super-cell N times
    for n in range(1, N + 1):
        # Load original structure
        poscar.atoms = poscar.atom_iterations[0]

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