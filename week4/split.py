# Import the classes.
from copy import deepcopy

from classes.POSCAR import POSCAR
from classes.Species import Species
from classes.IonPosition import IonPosition
from classes.Position import Position

# Import the list of structures to process.
from structures_list import *

# Split the species of the structure.
split = [
    ["Si", 50],
    ["Ge", 50]
]

# Number of times to expand the structure to a super-cell.
N = 1

# Process structures, generate files, and relax structures.
for st in structures:
    comment = st["comment"] + f"-split"

    if len(split):
        # split randomise into species and weights
        split_species = []
        split_weights = []
        for r in split:
            comment += f"-{r[0]}-{r[1]}"
            split_species.append(r[0])
            split_weights.append(r[1])

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

    # Expand the structure to a super-cell NxNxN.
    poscar.expand_to_super_cell(
        x=N,
        y=N,
        z=N
    );
    poscar.write(f"generated-expnded-{N}.vasp")

    if len(split):
        # Randomise the species of the structure.
        poscar.split_species_x(
            species=split_species,
            weights=split_weights
        )

        # Write the randomised structure to a VASP POSCAR file.
        poscar.write(f"generated-expnded-{N}-split.vasp")

    original_poscar = deepcopy(poscar)

    # Load randomised structure into ASE.
    poscar.load_into_ace()

    # Take snapshot of the structure, pre-relaxation.
    poscar.image(f"pre-relaxation-expnded-{N}-split.png")

    # Relax the structure.
    poscar.relax(
        write=True,
        filename=f"relaxed-expnded-{N}-split.vasp"
    )

    # Take snapshot of the structure, post-relaxation.
    poscar.image(f"post-relaxation-expnded-{N}-split.png")

    # Write the energy of the relaxed structure.
    poscar.write_energy(filename=f"energy-expnded-{N}-split.txt")

    # view all iterations of the structure
    for i, atoms in enumerate(poscar.atom_iterations):
        poscar.atoms = atoms
        poscar.view()