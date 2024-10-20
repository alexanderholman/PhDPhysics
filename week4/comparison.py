# ==================== setup start ===================
# ===================== setup end ====================

# imports
from copy import deepcopy

#import classes
from classes.POSCAR import POSCAR
from classes.Species import Species
from classes.IonPosition import IonPosition
from classes.Position import Position
from classes.Lattice import Lattice

# load initial structure
# Structures from Materials Project
# source https://next-gen.materialsproject.org/materials?_limit=75&chemsys=Si
mp_149 = {
    "comment": "mp-149",
    "lattice": Lattice( # trimmed 0.0...03 and 0.0...09 to 0.0
        a1=[5.4437023729394527, 0.0, 0.0],
        a2=[0.0, 5.4437023729394527, 0.0],
        a3=[0.0, 0.0, 5.4437023729394527]
    ),
    "species": [
        {
            "name": "Si",
            "positions": [ # rearranged to order [0,0,0] out to [1,1,1] in their pairs
                [0.0, 0.0, 0.0],
                [0.25, 0.25, 0.25],
                [0.0, 0.5, 0.5],
                [0.25, 0.75, 0.75],
                [0.5, 0.0, 0.5],
                [0.75, 0.25, 0.75],
                [0.5, 0.5, 0.0],
                [0.75, 0.75, 0.25],
            ]
        }
    ]
}

# load into POSCAR class
mp_149_poscar = POSCAR(
        comment=mp_149['comment'],
        lattice=mp_149['lattice'],
        species=list(map(
            lambda sp: Species(
                name=sp["name"],
                ion_positions=list(map(
                    lambda p: IonPosition(Position(p[0]), Position(p[1]), Position(p[2])),
                    sp["positions"]
                ))
            ),
            mp_149["species"])
        )
    )

# species and wights
species = ["Si", "Ge"]
weights = [1, 1]

# ==================== setup end ====================
# ---------------------------------------------------
# =================== build start ===================

# build comparisons
comparisons = {
    # 'split_no_alloy': [],
    # 'split_with_alloy': [], # not yet implemented todo implement
    # 'random_pre_expansion': [],
    'random_per_expansion': [], # not needed perhaps, energies will be random
}

# randomise before expansion
random_pre_expansion = deepcopy(mp_149_poscar)
random_pre_expansion.randomise_species(
    species=species,
    weights=weights,
)

while len(random_pre_expansion.species[0].ion_positions) is not len(random_pre_expansion.species[1].ion_positions):
    random_pre_expansion.randomise_species(
        species=species,
        weights=weights,
    )

# for each size of expansion to a super-cell 1 -> n
N = 5

for n in range(1, N + 1):
    # build split without alloy interface
    if "split_no_alloy" in comparisons:
        split_no_alloy = deepcopy(mp_149_poscar)
        split_no_alloy.expand_to_super_cell(
            x=n,
            y=n,
            z=n,
        )
        split_no_alloy.split_species_x(
            species=species,
            weights=weights,
        )
        comparisons['split_no_alloy'].append(split_no_alloy)

    # todo build split with alloy interface
    if "split_with_alloy" in comparisons:
        pass

    # build randomised at start and then expanded, force perfect ratio
    if "random_pre_expansion" in comparisons:
        random_pre_expansion_copy = deepcopy(random_pre_expansion)
        random_pre_expansion_copy.expand_to_super_cell(
            x=n,
            y=n,
            z=n,
        )
        comparisons['random_pre_expansion'].append(random_pre_expansion_copy)

    # build randomised at each expansion, force perfect radio
    if "random_per_expansion" in comparisons:
        random_per_expansion_copy = deepcopy(mp_149_poscar)
        random_per_expansion_copy.expand_to_super_cell(
            x=n,
            y=n,
            z=n,
        )
        random_per_expansion_copy.randomise_species(
            species=species,
            weights=weights,
        )
        while len(random_per_expansion_copy.species[0].ion_positions) is not len(random_per_expansion_copy.species[1].ion_positions):
            random_per_expansion_copy.randomise_species(
                species=species,
                weights=weights,
            )

# ==================== build end ====================

# for each comparison
for type, poscars in comparisons.items():
    dirname = f"./structures/comparison/{type}/"
    # for each build
    for i, poscar in enumerate(poscars):
        n = i + 1
        # generate initial VASP, and image
        poscar.write(
            filename=f"generated-{n}.vasp",
            dirname=dirname
        )
        poscar.load_into_ace()
        poscar.image(
            filename=f"generated-{n}.png",
            dirname=dirname
        )

        # generate relaxed VASP, energies, and image
        poscar.relax(
            write=True,
            filename=f"relaxed-{n}.vasp",
            dirname=dirname
        )
        poscar.write_energy(
            filename=f"energies-{n}.txt",
            dirname=dirname
        )
        poscar.image(
            filename=f"relaxed-{n}.png",
            dirname=dirname
        )

        # view
        poscar.view()

# todo: generate graph comparing energies of each comparison and how the energies change or not as the structure is expanded into super-cells