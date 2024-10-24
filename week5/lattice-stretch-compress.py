# ==================== setup start ===================
# ===================== setup end ====================

# imports
from copy import deepcopy
from numpy import arange
import os

#import classes
from classes.POSCAR import POSCAR
from classes.Species import Species
from classes.IonPosition import IonPosition
from classes.Position import Position
from classes.Lattice import Lattice

# load initial structure
# load split_no_alloy
# load random_post_expansion

stretch = 1.1
compress = 0.9
increments = 0.1
structures_list = []
hull_energies = {}

class StretchedStructure:
    structure: POSCAR
    group: str
    original: Lattice
    stretched: Lattice
    diff = Lattice
    energies_above_hull: float

    def __init__(self, structure: POSCAR, group: str, original: Lattice, stretched: Lattice) -> None:
        self.structure = structure
        self.group = group
        self.original = original
        self.stretched = stretched
        self.diff = Lattice(
            a1=[i - j for i, j in zip(stretched.a1, original.a1)],
            a2=[i - j for i, j in zip(stretched.a2, original.a2)],
            a3=[i - j for i, j in zip(stretched.a3, original.a3)]
        )

# ==================== setup end ====================
# ---------------------------------------------------
# =================== build start ===================

# for each split_no_alloy, random_post_expansion
for comparison in ['split_no_alloy', 'random_post_expansion']:
    # load structure
    poscar = POSCAR.from_file(filename="generated-1.vasp", dirname=f"./structures/comparison/{comparison}/")

    # stretch/compress x lattice
    for i in arange(compress, stretch, increments):
        for j in arange(compress, stretch, increments):
            for k in arange(compress, stretch, increments):
                new_lattice = Lattice(
                    a1=[i * a for a in poscar.lattice.a1],
                    a2=[j * a for a in poscar.lattice.a2],
                    a3=[k * a for a in poscar.lattice.a3]
                )
                stretch_structure = StretchedStructure(
                    structure=deepcopy(poscar),
                    group=comparison,
                    original=poscar.lattice,
                    stretched=new_lattice
                )
                stretch_structure.structure.lattice = new_lattice
                structures_list.append(stretch_structure)

    poscar.load_into_ace()
    poscar.relax()
    hull_energies[comparison] = poscar.atoms.get_potential_energy() / len(poscar.atoms)
    del poscar

# ==================== build end ====================

# open new csv file
csv_filename = "./structures/lattice-stretch-compress/energies.csv"
os.makedirs(os.path.dirname(csv_filename), exist_ok=True)
csv = open(csv_filename, "w")
csv.write("group, a1, a2, a3, energy_above_hull\n")

print(f"relaxing stretched structures: {len(structures_list)}")

# for each stretched_structure
for i, stretched_structure in enumerate(structures_list):
    n = i + 1
    print(f"relaxing structure {n}/{len(structures_list)}")
    dirname = f"./structures/lattice-stretch-compress/{stretched_structure.group}/"
    filename_postfix = f"{stretched_structure.diff.a1[0]}-{stretched_structure.diff.a2[1]}-{stretched_structure.diff.a3[2]}"
    poscar = stretched_structure.structure
    poscar.write(
        filename=f"loaded-{filename_postfix}.vasp",
        dirname=dirname
    )
    poscar.load_into_ace()
    poscar.image(
        filename=f"loaded-{filename_postfix}.png",
        dirname=dirname
    )
    poscar.relax(
        write=True,
        filename=f"relaxed-{filename_postfix}.vasp",
        dirname=dirname
    )
    poscar.image(
        filename=f"relaxed-{filename_postfix}.png",
        dirname=dirname
    )
    stretched_structure.energies_above_hull = poscar.atoms.get_potential_energy() / len(poscar.atoms) - hull_energies[stretched_structure.group]
    csv.write(f"{stretched_structure.group}, {stretched_structure.diff.a1[0]}, {stretched_structure.diff.a2[1]}, {stretched_structure.diff.a3[2]}, {stretched_structure.energies_above_hull}\n")

csv.close()