# ==================== setup start ===================
# ===================== setup end ====================

# imports
from copy import deepcopy
from datetime import datetime, date

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

stretch = [1.05, 1, 1]
compress = [0.95, 1, 1]
increments = [0.001, 0.01, 0.01]
structures_list = []
hull_energies = {}
timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

class StretchedStructure:
    structure: POSCAR
    group: str
    original: Lattice
    stretched: Lattice
    diff = Lattice
    energies_above_hull: float
    multipliers: list[float]

    def __init__(self, structure: POSCAR, group: str, original: Lattice, stretched: Lattice, multipliers: list[float]) -> None:
        self.structure = structure
        self.group = group
        self.original = original
        self.stretched = stretched
        self.diff = Lattice(
            a1=[i - j for i, j in zip(stretched.a1, original.a1)],
            a2=[i - j for i, j in zip(stretched.a2, original.a2)],
            a3=[i - j for i, j in zip(stretched.a3, original.a3)]
        )
        self.multipliers = multipliers

# ==================== setup end ====================
# ---------------------------------------------------
# =================== build start ===================

# for each split_no_alloy, random_post_expansion
for comparison in [
    'split_no_alloy',
    # 'random_post_expansion',
]:
    for n in range(3, 4):
        # load structure
        poscar = POSCAR.from_file(filename=f"generated-{n}.vasp", dirname=f"./structures/comparison/{comparison}/")

        # stretch/compress x lattice
        i_range = arange(compress[0], stretch[0], increments[0])
        if len(i_range) < 1:
            i_range = [1]
        for i in i_range:
            j_range = arange(compress[1], stretch[1], increments[1])
            if len(j_range) < 1:
                j_range = [1]
            for j in j_range:
                k_range = arange(compress[2], stretch[2], increments[2])
                if len(k_range) < 1:
                    k_range = [1]
                for k in k_range:
                    new_lattice = Lattice(
                        a1=[i * a for a in poscar.lattice.a1],
                        a2=[j * a for a in poscar.lattice.a2],
                        a3=[k * a for a in poscar.lattice.a3]
                    )
                    stretch_structure = StretchedStructure(
                        structure=deepcopy(poscar),
                        group=comparison,
                        original=poscar.lattice,
                        stretched=new_lattice,
                        multipliers=[i, j, k]
                    )
                    stretch_structure.structure.lattice = new_lattice
                    structures_list.append(stretch_structure)
    poscar.load_into_ace()
    poscar.relax()
    hull_energies[comparison] = poscar.atoms.get_potential_energy() / len(poscar.atoms)
    del poscar
# ==================== build end ====================

# open new csv file
csv_filename = f"./structures/lattice-stretch-compress/energies.csv"
os.makedirs(os.path.dirname(csv_filename), exist_ok=True)
csv = open(csv_filename, "w")
csv.write("group,multiplier_x,multiplier_y,multiplier_z,a1,a2,a3,energy_above_hull\n")

# for each stretched_structure
for i, stretched_structure in enumerate(structures_list):
    n = i + 1
    dirname = f"./structures/lattice-stretch-compress/{stretched_structure.group}/"
    filename_postfix = f"{stretched_structure.multipliers[0]}-{stretched_structure.multipliers[1]}-{stretched_structure.multipliers[2]}"
    poscar = stretched_structure.structure
    os.makedirs(os.path.dirname(dirname+"loaded/vasp/"), exist_ok=True)
    poscar.write(
        filename=f"loaded-{filename_postfix}.vasp",
        dirname=dirname+"loaded/vasp/"
    )
    poscar.load_into_ace()
    os.makedirs(os.path.dirname(dirname+"relaxed/vasp/"), exist_ok=True)
    poscar.relax(
        write=True,
        filename=f"relaxed-{filename_postfix}.vasp",
        dirname=dirname+"relaxed/vasp/"
    )
    stretched_structure.energies_above_hull = poscar.atoms.get_potential_energy() / len(poscar.atoms) - hull_energies[stretched_structure.group]
    csv.write(f"{stretched_structure.group},{stretched_structure.multipliers[0]},{stretched_structure.multipliers[1]},{stretched_structure.multipliers[2]},{stretched_structure.diff.a1[0]},{stretched_structure.diff.a2[1]},{stretched_structure.diff.a3[2]},{stretched_structure.energies_above_hull}\n")
    del stretched_structure
    structures_list[i] = None

csv.close()