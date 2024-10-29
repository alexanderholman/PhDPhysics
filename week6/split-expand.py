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
structures_list = []
bulk_structures: dict[str, POSCAR] = {
    "Si": POSCAR.from_file(filename="mp-149.poscar", dirname=f"./structures/"),
    "Ge": POSCAR.from_file(filename="mp-32.poscar", dirname=f"./structures/"),
}

# ==================== setup end ====================
# ---------------------------------------------------
# =================== build start ===================


# load structure
poscar = POSCAR.from_file(filename=f"generated-5.vasp", dirname=f"./structures/comparison/split_no_alloy/")
for i in range(1, 4):
    new_poscar = deepcopy(poscar)
    new_poscar.expand_to_super_cell(x=i)
    new_poscar.split_species_x(['Si', 'Ge'], [1, 1])
    structures_list.append(new_poscar)
del poscar
# ==================== build end ====================

# open new csv file
csv_filename = f"./structures/split-expand/energies.csv"
os.makedirs(os.path.dirname(csv_filename), exist_ok=True)
csv = open(csv_filename, "w")
csv.write("multiplier_x,multiplier_y,multiplier_z,energies\n")

# for each stretched_structure
for i, stretched_structure in enumerate(structures_list):
    n = i + 1
    dirname = f"./structures/split-expand/"
    filename_postfix = f"{n}-1-1"
    poscar = stretched_structure
    poscar.load_into_ace()
    os.makedirs(os.path.dirname(dirname+"loaded/vasp/"), exist_ok=True)
    poscar.write(
        filename=f"{filename_postfix}.vasp",
        dirname=dirname+"loaded/vasp/"
    )
    os.makedirs(os.path.dirname(dirname+"loaded/png/"), exist_ok=True)
    poscar.image(
        filename=f"{filename_postfix}.png",
        dirname=dirname+"loaded/png/"
    )
    os.makedirs(os.path.dirname(dirname+"relaxed/vasp/"), exist_ok=True)
    os.makedirs(os.path.dirname(dirname+"relaxed/png/"), exist_ok=True)
    formation_energy_per_atom = poscar.calculate_formation_energies(
        bulk_structures=bulk_structures,
        relax_write=True,
        relax_filename=f"{filename_postfix}.vasp",
        relax_dirname=dirname+"relaxed/vasp/",
        image=True,
        image_filename=f"{filename_postfix}.png",
        image_dirname=dirname+"relaxed/png/"
    )
    csv.write(f"{n},1,1,{formation_energy_per_atom}\n")
    del stretched_structure
    structures_list[i] = None

csv.close()