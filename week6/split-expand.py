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

structures_list = []
structure_names = []

# ==================== setup end ====================
# ---------------------------------------------------
# =================== build start ===================


# load structure
for i in range(1, 2):
    poscar = POSCAR.from_file(filename=f"generated-{i}.vasp", dirname=f"./structures/comparison/split_no_alloy/")
    for j in range(1, 4):
        new_poscar = deepcopy(poscar)
        new_poscar.expand_to_super_cell(x=1, y=j, z=j)
        new_poscar.split_species_x(['Si', 'Ge'], [1, 1])
        structures_list.append(new_poscar)
        structure_names.append(f"generated-{i}-1-{j}-{j}")
    del poscar
# ==================== build end ====================

# open new csv file
csv_filename = f"./structures/split-expand/energies.csv"
os.makedirs(os.path.dirname(csv_filename), exist_ok=True)
csv = open(csv_filename, "w")
csv.write("multiplier_x,multiplier_y,multiplier_z,energy_per_atom,formation_energy\n")

# for each stretched_structure
for i, stretched_structure in enumerate(structures_list):
    n = i + 1
    dirname = f"./structures/split-expand/"
    filename_postfix = structure_names[i]
    poscar = stretched_structure
    poscar.load_into_ace()
    os.makedirs(os.path.dirname(dirname+"loaded/vasp/"), exist_ok=True)
    os.makedirs(os.path.dirname(dirname+"loaded/png/"), exist_ok=True)
    os.makedirs(os.path.dirname(dirname+"relaxed/vasp/"), exist_ok=True)
    os.makedirs(os.path.dirname(dirname+"relaxed/png/"), exist_ok=True)
    poscar.write(
        filename=f"{filename_postfix}.vasp",
        dirname=dirname+"loaded/vasp/"
    )
    poscar.image(
        filename=f"{filename_postfix}.png",
        dirname=dirname+"loaded/png/"
    )
    formation_energy_per_atom = poscar.formation_energies(
        fmax=0.1,
        relax_write=True,
        relax_filename=f"{filename_postfix}.vasp",
        relax_dirname=dirname+"relaxed/vasp/"
    )
    poscar.image(
        filename=f"{filename_postfix}.png",
        dirname=dirname+"relaxed/png/"
    )
    csv.write(f"{structure_names[i]},1,1,{poscar.energy_per_atom()}, {formation_energy_per_atom}\n")
    del stretched_structure
    structures_list[i] = None

csv.close()