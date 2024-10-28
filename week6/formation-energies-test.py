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
    "Ga": POSCAR.from_file(filename="mp-142.poscar", dirname=f"./structures/"),
}
# expected_formation_energies = {
#     "Si7Ge": 0.003, # calculated -0.006545196732207437
#     "SiGe": 0.053, # calculated -0.0026943691455496577
#     "Ga3Si": 0.210, # calculated 0.18105908160402961
#     "mp-1095269": 0.111, #calculated 0.11059346380659463
# }

# ==================== setup end ====================
# ---------------------------------------------------
# =================== build start ===================


# load structure
poscar = POSCAR.from_file(filename="mp-149.poscar", dirname=f"./structures/")
poscar.expand_to_super_cell(
    x=3,
    y=3,
    z=3
)
poscar.expand_to_super_cell(x=2, y=1, z=1)
poscar.split_species_x(['Si', 'Ge'], [1, 1])
# ==================== build end ====================

dirname = f"./structures/formation-energies-test/"
filename = f"2-1-1"
poscar.load_into_ace()
os.makedirs(os.path.dirname(dirname+"loaded/vasp/"), exist_ok=True)
poscar.write(
    filename=f"{filename}.vasp",
    dirname=dirname+"loaded/vasp/"
)
os.makedirs(os.path.dirname(dirname+"loaded/png/"), exist_ok=True)
poscar.image(
    filename=f"{filename}.png",
    dirname=dirname+"loaded/png/"
)
os.makedirs(os.path.dirname(dirname+"relaxed/vasp/"), exist_ok=True)
os.makedirs(os.path.dirname(dirname+"relaxed/png/"), exist_ok=True)
formation_energy_per_atom = poscar.calculate_formation_energies(
    bulk_structures=bulk_structures,
    relax_write=True,
    relax_filename=f"{filename}.vasp",
    relax_dirname=f"{dirname}relaxed/vasp/",
    image=True,
    image_filename=f"{filename}.png",
    image_dirname=f"{dirname}relaxed/png/"
)
print(f"formation energies per atom: {formation_energy_per_atom}")