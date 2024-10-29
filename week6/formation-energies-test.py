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
bulk_structures: dict[str, POSCAR] = {
    "H": POSCAR.from_file(filename="mp-730101.poscar", dirname=f"./structures/"),
    "He": POSCAR.from_file(filename="mp-23158.poscar", dirname=f"./structures/"),
    "Li": POSCAR.from_file(filename="mp-1018134.poscar", dirname=f"./structures/"),
    "Be": POSCAR.from_file(filename="mp-87.poscar", dirname=f"./structures/"),
    "B": POSCAR.from_file(filename="mp-160.poscar", dirname=f"./structures/"),
    "C": POSCAR.from_file(filename="mp-2516584.poscar", dirname=f"./structures/"),
    "N": POSCAR.from_file(filename="mp-154.poscar", dirname=f"./structures/"),
    "O": POSCAR.from_file(filename="mp-12957.poscar", dirname=f"./structures/"),
    "F": POSCAR.from_file(filename="mp-561203.poscar", dirname=f"./structures/"),
    "Ne": POSCAR.from_file(filename="mp-111.poscar", dirname=f"./structures/"),
    "Na": POSCAR.from_file(filename="mp-10172.poscar", dirname=f"./structures/"),
    "Mg": POSCAR.from_file(filename="mp-153.poscar", dirname=f"./structures/"),
    "Al": POSCAR.from_file(filename="mp-134.poscar", dirname=f"./structures/"),
    "Si": POSCAR.from_file(filename="mp-149.poscar", dirname=f"./structures/"),
    "P": POSCAR.from_file(filename="mp-568348.poscar", dirname=f"./structures/"),
    "S": POSCAR.from_file(filename="mp-77.poscar", dirname=f"./structures/"),
    "Cl": POSCAR.from_file(filename="mp-22848.poscar", dirname=f"./structures/"),
    "Ar": POSCAR.from_file(filename="mp-23155.poscar", dirname=f"./structures/"),
    "K": POSCAR.from_file(filename="mp-1184804.poscar", dirname=f"./structures/"),
    "Ca": POSCAR.from_file(filename="mp-45.poscar", dirname=f"./structures/"),
    "Sc": POSCAR.from_file(filename="mp-67.poscar", dirname=f"./structures/"),
    "Ti": POSCAR.from_file(filename="mp-72.poscar", dirname=f"./structures/"),
    "V": POSCAR.from_file(filename="mp-146.poscar", dirname=f"./structures/"),
    "Cr": POSCAR.from_file(filename="mp-90.poscar", dirname=f"./structures/"),
    "Mn": POSCAR.from_file(filename="mp-35.poscar", dirname=f"./structures/"),
    "Fe": POSCAR.from_file(filename="mp-13.poscar", dirname=f"./structures/"),
    "Co": POSCAR.from_file(filename="mp-102.poscar", dirname=f"./structures/"),
    "Ni": POSCAR.from_file(filename="mp-23.poscar", dirname=f"./structures/"),
    "Cu": POSCAR.from_file(filename="mp-30.poscar", dirname=f"./structures/"),
    "Zn": POSCAR.from_file(filename="mp-79.poscar", dirname=f"./structures/"),
    "Ga": POSCAR.from_file(filename="mp-142.poscar", dirname=f"./structures/"),
    "Ge": POSCAR.from_file(filename="mp-32.poscar", dirname=f"./structures/"),
}
structures_list = []
expected = {
    "mp-1094056": {
        "energies": {
            "above_hull": 0.003,
            "predicted_formation": 0.003 # -0.006545196732207437
        },
        "experimentally_observed": False
    },
    "mp-1096549": {
        "energies": {
            "above_hull": 0.053,
            "predicted_formation": 0.053 # -0.0026943691455496577
        },
        "experimentally_observed": False
    },
    "mp-1183982": {
        "energies": {
            "above_hull": 0.210,
            "predicted_formation": 0.210 # 0.18105908160402961
        },
        "experimentally_observed": False
    },
    "mp-1095269": {
        "energies": {
            "above_hull": 0.111,
            "predicted_formation": 0.111 # 0.11059346380659463
        },
        "experimentally_observed": True
    },
    "mp-569128": {
        "energies": {
            "above_hull": 0.0,
            "predicted_formation": -0.044 #calculated -0.1742738611235124
        },
        "experimentally_observed": True
    },
    "mp-7700": {
        "energies": {
            "above_hull": 0.479,
            "predicted_formation": 0.454 #calculated 0.38623465313062055
        },
        "experimentally_observed": True
    },
}

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