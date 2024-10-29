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
    "As": POSCAR.from_file(filename="mp-158.poscar", dirname=f"./structures/"),
    "Se": POSCAR.from_file(filename="mp-570481.poscar", dirname=f"./structures/"),
    "Br": POSCAR.from_file(filename="mp-998864.poscar", dirname=f"./structures/"),
    "Kr": POSCAR.from_file(filename="mp-975590.poscar", dirname=f"./structures/"),
    "Rb": POSCAR.from_file(filename="mp-1179656.poscar", dirname=f"./structures/"),
    "Sr": POSCAR.from_file(filename="mp-139.poscar", dirname=f"./structures/"),
    "Y": POSCAR.from_file(filename="mp-112.poscar", dirname=f"./structures/"),
    "Zr": POSCAR.from_file(filename="mp-131.poscar", dirname=f"./structures/"),
    "Nb": POSCAR.from_file(filename="mp-75.poscar", dirname=f"./structures/"),
    "Mo": POSCAR.from_file(filename="mp-129.poscar", dirname=f"./structures/"),
    "Tc": POSCAR.from_file(filename="mp-113.poscar", dirname=f"./structures/"),
    "Ru": POSCAR.from_file(filename="mp-33.poscar", dirname=f"./structures/"),
    "Rh": POSCAR.from_file(filename="mp-74.poscar", dirname=f"./structures/"),
    "Pd": POSCAR.from_file(filename="mp-2.poscar", dirname=f"./structures/"),
    "Ag": POSCAR.from_file(filename="mp-8566.poscar", dirname=f"./structures/"),
    "Cd": POSCAR.from_file(filename="mp-94.poscar", dirname=f"./structures/"),
    "In": POSCAR.from_file(filename="mp-85.poscar", dirname=f"./structures/"),
    "Sn": POSCAR.from_file(filename="mp-117.poscar", dirname=f"./structures/"),
    "Sb": POSCAR.from_file(filename="mp-104.poscar", dirname=f"./structures/"),
    "Te": POSCAR.from_file(filename="mp-19.poscar", dirname=f"./structures/"),
    "I": POSCAR.from_file(filename="mp-639751.poscar", dirname=f"./structures/"),
    "Xe": POSCAR.from_file(filename="mp-972256.poscar", dirname=f"./structures/"),
    "Cs": POSCAR.from_file(filename="mp-1055940.poscar", dirname=f"./structures/"),
    "Ba": POSCAR.from_file(filename="mp-122.poscar", dirname=f"./structures/"),
    # "La-Lu": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    "Hf": POSCAR.from_file(filename="mp-103.poscar", dirname=f"./structures/"),
    "Ta": POSCAR.from_file(filename="mp-569794.poscar", dirname=f"./structures/"),
    "W": POSCAR.from_file(filename="mp-91.poscar", dirname=f"./structures/"),
    "Re": POSCAR.from_file(filename="mp-1186901.poscar", dirname=f"./structures/"),
    "Os": POSCAR.from_file(filename="mp-49.poscar", dirname=f"./structures/"),
    "Ir": POSCAR.from_file(filename="mp-101.poscar", dirname=f"./structures/"),
    "Pt": POSCAR.from_file(filename="mp-126.poscar", dirname=f"./structures/"),
    "Au": POSCAR.from_file(filename="mp-81.poscar", dirname=f"./structures/"),
    "Hg": POSCAR.from_file(filename="mp-1017981.poscar", dirname=f"./structures/"),
    "Tl": POSCAR.from_file(filename="mp-82.poscar", dirname=f"./structures/"),
    "Pb": POSCAR.from_file(filename="mp-20483.poscar", dirname=f"./structures/"),
    "Bi": POSCAR.from_file(filename="mp-567597.poscar", dirname=f"./structures/"),
    # "Po": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "At": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Rn": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Fr": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Ra": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Ac-Lr": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Rf": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Db": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Sg": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Bh": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Hs": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Mt": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Ds": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Rg": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Cn": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Nh": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Fl": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Mc": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Lv": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Ts": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Og": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    "La": POSCAR.from_file(filename="mp-26.poscar", dirname=f"./structures/"),
    "Ce": POSCAR.from_file(filename="mp-28.poscar", dirname=f"./structures/"),
    "Pr": POSCAR.from_file(filename="mp-97.poscar", dirname=f"./structures/"),
    "Nd": POSCAR.from_file(filename="mp-123.poscar", dirname=f"./structures/"),
    "Pm": POSCAR.from_file(filename="mp-867200.poscar", dirname=f"./structures/"),
    "Sm": POSCAR.from_file(filename="mp-69.poscar", dirname=f"./structures/"),
    "Eu": POSCAR.from_file(filename="mp-21462.poscar", dirname=f"./structures/"),
    "Gd": POSCAR.from_file(filename="mp-155.poscar", dirname=f"./structures/"),
    "Tb": POSCAR.from_file(filename="mp-18.poscar", dirname=f"./structures/"),
    "Dy": POSCAR.from_file(filename="mp-88.poscar", dirname=f"./structures/"),
    "Ho": POSCAR.from_file(filename="mp-144.poscar", dirname=f"./structures/"),
    "Er": POSCAR.from_file(filename="mp-1184115.poscar", dirname=f"./structures/"),
    "Tm": POSCAR.from_file(filename="mp-143.poscar", dirname=f"./structures/"),
    # "Yb": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    "Lu": POSCAR.from_file(filename="mp-973571.poscar", dirname=f"./structures/"),
    "Ac": POSCAR.from_file(filename="mp-862690.poscar", dirname=f"./structures/"),
    "Th": POSCAR.from_file(filename="mp-37.poscar", dirname=f"./structures/"),
    "Pa": POSCAR.from_file(filename="mp-62.poscar", dirname=f"./structures/"),
    "U": POSCAR.from_file(filename="mp-44.poscar", dirname=f"./structures/"),
    "Np": POSCAR.from_file(filename="mp-11534.poscar", dirname=f"./structures/"),
    "Pu": POSCAR.from_file(filename="mp-582819.poscar", dirname=f"./structures/"),
    # "Am": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Cm": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Bk": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Cf": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Es": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Fm": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Md": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "No": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
    # "Lr": POSCAR.from_file(filename=".poscar", dirname=f"./structures/"),
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