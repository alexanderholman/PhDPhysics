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
bulk_structures: dict[str, str] = {
    "H":  "mp-730101.poscar",
    "He": "mp-23158.poscar",
    "Li": "mp-1018134.poscar",
    "Be": "mp-87.poscar",
    "B":  "mp-160.poscar",
    "C":  "mp-2516584.poscar",
    "N":  "mp-154.poscar",
    "O":  "mp-12957.poscar",
    "F":  "mp-561203.poscar",
    "Ne": "mp-111.poscar",
    "Na": "mp-10172.poscar",
    "Mg": "mp-153.poscar",
    "Al": "mp-134.poscar",
    "Si": "mp-149.poscar",
    "P":  "mp-568348.poscar",
    "S":  "mp-77.poscar",
    "Cl": "mp-22848.poscar",
    "Ar": "mp-23155.poscar",
    "K":  "mp-1184804.poscar",
    "Ca": "mp-45.poscar",
    "Sc": "mp-67.poscar",
    "Ti": "mp-72.poscar",
    "V":  "mp-146.poscar",
    "Cr": "mp-90.poscar",
    "Mn": "mp-35.poscar",
    "Fe": "mp-13.poscar",
    "Co": "mp-102.poscar",
    "Ni": "mp-23.poscar",
    "Cu": "mp-30.poscar",
    "Zn": "mp-79.poscar",
    "Ga": "mp-142.poscar",
    "Ge": "mp-32.poscar",
    "As": "mp-158.poscar",
    "Se": "mp-570481.poscar",
    "Br": "mp-998864.poscar",
    "Kr": "mp-975590.poscar",
    "Rb": "mp-1179656.poscar",
    "Sr": "mp-139.poscar",
    "Y":  "mp-112.poscar",
    "Zr": "mp-131.poscar",
    "Nb": "mp-75.poscar",
    "Mo": "mp-129.poscar",
    "Tc": "mp-113.poscar",
    "Ru": "mp-33.poscar",
    "Rh": "mp-74.poscar",
    "Pd": "mp-2.poscar",
    "Ag": "mp-8566.poscar",
    "Cd": "mp-94.poscar",
    "In": "mp-85.poscar",
    "Sn": "mp-117.poscar",
    "Sb": "mp-104.poscar",
    "Te": "mp-19.poscar",
    "I":  "mp-639751.poscar",
    "Xe": "mp-972256.poscar",
    "Cs": "mp-1055940.poscar",
    "Ba": "mp-122.poscar",
    "La-Lu": None,
    "Hf": "mp-103.poscar",
    "Ta": "mp-569794.poscar",
    "W":  "mp-91.poscar",
    "Re": "mp-1186901.poscar",
    "Os": "mp-49.poscar",
    "Ir": "mp-101.poscar",
    "Pt": "mp-126.poscar",
    "Au": "mp-81.poscar",
    "Hg": "mp-1017981.poscar",
    "Tl": "mp-82.poscar",
    "Pb": "mp-20483.poscar",
    "Bi": "mp-567597.poscar",
    "Po": None,
    "At": None,
    "Rn": None,
    "Fr": None,
    "Ra": None,
    "Ac-Lr": None,
    "Rf": None,
    "Db": None,
    "Sg": None,
    "Bh": None,
    "Hs": None,
    "Mt": None,
    "Ds": None,
    "Rg": None,
    "Cn": None,
    "Nh": None,
    "Fl": None,
    "Mc": None,
    "Lv": None,
    "Ts": None,
    "Og": None,
    "La": "mp-26.poscar",
    "Ce": "mp-28.poscar",
    "Pr": "mp-97.poscar",
    "Nd": "mp-123.poscar",
    "Pm": "mp-867200.poscar",
    "Sm": "mp-69.poscar",
    "Eu": "mp-21462.poscar",
    "Gd": "mp-155.poscar",
    "Tb": "mp-18.poscar",
    "Dy": "mp-88.poscar",
    "Ho": "mp-144.poscar",
    "Er": "mp-1184115.poscar",
    "Tm": "mp-143.poscar",
    "Yb": None,
    "Lu": "mp-973571.poscar",
    "Ac": "mp-862690.poscar",
    "Th": "mp-37.poscar",
    "Pa": "mp-62.poscar",
    "U":  "mp-44.poscar",
    "Np": "mp-11534.poscar",
    "Pu": "mp-582819.poscar",
    "Am": None,
    "Cm": None,
    "Bk": None,
    "Cf": None,
    "Es": None,
    "Fm": None,
    "Md": None,
    "No": None,
    "Lr": None,
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