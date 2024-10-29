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
# ==================== build end ====================

for filename in expected.keys():
    # load structure
    poscar = POSCAR.from_file(filename=f"{expected[filename]}.poscar", dirname=f"./structures/formation-energy-known/")
    poscar.load_into_ace()
    energy_above_hull = poscar.energy_above_hull()
    formation_energy_per_atom = poscar.formation_energies()
    print(f"energy above hull per atom: {formation_energy_per_atom}, expected: {expected[filename]['energies']['above_hull']}")
    print(f"formation energies per atom: {formation_energy_per_atom}, expected: {expected[filename]['energies']['predicted_formation']}")