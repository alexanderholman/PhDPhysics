import os
from classes.POSCAR import POSCAR

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

dirname = "./structures/formation-energy-known/"

csv_filename = f"{dirname}energies.csv"
os.makedirs(os.path.dirname(csv_filename), exist_ok=True)
csv = open(csv_filename, "w")
csv.write("mp-id,species,energies_per_atom,energies_above_hull_expected,energies_above_hull_calculated,formation_energies_expected,formation_energies_calculated\n")

for filename in expected.keys():
    poscar = POSCAR.from_file(filename=f"{filename}.poscar", dirname=dirname, load_into_ace=True)
    species = ",".join([f"{s.name}({len(s.ion_positions)})" for s in poscar.species])
    per_species_energy_above_hull = ",".join([str(e) for e in poscar.energies_above_hull_multi_species()])
    csv_line = f"{filename},'{species}',{poscar.energy_per_atom()},{expected[filename]['energies']['above_hull']},'{per_species_energy_above_hull}',{expected[filename]['energies']['predicted_formation']},{poscar.formation_energies()}\n"
    csv.write(csv_line)
    del poscar

csv.close()