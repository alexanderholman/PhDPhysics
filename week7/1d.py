import os
from itertools import combinations
from ase.io import read, write
from ase.visualize import view
from ase.build import bulk, stack, make_supercell
from ase.optimize import FIRE
from mace.calculators import mace_mp
from copy import deepcopy
from datetime import datetime

TIMESTAMP = datetime.now().strftime("%Y%m%d%H%M%S")

CALC = mace_mp(model="large", dispersion=True, default_dtype="float64")

SPECIES = GROUP4 = [
    # "C",
    "Si",
    "Ge",
    # "Sn",
    # "Pb",
    # "Fl"
]
LATTICE_CONSTANTS = {
    "C": 3.5668,
    "Si": 5.4310,
    "Ge": 5.6578,
    "Sn": 6.4892,
    "Pb": 4.9503,
    "Fl": 4.8890,
}
ATOMIC_RADII = {
    "C": 0.77,
    "Si": 1.17,
    "Ge": 1.22,
    "Sn": 1.40,
    "Pb": 1.50,
    "Fl": 1.47,
}

BULKS = {s: bulk(s, crystalstructure='diamond', a=LATTICE_CONSTANTS[s], cubic=True) for s in SPECIES}

FLAVOURS = [s for s in combinations(SPECIES, 2)]

base_dirname = f"./structures/stacks/{TIMESTAMP}/"
bulk_dirname = f"{base_dirname}bulk/"
bulk_vasp_dirname = f"{base_dirname}bulk/vasp/"
bulk_images_dirname = f"{base_dirname}bulk/images/"
built_dirname = f"{base_dirname}built/"
built_vasp_dirname = f"{base_dirname}built/vasp/"
built_images_dirname = f"{base_dirname}built/images/"
relaxed_dirname = f"{base_dirname}relaxed/"
relaxed_vasp_dirname = f"{base_dirname}relaxed/vasp/"
relaxed_images_dirname = f"{base_dirname}relaxed/images/"

DIRS = [
    base_dirname,
    bulk_dirname,
    bulk_vasp_dirname,
    bulk_images_dirname,
    built_dirname,
    built_vasp_dirname,
    built_images_dirname,
    relaxed_dirname,
    relaxed_vasp_dirname,
    relaxed_images_dirname,
]
for dir in DIRS:
    os.makedirs(dir, exist_ok=True)

csv_filename = f"{base_dirname}energies.csv"
csv_handler = open(csv_filename, "w")
csv_handler.write("total,species_a,species_b,number_of_a_atoms,number_of_b_atoms,potential_energy,potential_energy_per_atom\n")
csv_handler.close()

MIN = 2
MAX = 5

for i in range(MIN, MAX + 1):
    for j in range(0, i + 1):
        n = j
        m = i - j
        t = n + m
        for f in FLAVOURS:
            bulk_a = deepcopy(BULKS[f[0]])
            bulk_b = deepcopy(BULKS[f[1]])
            lat_a = bulk_a.get_cell()[0][0]
            lat_b = bulk_b.get_cell()[0][0]
            lat_avg = (lat_a + lat_b) / 2
            n_matrix = [[n, 0, 0], [0, 1, 0], [0, 0, 1]]
            m_matrix = [[m, 0, 0], [0, 1, 0], [0, 0, 1]]
            if n == 0:
                # structure is all bulk_m
                atoms = make_supercell(bulk_b, m_matrix)
                atoms.set_cell([[lat_avg * m, 0, 0], [0, lat_avg, 0], [0, 0, lat_avg]], scale_atoms=True)
            elif m == 0:
                # structure is all bulk_n
                atoms = make_supercell(bulk_a, n_matrix)
                atoms.set_cell([[lat_avg * n, 0, 0], [0, lat_avg, 0], [0, 0, lat_avg]], scale_atoms=True)
            else:
                # structure is a stack of bulk_n and bulk_m
                bulk_a_supercell = make_supercell(bulk_a, n_matrix)
                bulk_b_supercell = make_supercell(bulk_b, m_matrix)
                atoms = stack(bulk_a_supercell, bulk_b_supercell, axis=0)

            # set calculator
            bulk_a.calc = bulk_b.calc = atoms.calc = CALC

            # pre-relaxation

            # write built atoms to file
            write(f"{built_vasp_dirname}{t}_{f[0]}_{n}_{f[1]}_{m}.vasp", atoms, 'vasp')

            # strain energies

            # tension energies


            # compression energies

            # interface energies

            # formation energies

            # write data to csv
            csv_handler = open(csv_filename, "a")
            csv_handler.write(f"{t},{f[0]},{f[1]},{n},{m},{atoms.get_potential_energy()},{atoms.get_potential_energy() / len(atoms)}\n")
            csv_handler.close()

            # relaxation
            bulk_optimizer = FIRE(bulk)
            bulk_optimizer.run(fmax=0.005)

            # post-relaxation