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
    # bulk_images_dirname,
    built_dirname,
    built_vasp_dirname,
    # built_images_dirname,
    relaxed_dirname,
    relaxed_vasp_dirname,
    # relaxed_images_dirname,
]
for dir in DIRS:
    os.makedirs(dir, exist_ok=True)

for b in BULKS:
    write(f"{bulk_vasp_dirname}{b}.vasp", BULKS[b], 'vasp', direct=True)

HEADER = "total,species_a,species_b,n,m,potential_energy,potential_energy_per_atom,formation_energy_per_atom,strain_tension_per_a_atom,strain_compression_per_b_atom,interface_energy\n"

csv_built_filename = f"{built_dirname}energies.csv"
csv_built_handler = open(csv_built_filename, "w")
csv_built_handler.write(HEADER)
csv_built_handler.close()

MIN = 2
MAX = 25

RELAX = True

if RELAX:
    csv_relaxed_filename = f"{relaxed_dirname}energies.csv"
    csv_relaxed_handler = open(csv_relaxed_filename, "w")
    csv_relaxed_handler.write(HEADER)
    csv_relaxed_handler.close()

for i in range(MIN, MAX + 1):
    for j in range(1, i):
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

            bulk_a_supercell = make_supercell(bulk_a, n_matrix)
            bulk_b_supercell = make_supercell(bulk_b, m_matrix)
            atoms = stack(bulk_a_supercell, bulk_b_supercell, axis=0)

            # set calculator
            bulk_a.calc = bulk_b.calc = bulk_a_supercell.calc = bulk_b_supercell.calc = atoms.calc = CALC

            # pre-relaxation

            # write built atoms to file
            write(f"{built_vasp_dirname}{t}_{f[0]}_{n}_{f[1]}_{m}.vasp", atoms, 'vasp', direct=True)

            # strain energies

            # tension energies
            strained_a_supercell = deepcopy(bulk_a_supercell)
            strained_a_supercell.set_cell([[lat_avg * n, 0, 0], [0, lat_avg, 0], [0, 0, lat_avg]], scale_atoms=True)
            strained_a_supercell_tension_energy = (strained_a_supercell.get_potential_energy() - bulk_a_supercell.get_potential_energy()) / len(bulk_a_supercell)

            # compression energies
            strained_b_supercell = deepcopy(bulk_b_supercell)
            strained_b_supercell.set_cell([[lat_avg * m, 0, 0], [0, lat_avg, 0], [0, 0, lat_avg]], scale_atoms=True)
            strained_b_supercell_compression_energy = (strained_b_supercell.get_potential_energy() - bulk_b_supercell.get_potential_energy()) / len(bulk_b_supercell)

            # interface energies
            atoms_interface_energy = atoms.get_potential_energy() - (strained_a_supercell.get_potential_energy() + strained_b_supercell.get_potential_energy())

            # formation energies
            E_super = atoms.get_potential_energy()
            n_Si_super = len([a for a in atoms if a.symbol == f[0]])
            E_Si_bulk = bulk_a.get_potential_energy()
            n_Si_bulk = len(bulk_a)
            m_Ge_super = len([a for a in atoms if a.symbol == f[1]])
            E_Ge_bulk = bulk_b.get_potential_energy()
            m_Ge_bulk = len(bulk_b)
            E_f =  (E_super - (n_Si_super * E_Si_bulk / n_Si_bulk) - (m_Ge_super * E_Ge_bulk / m_Ge_bulk))/(n_Si_super + m_Ge_super)
            formation_energy = E_f

            # write data to csv
            csv_built_handler = open(csv_built_filename, "a")
            csv_built_handler.write(f"{t},{f[0]},{f[1]},{n},{m},{E_super},{E_super / len(atoms)},{formation_energy},{strained_a_supercell_tension_energy},{strained_b_supercell_compression_energy},{atoms_interface_energy}\n")
            csv_built_handler.close()

            if RELAX:
                # relaxation
                bulk_a_optimizer = FIRE(bulk_a)
                bulk_a_optimizer.run(fmax=0.005)
                bulk_a_supercell_optimizer = FIRE(bulk_a_supercell)
                bulk_a_supercell_optimizer.run(fmax=0.005)
                strained_a_supercell_optimizer = FIRE(strained_a_supercell)
                strained_a_supercell_optimizer.run(fmax=0.005)
                bulk_b_optimizer = FIRE(bulk_b)
                bulk_b_optimizer.run(fmax=0.005)
                bulk_b_supercell_optimizer = FIRE(bulk_b_supercell)
                bulk_b_supercell_optimizer.run(fmax=0.005)
                strained_b_supercell_optimizer = FIRE(strained_b_supercell)
                strained_b_supercell_optimizer.run(fmax=0.005)
                atoms_optimizer = FIRE(atoms)
                atoms_optimizer.run(fmax=0.005)

                # write relaxed atoms to file
                write(f"{relaxed_vasp_dirname}{t}_{f[0]}_{n}_{f[1]}_{m}.vasp", atoms, 'vasp', direct=True)

                # strain energies

                # tension energies
                strained_a_supercell_tension_energy = (strained_a_supercell.get_potential_energy() - bulk_a_supercell.get_potential_energy()) / len(bulk_a_supercell)

                # compression energies
                strained_b_supercell_compression_energy = (strained_b_supercell.get_potential_energy() - bulk_b_supercell.get_potential_energy()) / len(bulk_b_supercell)

                # interface energies
                atoms_interface_energy = atoms.get_potential_energy() - (strained_a_supercell.get_potential_energy() + strained_b_supercell.get_potential_energy())

                # formation energies
                E_super = atoms.get_potential_energy()
                E_Si_bulk = bulk_a.get_potential_energy()
                E_Ge_bulk = bulk_b.get_potential_energy()
                E_f =  (E_super - (n_Si_super * E_Si_bulk / n_Si_bulk) - (m_Ge_super * E_Ge_bulk / m_Ge_bulk))/(n_Si_super + m_Ge_super)
                formation_energy = E_f

                # write data to csv
                csv_relaxed_handler = open(csv_relaxed_filename, "a")
                csv_relaxed_handler.write(f"{t},{f[0]},{f[1]},{n},{m},{E_super},{E_super / len(atoms)},{formation_energy},{strained_a_supercell_tension_energy},{strained_b_supercell_compression_energy},{atoms_interface_energy}\n")
                csv_relaxed_handler.close()