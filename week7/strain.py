from copy import deepcopy
from ase.build import stack
from ase.optimize import FIRE
from mace.calculators import mace_mp
from ase.io import write
from numpy.ma.extras import average

# here we will calculate the stain energies, formation energies, and interface energies for given structures
# we will as usual be using bulk Si and Ge structures
from classes.POSCAR import POSCAR

macemp = mace_mp(model="large", dispersion=True, default_dtype="float64")

atoms_Si = POSCAR.get_relaxed_bulk_structure('Si').atoms
atoms_Si.calc = macemp

a_Si = atoms_Si.cell[0][0]

write('structures/strain/Si_bulk.vasp', atoms_Si, 'vasp')

E_Si_a_Si = atoms_Si.get_potential_energy()
E_Si_a_Si_per_atom = E_Si_a_Si / len(atoms_Si)

results = []

for i in [1.1, 1.09, 1.08, 1.07, 1.06, 1.05, 1.04, 1.03, 1.02, 1.01, 1.0, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90]:
    atoms_Si_strained = deepcopy(atoms_Si)
    atoms_Si_strained.set_cell([[a_Si * i, 0, 0], [0, a_Si * i, 0], [0, 0, a_Si * i]], scale_atoms=True)
    atoms_Si_strained.calc = atoms_Si.calc
    atoms_Si_strained_optimizer = FIRE(atoms_Si_strained)
    atoms_Si_strained_optimizer.run(fmax=0.001)
    write(f"structures/energies/Si_a_{i}.vasp", atoms_Si_strained, 'vasp')
    E_Si_strained = atoms_Si_strained.get_potential_energy()
    E_Si_strained_per_atom = E_Si_strained / len(atoms_Si_strained)
    E_Si_strain = E_Si_strained - E_Si_a_Si
    E_Si_strain_per_atom = E_Si_strain / len(atoms_Si_strained)
    results.append([i, a_Si * i, E_Si_strained, E_Si_strained_per_atom, E_Si_strain, E_Si_strain_per_atom])

for result in results:
    print(','.join(map(str, result)))