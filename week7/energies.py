from copy import deepcopy
from ase.build import stack, make_supercell
from ase.optimize import FIRE, BFGS
from mace.calculators import mace_mp
from ase.io import write
from numpy.ma.extras import average

# here we will calculate the stain energies, formation energies, and interface energies for given structures
# we will as usual be using bulk Si and Ge structures
from classes.POSCAR import POSCAR

macemp = mace_mp(model="large", dispersion=True, default_dtype="float64")

atoms_Si, atoms_Ge = POSCAR.get_relaxed_bulk_structure('Si').atoms, POSCAR.get_relaxed_bulk_structure('Ge').atoms
atoms_Si.calc = atoms_Ge.calc = macemp

atoms_SiGe_8_poscar = POSCAR.get_bulk_structure('Si')
atoms_SiGe_8_poscar.split_species_x(['Si', 'Ge'], [1, 1])
atoms_SiGe_8_poscar.load_into_ace()

a_Si, a_Ge = atoms_Si.cell[0][0], atoms_Ge.cell[0][0]

a_SiGe = average((a_Si, a_Ge))

atoms_SiGe_8 = atoms_SiGe_8_poscar.atoms
atoms_SiGe_8.cell = [[a_SiGe, 0, 0], [0, a_SiGe, 0], [0, 0, a_SiGe]]
atoms_SiGe_8.calc = atoms_Si.calc
atoms_Si_Ge_8_optimizer = FIRE(atoms_SiGe_8)
atoms_Si_Ge_8_optimizer.run(fmax=0.001)

E_SiGe_8_a_SiGe = atoms_SiGe_8.get_potential_energy()
E_SiGe_8_a_SiGe_per_atom = E_SiGe_8_a_SiGe / len(atoms_SiGe_8)

write('structures/energies-unrelaxed/Si_a_Si.vasp', atoms_Si, 'vasp')
write('structures/energies-unrelaxed/Ge_a_Ge.vasp', atoms_Ge, 'vasp')
write('structures/energies-unrelaxed/SiGe_8_a_SiGe.vasp', atoms_SiGe_8, 'vasp')

E_Si_a_Si = atoms_Si.get_potential_energy()
E_Si_a_Si_per_atom = E_Si_a_Si / len(atoms_Si)
E_Ge_a_Ge = atoms_Ge.get_potential_energy()
E_Ge_a_Ge_per_atom = E_Ge_a_Ge / len(atoms_Ge)

atoms_Si_strained = deepcopy(atoms_Si)
atoms_Si_strained.cell = [[a_SiGe, 0, 0], [0, a_SiGe, 0], [0, 0, a_SiGe]]
atoms_Si_strained.calc = atoms_Si.calc
atoms_Si_strained_optimizer = FIRE(atoms_Si_strained)
atoms_Si_strained_optimizer.run(fmax=0.001)

write('structures/energies-unrelaxed/Si_a_SiGe.vasp', atoms_Si_strained, 'vasp')

E_Si_a_SiGe = atoms_Si_strained.get_potential_energy()
E_Si_a_SiGe_per_atom = E_Si_a_SiGe / len(atoms_Si_strained)

E_Si_strain = E_Si_a_SiGe - E_Si_a_Si

atoms_Ge_strained = deepcopy(atoms_Ge)
atoms_Ge_strained.cell = [[a_SiGe, 0, 0], [0, a_SiGe, 0], [0, 0, a_SiGe]]
atoms_Ge_strained.calc = atoms_Ge.calc
atoms_Ge_strained_optimizer = FIRE(atoms_Ge_strained)
atoms_Ge_strained_optimizer.run(fmax=0.001)

write('structures/energies-unrelaxed/Ge_a_SiGe.vasp', atoms_Ge_strained, 'vasp')

E_Ge_a_SiGe = atoms_Ge_strained.get_potential_energy()
E_Ge_a_SiGe_per_atom = E_Ge_a_SiGe / len(atoms_Ge_strained)

E_Ge_strain = E_Ge_a_SiGe - E_Ge_a_Ge

atoms_SiGe = stack(atoms_Si, atoms_Ge)
atoms_SiGe.calc = atoms_Si.calc
atoms_SiGe_optimizer = FIRE(atoms_SiGe)
atoms_SiGe_optimizer.run(fmax=0.001)

write('structures/energies-unrelaxed/SiGe_a_SiGe.vasp', atoms_SiGe, 'vasp')

E_SiGe_a_SiGe = atoms_SiGe.get_potential_energy()
E_SiGe_a_SiGe_per_atom = E_SiGe_a_SiGe / len(atoms_SiGe)

E_interface = E_SiGe_a_SiGe - (E_Si_a_SiGe + E_Ge_a_SiGe)

print(f"Si energy: {E_Si_a_Si} ({a_Si}) vs {E_Si_a_SiGe} ({a_SiGe}), per atom: {E_Si_a_Si_per_atom} ({a_Si}) vs {E_Si_a_SiGe_per_atom} ({a_SiGe})")
print(f"Ge energy: {E_Ge_a_Ge} ({a_Ge}) vs {E_Ge_a_SiGe} ({a_SiGe}), per atom: {E_Ge_a_Ge_per_atom} ({a_Ge}) vs {E_Ge_a_SiGe_per_atom} ({a_SiGe})")
print(f"SiGe8 energy: {E_SiGe_8_a_SiGe}, per atom: {E_SiGe_8_a_SiGe_per_atom}")
print(f"SiGe energy: {E_SiGe_a_SiGe}, per atom: {E_SiGe_a_SiGe_per_atom}")
print("----------")
print(f"Expected energy of E_SiGe_8_a_SiGe should be greater than E_Si_a_Si({a_Si}): {E_Si_a_Si} < {E_SiGe_8_a_SiGe}: {E_Si_a_Si < E_SiGe_8_a_SiGe} ({E_SiGe_8_a_SiGe - E_Si_a_Si}), and E_Ge_a_Ge({a_Ge}) {E_Ge_a_Ge} < {E_SiGe_8_a_SiGe}: {E_Ge_a_Ge < E_SiGe_8_a_SiGe} ({E_SiGe_8_a_SiGe - E_Ge_a_Ge})")
print(f"Expected energy of E_SiGe_a_SiGe should be greater than E_Si_a_Si({a_Si}) + E_Ge_a_Ge({a_Ge}): {E_Si_a_Si + E_Ge_a_Ge} < {E_SiGe_a_SiGe}: {E_Si_a_Si + E_Ge_a_Ge < E_SiGe_a_SiGe} ({E_SiGe_a_SiGe - (E_Si_a_Si + E_Ge_a_Ge)})")
print(f"Expected energy of E_SiGe_a_SiGe should be less than E_Si_a_SiGe({a_SiGe}) + E_Ge_a_SiGe({a_SiGe}): {E_Si_a_SiGe + E_Ge_a_SiGe} > {E_SiGe_a_SiGe}: {E_Si_a_SiGe + E_Ge_a_SiGe > E_SiGe_a_SiGe} ({E_SiGe_a_SiGe - (E_Si_a_SiGe + E_Ge_a_SiGe)})")
print("----------")
print(f"Strain energy from from application of tension to Si bulk ({a_Si} to {a_SiGe}): {E_Si_strain}")
print(f"Strain energy from from application of compression to Ge bulk ({a_Ge} to {a_SiGe}): {E_Ge_strain}")
print(f"Interface energy: {E_interface}")