from ase.build import stack
from ase.optimize import FIRE
from ase.visualize import view

# here we will calculate the stain energies, formation energies, and interface energies for given structures
# we will as usual be using bulk Si and Ge structures
from classes.POSCAR import POSCAR

si_bulk, ge_bulk = POSCAR.get_relaxed_bulk_structure('Si'), POSCAR.get_relaxed_bulk_structure('Ge')
si_bulk_energy, ge_bulk_energy = si_bulk.atoms.get_potential_energy(), ge_bulk.atoms.get_potential_energy()

# stack the structures
stacked = stack(si_bulk.atoms, ge_bulk.atoms)
stacked.calc = si_bulk.atoms.calc

optimizer = FIRE(stacked)
optimizer.run(fmax=0.005)

stacked_energy = stacked.get_potential_energy()

print(f"Si bulk energy: {si_bulk_energy}")
print(f"Ge bulk energy: {ge_bulk_energy}")
print(f"Stacked energy: {stacked_energy}")
print(f"Interface energy: {stacked_energy - (si_bulk_energy + ge_bulk_energy)}")

view(stacked)