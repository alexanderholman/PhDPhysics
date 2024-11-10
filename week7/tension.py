from copy import deepcopy

from ase.build import bulk, stack
from ase.visualize import view
from ase.optimize import FIRE
from mace.calculators import mace_mp

si_bulk = bulk("Si", crystalstructure='diamond', a=5.4310, cubic=True)
ge_bulk = bulk("Ge", crystalstructure='diamond', a=5.6578, cubic=True)
si_bulk.calc = mace_mp(model="large", dispersion=True, default_dtype="float64")

stack_bulk = stack(si_bulk, ge_bulk, axis=0)

lattice_constant = stack_bulk.cell[1][1]

si_tension = deepcopy(si_bulk)
si_tension.set_cell([[lattice_constant, 0, 0], [0, lattice_constant, 0], [0, 0, lattice_constant]], scale_atoms=True)

si_bulk_optimizer = FIRE(si_bulk)
si_bulk_optimizer.run(fmax=0.001)

si_tension_optimizer = FIRE(si_tension)
si_tension_optimizer.run(fmax=0.001)

strain_energy = tension_energy = si_tension.get_potential_energy() - si_bulk.get_potential_energy()

print(f"Strain/tension energy: {strain_energy}")

view(si_bulk)
view(si_tension)