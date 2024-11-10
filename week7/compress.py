from copy import deepcopy

from ase.build import bulk, stack
from ase.visualize import view
from ase.optimize import FIRE
from mace.calculators import mace_mp

si_bulk = bulk("Si", crystalstructure='diamond', a=5.4310, cubic=True)
ge_bulk = bulk("Ge", crystalstructure='diamond', a=5.6578, cubic=True)
ge_bulk.calc = mace_mp(model="large", dispersion=True, default_dtype="float64")

stack_bulk = stack(si_bulk, ge_bulk, axis=0)

lattice_constant = stack_bulk.cell[1][1]

ge_compression = deepcopy(ge_bulk)
ge_compression.set_cell([[lattice_constant, 0, 0], [0, lattice_constant, 0], [0, 0, lattice_constant]], scale_atoms=True)

ge_bulk_optimizer = FIRE(ge_bulk)
ge_bulk_optimizer.run(fmax=0.001)

ge_compression_optimizer = FIRE(ge_compression)
ge_compression_optimizer.run(fmax=0.001)

strain_energy = compression_energy = ge_compression.get_potential_energy() - ge_bulk.get_potential_energy()

print(f"Strain/compression energy: {strain_energy}")

view(ge_bulk)
view(ge_compression)