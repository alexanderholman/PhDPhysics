from copy import deepcopy
from ase.build import bulk, stack
from ase.optimize import FIRE
from mace.calculators import mace_mp

SiLATCONST = 5.4310
GeLATCONST = 5.6578

si_bulk = bulk("Si", crystalstructure='diamond', a=SiLATCONST, cubic=True)
ge_bulk = bulk("Ge", crystalstructure='diamond', a=GeLATCONST, cubic=True)
ge_bulk.calc = mace_mp(model="large", dispersion=True, default_dtype="float64")

stack_bulk = stack(si_bulk, ge_bulk, axis=0)

SiGeLATCONST = stack_bulk.cell[1][1]

ge_compression = deepcopy(ge_bulk)
ge_compression.set_cell([[SiGeLATCONST, 0, 0], [0, SiGeLATCONST, 0], [0, 0, SiGeLATCONST]], scale_atoms=True)

ge_bulk_optimizer = FIRE(ge_bulk)
ge_bulk_optimizer.run(fmax=0.001)

ge_compression_optimizer = FIRE(ge_compression)
ge_compression_optimizer.run(fmax=0.001)

strain_energy = compression_energy = ge_compression.get_potential_energy() - ge_bulk.get_potential_energy()

print(f"Ge Strain/compression energy ({GeLATCONST} to {SiGeLATCONST}): {strain_energy}, per atom: {strain_energy / len(ge_compression)}")

# SI equivalent
si_bulk.calc = ge_bulk.calc

EqSiGeLATCONST = SiGeLATCONST / GeLATCONST * SiLATCONST

si_compression = deepcopy(si_bulk)
si_compression.set_cell([[EqSiGeLATCONST, 0, 0], [0, EqSiGeLATCONST, 0], [0, 0, EqSiGeLATCONST]], scale_atoms=True)

si_bulk_optimizer = FIRE(si_bulk)
si_bulk_optimizer.run(fmax=0.001)

si_compression_optimizer = FIRE(si_compression)
si_compression_optimizer.run(fmax=0.001)

strain_energy = compression_energy = si_compression.get_potential_energy() - si_bulk.get_potential_energy()

print(f"Si Strain/compression energy ({SiLATCONST} to {EqSiGeLATCONST}): {strain_energy}, per atom: {strain_energy / len(si_compression)}")
