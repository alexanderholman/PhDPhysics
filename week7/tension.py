from copy import deepcopy
from ase.build import bulk, stack
from ase.optimize import FIRE
from mace.calculators import mace_mp

SiLATCONST = 5.4310
GeLATCONST = 5.6578

si_bulk = bulk("Si", crystalstructure='diamond', a=SiLATCONST, cubic=True)
ge_bulk = bulk("Ge", crystalstructure='diamond', a=GeLATCONST, cubic=True)
si_bulk.calc = mace_mp(model="large", dispersion=True, default_dtype="float64")

stack_bulk = stack(si_bulk, ge_bulk, axis=0)

SiGeLATCONST = stack_bulk.cell[1][1]

si_tension = deepcopy(si_bulk)
si_tension.set_cell([[SiGeLATCONST, 0, 0], [0, SiGeLATCONST, 0], [0, 0, SiGeLATCONST]], scale_atoms=True)

si_bulk_optimizer = FIRE(si_bulk)
si_bulk_optimizer.run(fmax=0.001)

si_tension_optimizer = FIRE(si_tension)
si_tension_optimizer.run(fmax=0.001)

strain_energy = tension_energy = si_tension.get_potential_energy() - si_bulk.get_potential_energy()

print(f"Si Strain/tension energy ({SiLATCONST} to {SiGeLATCONST}): {strain_energy}, per atom: {strain_energy / len(si_tension)}")

# GE equivalent
ge_bulk.calc = si_bulk.calc

EqGeLATCONST = SiGeLATCONST / SiLATCONST * GeLATCONST

ge_tension = deepcopy(ge_bulk)

ge_tension.set_cell([[EqGeLATCONST, 0, 0], [0, EqGeLATCONST, 0], [0, 0, EqGeLATCONST]], scale_atoms=True)

ge_bulk_optimizer = FIRE(ge_bulk)
ge_bulk_optimizer.run(fmax=0.001)

ge_tension_optimizer = FIRE(ge_tension)
ge_tension_optimizer.run(fmax=0.001)

strain_energy = tension_energy = ge_tension.get_potential_energy() - ge_bulk.get_potential_energy()

print(f"Ge Strain/tension energy ({GeLATCONST} to {EqGeLATCONST}): {strain_energy}, per atom: {strain_energy / len(ge_tension)}")
