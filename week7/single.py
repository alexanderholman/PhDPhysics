from ase.build import bulk, stack, make_supercell
from ase.visualize import view
from ase.optimize import FIRE
from mace.calculators import mace_mp

si_bulk = bulk("Si", crystalstructure='diamond', a=5.4310, cubic=True)
ge_bulk = bulk("Ge", crystalstructure='diamond', a=5.6578, cubic=True)

si_supercell_2x1x1 = make_supercell(si_bulk, [[2, 0, 0], [0, 1, 0], [0, 0, 1]])
stacked = stack(si_supercell_2x1x1, ge_bulk, axis=0)

view(stacked)