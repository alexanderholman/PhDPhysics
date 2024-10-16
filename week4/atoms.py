import os
from ase.io import read, write
from ase.optimize import FIRE
from ase.visualize import view
from mace.calculators import mace_mp

# Import the classes.
from classes.POSCAR import POSCAR
from classes.Species import Species
from classes.IonPosition import IonPosition
from classes.Position import Position
from classes.Helper import Helper

# Import the list of structures to process.
from structures_list import *

# add ordered list of species to basis, e.g. could be exclusively Si, or Ge, or mixture

# lattice from basis

# Process structures, generate files, and relax structures.
for st in structures:
    poscar = POSCAR(
        comment=st["comment"],
        species=list(map(
            lambda sp: Species(
                name=sp["name"], 
                ion_positions=list(map(
                    lambda p: IonPosition(Position(p[0]), Position(p[1]), Position(p[2])), 
                    sp["positions"]
                ))
            ), 
            st["species"])
        )

    )

    #if lattice is provided, use it
    if "lattice" in st:
        poscar.lattice = st["lattice"]

    generated_poscar_file = poscar.write("generated.vasp")
    
    # Read the generated structure from a VASP POSCAR file.
    poscar.load_into_ace()

    # Randomise the species of the structure.
    poscar.randomise_species(
        species=["Si", "Ge"],
        weights=[50, 50]
    )

    randomised_poscar_file = poscar.write("generated-random.vasp")

    # Read the generated structure from a VASP POSCAR file.
    poscar.load_into_ace()

    # Expand the structure to a super-cell.
    poscar.expand_to_super_cell(
        x = 2,
        y = 2,
        z = 2
    )

    expanded_poscar_file = poscar.write("generated-expanded.vasp")

    # Read the generated structure from a VASP POSCAR file.
    poscar.load_into_ace()

    poscar.relax(True)

    #view all iterations of the structure
    for i, atoms in enumerate(poscar.atom_iterations):
        poscar.atoms = atoms
        poscar.view()
        poscar.image(f"image-{i}.png")
