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

    # generate filenames
    dirname = "./structures/" + Helper.slugify(st["comment"]) + "/"
    generated_file = dirname + "generated.vasp"
    relaxed_file = dirname + "relaxed.vasp"
    image_file = dirname + "image.png"
    energy_file = dirname + "potential_energy.txt"

    # Create the directory if it does not exist.
    os.makedirs(os.path.dirname(generated_file), exist_ok=True)

    # Write the generated structure to a VASP POSCAR file.
    file = open(generated_file, "w")
    file.write(str(poscar))
    file.close()
    
    # Read the generated structure from a VASP POSCAR file.
    atoms = read(generated_file)

    # Set the periodic boundary conditions to True.
    atoms.set_pbc(True)

    # Initialize the MACE machine learning potential.
    macemp = mace_mp(model="large", dispersion=True, default_dtype="float64")

    # Assign the MACE potential to the atoms object for subsequent calculations.
    atoms.calc = macemp

    # Set up and run the optimization.
    optimizer = FIRE(atoms)

    # Run the optimization until the maximum force on each atom is less than 0.005 eV/Å.
    optimizer.run(fmax=0.005)

    # Write the relaxed structure to a VASP POSCAR file.
    atoms.write(relaxed_file, format='vasp', direct=True)

    # Write potential energy to file
    file = open(energy_file, "w")
    file.write("total: " + str(atoms.get_potential_energy()) + "\nper atom: " + str(atoms.get_potential_energy() / len(atoms)))
    file.close()

    # Write the image of the relaxed structure.
    write(image_file, atoms, rotation="45x,45y,45z", scale=150)

    # View the relaxed structure.
    view(atoms)