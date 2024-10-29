import os
from ase.io import read, write
from ase.optimize import FIRE
from ase.visualize import view
from mace.calculators import mace_mp
from copy import deepcopy
from random import randint

from .Helper import Helper
from .Lattice import Lattice
from .Species import Species
from .IonPosition import IonPosition
from .Position import Position

class POSCAR:
    COORDINATES: list[str] = [
        "Direct",
        "Cartesian"
    ]
    BULKS: dict[str, str] = {
        "H":  "mp-730101.poscar",
        "He": "mp-23158.poscar",
        "Li": "mp-1018134.poscar",
        "Be": "mp-87.poscar",
        "B":  "mp-160.poscar",
        "C":  "mp-2516584.poscar",
        "N":  "mp-154.poscar",
        "O":  "mp-12957.poscar",
        "F":  "mp-561203.poscar",
        "Ne": "mp-111.poscar",
        "Na": "mp-10172.poscar",
        "Mg": "mp-153.poscar",
        "Al": "mp-134.poscar",
        "Si": "mp-149.poscar",
        "P":  "mp-568348.poscar",
        "S":  "mp-77.poscar",
        "Cl": "mp-22848.poscar",
        "Ar": "mp-23155.poscar",
        "K":  "mp-1184804.poscar",
        "Ca": "mp-45.poscar",
        "Sc": "mp-67.poscar",
        "Ti": "mp-72.poscar",
        "V":  "mp-146.poscar",
        "Cr": "mp-90.poscar",
        "Mn": "mp-35.poscar",
        "Fe": "mp-13.poscar",
        "Co": "mp-102.poscar",
        "Ni": "mp-23.poscar",
        "Cu": "mp-30.poscar",
        "Zn": "mp-79.poscar",
        "Ga": "mp-142.poscar",
        "Ge": "mp-32.poscar",
        "As": "mp-158.poscar",
        "Se": "mp-570481.poscar",
        "Br": "mp-998864.poscar",
        "Kr": "mp-975590.poscar",
        "Rb": "mp-1179656.poscar",
        "Sr": "mp-139.poscar",
        "Y":  "mp-112.poscar",
        "Zr": "mp-131.poscar",
        "Nb": "mp-75.poscar",
        "Mo": "mp-129.poscar",
        "Tc": "mp-113.poscar",
        "Ru": "mp-33.poscar",
        "Rh": "mp-74.poscar",
        "Pd": "mp-2.poscar",
        "Ag": "mp-8566.poscar",
        "Cd": "mp-94.poscar",
        "In": "mp-85.poscar",
        "Sn": "mp-117.poscar",
        "Sb": "mp-104.poscar",
        "Te": "mp-19.poscar",
        "I":  "mp-639751.poscar",
        "Xe": "mp-972256.poscar",
        "Cs": "mp-1055940.poscar",
        "Ba": "mp-122.poscar",
        "La-Lu": None,
        "Hf": "mp-103.poscar",
        "Ta": "mp-569794.poscar",
        "W":  "mp-91.poscar",
        "Re": "mp-1186901.poscar",
        "Os": "mp-49.poscar",
        "Ir": "mp-101.poscar",
        "Pt": "mp-126.poscar",
        "Au": "mp-81.poscar",
        "Hg": "mp-1017981.poscar",
        "Tl": "mp-82.poscar",
        "Pb": "mp-20483.poscar",
        "Bi": "mp-567597.poscar",
        "Po": None,
        "At": None,
        "Rn": None,
        "Fr": None,
        "Ra": None,
        "Ac-Lr": None,
        "Rf": None,
        "Db": None,
        "Sg": None,
        "Bh": None,
        "Hs": None,
        "Mt": None,
        "Ds": None,
        "Rg": None,
        "Cn": None,
        "Nh": None,
        "Fl": None,
        "Mc": None,
        "Lv": None,
        "Ts": None,
        "Og": None,
        "La": "mp-26.poscar",
        "Ce": "mp-28.poscar",
        "Pr": "mp-97.poscar",
        "Nd": "mp-123.poscar",
        "Pm": "mp-867200.poscar",
        "Sm": "mp-69.poscar",
        "Eu": "mp-21462.poscar",
        "Gd": "mp-155.poscar",
        "Tb": "mp-18.poscar",
        "Dy": "mp-88.poscar",
        "Ho": "mp-144.poscar",
        "Er": "mp-1184115.poscar",
        "Tm": "mp-143.poscar",
        "Yb": None,
        "Lu": "mp-973571.poscar",
        "Ac": "mp-862690.poscar",
        "Th": "mp-37.poscar",
        "Pa": "mp-62.poscar",
        "U":  "mp-44.poscar",
        "Np": "mp-11534.poscar",
        "Pu": "mp-582819.poscar",
        "Am": None,
        "Cm": None,
        "Bk": None,
        "Cf": None,
        "Es": None,
        "Fm": None,
        "Md": None,
        "No": None,
        "Lr": None,
    }
    RELAXED_BULK_STRUCTURES: dict[str, "POSCAR"] = {}
    comment: str = ""
    scaling_factor: float = 1.0
    lattice: Lattice
    species: list[Species] = []
    selective_dynamics: bool = True
    coordinate_mode: str = "Direct"
    atoms: None
    atom_iterations = []
    written_to: str = None
    is_relaxed: bool = False

    def __init__(
            self,
            comment: str,
            species: list[Species],
            scaling_factor: float = 1.0,
            lattice: Lattice = Lattice(
                a1=[5.43, 0.0, 0.0],
                a2=[0.0, 5.43, 0.0],
                a3=[0.0, 0.0, 5.43]
            ),
            selective_dynamics: bool = True,
            coordinate_mode: str = "Direct"
    ) -> None:
        self.comment = comment
        self.scaling_factor = scaling_factor
        self.lattice = lattice
        self.species = species
        self.selective_dynamics = selective_dynamics
        self.coordinate_mode = coordinate_mode

    @staticmethod
    def get_bulk_structure(species: str) -> "POSCAR":
        if species not in POSCAR.BULKS:
            raise ValueError(f"Unknown species {species}")
        if POSCAR.BULKS[species] is None:
            raise ValueError(f"Missing bulk structure for species {species}")
        return POSCAR.from_file(filename=POSCAR.BULKS[species], dirname=f"./structures/bulk/", load_into_ace=True)

    @staticmethod
    def get_relaxed_bulk_structure(species: str) -> "POSCAR":
        if species not in POSCAR.RELAXED_BULK_STRUCTURES:
            POSCAR.RELAXED_BULK_STRUCTURES[species] = POSCAR.get_bulk_structure(species)
            POSCAR.RELAXED_BULK_STRUCTURES[species].relax()
        return POSCAR.RELAXED_BULK_STRUCTURES[species]

    @staticmethod
    def get_hull_energy(species: str):
        return POSCAR.get_relaxed_bulk_structure(species).energy_per_atom()

    def energy_per_atom(self, relax: bool = False, relax_write: bool = False, relax_filename: str = None, relax_dirname: str = None) -> float:
        if relax or self.is_relaxed is False:
            self.relax(
                write=relax_write,
                filename=relax_filename,
                dirname=relax_dirname
            )
        return self.atoms.get_potential_energy() / len(self.atoms)

    def energies_above_hull(self, relax: bool = False, relax_write: bool = False, relax_filename: str = None, relax_dirname: str = None) -> float:
        if relax or self.is_relaxed is False:
            self.relax(
                write=relax_write,
                filename=relax_filename,
                dirname=relax_dirname
            )
        if len(self.species) > 1:
            raise ValueError("The structure must contain only one species to calculate the energy above hull")
        species = self.species[0]
        if species.name not in self.BULKS:
            raise ValueError(f"Unknown species {species.name}")
        if self.BULKS[species.name] is None:
            raise ValueError(f"Missing bulk structure for species {species.name}")
        return self.atoms.get_potential_energy() / len(self.atoms) - POSCAR.get_hull_energy(species.name)

    def energies_above_hull_multi_species(self, relax: bool = False, relax_write: bool = False, relax_filename: str = None, relax_dirname: str = None) -> list[float]:
        if relax or self.is_relaxed is False:
            self.relax(
                write=relax_write,
                filename=relax_filename,
                dirname=relax_dirname
            )
        energies = []
        for species in self.species:
            if species.name not in self.BULKS:
                raise ValueError(f"Unknown species {species.name}")
            if self.BULKS[species.name] is None:
                raise ValueError(f"Missing bulk structure for species {species.name}")
            species_energies = [map(lambda atom: atom.get_potential_energy(), filter(lambda atom: atom.symbol is species.name, self.atoms))]
            energies.append(sum(species_energies) / len(species_energies) - POSCAR.get_hull_energy(species.name))
        return energies

    def formation_energies(self, relax: bool = False, relax_write: bool = False, relax_filename: str = None, relax_dirname: str = None) -> float:
        if relax or self.is_relaxed is False:
            self.relax(
                write=relax_write,
                filename=relax_filename,
                dirname=relax_dirname
            )
        Esuper = self.atoms.get_potential_energy()
        Especies = 0
        Nspecies = 0
        for species in self.species:
            if species.name not in self.BULKS:
                raise ValueError(f"Unknown species {species.name}")
            if self.BULKS[species.name] is None:
                raise ValueError(f"Missing bulk structure for species {species.name}")
            bulk_structure = POSCAR.get_relaxed_bulk_structure(species.name)
            Especies += len(species.ion_positions) * (bulk_structure.atoms.get_potential_energy() / len(bulk_structure.atoms))
            Nspecies += len(species.ion_positions)
        Ef = (Esuper - Especies) / Nspecies
        return Ef

    def append_to_x(self, poscar: "POSCAR") -> None:
        pass

    def append_to_x_negative(self, poscar: "POSCAR") -> None:
        poscar.append_to_x(self)

    def append_to_y(self, poscar: "POSCAR") -> None:
        pass

    def append_to_y_negative(self, poscar: "POSCAR") -> None:
        poscar.append_to_y(self)

    def append_to_z(self, poscar: "POSCAR") -> None:
        pass

    def append_to_z_negative(self, poscar: "POSCAR") -> None:
        poscar.append_to_z(self)

    def expand_to_super_cell(self, x: int = 1, y:int = 1, z: int = 1) -> None:
        if x == 1 and y == 1 and z == 1:
            return
        if x < 1 or y < 1 or z < 1:
            raise ValueError("The super cell dimensions must be greater than 0.")

        # self.comment += f"-expanded_by_{x}x{y}x{z}"

        # expand lattice
        self.lattice.a1 = [x * a for a in self.lattice.a1]
        self.lattice.a2 = [y * a for a in self.lattice.a2]
        self.lattice.a3 = [z * a for a in self.lattice.a3]

        # expand species
        result_species = []
        i = 0
        while i < x:
            j = 0
            while j < y:
                k = 0
                while k < z:
                    for l, species in enumerate(self.species):
                        if (len(result_species) < l + 1):
                            result_species.append(
                                Species(
                                    name=species.name,
                                    ion_positions=[]
                                )
                            )
                        for ion_position in species.ion_positions:
                            old_x, old_y, old_z = ion_position.x.value, ion_position.y.value, ion_position.z.value
                            new_x, new_y, new_z = old_x + i, old_y + j, old_z + k
                            result_species[l].ion_positions.append(
                                IonPosition(
                                    x=Position(
                                        value=new_x
                                    ),
                                    y=Position(
                                        value=new_y
                                    ),
                                    z=Position(
                                        value=new_z
                                    )
                                )
                            )
                    k += 1
                j += 1
            i += 1
        for species in result_species:
            for ion_position in species.ion_positions:
                ion_position.x.value = ion_position.x.value / x
                ion_position.y.value = ion_position.y.value / y
                ion_position.z.value = ion_position.z.value / z
        self.species = result_species

        # move layers of species to fit the new lattice

    def randomise_species(self, species: list[str] = None, weights: list[int] = None) -> None:
        result_species = []
        # self.comment += "-randomised"
        if species is None:
            species = [s.name for s in self.species]
        if weights is None:
            weights = [1 for _ in species]
        if len(weights) != len(species):
            raise ValueError("The number of species and weights must be the same.")
        weighted_species = []
        for i, s in enumerate(species):
            for _ in range(weights[i]):
                weighted_species.append(s)

        # self.comment += "_to_" + "_".join(species)
        for s in species:
            result_species.append(
                Species(
                    name=s,
                    ion_positions=[]
                )
            )
        for s in self.species:
            for i, ion in enumerate(s.ion_positions):
                # randomly assign a new species
                result_species[species.index(weighted_species[randint(0, len(weighted_species) - 1)])].ion_positions.append(
                    IonPosition(
                        x=ion.x,
                        y=ion.y,
                        z=ion.z
                    )
                )
        self.species = result_species

    def split_species_x(self, species: list[str] = None, weights: list[int] = None) -> None:
        result_species = []
        # self.comment += "-split"
        if species is None:
            species = [s.name for s in self.species]
        if weights is None:
            weights = [1 for _ in species]
        if len(weights) != len(species):
            raise ValueError("The number of species and weights must be the same.")
        species_distances = []
        total_distance = sum(weights)
        distance_travelled = 0
        for i, s in enumerate(species):
            travel_from = distance_travelled
            travel_distance = weights[i]/total_distance
            travel_to = travel_from + travel_distance
            species_distances.append([travel_from, travel_to, s])
            distance_travelled = travel_to

        # self.comment += "_to_" + "_".join(species)
        species_to_index = {}
        for i, s in enumerate(species):
            result_species.append(
                Species(
                    name=s,
                    ion_positions=[]
                )
            )
            species_to_index[s] = i
        for s in self.species:
            for i, ion in enumerate(s.ion_positions):
                # assign new species. ion.x is passed the "relative weight" of the species
                for j, distances in enumerate(species_distances):
                    travel_from, travel_to, new_s = distances
                    if travel_from <= ion.x.value < travel_to:
                        result_species[species_to_index.get(new_s)].ion_positions.append(
                            IonPosition(
                                x=ion.x,
                                y=ion.y,
                                z=ion.z
                            )
                        )
                        break
        self.species = result_species

    def split_species_y(self, species: list[str] = None, weights: list[int] = None) -> None:
        result_species = []
        # self.comment += "-split"
        if species is None:
            species = [s.name for s in self.species]
        if weights is None:
            weights = [1 for _ in species]
        if len(weights) != len(species):
            raise ValueError("The number of species and weights must be the same.")
        species_distances = []
        total_distance = sum(weights)
        distance_travelled = 0
        for i, s in enumerate(species):
            travel_from = distance_travelled
            travel_distance = weights[i]/total_distance
            travel_to = travel_from + travel_distance
            species_distances.append([travel_from, travel_to, s])
            distance_travelled = travel_to

        # self.comment += "_to_" + "_".join(species)
        species_to_index = {}
        for i, s in enumerate(species):
            result_species.append(
                Species(
                    name=s,
                    ion_positions=[]
                )
            )
            species_to_index[s] = i
        for s in self.species:
            for i, ion in enumerate(s.ion_positions):
                # assign new species. ion.y is passed the "relative weight" of the species
                for j, distances in enumerate(species_distances):
                    travel_from, travel_to, new_s = distances
                    if travel_from <= ion.y.value < travel_to:
                        result_species[species_to_index.get(new_s)].ion_positions.append(
                            IonPosition(
                                x=ion.x,
                                y=ion.y,
                                z=ion.z
                            )
                        )
                        break
        self.species = result_species

    def split_species_z(self, species: list[str] = None, weights: list[int] = None) -> None:
        result_species = []
        # self.comment += "-split"
        if species is None:
            species = [s.name for s in self.species]
        if weights is None:
            weights = [1 for _ in species]
        if len(weights) != len(species):
            raise ValueError("The number of species and weights must be the same.")
        species_distances = []
        total_distance = sum(weights)
        distance_travelled = 0
        for i, s in enumerate(species):
            travel_from = distance_travelled
            travel_distance = weights[i]/total_distance
            travel_to = travel_from + travel_distance
            species_distances.append([travel_from, travel_to, s])
            distance_travelled = travel_to

        # self.comment += "_to_" + "_".join(species)
        species_to_index = {}
        for i, s in enumerate(species):
            result_species.append(
                Species(
                    name=s,
                    ion_positions=[]
                )
            )
            species_to_index[s] = i
        for s in self.species:
            for i, ion in enumerate(s.ion_positions):
                # assign new species. ion.z is passed the "relative weight" of the species
                for j, distances in enumerate(species_distances):
                    travel_from, travel_to, new_s = distances
                    if travel_from <= ion.z.value < travel_to:
                        result_species[species_to_index.get(new_s)].ion_positions.append(
                            IonPosition(
                                x=ion.x,
                                y=ion.y,
                                z=ion.z
                            )
                        )
                        break
        self.species = result_species

    def write(self, filename: str, dirname: str = None) -> str:
        if dirname is None:
            dirname = "./structures/" + Helper.slugify(self.comment) + "/"
        filename = dirname + filename

        # Create the directory if it does not exist.
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        # Write the generated structure to a VASP POSCAR file.
        file = open(filename, "w")
        file.write(str(self))
        file.close()

        self.written_to = filename

        return filename

    def load_into_ace(self):
        temp = False
        if self.written_to is None:
            temp = True
            dirname = "./temp/"
            filename = dirname + "POSCAR"
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            file = open(filename, "w")
            file.write(str(self))
            file.close()
            self.written_to = filename
        self.atoms = read(self.written_to)
        self.atom_iterations.append(deepcopy(self.atoms))
        if temp:
            os.remove(self.written_to)
            self.written_to = None
        return self.atoms

    def relax(self, fmax: float = 0.005, write: bool = False, filename: str = None, dirname: str = None) -> None:
        # Set the periodic boundary conditions to True.
        self.atoms.set_pbc(True)

        # Initialize the MACE machine learning potential.
        macemp = mace_mp(model="large", dispersion=True, default_dtype="float64")

        # Assign the MACE potential to the atoms object for subsequent calculations.
        self.atoms.calc = macemp

        # Set up and run the optimization.
        optimizer = FIRE(self.atoms)

        # Run the optimization until the maximum force on each atom is less than 0.005 eV/Å.
        optimizer.run(fmax=fmax)

        # Set the relaxed flag to True.
        self.is_relaxed = True

        if write:
            if dirname is None:
                dirname = "./structures/" + Helper.slugify(self.comment) + "/"
            if filename is None:
                filename = "relaxed.vasp"
            filename = dirname + filename
            # Write the relaxed structure to a VASP POSCAR file.
            self.atoms.write(filename, format='vasp', direct=True)

        self.atom_iterations.append(deepcopy(self.atoms))

    def write_energy(self, filename: str = None, dirname: str = None) -> None:
        if dirname is None:
            dirname = "./structures/" + Helper.slugify(self.comment) + "/"
        if filename is None:
            filename = "energy.txt"
        filename = dirname + filename

        # Write potential energy to file
        file = open(filename, "w")
        file.write("total: " + str(self.atoms.get_potential_energy()) + "\nper atom: " + str(self.atoms.get_potential_energy() / len(self.atoms)))
        file.close()

    def image(self, filename: str = None, dirname: str = None) -> None:
        if dirname is None:
            dirname = "./structures/" + Helper.slugify(self.comment) + "/"
        if filename is None:
            filename = "image.png"
        filename = dirname + filename
        write(filename, self.atoms, rotation="45x,45y,45z", scale=150)

    def view(self) -> None:
        view(self.atoms)

    @staticmethod
    def from_file(filename: str, dirname: str = None, load_into_ace: bool = False) -> "POSCAR":
        filename = dirname + filename
        ase_atoms = read(filename)
        return_poscar = POSCAR(
            comment=f"loaded from {filename}",
            species=[], # species will be added in the next loop
            lattice=Lattice(
                a1=ase_atoms.cell[0],
                a2=ase_atoms.cell[1],
                a3=ase_atoms.cell[2]
            ),
        )
        for s in ase_atoms.get_chemical_symbols():
            if s not in [species.name for species in return_poscar.species]:
                return_poscar.species.append(
                    Species(
                        name=s,
                        ion_positions=[]
                    )
                )
        for i, s in enumerate(ase_atoms.get_chemical_symbols()):
            return_poscar.species[[species.name for species in return_poscar.species].index(s)].ion_positions.append(
                IonPosition(
                    x=Position(
                        value=ase_atoms.get_scaled_positions()[i][0]
                    ),
                    y=Position(
                        value=ase_atoms.get_scaled_positions()[i][1]
                    ),
                    z=Position(
                        value=ase_atoms.get_scaled_positions()[i][2]
                    )
                )
            )
        if (load_into_ace):
            return_poscar.load_into_ace()
        return return_poscar

    def __str__(self) -> str:
        # generate a VASP POSCAR string

        # add comment to the POSCAR string
        poscar = self.comment + "\n"

        # add scaling factor to the POSCAR string
        poscar += str(self.scaling_factor) + "\n"

        # add lattice vectors to the POSCAR string
        poscar += " ".join(map(str, self.lattice.a1)) + "\n"
        poscar += " ".join(map(str, self.lattice.a2)) + "\n"
        poscar += " ".join(map(str, self.lattice.a3)) + "\n"

        # add species names and quantities to the POSCAR string
        for s in self.species:
            if (len(s.ion_positions) > 0):
                poscar += s.name + " "
        poscar += "\n"
        for s in self.species:
            if (len(s.ion_positions) > 0):
                poscar += str(len(s.ion_positions)) + " "
        poscar += "\n"

        # add selective dynamics and coordinate mode to the POSCAR string
        if self.selective_dynamics:
            poscar += "Selective Dynamics\n"
        poscar += self.coordinate_mode + "\n"

        # add ion positions to the POSCAR string
        for s in self.species:
            for ion_position in s.ion_positions:
                poscar += str(ion_position.x) + " "
                poscar += str(ion_position.y) + " "
                poscar += str(ion_position.z) + " "
                if self.selective_dynamics:
                    poscar += ("T" if ion_position.x.selective_dynamics else "F") + " "
                    poscar += ("T" if ion_position.y.selective_dynamics else "F") + " "
                    poscar += ("T" if ion_position.z.selective_dynamics else "F") + " "
                poscar += "\n"
        return poscar.rstrip()