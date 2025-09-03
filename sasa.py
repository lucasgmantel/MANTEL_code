import Bio.PDB
import sys
from pathlib import Path

def check_arguments():
    if len(sys.argv) != 2:
        sys.exit("ERROR : expecting exactly 1 argument.")

def verify_file_type(file_name):
    if not file_name.endswith(".pdb"):
        sys.exit("ERROR : file must be a .pdb file.")

def load_structure(file_name):
    if Path(file_name).exists():
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", file_name)
        return structure
    else:
        sys.exit(f"ERROR : file {file_name} does not exist.")

def get_atom_coordinates(structure):
    all_atoms = list(structure.get_atoms())
    coords = [atom.coord for atom in all_atoms]
    return [(all_atoms[i], coords[i]) for i in range(len(all_atoms))]

def main():
    check_arguments()
    pdb_file = sys.argv[1]
    verify_file_type(pdb_file)
    protein = load_structure(pdb_file)
    atom_coords = get_atom_coordinates(protein)
    for atom, coord in atom_coords:
        print(atom, coord)

if __name__ == "__main__":
    main()
