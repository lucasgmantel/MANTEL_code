import Bio.PDB
import sys
import freesasa
from pathlib import Path
import matplotlib.pyplot as plt

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

def verify_results(structure):
    print("\n### Computing total solvent-accessible surface area (SASA) using FreeSASA... ###\n")
    structure_conv = freesasa.structureFromBioPDB(structure)
    result = freesasa.calc(structure_conv)
    print(f"Total SASA: {result.totalArea():.1f} Å²")
    print("\n### Computing relative solvent-accessible surface area (rSASA) using FreeSASA... ###\n")
    max_asa = { "ALA": 121, "ARG": 265, "ASN": 187, "ASP": 187, "CYS": 148,
               "GLN": 214, "GLU": 214, "GLY": 97, "HIS": 216, "ILE": 195,
               "LEU": 191, "LYS": 230, "MET": 203, "PHE": 228, "PRO": 154,
               "SER": 143, "THR": 163, "TRP": 264, "TYR": 255, "VAL": 165 }
    res_areas = result.residueAreas()
    residues_labels = []
    relative_sasa_values = []
    for chain, residues in res_areas.items():
        for resnum, resinfo in residues.items():
            resname = resinfo.residueType
            sasa_abs = resinfo.total
            # Calcul du relatif si dispo
            if resname in max_asa:
                sasa_rel = (sasa_abs / max_asa[resname]) * 100
                residues_labels.append(f"{resname}{resnum}")
                relative_sasa_values.append(sasa_rel)
            else:
                print(f"{resname}{resnum} (chaîne {chain}): {sasa_abs:.1f} Å², (référence manquante)")
    plt.figure(figsize=(10, 6))
    plt.bar(residues_labels, relative_sasa_values, color='skyblue')
    plt.xlabel("Residues")
    plt.ylabel("Relative solvent accessibility (RSA) (%)")
    plt.xticks(rotation=90, fontsize=4)
    plt.tight_layout()
    plt.show()

def main():
    check_arguments()
    pdb_file = sys.argv[1]
    verify_file_type(pdb_file)
    protein = load_structure(pdb_file)
    verify_results(protein)

if __name__ == "__main__":
    main()
