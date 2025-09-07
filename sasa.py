import Bio.PDB
import sys
import freesasa
from pathlib import Path
import matplotlib.pyplot as plt
import math

# Dictionnaire des rayons de Van der Waals (en Å)
vdw_radii = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    # Ajoute d'autres éléments si nécessaire
}

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
    all_atoms = [atom for atom in structure.get_atoms()
                 if atom.get_parent().get_resname() != "HOH"]
    coords = [atom.coord for atom in all_atoms]
    return [(all_atoms[i], coords[i]) for i in range(len(all_atoms))]

def load_sphere_points(file_path="sphere.txt"):
    """Charge les 92 points de la sphère de Saff et Kuijlaars."""
    with open(file_path, "r") as f:
        lines = f.readlines()[2:]  # On saute les 2 premières lignes
    points = []
    for line in lines:
        x, y, z = map(float, line.strip().split())
        points.append((x, y, z))
    return points

def generate_atom_sphere(atom_coord, vdw_radius, sphere_points):
    """Génère une sphère centrée sur un atome, avec le bon rayon."""
    total_radius = vdw_radius + 1.4  # Rayon de Van der Waals + rayon de la sonde
    scaled_points = []
    for (x, y, z) in sphere_points:
        new_x = atom_coord[0] + x * total_radius
        new_y = atom_coord[1] + y * total_radius
        new_z = atom_coord[2] + z * total_radius
        scaled_points.append((new_x, new_y, new_z))
    return scaled_points

def calculate_sasa(all_atoms, sphere_points, vdw_radii):
    # Dictionnaire pour stocker les points de chaque sphère par atome
    atom_spheres = {}

    # Génération des sphères pour chaque atome
    for atom, coord in all_atoms:
        element = atom.element.upper()
        if element not in vdw_radii:
            continue  # Ignore les atomes sans rayon défini

        # Génère la sphère pour cet atome
        R = vdw_radii[element] + 1.4
        sphere = generate_atom_sphere(coord, vdw_radii[element], sphere_points)
        atom_spheres[atom] = sphere  # Stocke la sphère pour cet atome

    # Calcul des points exposés pour chaque atome
    sasa_per_residue = {}
    for atom, coord in all_atoms:
        element = atom.element.upper()
        if element not in vdw_radii:
            continue

        # Récupère la sphère de cet atome
        sphere = atom_spheres[atom]
        exposed_points = 0

        # Pour chaque point de la sphère de cet atome
        for point in sphere:
            # Vérifie les collisions avec les sphères des AUTRES atomes
            if is_point_exposed(point, all_atoms, atom, vdw_radii):
                exposed_points += 1

        # Surface par point pour une sphère de rayon R
        R = vdw_radii[element] + 1.4
        point_area = (4 * math.pi * (R**2)) / len(sphere_points)
        atom_sasa = point_area * exposed_points

        # Cumule par résidu
        residue_id = atom.get_parent().get_id()[1]
        residue_name = atom.get_parent().get_resname()
        residue_key = f"{residue_name}{residue_id}"
        if residue_key not in sasa_per_residue:
            sasa_per_residue[residue_key] = 0
        sasa_per_residue[residue_key] += atom_sasa

    return sasa_per_residue

def is_point_exposed(point, all_atoms, current_atom, vdw_radii):
    """Vérifie si un point est exposé au solvant en excluant l'atome courant."""
    for atom, coord in all_atoms:
        if atom == current_atom:
            continue  # Ignore l'atome courant

        element = atom.element.upper()
        if element not in vdw_radii:
            continue

        # Distance entre le point et l'atome
        distance = math.sqrt(
            (point[0] - coord[0])**2 +
            (point[1] - coord[1])**2 +
            (point[2] - coord[2])**2
        )

        # Rayon de l'atome + rayon de la sonde (1.4 Å)
        if distance < vdw_radii[element] + 1.4:
            return False  # Point en collision avec un autre atome
    return True  # Point exposé

def verify_results(structure):
    print("\n### Computing total solvent-accessible surface area (SASA) using FreeSASA... ###\n")
    structure_conv = freesasa.structureFromBioPDB(structure)
    result = freesasa.calc(structure_conv)
    print(f"Total SASA (FreeSASA): {result.totalArea():.1f} Å²")

    # Calcul de la SASA "à la main"
    all_atoms = get_atom_coordinates(structure)
    sphere_points = load_sphere_points()
    sasa_per_residue = calculate_sasa(all_atoms, sphere_points, vdw_radii)
    total_sasa = sum(sasa_per_residue.values())
    print(f"Total SASA (méthode Saff & Kuijlaars) : {total_sasa:.1f} Å²")

    # Dictionnaire des valeurs maximales de SASA par résidu (référence)
    max_asa = {
        "ALA": 121, "ARG": 265, "ASN": 187, "ASP": 187, "CYS": 148,
        "GLN": 214, "GLU": 214, "GLY": 97, "HIS": 216, "ILE": 195,
        "LEU": 191, "LYS": 230, "MET": 203, "PHE": 228, "PRO": 154,
        "SER": 143, "THR": 163, "TRP": 264, "TYR": 255, "VAL": 165
    }

    # Extraction des SASA par résidu depuis FreeSASA
    res_areas_freesasa = result.residueAreas()

    # Préparation des données pour le barplot
    residues_labels = []
    rsa_freesasa_values = []
    rsa_manual_values = []

    for chain, residues in res_areas_freesasa.items():
        for resnum, resinfo in residues.items():
            resname = resinfo.residueType
            residue_key = f"{resname}{resnum}"
            sasa_freesasa = resinfo.total
            sasa_manual = sasa_per_residue.get(residue_key, 0.0)

            # Calcul des RSA (Relative Solvent Accessibility)
            if resname in max_asa:
                rsa_freesasa = (sasa_freesasa / max_asa[resname]) * 100
                rsa_manual = (sasa_manual / max_asa[resname]) * 100

                residues_labels.append(residue_key)
                rsa_freesasa_values.append(rsa_freesasa)
                rsa_manual_values.append(rsa_manual)

    # Tracé du barplot
    x = range(len(residues_labels))
    width = 0.35  # Largeur des barres

    fig, ax = plt.subplots(figsize=(15, 8))
    rects1 = ax.bar([i - width/2 for i in x], rsa_freesasa_values, width, label='RSA FreeSASA', color='blue')
    rects2 = ax.bar([i + width/2 for i in x], rsa_manual_values, width, label='RSA Ton Algo', color='orange')

    ax.set_ylabel('RSA (%)')
    ax.set_title('Comparaison des RSA par résidu')
    ax.set_xticks(x)
    ax.set_xticklabels(residues_labels, rotation=90, fontsize=6)
    ax.legend()

    plt.tight_layout()
    plt.savefig('rsa_comparison.png', dpi=300, bbox_inches='tight')  # Sauvegarde au format PNG
    plt.close()  # Ferme la figure pour éviter qu'elle ne s'affiche

def main():
    check_arguments()
    pdb_file = sys.argv[1]
    verify_file_type(pdb_file)
    protein = load_structure(pdb_file)
    verify_results(protein)

if __name__ == "__main__":
    main()
