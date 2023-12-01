# docking_utils.py
import os
import subprocess

def run_fpocket(protein):
    command = f"fpocket -f {protein}.pdb -o {protein}_out"
    subprocess.run(command, shell=True)

def get_pocket_center(protein, pocket_number):
    run_fpocket(protein)
    time.sleep(10)

    file_path = f"{script_dir}/Protein/{protein}_out/{protein}_pockets.pqr"
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{file_path} does not exist. Check if fpocket ran correctly and the protein file is correctly formatted.")
    with open(file_path, "r") as file:
        lines = file.readlines()
    lines = [line for line in lines if line.startswith("ATOM") and int(line[22:26].strip()) == int(pocket_number)]
    coords = [list(map(float, re.findall(r"[\d\.-]+", line[30:54]))) for line in lines]
    x_coords = [coord[0] for coord in coords if len(coord) >= 3]
    y_coords = [coord[1] for coord in coords if len(coord) >= 3]
    z_coords = [coord[2] for coord in coords if len(coord) >= 3]

    center_x = sum(x_coords) / len(x_coords) if x_coords else None
    center_y = sum(y_coords) / len(y_coords) if y_coords else None
    center_z = sum(z_coords) / len(z_coords) if z_coords else None
    return center_x, center_y, center_z

def get_ligand_center(ligand):
    file_path = f"{ligand}"
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{file_path} does not exist.")
    with open(file_path, "r") as file:
        lines = file.readlines()
    lines = [line for line in lines if line.startswith("HETATM")]
    coords = [[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in lines]
    x_coords = [coord[0] for coord in coords if len(coord) >= 3]
    y_coords = [coord[1] for coord in coords if len(coord) >= 3]
    z_coords = [coord[2] for coord in coords if len(coord) >= 3]

    center_x = sum(x_coords) / len(x_coords) if x_coords else None
    center_y = sum(y_coords) / len(y_coords) if y_coords else None
    center_z = sum(z_coords) / len(z_coords) if z_coords else None
    return center_x, center_y, center_z
