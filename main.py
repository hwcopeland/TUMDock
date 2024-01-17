# main.py

import os
import sys
from utils.protein_utils import *
from utils.ligand_utils import *
from utils.docking_utils import *
from utils.rcsb_utils import create_directories, download_protein

# ENV

# Vina path
vina_path = "/home/hwcopeland/Sandbox/autodock_vina_1_1_2_linux_x86/bin/vina"
num_modes = 9


## Set up

# Get the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))
# Working directory
os.chdir(script_dir)
# Define directories
dirs = [f"{script_dir}/Ligands/PDB", f"{script_dir}/Ligands/PDBQT", f"{script_dir}/Ligands/Ions",f"{script_dir}/Ligands/Pharmacophores" , "Protein", "Results"]
# Create directories
create_directories(dirs)

## RCSB & Preprocessing

# Download protein
protein = input("Please enter the protein RCSB code: ")
download_protein(protein)
### Process protein ###
parser = PDBParser()
# Assuming that the protien name is the same as the RCSB code, we strip the protien stucture from the PDB file.
structure = parser.get_structure(protein, f"{script_dir}/Protein/{protein}.pdb")
remove_metals_and_ions(structure)
ligand_residues = get_ligand_residues(structure)
write_ligands_to_files(ligand_residues, f'{script_dir}/Ligands/PDB', f'{script_dir}/Ligands/Ions', metals_and_ions)
water_residues = get_water_residues(structure)
clean_protein_structure(structure, ligand_residues, water_residues, protein)
prepare_receptor(protein)
# Process ligands
smiles_list = pdb_to_smiles(script_dir)
pharmacophore_models = generate_pharmacophore_models(smiles_list, f'{script_dir}/Ligands/Pharmacophores/')
ligand_dir = f'{script_dir}/Ligands/PDB/'
ligand_files = [os.path.join(ligand_dir, file) for file in os.listdir(ligand_dir) if file.endswith('.pdb')]

for ligand in ligand_files:
    ligand_p = f"{ligand}"
    output_path = f"{script_dir}/Ligands/PDBQT/{os.path.basename(ligand)[:-4]}.pdbqt"
    convert_ligand_to_pdbqt(ligand_p, output_path)
## Docking
if input("Would you like to view the pockets? (y/n): ") == "y":
    run_fpocket(protein)
    time.sleep(10)
    subprocess.run(f"{script_dir}/Protein/{protein}_out/{protein}_PYMOL.sh", shell=True)
    pocket_number = input("Please identify the pocket you would like to target: ")
    center_x, center_y, center_z = get_pocket_center(protein, pocket_number)
else:
    if input("Would you like to select a pocket based on ligand center? (y/n): ") == "y":
        ligpocket = input("What Ligand would you like to use?: ")
        center_x, center_y, center_z = get_ligand_center(f"{script_dir}/Ligands/PDB/ligand_{ligpocket}.pdb")
    else:
        print("No pocket selection method chosen. Exiting...")
        sys.exit(1)
size_x = 20
size_y = 20
size_z = 20
if input("Would you like to use a custom box size? (y/n): ") == "y":
    size_x = input("Please enter the size of the box in the x direction: ")
    size_y = input("Please enter the size of the box in the y direction: ")
    size_z = input("Please enter the size of the box in the z direction: ") 
if input("Would you like to use a custom number of modes? (y/n): ") == "y":
    num_modes = input("Please enter the number of modes: ")

receptor_path = f"{script_dir}/Protein/{protein}_clean.pdbqt"
ligand_dir = f'{script_dir}/Ligands/PDBQT/'
ligand_files = [os.path.join(ligand_dir, file) for file in os.listdir(ligand_dir) if file.endswith('.pdbqt')]
for ligand in ligand_files:
    ligandqt_path = f"{ligand}"
    output_path = f"{script_dir}/Results/{os.path.basename(ligand)[:-6]}_output.pdbqt"
    command = f"{vina_path} --receptor {receptor_path} --ligand {ligandqt_path} --out {output_path} --center_x {center_x} --center_y {center_y} --center_z {center_z} --size_x {size_x} --size_y {size_y} --size_z {size_z} --num_modes {num_modes}"
    print(f"Running command: {command}")
    subprocess.run(command, shell=True)

## Analysis

# 


## MegaMolBart Iteration


## Iteration

## Analysis