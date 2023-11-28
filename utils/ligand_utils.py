# ligand_utils.py

import subprocess
import re

def pdb_to_smiles(script_dir):
    ligand_dir = os.path.join(script_dir, 'Ligands', 'PDB')
    ligand_files = [os.path.join(ligand_dir, f) for f in os.listdir(ligand_dir) if f.endswith('.pdb')]

    smiles_list = []

    for ligand_file in ligand_files:
        command = f'obabel -ipdb "{ligand_file}" -osmi'
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        output_parts = result.stdout.strip().split()
        if output_parts:
            smiles = output_parts[0]
            smiles_list.append(smiles)

    return smiles_list

def generate_pharmacophore_models(smiles_list, output_dir):
    pharmacophore_models = []
    for i, smiles in enumerate(smiles_list):
        if smiles:
            m = Chem.MolFromSmiles(smiles)
            if m is not None:
                AllChem.EmbedMolecule(m)
                pharmacophore_models.append(m)
                Chem.rdmolfiles.MolToMolFile(m, os.path.join(output_dir, f'molecule_{i}.sdf'))
    return pharmacophore_models

def convert_ligand_to_pdbqt(ligand_path, output_path):
    flags = '--minimize --steps 1500 --sd'
    command = f'obabel {ligand_path} -O {output_path} -h {flags}'
    subprocess.run(command, shell=True)
