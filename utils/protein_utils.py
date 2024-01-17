# protein_utils.py
import subprocess
import os
import shutil
from Bio.PDB import PDBList, PDBParser, PDBIO, Structure, Model, Chain
from rdkit.Chem import rdmolfiles, AllChem

metals_and_ions = ['NA', 'MG', 'K', 'CA', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'CD', 'W', 'AU', 'HG', 'CL', 'BR', 'F', 'I', 'SO4', 'PO4', 'NO3', 'CO3']
standard_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

def create_directories(dirs):
    for dir_name in dirs:
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        elif dir_name in ["Protein", "Results"]:
            shutil.rmtree(dir_name)
            os.makedirs(dir_name)

def remove_metals_and_ions(structure):
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in metals_and_ions:
                    chain.detach_child(residue.get_id())

def get_ligand_residues(structure):
    ligand_residues = []
    for model in structure:
        for chain in model:
            ligand_residues.extend((chain, res) for res in chain if res.get_resname() not in standard_aa and res.get_resname() != "HOH")
    return ligand_residues

def write_ligands_to_files(ligand_residues, output_dir, metals_and_ions_dir, metals_and_ions):
    for chain, ligand in ligand_residues:
        io = PDBIO()
        s = Structure.Structure('Ligand')
        m = Model.Model(0)
        s.add(m)
        c = Chain.Chain(chain.id)
        m.add(c)
        c.add(ligand.copy())
        io.set_structure(s)
        filename = f'ligand_{ligand.get_resname()}.pdb'
        if ligand.get_resname() in metals_and_ions:
            io.save(os.path.join(metals_and_ions_dir, filename))
        else:
            io.save(os.path.join(output_dir, filename))

def get_water_residues(structure):
    water_residues = []
    for model in structure:
        for chain in model:
            water_residues.extend((chain, res.get_id()) for res in chain if res.get_resname() == "HOH")
    return water_residues


def clean_protein_structure(structure, ligand_residues, water_residues, protein):
    for chain, ligand in ligand_residues:
        chain.detach_child(ligand.id)
    for chain, res_id in water_residues:
        chain.detach_child(res_id)
    io = PDBIO()
    io.set_structure(structure)
    io.save(f"{protein}_clean.pdb")

def prepare_receptor(protein):
    os.system(f"prepare_receptor -r {protein}_clean.pdb -o {protein}_clean.pdbqt -v -A hydrogens")
    