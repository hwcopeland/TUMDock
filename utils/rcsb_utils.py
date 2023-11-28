# RCSB_utils.py

import os
import shutil
from Bio.PDB import PDBList

def create_directories(dirs):
    for dir_name in dirs:
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        elif dir_name in ["Protein", "Results","Ions"]:
            shutil.rmtree(dir_name)
            os.makedirs(dir_name)

def download_protein(protein):
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(protein, pdir='./Protein/', file_format='pdb')
    os.chdir('./Protein/')
    os.rename(f'pdb{protein.lower()}.ent', f'{protein}.pdb')