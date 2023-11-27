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

def download_protein(protein, download_similar_ligands=False):
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(protein, pdir='./Protein/', file_format='pdb')
    os.chdir('./Protein/')
    os.rename(f'pdb{protein.lower()}.ent', f'{protein}.pdb')

    if download_similar_ligands:
        query_json = json.dumps(query)
        response = query_rcsb(query_json)

        # Process the response to get the similar ligands
        similar_ligands = process_response(response)

        # Download the similar ligands
        for ligand in similar_ligands:
            download_ligand(ligand)

def query_rcsb(query_json):
    # Create the JSON query
    query = {
        "query": {
            "type": "terminal",
            "service": "structure",
            "parameters": {
                "value": {
                    "entry_id": "",
                    "assembly_id": "1"
                },
                "operator": "strict_shape_match"
            }
        },
        "return_type": "assembly",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 25
            },
            "results_content_type": [
                "experimental"
            ],
            "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc"
                }
            ],
            "scoring_strategy": "combined"
        }
    }
    
    url = "https://search.rcsb.org/rcsbsearch/v1/query"
    headers = {"Content-Type": "application/json"}



    return response_json()
