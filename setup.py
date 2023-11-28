import os

# Get the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))
# Working directory
os.chdir(script_dir)

# User input for protein
protein = input("Please enter the protein name: ")
# Define directories
dirs = [f"{script_dir}/Ligands/PDB", f"{script_dir}/Ligands/PDBQT", f"{script_dir}/Ligands/Ions", "Protein", "Results"]

# Ensure directories exist and are clean
for dir_name in dirs:
    # Create directory if it doesn't exist
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    # If it's Protein or Results directory, clean it
    elif dir_name in ["Protein", "Results"]:
        for file_name in os.listdir(dir_name):
            file_path = os.path.join(dir_name, file_name)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)