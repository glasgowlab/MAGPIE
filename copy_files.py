# Script that copies all rank_001 pdb files to a folder called AF_output in the current directory.

import os
import shutil

def copy_pdb_files(input_folder, output_folder):
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Iterate through each directory in the input folder
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            # Check if the file is a .pdb file containing "rank_001"
            if file.endswith('.pdb') in file:
                # Build the source and destination paths
                source_path = os.path.join(root, file)
                destination_path = os.path.join(output_folder, file)

                # Copy the file to the output folder
                shutil.copy(source_path, destination_path)
                print(f"Copied {file} to {output_folder}")

# Specify the input folder and output folder
input_folder = "/Users/belensundberg/Desktop/Glasgow_Lab/github_projects/MAGPIE/with_PEP_UnRx/"
output_folder = "./BS_PDBs"

# Call the function to copy .pdb files
copy_pdb_files(input_folder, output_folder)
