# A script that moves files for each protein into a separate folder:

import os
import shutil

# Define the path to the directory containing the ColabFold output:
directory_path = '/Users/belensundberg/Desktop/Glasgow_Lab/github_projects/MAGPIE/with_PEP_UnRx'


def move_files_to_folders(directory_path):

    # Get a list of all files in the specified directory
    files = [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))]

    # Create a dictionary to store file keys and corresponding files
    file_dictionary = {}

    # Iterate through the files
    for file_name in files:
        # Extract the first word (you may need to adjust this based on your file naming convention)
        first_word = file_name.split('_')[0]

        # Append the file to the dictionary based on the first word
        file_dictionary.setdefault(first_word, []).append(file_name)

    # Move files into folders based on the first word
    for first_word, files in file_dictionary.items():
        # Create a new folder for each set of files with the same first word
        new_folder = os.path.join(directory_path, first_word)
        os.makedirs(new_folder, exist_ok=True)

        # Move files into the new folder
        for file_name in files:
            source_path = os.path.join(directory_path, file_name)
            destination_path = os.path.join(new_folder, file_name)
            shutil.move(source_path, destination_path)
            print(f"Moved '{file_name}' to '{new_folder}'")

if __name__ == "__main__":
    # Replace 'your_directory_path' with the path to your target directory
    target_directory = '/Users/belensundberg/Desktop/Glasgow_Lab/github_projects/MAGPIE/with_PEP_UnRx/PDBs'
    
    # Move files to folders based on the first word
    move_files_to_folders(target_directory)
