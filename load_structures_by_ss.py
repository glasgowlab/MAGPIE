import os
import pymol
from pymol import stored

from pymol import cmd

# Change the directory to the current working directory (should have your PDB files)
pdb_dir = os.getcwd()
print('pdb directory: ', pdb_dir)

#create directory: 
output_dir_path = os.getcwd() 
print('output_dir_path: ', output_dir_path)

# Change to the PDB directory
os.chdir(pdb_dir)

# Get a list of all PDB files in the directory
pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]

extra_helices = ['A0A0P6YHQ1_PEP_rx', 'A0A357ARF4_PEP_rx', 'A0A1P8UHZ0_PEP_rx', 'A0A845QDK5_PEP_rx', 'A0A940J4H3_PEP_rx', 'A0A6P2BVL8_PEP_rx', 'A0A957EMT2_PEP_rx', 'A0A7W2T3Z8_PEP_rx', 'A0A5C6B571_PEP_rx', 'A0A356DVB6_PEP_rx', 'A0A4Y8ZXG8_PEP_rx', 'I8HY23_PEP_rx', 'A0A0P6WUL9_PEP_rx', 'A0A7W1TMQ5_PEP_rx', 'A0A9D1CHV2_PEP_rx']

extra_sheets = ['A0A953BG47_PEP_rx', 'A0A357KWS7_PEP_rx', 'A0A6N7EXS0_PEP_rx', 'A0A9D5T0S7_PEP_rx', 'A0A143DEI8_PEP_rx', 'A0A9D9EWD2_PEP_rx', 'A0A1H8ZQQ9_PEP_rx', 'A0A077AYE3_PEP_rx', 'A0A0F6YL10_PEP_rx', 'A0A960DEI3_PEP_rx', 'A0A0M8K9I8_PEP_rx', 'A0A1G9V4P9_PEP_rx', 'A0A2N3Q0X7_PEP_rx', 'A0A0P9EST4_PEP_rx', 'A0A350LL64_PEP_rx', 'A0A2D9T5C8_PEP_rx', 'A0A4V2AJ23_PEP_rx', 'A0A7L9RTV1_PEP_rx', 'A0A931J453_PEP_rx']

extra_loops = ['A0A9D6WYZ8_PEP_rx', 'A0A949N606_PEP_rx', 'A0A085W4P3_PEP_rx', 'A0A6B1B3N7_PEP_rx', 'A0A6J4Q7L3_PEP_rx', 'A0A1G9NID4_PEP_rx', 'A0A959ZK87_PEP_rx', 'A0A117S2W8_PEP_rx', 'A0A1I5YU75_PEP_rx', 'A0A3L7Y257_PEP_rx', 'A0A255E815_PEP_rx', 'A0A954S093_PEP_rx', 'A0A1I5PL87_PEP_rx']

reference_pdb = '1PFK'
print(str(reference_pdb))
cmd.fetch(reference_pdb, '1PFK')
cmd.remove('resname HOH')
os.remove('1PFK.cif')

found_extra_helices = []
found_extra_sheets = []
found_extra_loops = []

# Iterate through the remaining PDB files and align them to the reference
for pdb_file in pdb_files[0:]:
    pdb_file_name = pdb_file.split('.')[0].split('_relaxed')[0]
    pdb_file_name = pdb_file_name+'_rx'

    if pdb_file_name in extra_helices:
        name = pdb_file
        found_extra_helices.append(name)
    if pdb_file_name in extra_sheets:
        name = pdb_file
        found_extra_sheets.append(name)
    if pdb_file_name in extra_loops:
        name = pdb_file
        found_extra_loops.append(name)

def process_and_print(file_list, feature_name, save_path):
    if file_list:
        print(f'{len(file_list)} structures with extra {feature_name} found:')
        for name in file_list:
            print(name)
            object_name = os.path.splitext(os.path.basename(name))[0]  # Extract base name without extension
            cmd.load(name, object_name)
            # Align each structure to the reference
            cmd.align(object_name, reference_pdb + " and chain A", cutoff=2.0, cycles=5, gap=-10.0,
                      extend=-0.5, max_gap=50, object=None,
                      matrix="BLOSUM62", mobile_state=0, target_state=0,
                      quiet=1, max_skip=0, transform=1, reset=0)

        save_path = os.path.join(output_dir_path, f'extra_{feature_name}.pse')
        cmd.save(save_path)
    else:
        print(f'No structures with extra {feature_name} found')

# Process each category
process_and_print(found_extra_helices, 'helices', output_dir_path)
process_and_print(found_extra_loops, 'loops', output_dir_path)
process_and_print(found_extra_sheets, 'sheets', output_dir_path)


