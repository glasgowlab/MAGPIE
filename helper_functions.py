import pdb_parser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
def extract_info_pdb(pdb_file, chain_id) -> list:
    all_Ca = []
    shaply = {
        'ALA': '#8CFF8C',
        'GLY': '#FFFFFF',
        'LEU': '#455E45',
        'SER': '#FF7042',
        'VAL': '#FF8CFF',
        'THR': '#B84C00',
        'LYS': '#4747B8',
        'ASP': '#A00042',
        'ILE': '#004C00',
        'ASN': '#FF7C70',
        'GLU': '#660000',
        'PRO': '#525252',
        'ARG': '#00007C',
        'PHE': '#534C42',
        'GLN': '#FF4C4C',
        'TYR': '#8C704C',
        'HIS': '#7070FF',
        'CYS': '#FFFF70',
        'MET': '#B8A042',
        'TRP': '#4F4600'
    }
    polar = {
        'ASP': '#E60A0A',
        'GLU': '#E60A0A',
        'CYS': '#E6E600',
        'MET': '#E6E600',
        'LYS': '#145AFF',
        'ARG': '#145AFF',
        'SER': '#FA9600',
        'THR': '#FA9600',
        'PHE': '#3232AA',
        'TYR': '#3232AA',
        'ASN': '#00DCDC',
        'GLN': '#00DCDC',
        'GLY': '#EBEBEB',
        'LEU': '#0F820F',
        'VAL': '#0F820F',
        'ILE': '#0F820F',
        'ALA': '#C8C8C8',
        'TRP': '#B45AB4',
        'HIS': '#8282D2',
        'PRO': '#DC9682',
        'Others': '#BEA06E'
    }
    with open(pdb_file, "r") as file:
        found_chain = False
        for line in file:
            dict_of_atoms = {}
            if "ATOM " in line or "HETATM"  in line:
                if "TER" in line and found_chain:
                    break
                parsed_line = pdb_parser.PDBLineParser(line)
                parsed_line.parse_line()
                if parsed_line.chain_identifier == chain_id and parsed_line.atom_name == 'CA':
                    dict_of_atoms['X'] = parsed_line.x_cord
                    dict_of_atoms['Y'] = parsed_line.y_cord
                    dict_of_atoms['Z'] = parsed_line.z_cord
                    dict_of_atoms['chain'] = parsed_line.chain_identifier
                    dict_of_atoms['residue_index'] = parsed_line.residue_sequence_number
                    dict_of_atoms['AA'] = parsed_line.residue_name
                    dict_of_atoms['atom_name'] = parsed_line.atom_name
                    dict_of_atoms['file'] = pdb_file
                    dict_of_atoms['shapely'] = shaply[parsed_line.residue_name]
                    dict_of_atoms['polar'] = polar[parsed_line.residue_name]
                    found_chain = True
                    all_Ca.append(dict_of_atoms)
    
    return all_Ca
def extract_list_pdb(pdb_list, chain_id) ->list:
    all_Ca = []
    for pdb_file in pdb_list:
        shaply = {
            'ALA': '#8CFF8C',
            'GLY': '#FFFFFF',
            'LEU': '#455E45',
            'SER': '#FF7042',
            'VAL': '#FF8CFF',
            'THR': '#B84C00',
            'LYS': '#4747B8',
            'ASP': '#A00042',
            'ILE': '#004C00',
            'ASN': '#FF7C70',
            'GLU': '#660000',
            'PRO': '#525252',
            'ARG': '#00007C',
            'PHE': '#534C42',
            'GLN': '#FF4C4C',
            'TYR': '#8C704C',
            'HIS': '#7070FF',
            'CYS': '#FFFF70',
            'MET': '#B8A042',
            'TRP': '#4F4600'
        }
        polar = {
            'ASP': '#E60A0A',
            'GLU': '#E60A0A',
            'CYS': '#E6E600',
            'MET': '#E6E600',
            'LYS': '#145AFF',
            'ARG': '#145AFF',
            'SER': '#FA9600',
            'THR': '#FA9600',
            'PHE': '#3232AA',
            'TYR': '#3232AA',
            'ASN': '#00DCDC',
            'GLN': '#00DCDC',
            'GLY': '#EBEBEB',
            'LEU': '#0F820F',
            'VAL': '#0F820F',
            'ILE': '#0F820F',
            'ALA': '#C8C8C8',
            'TRP': '#B45AB4',
            'HIS': '#8282D2',
            'PRO': '#DC9682',
            'Others': '#BEA06E'
        }
        with open(pdb_file, "r") as file:
            found_chain = False
            for line in file:
                dict_of_atoms = {}
                if "ATOM " in line or "HETATM"  in line:
                    if "TER" in line and found_chain:
                        break
                    parsed_line = pdb_parser.PDBLineParser(line)
                    parsed_line.parse_line()
                    if parsed_line.chain_identifier == chain_id and parsed_line.atom_name == 'CA':
                        dict_of_atoms['X'] = parsed_line.x_cord
                        dict_of_atoms['Y'] = parsed_line.y_cord
                        dict_of_atoms['Z'] = parsed_line.z_cord
                        dict_of_atoms['chain'] = parsed_line.chain_identifier
                        dict_of_atoms['residue_index'] = parsed_line.residue_sequence_number
                        dict_of_atoms['AA'] = parsed_line.residue_name
                        dict_of_atoms['atom_name'] = parsed_line.atom_name
                        dict_of_atoms['file'] = pdb_file
                        dict_of_atoms['shapely'] = shaply[parsed_line.residue_name]
                        dict_of_atoms['polar'] = polar[parsed_line.residue_name]
                        found_chain = True
                        all_Ca.append(dict_of_atoms)

    return all_Ca
def extract_connections(pdb_file : str) -> list:
    with open(pdb_file, "r") as file:
        dict_of_atoms = {}
        for line in file:
            if "CONECT" in line:
                positions = []
                if len(line.split("  ")) != 1:
                  line_list = line.split("  ")
                else:
                  line_list = line.split(" ")
                for i in range(len(line_list)):
                    if i == 0 or i == 1:
                        continue
                    if i == len(line_list) - 1:  
                        positions.append(line_list[i][0:-1])  
                    else:
                        positions.append(line_list[i])
                dict_of_atoms[line_list[1]] = positions
    return dict_of_atoms
def extract_info_ligand (pdb_file, chain_id) -> list:
    all_atoms = []

    with open(pdb_file, "r") as file:
        found_chain =  False
        for line in file:
            dict_of_atoms = {}
            if "ATOM " in line or "HETATM" in line:
                if "TER" in line and  found_chain:
                    continue
                parsed_line = pdb_parser.PDBLineParser(line)
                parsed_line.parse_line()
                if parsed_line.chain_identifier == chain_id and not "H" in  parsed_line.atom_name:
                    dict_of_atoms['X'] = parsed_line.x_cord
                    dict_of_atoms['Y'] = parsed_line.y_cord
                    dict_of_atoms['Z'] = parsed_line.z_cord
                    dict_of_atoms['chain'] = parsed_line.chain_identifier
                    dict_of_atoms['atom_serial_number'] = str(parsed_line.atom_serial_number)
                    dict_of_atoms['atom_name'] = parsed_line.atom_name
                    dict_of_atoms['file'] = pdb_file
                    dict_of_atoms['color'] = get_color_code(parsed_line.atom_name[0])
                    all_atoms.append(dict_of_atoms)
                    found_chain =True
    return all_atoms


def get_color_code(atom_name):
    color = {
        'O': '#FF0000',  # Oxygen - Red
        'N': '#0000FF',  # Nitrogen - Blue
        'C': '#808080',  # Carbon - Grey
        'S': '#FFFF00',  # Sulfur - Yellow
        'P': '#FFA500'   # Phosphorus - Orange 
       }
    return color.get(atom_name, "#000000") 


def make_chunks(data: list, thread_count) -> dict:
    """
    Takes a list and splits it into parts based on the thread count
    :param data: a list that needs to be split up amongst the threads
    :param thread_count: the number of threads to use
    :return: None
    """
    threads = {}

    for x in range(0, thread_count):
        threads[x] = []

    thread = 0
    for x in range(0, len(data)):
        threads[thread].append(data[x])
        thread += 1
        if thread == thread_count:
            thread = 0
    return threads


def create_directory(directory_path):
    if os.path.exists(directory_path):
        return f"The directory '{directory_path}' already exists."
    else:
        os.makedirs(directory_path)