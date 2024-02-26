import glob

import Bio.PDB.Chain

import pdb_parser
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import os
from Bio import *
import sys
def extract_atoms_from_chain(chain: Bio.PDB.Chain.Chain, atom_index: str = "" , file = "" ) -> list:

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

    all_atoms = []
    for residue in chain:
        dict_of_atoms = {}
        for atom in residue:
            atom_coordinates = atom.get_coord()
            atom_flag = False
            if atom_index == "" or atom_index == atom.get_name():
                dict_of_atoms['X'] = atom_coordinates[0]
                dict_of_atoms['Y'] = atom_coordinates[1]
                dict_of_atoms['Z'] = atom_coordinates[2]
                dict_of_atoms['chain'] = chain.get_id()
                dict_of_atoms['residue_index'] = int(residue.__repr__().split("resseq=")[1].split(" icode")[0])
                dict_of_atoms['AA'] = residue.get_resname()
                dict_of_atoms['atom_name'] = atom.get_name()
                dict_of_atoms['shapely'] = shaply.get(residue.get_resname(), "#BEA06E")
                dict_of_atoms['polar'] = polar.get(residue.get_resname(), "#BEA06E")
                dict_of_atoms['file'] = file.split("/")[-1]

                if  atom_index == atom.get_name():
                    atom_flag = True
            if atom_flag:
                break
        all_atoms.append(dict_of_atoms)



    return all_atoms


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
                        if residue.get_resname() == "UNK":
                            continue
                        dict_of_atoms['X'] = parsed_line.x_cord
                        dict_of_atoms['Y'] = parsed_line.y_cord
                        dict_of_atoms['Z'] = parsed_line.z_cord
                        dict_of_atoms['chain'] = parsed_line.chain_identifier
                        dict_of_atoms['residue_index'] = parsed_line.residue_sequence_number
                        dict_of_atoms['AA'] = parsed_line.residue_name
                        dict_of_atoms['atom_name'] = parsed_line.atom_name
                        dict_of_atoms['file'] = file.split("/")[-1]
                        dict_of_atoms['shapely'] = shaply.get(residue.get_resname(), "#BEA06E")
                        dict_of_atoms['polar'] = polar.get(residue.get_resname(), "#BEA06E")
                        found_chain = True
                        all_Ca.append(dict_of_atoms)

    return all_Ca
def filter_list(original_list, filter_list):
    # Create a set for faster membership testing
    filter_set = set(filter_list)

    # Use list comprehension to filter elements
    filtered_list = [x for x in original_list if x in filter_set]

    return filtered_list
def extract_info_ligand (pdb_file, chain_id) -> [list,dict]:
    all_atoms = []
    conect_list = []
    with open(pdb_file, "r") as file:

        smaller_id = 1000000000
        larger_id = 0
        smaller_index = 10000000
        larger_index = 0

        for i,line in enumerate(file):
            bonds = {}
            if "HETATM" in line or "CONECT" in line:
                if "TER" in line:
                    continue
                if "CONECT" in line:
                    conect_list.append(line)
                    continue
                parsed_line = pdb_parser.PDBLineParser(line)
                parsed_line.parse_line()
                if parsed_line.chain_identifier == chain_id:
                    if parsed_line.atom_name == "UNK":
                        continue
                    bonds['X'] = parsed_line.x_cord
                    bonds['Y'] = parsed_line.y_cord
                    bonds['Z'] = parsed_line.z_cord
                    bonds['chain'] = parsed_line.chain_identifier
                    bonds['atom_serial_number'] = str(parsed_line.atom_serial_number)
                    bonds['atom_name'] = parsed_line.atom_name
                    bonds['file'] = pdb_file
                    bonds['color'] = get_color_code(parsed_line.atom_name[0])
                    all_atoms.append(bonds)

                    if len(str(parsed_line.atom_serial_number)) < smaller_id:
                        smaller_id = len(str(parsed_line.atom_serial_number))
                    if len(str(parsed_line.atom_serial_number)) > larger_id:
                        larger_id =  len(str(parsed_line.atom_serial_number))
                    if int(parsed_line.atom_serial_number) > larger_index:
                        larger_index = int(parsed_line.atom_serial_number)
                    if int(parsed_line.atom_serial_number) < smaller_index:
                        smaller_index  = int(parsed_line.atom_serial_number)
        bonds = {}
        pairs_created = []
        list_range = list(range(smaller_index,larger_index+1))
        list_range_str = []
        for x in list_range:
            list_range_str.append(str(x))

        if (smaller_id >= 5 or  larger_id >=5) and smaller_id != larger_id:
            print("Unable to solve small-molecule ligand bonds.")
            return [all_atoms,bonds]
        for line in conect_list:
            if "CONECT" in line:
                positions = []
                line_list = line.split(" ")
                if len(line_list) == 1:
                    line_list = line.split("CONECT")
                    line_string = line_list[1].strip()
                    line_list =[]
                    for i in range(larger_id, len(line_string)+larger_id, larger_id):
                        line_list.append(line_string[i-larger_id: i])
                else:
                    line_list = [x.strip() for x in line_list if (not x == "") and not x == 'CONECT']
                line_list = filter_list(line_list,list_range_str)
                unique_pairs = []
                for id in  range(len(line_list)):
                    if id  == 0:
                        continue
                    current_pair = set((line_list[0],line_list[id]))
                    if current_pair not in pairs_created:
                        unique_pairs.append(line_list[id])
                        pairs_created.append(current_pair)
                    if len(unique_pairs) >0:
                        bonds[line_list[0]] = unique_pairs

    return [all_atoms,bonds]



def get_color_code(atom_name):
    color = {
        'O': '#FF0000',  # Oxygen - Red
        'N': '#0000FF',  # Nitrogen - Blue
        'C': '#808080',  # Carbon - Grey
        'S': '#FFFF00',  # Sulfur - Yellow
        'P': '#FFA500',   # Phosphorus - Orange
        'H': '#33FFFF',
        'F': "#FF00FF"
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


def process_residues_to_graph(sequence_logo_targets, is_ligand):

    if is_ligand:
        if "-" in sequence_logo_targets:
            print("Range not allowed in small-molecule ligand-protein interactions!")
            sys.exit()
        list_to_plot = [str(x) for x in sequence_logo_targets.split(",")]

    else:
        plot_list = sequence_logo_targets.split(",")
        list_to_plot = []
        for item in plot_list:
            if "-" in item:
                if int(item.split("-")[0]) > int(item.split("-")[1]):
                    print(f"Invalid Selection {item}. Selections must be in increasing order")
                    continue
                range_plot = [x for x in range(int(item.split("-")[0]), int(item.split("-")[1]) + 1)]
                list_to_plot += range_plot
            else:
                list_to_plot.append(int(item))
    return  list_to_plot
def create_directory(directory_path):
    if os.path.exists(directory_path):
        return f"The directory '{directory_path}' already exists."
    else:
        os.makedirs(directory_path)

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
