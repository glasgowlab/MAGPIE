import math
from difflib import SequenceMatcher
import argparse
import os
import pathlib

from Bio import pairwise2
from Bio.pairwise2 import format_alignment


import os


VALID_3_LINE_STARTERS = {"ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I", "LYS":"K", "LEU":"L", "MET":"M",
                         "ASN":"N", "PHE":"F", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S", "THR":"T", "VAL":"V", "TRP":"W", "TYR":"T"}
REF_SEQ = ""
def similar(string_a : str, string_b: str) -> int:
    """
    This function is used to get the similarity of two strings and return a score in percent
    :param string_a: string A
    :param string_b: string B
    :return: sim score
    """

    if not string_b or not string_a:
        return 0
    elif string_b == "\n" or string_a == "\n":
        return 0

    if len(string_a) < len(string_b):
        string_c = string_a
        string_a = string_b
        string_b = string_c

    best_score = 0
    for char_index in range(0,len(string_a)+1-len(string_b)):
        raw_score = SequenceMatcher(None, string_a[char_index:char_index+len(string_b)], string_b).ratio()
        if raw_score > best_score:
            best_score = raw_score
    return int(best_score*100)



def get_sequence(list_pdb: list) -> str:
    pre_tested_seq = ""
    last_res_index = 0
    last_chain = "5"
    for line in list_pdb:
        if "TER" in line[0:6] or "END" in line[0:6]:
            continue
        elif line[17:21].strip() in VALID_3_LINE_STARTERS:
            chain = line[21]
            res_index = int(line[22:26].strip())
            if res_index != last_res_index or chain != last_chain:
                last_chain = chain
                last_res_index = res_index
                pre_tested_seq += VALID_3_LINE_STARTERS[line[17:21].strip()]
    return pre_tested_seq


def align_proteins(protein_seq_1 : str, protein_seq_2: str) -> list:
    alignments = format_alignment(*pairwise2.align.globalxx(protein_seq_1, protein_seq_2)[0])
    return alignments



def find_similar_seqs(query_sequences : list, list_pdb : list, seq_identity : int, first_only: bool = True) -> list:
    """
    This function is used to scan the query sequence against the pdb file via a sliding window of 1 and compares using Python's SequenceMatcher
    :param list_pdb:
    :param seq_identity: threshold for the seq identity
    :return: list of filtered pdb lines
    """
    filtered_list_pdb = []
    pre_tested_pdb = []
    pre_tested_seq = ""
    last_res_index = 0
    last_chain = "5"
    for line in list_pdb:
        if "TER" in line[0:6] or "END" in line[0:6]:
            for query_sequence in query_sequences:
                if (similar(pre_tested_seq.lower(), query_sequence.lower())) > seq_identity:
                    filtered_list_pdb += pre_tested_pdb
                    filtered_list_pdb.append(line)
                    if first_only:
                        return filtered_list_pdb
            pre_tested_seq = ""
            pre_tested_pdb = []
        elif line[17:21].strip() in VALID_3_LINE_STARTERS:
            chain = line[21]
            res_index = int(line[22:26].strip())
            if res_index != last_res_index or chain != last_chain:
                last_chain = chain
                last_res_index = res_index
                pre_tested_seq += VALID_3_LINE_STARTERS[line[17:21].strip()]
            pre_tested_pdb.append(line)
    return filtered_list_pdb


def parse_full_pdb(input_pdb_path: str) -> list:
    """
    This function opens and parses the pdb file to look for atoms
    :param input_pdb_path: the path of the pdb
    :return: list form of the filtered pdb
    """
    parsed_pdb = []
    with open(input_pdb_path, "r") as inputPDB:
        for line in inputPDB:
            if "ATOM" in line[0:6] or "TER" in line[0:6] or "END" in line[0:6] or "HETATM" in line[0:6] and not "#" in line:
                parsed_pdb.append(line)
    return parsed_pdb


def parse_pdb(input_pdb_path: str) -> list:
    """
    This function opens and parses the pdb file to look for atoms
    :param input_pdb_path: the path of the pdb
    :return: list form of the filtered pdb
    """
    parsed_pdb = []
    with open(input_pdb_path, "r") as inputPDB:
        for line in inputPDB:
            if "ATOM" in line[0:6] or "TER" in line[0:6] or "END" in line[0:6] and not "#" in line:
                parsed_pdb.append(line)
    return parsed_pdb

def parse_small_molecule_pdb(input_pdb_path: str = None, parsed_pdb:list = None) -> list:
    """
    This function opens and parses the pdb file to look for hetatoms
    :param input_pdb_path: the path of the pdb
    :return: list form of the filtered pdb
    """
    parsed_sm_pdb = []
    if parsed_pdb:
        for line in parsed_pdb:
            if "HETATM" in line[0:6] or "TER" in line[0:6] or "END" in line[0:6] and not "#" in line:
                parsed_sm_pdb.append(line)
    else:
        with open(input_pdb_path, "r") as inputPDB:
            for line in inputPDB:
                if "HETATM" in line[0:6] or "TER" in line[0:6] or "END" in line[0:6] and not "#" in line:
                    parsed_sm_pdb.append(line)
    return parsed_sm_pdb


def arg_spliter(arg_str : str) -> list:
    """
    this function just splits a string into a list
    :param arg_str: the string to be split
    :return: the split string -> list
    """
    if "," in arg_str:
        return arg_str.split(",")
    return [arg_str]


def sm_search(sm_names : str = None, target_sm_chains : str = None, order: str = "chains", take_first_ligand_only: bool = True, input_pdb_path: str = None, filtered_pdb: list = None) -> list:
    """
    This function filtered the pdb file by your sm names and chains
    :param input_pdb_path: the pdb input path
    :param sm_names: the small molecule name you want to filter by
    :param target_sm_chains: the small molecule chains you want to filter by
    :param order: Should you search by name or by chains first [chains,sequence]
    :return: list of filtered pdb lines
    """

    parsed = parse_small_molecule_pdb(input_pdb_path,filtered_pdb)
    sm_filtered = []

    for line in parsed:
        if "TER" in line[0:6] or "END" in line[0:6]:
            sm_filtered.append(line)
            continue
        chain = line[21]
        if order == "chains" and target_sm_chains:
            if chain in target_sm_chains:
                if sm_names:
                    if line[17:20] in sm_names:
                        sm_filtered.append(line)
                else:
                    sm_filtered.append(line)
        elif order == "chains" and sm_names:
            if line[17:20] in sm_names:
                sm_filtered.append(line)
        elif order == "name" and sm_names:
            if line[17:20] in sm_names:
                if target_sm_chains:
                    if chain in target_sm_chains:
                        sm_filtered.append(line)
                else:
                    sm_filtered.append(line)
        elif order == "name" and target_sm_chains:
            if chain in target_sm_chains:
                sm_filtered.append(line)

    #only returning the first ligand if the option is set
    if take_first_ligand_only:
        kept_data = []
        current_chain = ""
        current_res_num = ""
        first_chain = ""
        first_res_num = ""
        for line in sm_filtered:
            if "TER" in line[0:6] or "END" in line[0:6]:
                kept_data.append(line)
            else:
                if current_chain == "":
                    first_chain = line[21]
                    first_res_num = int(line[22:26].strip())
                current_chain = line[21]
                current_res_num= int(line[22:26].strip())
                if current_chain != first_chain or current_res_num != first_res_num:
                    break
                else:
                    kept_data.append(line)
        return kept_data

    #returing all ligands
    return sm_filtered


class Mesh:
    def __init__(self, filtered_atoms):
        import numpy as np
        from sklearn.neighbors import KDTree
        self.points = []
        for line in filtered_atoms:
            if "ATOM" in line[0:6] or "HETATM" in line[0:6]:
                xyz = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
                self.points.append(xyz)
        self.array = np.array(self.points)
        self.tree = KDTree(self.array, leaf_size=max(len(self.array) / 1000, 40))

    def get_close_atoms(self, parsed_pdb, search_radius):
        kept_atoms = []
        for line in parsed_pdb:
            if "TER" in line[0:6] or "END" in line[0:6]:
                continue
            if not ("ATOM" in line[0:6] or "HETATM" in line[0:6]):
                continue
            xyz = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
            distances, ind = self.tree.query([xyz], k=1)
            if distances[0][0] < search_radius:
                kept_atoms.append(line)
        return kept_atoms


def get_sub_struct_lines(hits:list, parsed_pdb: list):

    last_chain = 99999
    sub_struct_count = 0
    last_index = 99999
    last_name = ""
    structs = {}
    accounted_for_ter = False
    for x, line in enumerate(parsed_pdb):
        if "TER" in line[0:6]:
            sub_struct_count +=1
            accounted_for_ter = True
            continue
        if "ATOM" in line[0:6] or "HETATM" in line[0:6]:
            chain = line[21]
            res_name = line[17:20].strip()
            if last_name == "":
                last_name = res_name
            res_index = int(line[22:26].strip())
            if chain != last_chain:
                accounted_for_ter = True
                last_chain = chain
                if chain not in structs:
                    structs[chain] = {}
                    sub_struct_count = 0
                else:
                    sub_struct_count = len(structs[chain])
                structs[chain][sub_struct_count] = []

            if res_name not in VALID_3_LINE_STARTERS:
                    if res_index != last_index:
                        last_index = res_index
                        last_name = res_name
                        if accounted_for_ter:
                            accounted_for_ter = False
                        else:
                            sub_struct_count +=1
                        if sub_struct_count not in structs[chain]:
                            structs[chain][sub_struct_count] = []
            structs[chain][sub_struct_count].append(line)
            accounted_for_ter = False
    final_pdb = []
    used = []
    for line in hits:
        for chain in structs:
            for sub_struct in structs[chain]:
                if [chain,sub_struct] in used:
                    continue

                if line in structs[chain][sub_struct]:
                    final_pdb+=structs[chain][sub_struct]
                    final_pdb+=["TER\n"]
                    used.append([chain,sub_struct])
    return final_pdb


def mesh_search(sm_3_names: str,chains: str, seqs: str,input_pdb_path: str = None, search_radius: float = 8.0,seq_id :int = 95,take_first_only: bool = True) -> list:
    mesh_atoms = []
    mesh_atoms_filtered = []
    parsed_pdb = parse_full_pdb(input_pdb_path)
    if seqs:
        mesh_atoms += seq_search(seqs, parsed_pdb, None, seq_identity=seq_id)
    if chains:
        mesh_atoms += chain_search(chains, parsed_pdb, None)
    if sm_3_names:
        mesh_atoms += sm_search(sm_names=sm_3_names,take_first_ligand_only=take_first_only,input_pdb_path=input_pdb_path)

    for line in mesh_atoms:
        if "TER" in line[0:6] or "END" in line[0:6]:
            continue
        mesh_atoms_filtered.append(line)
    if mesh_atoms_filtered == []:
        return []

    mesh = Mesh(mesh_atoms_filtered)

    hits = mesh.get_close_atoms(parsed_pdb,search_radius)
    filtered = get_sub_struct_lines(hits,parsed_pdb)

    return filtered



def chain_search(chains : str, parsed_pdb : list = None,input_pdb_path: str = None ) -> list:
    """
    This function checks to see if the pdb lines contain a certain chain. If they do they get added to a new list
    :param chains: the query chains
    :param parsed_pdb: the file to list pdb
    :param input_pdb_path: the pdb input path
    :return: list of filtered pdb lines
    """
    if input_pdb_path is not None:
        parsed_pdb = parse_pdb(input_pdb_path)
    if chains:
        binder_chains = arg_spliter(chains)
    else:
        binder_chains = []

    chain_filtered_pdb = []

    for line in parsed_pdb:
        if "HETATM" in line[0:6]:
            continue
        if "ATOM" not in line[0:6]:
            chain_filtered_pdb.append(line)
            continue
        chain = line[21]
        if chain in binder_chains:
            chain_filtered_pdb.append(line)

    return chain_filtered_pdb

def seq_search(seq : str, parsed_pdb : list = None,input_pdb_path: str = None, seq_identity: int = 95, search_first: bool = True) -> list:
    """
    This function parses the pdb and sends it to the similarity function for sequence identity
    :param seq: the query sequence to look for in the PDB
    :param parsed_pdb: the file to list pdb
    :param input_pdb_path: the pdb input path
    :param seq_identity: threshold for the seq identity
    :return: list of filtered pdb lines
    """

    if input_pdb_path is not None:
        parsed_pdb = parse_pdb(input_pdb_path)

    seq = arg_spliter(seq)
    pdb_seq_hits = find_similar_seqs(seq, parsed_pdb, seq_identity,search_first)

    return pdb_seq_hits




def seq_and_chain_search(input_pdb_path: str, binder_seq : str = None, target_protein_seq : str = None, binder_chains : str = None, target_protein_chains : str = None, order : str = "chains", seq_id : int = 95, sm_name : str = None, target_sm_chains : str = None, search_sm : str = "chains", take_first_ligand_only: bool = True, mesh_search_results:list = []) -> [dict,int]:
    """
    This function is used to determine the order of filtering and what to filter out using the user input. The binder/target_protein are identical calls, just different inputs
    All the params from main are passed into here
    :return: dict of the different filtered binder/target_protein/sm
    """

    pdb_hits = {"binder":[],"target_protein":[],"sm":[]}

    if (binder_seq and binder_chains):
        if order == "seq":
            if mesh_search_results:
                pdb_binder_seq_hits = seq_search(binder_seq, mesh_search_results, None, seq_identity=seq_id)
            else:
                pdb_binder_seq_hits = seq_search(binder_seq, None, input_pdb_path, seq_identity=seq_id)
            pdb_hits["binder"] = chain_search(binder_chains, pdb_binder_seq_hits)
        else:
            if mesh_search_results:
                pdb_binder_chain_hits = chain_search(binder_chains, mesh_search_results,None)
            else:
                pdb_binder_chain_hits = chain_search(binder_chains, None, input_pdb_path)
            pdb_hits["binder"] = seq_search(binder_seq, pdb_binder_chain_hits, None, seq_identity=seq_id)
    elif binder_seq:
        if mesh_search_results:
            pdb_hits["binder"] = seq_search(binder_seq,mesh_search_results,None, seq_identity=seq_id)
        else:
            pdb_hits["binder"] = seq_search(binder_seq, None, input_pdb_path, seq_identity=seq_id)
    elif binder_chains:
        if mesh_search_results:
            pdb_hits["binder"] = chain_search(binder_chains, mesh_search_results,None)
        else:
            pdb_hits["binder"] = chain_search(binder_chains, None, input_pdb_path)
    else:
        if(target_protein_seq or target_protein_chains):
            if(target_protein_seq):
                if mesh_search_results:
                    holder = seq_search(target_protein_seq, mesh_search_results, None,
                                                      seq_identity=seq_id, search_first=False)
                    for line in mesh_search_results:
                        if "TER" in line[0:6] or "END" in line[0:6]:
                            continue
                        if "HETATM" in line[0:6]:
                            continue
                        if line not in holder:
                            pdb_hits["binder"].append(line)
                else:
                    holder = seq_search(target_protein_seq, None, input_pdb_path,
                                                            seq_identity=seq_id)
                    parsed_pdb = parse_pdb(input_pdb_path)
                    for line in parsed_pdb:
                        if "TER" or "END" in line:
                            continue
                        if line not in holder:
                            pdb_hits["binder"].append(line)
            else:
                if mesh_search_results:
                    holder = chain_search(target_protein_chains, mesh_search_results)
                    for line in mesh_search_results:
                        if "TER" in line[0:6] or "END" in line[0:6]:
                            continue
                        if "HETATM" in line[0:6]:
                            continue
                        if line not in holder:
                            pdb_hits["binder"].append(line)
                else:
                    holder = chain_search(target_protein_chains, None, input_pdb_path)
                    parsed_pdb = parse_pdb(input_pdb_path)
                    for line in parsed_pdb:
                        if "TER" or "END" in line:
                            continue
                        if line not in holder:
                            pdb_hits["binder"].append(line)


        elif(target_sm_chains or sm_name):
            if mesh_search_results:
                for line in mesh_search_results:
                    if "HETATM" in line[0:6]:
                        continue
                    pdb_hits["binder"].append(line)
            else:
                print("Error, not sure what to do here")
                exit(1)
        else:
            print("ERROR no chain or seq selected for binder")
            exit(1)

    if (target_protein_seq and target_protein_chains):
        if order == "seq":
            if mesh_search_results:
                pdb_target_protein_seq_hits = seq_search(target_protein_seq,mesh_search_results, None, seq_identity=seq_id)
            else:
                pdb_target_protein_seq_hits = seq_search(target_protein_seq, None, input_pdb_path, seq_identity=seq_id)
            pdb_hits["target_protein"] = chain_search(target_protein_chains, pdb_target_protein_seq_hits)
        else:
            if mesh_search_results:
                pdb_target_protein_chain_hits = chain_search(target_protein_chains, mesh_search_results,None)
            else:
                pdb_target_protein_chain_hits = chain_search(target_protein_chains, None, input_pdb_path)
            pdb_hits["target_protein"] = seq_search(target_protein_seq, pdb_target_protein_chain_hits, None, seq_identity=seq_id)
    elif target_protein_seq:
        if mesh_search_results:
            pdb_hits["target_protein"] = seq_search(target_protein_seq,mesh_search_results,None, seq_identity=seq_id)
        else:
            pdb_hits["target_protein"] = seq_search(target_protein_seq, None, input_pdb_path, seq_identity=seq_id)
    elif target_protein_chains:
        if mesh_search_results:
            pdb_hits["target_protein"] = chain_search(target_protein_chains, mesh_search_results,None)
        else:
            pdb_hits["target_protein"] = chain_search(target_protein_chains, None, input_pdb_path)
    elif not sm_name and not target_sm_chains:
        print("ERROR no chain or seq selected for target_protein")
        exit(1)



    if target_sm_chains or sm_name:

        if target_sm_chains:
            target_sm_chains = arg_spliter(target_sm_chains)
        if sm_name:
            sm_name = arg_spliter(sm_name)
        if mesh_search_results:
            pdb_hits["sm"] = sm_search(sm_name,target_sm_chains,search_sm,take_first_ligand_only,input_pdb_path,mesh_search_results)
        else:
            pdb_hits["sm"] = sm_search(sm_name,target_sm_chains,search_sm,take_first_ligand_only,input_pdb_path)

    if (target_protein_seq or target_protein_chains):
        global REF_SEQ
        if REF_SEQ == "" or REF_SEQ is None:
            REF_SEQ = get_sequence(pdb_hits["target_protein"])
        else:
            target_seq = get_sequence(pdb_hits["target_protein"])
            alignment = align_proteins(protein_seq_2 = target_seq,protein_seq_1 = REF_SEQ)
            print(alignment)
            largest_match_ref = max(alignment.split("\n")[0].split("-"),key=len)
            largest_match_align = max(alignment.split("\n")[2].split("-"), key=len)

            match = SequenceMatcher(None, largest_match_ref, largest_match_align).find_longest_match(0, len(largest_match_ref), 0, len(largest_match_align))


            largest_match = largest_match_ref[match.a:match.a + match.size]

            largest_index_ref = REF_SEQ.find(largest_match)
            largest_index_align = target_seq.find(largest_match)


            print(largest_match,largest_index_align,largest_index_ref)
            print(input_pdb_path)

            if largest_index_align >= largest_index_ref:
                print(-1*(largest_index_align - largest_index_ref),"longer")
                return [pdb_hits,-1*(largest_index_align - largest_index_ref)]
            else:
                print(largest_index_ref - largest_index_align,"shorter")
                return [pdb_hits,  largest_index_ref - largest_index_align]

    return [pdb_hits, None]


def write_pdbs(pdb_name : str, pdb_hits: dict, output_path : str, name_conversion: dict = None,binder_chain_rename = "C",target_protein_chain_rename = "A",target_sm_chain_rename = "B",target_ref_start=None) -> None:
    """
    used to write the output pdb file. Takes the info from the filtering step and outputs it into a pdb
    :param pdb_name: the name of the output pdb file
    :param pdb_hits: the dict of lines for the binder, target_protein, and sm use to make the output file
    :param output_path: the output file path directory
    :return: None
    """
    out_atom_index = 0
    with open(output_path+"/"+pdb_name, "w") as output_pdb:
        ter_used_last = False
        for key in pdb_hits:
            out_res_index = 0
            last_res_index = 0
            last_chain = 0
            if key == "target_protein" and target_ref_start != None:
                out_res_index = target_ref_start
            for x, line in enumerate(pdb_hits[key]):
                if "TER" in line[0:6] and not ter_used_last:
                    ter_used_last = True
                    if x == 0:
                        continue
                    output_pdb.write(line)
                if "ATOM" in line[0:6] or "HETATM" in line[0:6]:
                    ter_used_last = False
                    chain = line[21]
                    res_index = int(line[22:26].strip())
                    if res_index != last_res_index or chain != last_chain:
                        out_res_index += 1
                        last_res_index = res_index
                        last_chain = chain
                    atom_name = line[12:16].strip()
                    if name_conversion and "HETATM" in line[0:6]:
                        if atom_name in name_conversion:
                            atom_name = str(name_conversion[atom_name])
                    atom_name += max(0,4 - len(atom_name))*" "
                    out_atom_index += 1
                    if key == "binder":
                        out_chain = ord(binder_chain_rename)
                    elif key == "target_protein":
                        out_chain = ord(target_protein_chain_rename)
                    else:
                        out_chain = ord(target_sm_chain_rename)
                    output_pdb.write(line[0:6] + ((5 - len(str(out_atom_index))) * " ") + str(out_atom_index) + line[11:12] + atom_name + line[16:21] + chr(out_chain) + ((4 - len(str(out_res_index))) * " ") + str(out_res_index) + line[26:])



def parse_FA_file(file_path: str) -> str:
    """
    Used to parse fasta files with comments inside
    :param file_path: path of the fa file
    :return:
    """
    output_data = ""
    with open(file_path, "r") as input_file:
        current_string = ""
        for index,line in enumerate(input_file):
            line = line.split("#")[0]
            if line:
                if ">" in line:
                    if current_string:
                        output_data+= current_string + ","
                    current_string = ""
                    continue
                else:
                    current_string += line
    return output_data


def get_connections(atom_name_o: str, atom_element_o: str, sm_pdb_lines: list,bond_length: float = 2.1)-> []:

    connections = []
    for line1 in sm_pdb_lines:
        if "TER" in line1 or "END" in line1:
            continue
        atom_element = line1[76:78].strip()
        if atom_element == "H":
            continue
        atom_name = line1[12:16].strip()
        if atom_name == atom_name_o:
            atom1_coords = [float(line1[30:38].strip()), float(line1[38:46].strip()), float(line1[46:54].strip())]
            for line2 in sm_pdb_lines:
                if "TER" in line2 or "END" in line2:
                    continue
                atom_element = line2[76:78].strip()
                if atom_element == "H":
                    continue
                atom_name = line2[12:16].strip()
                if atom_name == atom_name_o:
                    continue
                atom2_coords = [float(line2[30:38].strip()), float(line2[38:46].strip()), float(line2[46:54].strip())]
                distence = math.sqrt(math.pow(atom1_coords[0] - atom2_coords[0], 2) + math.pow(atom1_coords[1] - atom2_coords[1],2) + math.pow(atom1_coords[2] - atom2_coords[2], 2))
                if distence < bond_length:
                    connections.append([atom_name,atom_element])
    return connections



def compare_tree_get_rename(tree_o : dict, tree_q: dict) -> dict:
    name_conversion = {}
    used = []
    for atom_o in tree_o:
        for atom_q in tree_q:
            if atom_q[0] in used:
                continue
            if atom_o[1] != atom_q[1]:
                continue
            atom_o_connections = [x[1] for x in tree_o[atom_o]]
            atom_o_connections.sort()
            atom_q_connections = [x[1] for x in tree_q[atom_q]]
            atom_q_connections.sort()
            if atom_o_connections != atom_q_connections:
                continue
            name_conversion[list(atom_q)[0]] = list(atom_o)[0]
            used.append(atom_q[0])
            break

    return name_conversion


def get_sm_atom_names_and_connections(ligand_atoms: list,bond_length: float = 2.1, level_depth: int = 1) -> dict:
    tree = {}
    for line1 in ligand_atoms:
        if "TER" in line1 or "END" in line1:
            continue
        atom_element = line1[76:78].strip()
        atom_name = line1[12:16].strip()
        if atom_element == "H":
            continue
        # tree[atom_name,atom_element] = []
        # connections = get_connections(atom_name,atom_element,ligand_atoms)
        # for connection in connections:
        #     tree[atom_name, atom_element].append(connection)
        #     sec_connections = get_connections(connection[0],connection[1],ligand_atoms)
        #     for sec_connection in sec_connections:
        #         tree[atom_name, atom_element].append(sec_connection)
        #         if sec_connection[0] == atom_name:
        #             continue
        #         third_connections = get_connections(sec_connection[0],sec_connection[1],ligand_atoms)
        #         for third_connection in third_connections:
        #             tree[atom_name, atom_element].append(third_connection)

        tree[atom_name, atom_element] = []
        connections = get_connections(atom_name, atom_element, ligand_atoms)
        for connection in connections:
            tree[atom_name, atom_element].append(connection)
        for x in range(0,level_depth):
            for connection in tree:
                connections = get_connections(connection[0], connection[1], ligand_atoms)
                for connection_sub in connections:
                    tree[atom_name, atom_element].append(connection_sub)
    return tree




def main(input_pdb_path: str, output_path: str, binder_seq : str = None, target_protein_seq : str = None, binder_chains : str = None, anitgen_chains : str = None, order : str = None, seq_identity: int = 95, sm_name : str = None, target_sm_chains : str = None, search_sm : str = "chains", take_first_ligand: bool = True, name_conversion: dict = None, target_rename: str = "A", sm_rename: str = "B", binder_rename : str = "C",mesh_radius: float = 8.0, mesh_search_vars:str = ";;") -> None:

    """
    Main function that calls the parsing/filtering and writing of the new pdb file(s)
    :return: None
    """
    if len(mesh_search_vars) > 2:
        chains = mesh_search_vars.split(";")[0]
        seqs = mesh_search_vars.split(";")[1]
        sm = mesh_search_vars.split(";")[2]
        mesh_results = mesh_search(sm,chains,seqs,input_pdb_path,mesh_radius,seq_identity,take_first_ligand)
        pdb_to_write = seq_and_chain_search(input_pdb_path,binder_seq,target_protein_seq,binder_chains,anitgen_chains,order,seq_identity,sm_name,target_sm_chains,search_sm,take_first_ligand,mesh_results)
    else:
        pdb_to_write = seq_and_chain_search(input_pdb_path,binder_seq,target_protein_seq,binder_chains,anitgen_chains,order,seq_identity,sm_name,target_sm_chains,search_sm,take_first_ligand)

    pdb_name = os.path.basename(input_pdb_path).replace(".pdb","_cleaned.pdb")
    write_pdbs(pdb_name,pdb_to_write[0],output_path, name_conversion,binder_rename,target_rename,sm_rename,pdb_to_write[1])


if __name__ == "__main__":

    """
    add arguments to arg parse
    see below for what each term does
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-S", "--binder_seqs", type=str, help="the sequence(s) of the binder used for identification", default=None)
    parser.add_argument("-s", "--target_protein_seqs", type=str, help="the sequence(s) of the target protein used for identification", default=None)
    parser.add_argument("-F", "--binder_seq_fa", type=str,help="the sequence file of the binder used for identification", default=None)
    parser.add_argument("-f", "--target_protein_seq_fa", type=str,help="the sequence file of the target protein used for identification", default=None)
    parser.add_argument("-L", "--target_sm_3_names", type=str, help="the three letter code of the small molecule(s)", default=None)
    parser.add_argument("-C", "--binder_chains", type=str, help="the chain(s) of the binder used for identification", default=None)
    parser.add_argument("-c", "--target_protein_chains", type=str, help="the chains of the target protein used for identification", default=None)
    parser.add_argument("-l", "--target_sm_chains", type=str, help="the chain(s) the small molecule is/are on", default=None)
    parser.add_argument("-U", "--search_first_protein", type=str, help="if using both chains and sequences search, what should be used to be filtered first?",default="chains",choices=["chains","sequences"])
    parser.add_argument("-u", "--search_first_sm", type=str,help="if using both chains and name search for a small molecule, what should be used to be filtered first?", default="chains", choices=["chains", "name"])
    parser.add_argument("-d", "--seq_identity", type=int,help="sequence identity cutoff for finding similar chains",default=95)
    parser.add_argument("-i", "--input_path", type=str, help="path of the input pdb (file or directory)", required=True)
    parser.add_argument("-o", "--output_path", type=str, help="path of the output directory", required=True)
    parser.add_argument("-N", "--target_protein_chain_rename", type=str, help="what the target ligand chain is supposed to be named",default="A")
    parser.add_argument("-n", "--target_sm_chain_rename", type=str,help="what the target ligand chain is supposed to be named", default="B")
    parser.add_argument("-b", "--binder_chain_rename", type=str,help="what the binder chain is supposed to be named", default="C")
    parser.add_argument("-t", "--take_first_sm_only", type=bool, help="should we take the first instence of the ligand only?",default=True)
    parser.add_argument("-T", "--name_sm_atoms_same", type=bool,help="should rename all matching ligands with the same atom name. Uses the first file as a refrence", default=False)
    parser.add_argument("-r", "--sm_ligand_reference_path", type=str, help="path of the reference file for renaming sm ligands file default is first file in the input dir", default=None)
    parser.add_argument("-B", "--bond_length", type=float,help="distance that is considered a bond between 2 atoms for chemical graphs", default=2.1)
    parser.add_argument("-m", "--search_radius", type=float,help="distance that is considered for finding close things in mesh search", default=8.0)
    parser.add_argument("-M", "--mesh_search", type=str, help="add the chains, seqs, and sm for the mesh filter example: 'A,B;AWTRWARE,AWAWAWAW;TPA,ATP'", default=";;")
    parser.add_argument("-A", "--seq_target_align", type=bool,help="Align the target protein in sequence space and results in pdb numbering aligning. Do not use for small molecules",default=False)
    parser.add_argument("-a", "--seq_target_ref_pdb", type=str,help="Reference structure for target protein in seq_target_align",default=None,required=False)
    args = parser.parse_args()

    # checking to see if a sm or target_protein is present
    if not args.target_protein_seqs and not args.target_protein_chains and not args.target_sm_3_names and not args.target_sm_chains and not args.target_seq_fa:
        print("You must enter an target protein chain and/or sequence or small molecule name and/or chain or both!")
        exit(1)

    # checking to see if the input path is real
    if not os.path.exists(args.input_path):
        print("the input path you entered does not exist!")
        exit(1)

    # checking to see if the output path is real
    if not os.path.exists(args.output_path):
        print("the output path you entered does not exist!")
        exit(1)

    # checking to see if the fa files path is real and if it is a file or dir
    if args.binder_seq_fa:
        if (not os.path.exists(args.binder_seq_fa)):
            print("the binder or target seq fa file does not exist")
            exit(1)
        if not os.path.isfile(args.binder_seq_fa):
            print("the binder or target seq fa file is a directroy and not a file")
    if args.target_protein_seq_fa:
        if not os.path.exists(args.target_protein_seq_fa):
            print("the binder or target seq fa file does not exist")
            exit(1)
        if not os.path.isfile(args.target_protein_seq_fa):
            print("the binder or target seq fa file does not exist")
            exit(1)

    #combining the fa file with the string seqs for the binder
    if args.binder_seq_fa:
        if args.binder_seqs:
            args.binder_seqs += ","+parse_FA_file(args.binder_seq_fa)
        else:

            args.binder_seqs = parse_FA_file(args.binder_seq_fa)
    # combining the fa file with the string seqs for the target
    if args.target_protein_seq_fa:
        if args.target_protein_seqs:
            args.target_protein_seqs += ","+parse_FA_file(args.target_protein_seq_fa)
        else:
            args.target_protein_seqs = parse_FA_file(args.target_protein_seq_fa)

    if args.sm_ligand_reference_path:
        if not os.path.exists(args.sm_ligand_reference_path):
            print("the reference sm file does not exist!")
            exit(1)
        if not os.path.isfile(args.sm_ligand_reference_path):
            print("the reference sm file given was a directory not a file!")
            exit(1)


    if args.seq_target_ref_pdb:
        if not os.path.exists(args.seq_target_ref_pdb):
            print("the reference sequence alignment file does not exist!")
            exit(1)
        if not os.path.isfile(args.seq_target_ref_pdb):
            print("the reference sequence alignment file given was a directory not a file!")
            exit(1)


    ref_tree = {}
    if args.sm_ligand_reference_path:
        sm_name = ""
        target_sm_chains = ""
        if args.target_sm_chains:
            target_sm_chains = arg_spliter(args.target_sm_chains)
        if args.target_sm_3_names:
            sm_name = arg_spliter(args.target_sm_3_names)
        if args.name_sm_atoms_same:
            first_ligand = sm_search(sm_name, target_sm_chains, args.search_first_sm, True,args.sm_ligand_reference_path)
            ref_tree = get_sm_atom_names_and_connections(first_ligand,args.bond_length)




    # running the main method for a single file that is a .pdb
    if os.path.isfile(args.input_path):
        if args.seq_target_ref_pdb:
            print("seq_target_ref_pdb cannot be used with only one pdb file")
            exit(1)

        if pathlib.Path(args.input_path).suffix == ".pdb":
            if args.name_sm_atoms_same:
                sm_name = ""
                target_sm_chains = ""
                if args.target_sm_chains:
                    target_sm_chains = arg_spliter(args.target_sm_chains)
                if args.target_sm_3_names:
                    sm_name = arg_spliter(args.target_sm_3_names)
                first_ligand = sm_search(sm_name, target_sm_chains,args.search_first_sm, True,args.input_path)
                q_tree = get_sm_atom_names_and_connections(first_ligand,args.bond_length)
                name_conversion_key = compare_tree_get_rename(ref_tree, q_tree)


                main(args.input_path, args.output_path, args.binder_seqs, args.target_protein_seqs, args.binder_chains,args.target_protein_chains, args.search_first_protein, args.seq_identity, args.target_sm_3_names, args.target_sm_chains, args.search_first_sm, args.take_first_sm_only, name_conversion_key, args.target_protein_chain_rename, args.target_sm_chain_rename, args.binder_chain_rename, args.search_radius, args.mesh_search)
            else:
                main(args.input_path, args.output_path, args.binder_seqs, args.target_protein_seqs, args.binder_chains,args.target_protein_chains, args.search_first_protein, args.seq_identity, args.target_sm_3_names, args.target_sm_chains, args.search_first_sm, args.take_first_sm_only, {}, args.target_protein_chain_rename, args.target_sm_chain_rename, args.binder_chain_rename, args.search_radius, args.mesh_search)

        else:
            print("non .pdb file types are not supported!")
            exit(1)

    # running the main method for all pdbs in a directory
    else:
        file_counter = 0
        for file in os.listdir(args.input_path):
            if ".pdb" in file:
                file_counter += 1
                sm_name = ""
                target_sm_chains = ""
                if args.target_sm_chains:
                    target_sm_chains = arg_spliter(args.target_sm_chains)
                if args.target_sm_3_names:
                    sm_name = arg_spliter(args.target_sm_3_names)
                if args.name_sm_atoms_same and file_counter == 1 and not args.sm_ligand_reference_path:
                    first_ligand = sm_search( sm_name, target_sm_chains, args.search_first_sm,True,os.path.join(args.input_path,file))
                    ref_tree = get_sm_atom_names_and_connections(first_ligand,args.bond_length)
                    if args.seq_target_ref_pdb and args.seq_target_ref_pdb.split("/")[-1] != file:
                        main(args.seq_target_ref_pdb, args.output_path, args.binder_seqs,
                             args.target_protein_seqs, args.binder_chains, args.target_protein_chains,
                             args.search_first_protein, args.seq_identity, args.target_sm_3_names,
                             args.target_sm_chains, args.search_first_sm, args.take_first_sm_only, {},
                             args.target_protein_chain_rename, args.target_sm_chain_rename, args.binder_chain_rename,
                             args.search_radius, args.mesh_search)

                    main(os.path.join(args.input_path,file), args.output_path, args.binder_seqs, args.target_protein_seqs, args.binder_chains, args.target_protein_chains, args.search_first_protein, args.seq_identity, args.target_sm_3_names, args.target_sm_chains, args.search_first_sm,args.take_first_sm_only,{},args.target_protein_chain_rename, args.target_sm_chain_rename, args.binder_chain_rename, args.search_radius, args.mesh_search)
                else:
                    name_conversion_key = {}
                    first_ligand = sm_search(sm_name, target_sm_chains,args.search_first_sm, True,os.path.join(args.input_path, file))
                    q_tree = get_sm_atom_names_and_connections(first_ligand,args.bond_length)
                    if args.name_sm_atoms_same:
                        name_conversion_key = compare_tree_get_rename(ref_tree,q_tree)
                    if args.seq_target_ref_pdb and args.seq_target_ref_pdb.split("/")[-1] != file:
                        main(args.seq_target_ref_pdb, args.output_path, args.binder_seqs,
                             args.target_protein_seqs, args.binder_chains, args.target_protein_chains,
                             args.search_first_protein, args.seq_identity, args.target_sm_3_names,
                             args.target_sm_chains, args.search_first_sm, args.take_first_sm_only, {},
                             args.target_protein_chain_rename, args.target_sm_chain_rename, args.binder_chain_rename,
                            args.search_radius, args.mesh_search)
                    main(os.path.join(args.input_path,file), args.output_path, args.binder_seqs, args.target_protein_seqs, args.binder_chains, args.target_protein_chains, args.search_first_protein, args.seq_identity, args.target_sm_3_names, args.target_sm_chains, args.search_first_sm,args.take_first_sm_only,name_conversion_key,args.target_protein_chain_rename, args.target_sm_chain_rename, args.binder_chain_rename, args.search_radius, args.mesh_search)

