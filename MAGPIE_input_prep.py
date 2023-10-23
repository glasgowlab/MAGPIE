
from difflib import SequenceMatcher
import argparse
import os
import pathlib


VALID_3_LINE_STARTERS = {"ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I", "LYS":"K", "LEU":"L", "MET":"M",
                         "ASN":"N", "PHE":"F", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S", "THR":"T", "VAL":"V", "TRP":"W", "TYR":"T"}

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

def find_simular_chains(query_sequences : list, list_pdb : list, seq_identity : int) -> list:
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
                    break
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

def parse_small_molecule_pdb(input_pdb_path: str) -> list:
    """
    This function opens and parses the pdb file to look for hetatoms
    :param input_pdb_path: the path of the pdb
    :return: list form of the filtered pdb
    """
    parsed_pdb = []
    with open(input_pdb_path, "r") as inputPDB:
        for line in inputPDB:
            if "HETATM" in line[0:6] or "TER" in line[0:6] or "END" in line[0:6] and not "#" in line:
                parsed_pdb.append(line)
    return parsed_pdb


def arg_spliter(arg_str : str) -> list:
    """
    this function just splits a string into a list
    :param arg_str: the string to be split
    :return: the split string -> list
    """
    if "," in arg_str:
        return arg_str.split(",")
    return [arg_str]


def sm_search(input_pdb_path: str, sm_names : str = None, target_sm_chains : str = None, order: str = "chains") -> list:
    """
    This function filtered the pdb file by your sm names and chains
    :param input_pdb_path: the pdb input path
    :param sm_names: the small molecule name you want to filter by
    :param target_sm_chains: the small molecule chains you want to filter by
    :param order: Should you search by name or by chains first [chains,sequence]
    :return: list of filtered pdb lines
    """

    parsed = parse_small_molecule_pdb(input_pdb_path)
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
    return sm_filtered

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
        if "ATOM" not in line[0:6]:
            chain_filtered_pdb.append(line)
            continue
        chain = line[21]
        if chain in binder_chains:
            chain_filtered_pdb.append(line)

    return chain_filtered_pdb

def seq_search(seq : str, parsed_pdb : list = None,input_pdb_path: str = None, seq_identity: int = 95) -> list:
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
    pdb_seq_hits = find_simular_chains(seq, parsed_pdb, seq_identity)

    return pdb_seq_hits

def seq_and_chain_search(input_pdb_path: str, binder_seq : str = None, target_protein_seq : str = None, binder_chains : str = None, target_protein_chains : str = None, order : str = "chains", seq_id : int = 95, sm_name : str = None, target_sm_chains : str = None, search_sm : str = "chains") -> dict:
    """
    This function is used to determine the order of filtering and what to filter out using the user input. The binder/target_protein are identical calls, just different inputs
    All the params from main are passed into here
    :return: dict of the different filtered binder/target_protein/sm
    """

    pdb_hits = {"binder":[],"target_protein":[],"sm":[]}

    if (binder_seq and binder_chains):
        if order == "seq":
            pdb_binder_seq_hits = seq_search(binder_seq, None, input_pdb_path, seq_identity=seq_id)
            pdb_hits["binder"] = chain_search(binder_chains, pdb_binder_seq_hits)
        else:
            pdb_binder_chain_hits = chain_search(binder_chains, None, input_pdb_path)
            pdb_hits["binder"] = seq_search(binder_seq, pdb_binder_chain_hits, None, seq_identity=seq_id)
    elif binder_seq:
        pdb_hits["binder"] = seq_search(binder_seq,None,input_pdb_path, seq_identity=seq_id)
    elif binder_chains:
        pdb_hits["binder"] = chain_search(binder_chains, None, input_pdb_path)
    else:
        print("ERROR no chain or seq selected for binder")
        exit(1)

    if (target_protein_seq and target_protein_chains):
        if order == "seq":
            pdb_target_protein_seq_hits = seq_search(target_protein_seq, None, input_pdb_path, seq_identity=seq_id)
            pdb_hits["target_protein"] = chain_search(target_protein_chains, pdb_target_protein_seq_hits)
        else:
            pdb_target_protein_chain_hits = chain_search(target_protein_chains, None, input_pdb_path)
            pdb_hits["target_protein"] = seq_search(target_protein_seq, pdb_target_protein_chain_hits, None, seq_identity=seq_id)
    elif target_protein_seq:

        pdb_hits["target_protein"] = seq_search(target_protein_seq,None,input_pdb_path, seq_identity=seq_id)
    elif target_protein_chains:

        pdb_hits["target_protein"] = chain_search(target_protein_chains, None, input_pdb_path)
    elif not sm_name and not target_sm_chains:
        print("ERROR no chain or seq selected for target_protein")
        exit(1)

    if target_sm_chains or sm_name:

        if target_sm_chains:
            target_sm_chains = arg_spliter(target_sm_chains)
        if sm_name:
            sm_name = arg_spliter(sm_name)
        pdb_hits["sm"] = sm_search(input_pdb_path,sm_name,target_sm_chains,search_sm)


    return pdb_hits


def write_pdbs(pdb_name : str, pdb_hits: dict, output_path : str) -> None:
    """
    used to write the output pdb file. Takes the info from the filtering step and outputs it into a pdb
    :param pdb_name: the name of the output pdb file
    :param pdb_hits: the dict of lines for the binder, target_protein, and sm use to make the output file
    :param output_path: the output file path directory
    :return: None
    """
    out_atom_index = 0
    out_chain = 64
    with open(output_path+"/"+pdb_name, "w") as output_pdb:
        ter_used_last = False
        for key in pdb_hits:
            out_chain+=1
            out_res_index = 0
            last_res_index = 0
            last_chain = "5"
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

                    out_atom_index += 1
                    output_pdb.write(line[0:6] + ((5 - len(str(out_atom_index))) * " ") + str(out_atom_index) + line[11:21] + chr(out_chain) + ((4 - len(str(out_res_index))) * " ") + str(out_res_index) + line[26:])


def main(input_pdb_path: str, output_path: str, binder_seq : str = None, target_protein_seq : str = None, binder_chains : str = None, anitgen_chains : str = None, order : str = None, seq_identity: int = 95, sm_name : str = None, target_sm_chains : str = None, search_sm : str = "chains") -> None:

    """
    Main function that calls the parsing/filtering and writing of the new pdb file(s)
    :return: None
    """
    pdb_to_write = seq_and_chain_search(input_pdb_path,binder_seq,target_protein_seq,binder_chains,anitgen_chains,order,seq_identity,sm_name,target_sm_chains,search_sm)
    pdb_name = os.path.basename(input_pdb_path).replace(".pdb","_cleaned.pdb")
    write_pdbs(pdb_name,pdb_to_write,output_path)


if __name__ == "__main__":

    """
    add arguments to arg parse
    see below for what each term does
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-S", "--binder_seqs", type=str, help="the sequence(s) of the binder used for identification", default=None)
    parser.add_argument("-s", "--target_protein_seqs", type=str, help="the sequence(s) of the target protein used for identification", default="")
    parser.add_argument("-L", "--target_sm_3_names", type=str, help="the three letter code of the small molecule(s)", default=None)
    parser.add_argument("-C", "--binder_chains", type=str, help="the chain(s) of the binder used for identification", default=None)
    parser.add_argument("-c", "--target_protein_chains", type=str, help="the chains of the target protein used for identification", default="")
    parser.add_argument("-l", "--target_sm_chains", type=str, help="the chain(s) the small molecule is/are on", default=None)
    parser.add_argument("-F", "--search_first_protein", type=str, help="if using both chains and sequences search, what should be used to be filtered first?",default="chains",choices=["chains","sequences"])
    parser.add_argument("-f", "--search_first_sm", type=str,help="if using both chains and name search for a small molecule, what should be used to be filtered first?", default="chains", choices=["chains", "name"])
    parser.add_argument("-d", "--seq_identity", type=int,help="sequence identity cutoff for finding similar chains",default=95)
    parser.add_argument("-i", "--input_path", type=str, help="path of the input pdb (file or directory)", required=True)
    parser.add_argument("-o", "--output_path", type=str, help="path of the output directory", required=True)
    args = parser.parse_args()


    #checking to see if a binder is present

    if not args.binder_seqs and not args.binder_chains:
        print("You must enter a binder chain or sequence or both!")
        exit(1)

    # checking to see if a sm or target_protein is present
    if not args.target_protein_seqs and not args.target_protein_chains and not args.target_sm_3_names and not args.target_sm_chains:
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

    # running the main method for a single file that is a .pdb
    if os.path.isfile(args.input_path):
        if pathlib.Path(args.input_path).suffix == ".pdb":
            main(args.input_path, args.output_path, args.binder_seqs, args.target_protein_seqs, args.binder_chains,args.target_protein_chains, args.search_first_protein, args.seq_identity, args.target_sm_3_names, args.target_sm_chains, args.search_first_sm)
        else:
            print("non .pdb file types are not supported!")
            exit(1)

    # running the main method for all pdbs in a directory
    else:
        for file in os.listdir(args.input_path):
            if ".pdb" in file:
                main(os.path.join(args.input_path,file), args.output_path, args.binder_seqs, args.target_protein_seqs, args.binder_chains, args.target_protein_chains, args.search_first_protein, args.seq_identity, args.target_sm_3_names, args.target_sm_chains, args.search_first_sm)

