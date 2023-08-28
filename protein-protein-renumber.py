from difflib import SequenceMatcher
import argparse
import os
import pathlib

VALID_3_LINE_STARTERS = {"ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I", "LYS":"K", "LEU":"L", "MET":"M",
                         "ASN":"N", "PHE":"F", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S", "THR":"T", "VAL":"V", "TRP":"W", "TYR":"T"}

def similar(string_a : str, string_b: str) -> int:

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
    parsed_pdb = []
    with open(input_pdb_path, "r") as inputPDB:
        for line in inputPDB:
            if "ATOM" in line[0:6] or "TER" in line[0:6] or "END" in line[0:6] and not "#" in line:
                parsed_pdb.append(line)
    return parsed_pdb


def arg_spliter(arg_str : str) -> list:
    if "," in arg_str:
        return arg_str.split(",")
    return [arg_str]


def chain_search(chains : str, parsed_pdb : list = None,input_pdb_path: str = None ) -> list:

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

    if input_pdb_path is not None:
        parsed_pdb = parse_pdb(input_pdb_path)

    seq = arg_spliter(seq)
    pdb_seq_hits = find_simular_chains(seq, parsed_pdb, seq_identity)

    return pdb_seq_hits

def seq_and_chain_search(input_pdb_path: str, binder_seq : str = None, antigen_seq : str = None, binder_chains : str = None, antigen_chains : str = None, order : str = "chain", seq_id : int = 95) -> dict:
    pdb_hits = {"binder":[],"antigen":[]}
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

    if (antigen_seq and antigen_chains):
        if order == "seq":
            pdb_antigen_seq_hits = seq_search(antigen_seq, None, input_pdb_path, seq_identity=seq_id)
            pdb_hits["antigen"] = chain_search(antigen_chains, pdb_antigen_seq_hits)
        else:
            pdb_antigen_chain_hits = chain_search(antigen_chains, None, input_pdb_path)
            pdb_hits["antigen"] = seq_search(antigen_seq, pdb_antigen_chain_hits, None, seq_identity=seq_id)
    elif antigen_seq:

        pdb_hits["antigen"] = seq_search(antigen_seq,None,input_pdb_path, seq_identity=seq_id)
    elif antigen_chains:

        pdb_hits["antigen"] = chain_search(antigen_chains, None, input_pdb_path)
    else:
        print("ERROR no chain or seq selected for antigen")
        exit(1)


    return pdb_hits



def write_pdbs(pdb_name : str, pdb_hits: dict, output_path : str) -> None:

    out_atom_index = 0
    out_chain = 64

    with open(output_path+"/"+pdb_name, "w") as output_pdb:
        ter_used_last = False
        for key in pdb_hits:
            out_chain+=1
            out_res_index = 0
            last_res_index = 0
            last_chain = "5"
            for line in pdb_hits[key]:
                if "TER" in line[0:6] and not ter_used_last:
                    ter_used_last = True
                    output_pdb.write(line)
                if "ATOM" in line[0:6]:
                    ter_used_last = False
                    chain = line[21]
                    res_index = int(line[22:26].strip())
                    if res_index != last_res_index or chain != last_chain:
                        out_res_index += 1
                        last_res_index = res_index
                        last_chain = chain

                    out_atom_index += 1
                    output_pdb.write("ATOM  " + ((5 - len(str(out_atom_index))) * " ") + str(out_atom_index) + line[11:21] + chr(out_chain) + ((4 - len(str(out_res_index))) * " ") + str(out_res_index) + line[26:])


def main(input_pdb_path: str, output_path: str, binder_seq : str = None, antigen_seq : str = None, binder_chains : str = None, anitgen_chains : str = None, order : str = None, seq_identity: int = 95) -> None:


    pdb_to_write = seq_and_chain_search(input_pdb_path,binder_seq,antigen_seq,binder_chains,anitgen_chains,order,seq_identity)
    pdb_name = os.path.basename(input_pdb_path).replace(".pdb","_cleaned.pdb")
    write_pdbs(pdb_name,pdb_to_write,output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-S", "--binder_seqs", type=str, help="the sequence of the binder used for identification", default=None)
    parser.add_argument("-s", "--antigen_seqs", type=str, help="the sequence of the antigen used for identification", default=None)
    parser.add_argument("-C", "--binder_chains", type=str, help="the chains of the binder used for identification", default=None)
    parser.add_argument("-c", "--antigen_chains", type=str, help="the chains of the antigen used for identification", default=None)
    parser.add_argument("-f", "--search_first", type=str, help="if using both chains and sequences search, what should be used to be filtered first?",default="chain",choices=["chains","sequences"])
    parser.add_argument("-d", "--seq_identity", type=int,help="if using both chains and sequences search, what should be used to be filtered first?",default=95)
    parser.add_argument("-i", "--input_path", type=str, help="path of the input pdb (file or directory)", required=True)
    parser.add_argument("-o", "--output_path", type=str, help="path of the output directory", required=True)
    args = parser.parse_args()
    if not args.binder_seqs and not args.binder_chains:
        print("You must enter a binder chain or sequence or both!")
        exit(1)
    if not args.antigen_seqs and not args.antigen_chains:
        print("You must enter an antigen chain or sequence or both!")
        exit(1)

    if not os.path.exists(args.input_path):
        print("the input path you entered does not exist!")
        exit(1)

    if not os.path.exists(args.output_path):
        print("the output path you entered does not exist!")
        exit(1)

    if os.path.isfile(args.input_path):
        if pathlib.Path(args.input_path).suffix == ".pdb":
            main(args.input_path, args.output_path, args.binder_seqs, args.antigen_seqs, args.binder_chains,args.antigen_chains, args.search_first, args.seq_identity)
        else:
            print("non .pdb file types are not supported!")
            exit(1)
    else:
        for file in os.listdir(args.input_path):
            if ".pdb" in file:
                main(os.path.join(args.input_path,file), args.output_path, args.binder_seqs, args.antigen_seqs, args.binder_chains, args.antigen_chains, args.search_first, args.seq_identity)




