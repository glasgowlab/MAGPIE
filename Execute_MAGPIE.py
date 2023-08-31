import argparse

import time

import helper_functions
import pandas as pd
import glob
import multiprocessing

import sequence_logo_main


def process_files(list_of_paths, binder, thread_num):
    binder_chain_cordinates = []

    for file in list_of_paths:
        binder_chain_cordinates += helper_functions.extract_info_pdb(file, binder)

    data_frame_binders = pd.DataFrame(binder_chain_cordinates)

    data_frame_binders.to_json(f"temp/binders_{thread_num}.json")


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




def merge_json_files_to_dataframe(file_pattern):
    dataframes = []

    # Get a list of JSON files matching the file pattern
    json_files = glob.glob(file_pattern)

    # Iterate through each JSON file and read it into a DataFrame
    for json_file in json_files:
        with open(json_file, 'r') as f:
            data = pd.read_json(f)
            dataframes.append(data)

    # Concatenate the DataFrames into a single DataFrame
    merged_dataframe = pd.concat(dataframes, ignore_index=True)

    return merged_dataframe

if __name__ == "__main__":
    start_time = time.time()
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    argparser.add_argument("--path_to_pdbs", type=str, default='Ligand_Example',
                           help="path to pdbs, examples: 'Ligand_Example' or 'Protein Example'")

    argparser.add_argument("--target_chain", type=str, default="X",
                           help="Chain id to graph around, target chain")
    argparser.add_argument("--binder_chain", type=str, default="A", help="Chain id that binds to target chain.")
    argparser.add_argument("--is_ligand",  type=str, default= "False",
                           help="Is the target chain a ligand? ")
    argparser.add_argument("--distance", type=int, default=8,
                           help="distance in A to graph around target chain")
    argparser.add_argument("--to_plot", type=str, default="",
                           help="Amino acids to generate MAGPIE graphs.")
    argparser.add_argument("--combined_flag",   type=str, default= "False",
                           help="Show individual plots for residues in to_plot or only the combined plot? (True or False)")
    argparser.add_argument("--threads", type=int, default=1,
                           help="How many threads to use?")
    args = argparser.parse_args()
    input_folder = args.path_to_pdbs

    list_of_binders = glob.glob(input_folder + '/*.pdb')
    chunks = make_chunks(list_of_binders, args.threads)
    threads = []
    for thread_num in range(args.threads):
        current_thread = multiprocessing.Process(target=process_files, args=(
            chunks[thread_num],  # PDB Files
            args.binder_chain,
            thread_num
        ))
        threads.append(current_thread)
        current_thread.start()

    for t in threads:
        t.join()

    if args.is_ligand == "True":
        is_ligand = True
    else:
        is_ligand = False
    if args.combined_flag == "True":
        combined_flag = True
    else:
        combined_flag= False
    binder_df = merge_json_files_to_dataframe("temp/*.json")
    to_plot = args.to_plot
    if is_ligand:
        plot_list = [str(x) for x in to_plot.split(",")]
    else:
        plot_list = [int(x) for x in to_plot.split(",")]
    sequence_logo_main.sequence_logos(list_of_binders[0], args.target_chain, binder_df, plot_list,is_ligand,combined_flag,args.distance)

