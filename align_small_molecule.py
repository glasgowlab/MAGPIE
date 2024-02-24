import argparse
import pymol
import os
import glob
from ast import literal_eval

def create_directory(directory_path):
    if os.path.exists(directory_path):
        return f"The directory '{directory_path}' already exists."
    else:
        os.makedirs(directory_path)
def align_small(path, output, chain, rmsd_threshold,pair_fit):
    pymol.finish_launching(['pymol', '-qc'])  # Launch PyMOL in quiet and command-line mode
    pymol.cmd.feedback("disable", "all", "everything")
    pdb_files = glob.glob(path + "/*.pdb")
    references = []
    command = ""
    sm_atom_names = {"ref":{},"cur":{}}
    for i, file in enumerate(pdb_files):
        print(f"Processing file number {i + 1} out of {len(pdb_files)}")
        pymol.cmd.load(file, str(i))
        pymol.cmd.h_add("all")
        atom_count = pymol.cmd.count_atoms(f"chain {chain} and model {str(i)}")

        if atom_count == 0:
            print(f"Skipping {file} as it does not contain chain {chain}")
            continue
        if i == 0:
            ref_name = "reference_1"
            create_directory(os.path.join(output,ref_name))
            pymol.cmd.select("reference_1", f"chain {chain} and model {str(i)}")
            references.append(ref_name)
            pdb_name = os.path.basename(file).replace(".pdb", f"{ref_name}.pdb")
            pymol.cmd.save(os.path.join(output, ref_name, pdb_name), f"model {str(i)}")
            continue
        pymol.cmd.select("current", f"chain {chain} and model {str(i)}")

        found_cluster = False
        lowest_rmsd = 99999
        lowest_reference= ''
        if pair_fit:
            sm_atom_names["cur"] = {}
            pymol.cmd.iterate('current and not element H', 'cur[name]=index', space=sm_atom_names)
        for reference in references:
            if pair_fit:
                sm_atom_names["ref"] = {}
                pymol.cmd.iterate(f'{reference} and not element H', 'ref[name]=index', space=sm_atom_names)
                command = ""
                count = 0
                for atom_name in sm_atom_names["ref"]:
                    if atom_name in sm_atom_names["cur"]:
                        count+=1
                        command += f"'current and index {sm_atom_names['cur'][atom_name]}', '{reference} and index {sm_atom_names['ref'][atom_name]}', "
                if count < 3:
                    print(f"Using normal align as the atom names or number differ to much for {file}")
                    rmsd = pymol.cmd.align("current", str(reference))
                    rmsd_value = rmsd[0]
                else:
                    rmsd_value = float(pymol.cmd.pair_fit(*eval(command)))
            else:
                rmsd = pymol.cmd.align("current", str(reference))
                rmsd_value = rmsd[0]
            if rmsd_value <= rmsd_threshold and rmsd_value < lowest_rmsd:
                    found_cluster = True
                    lowest_rmsd = rmsd_value
                    lowest_reference = reference
        if found_cluster:
            pymol.cmd.align("current",str(lowest_reference))
            pdb_name = os.path.basename(file).replace(".pdb", f"{lowest_reference}.pdb")
            pymol.cmd.save(os.path.join(output, lowest_reference, pdb_name), f"model {str(i)}")
            pymol.cmd.remove(str(i))
            pymol.cmd.delete(str(i))
        else:
            num_cluster = len(references) + 1
            ref_name_new = f"reference_{num_cluster}"
            pymol.cmd.select(ref_name_new,  f"chain {chain} and model {str(i)}")
            references.append(ref_name_new)
            create_directory(os.path.join(output,ref_name_new))
            pdb_name = os.path.basename(file).replace(".pdb", f"{ref_name_new}.pdb")
            pymol.cmd.save(os.path.join(output, ref_name_new, pdb_name), f"model {str(i)}")

    pymol.cmd.reinitialize()

def create_directory(directory_path):
    if os.path.exists(directory_path):
        return f"The directory '{directory_path}' already exists."
    else:
        os.makedirs(directory_path)



if __name__ == "__main__":

    # Remove structures that were aligned already

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--chain_to_align", type=str, help="chain identifier for chains to align", default="A", required=True)
    parser.add_argument("-T", "--rmsd_threshold", type=float, help="RMSD Threshold Difference", default=0.3)
    parser.add_argument("-i", "--input_path", type=str, help="path of the input directory", required=True)
    parser.add_argument("-o", "--output_path", type=str, help="path of the output directory", required=True)
    parser.add_argument("-p", "--pair_fit", type=str, help="use pair_fit instead", required=False,choices=["True","False"],default="False")

    args = parser.parse_args()
    create_directory(args.output_path)
    align_small(args.input_path, args.output_path, args.chain_to_align, args.rmsd_threshold,args.pair_fit == "True")
