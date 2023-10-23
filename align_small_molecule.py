import argparse
import pymol
import os
import glob
from helper_functions import create_directory
def align_small(path, output, chain, rmsd_threshold):
    pymol.finish_launching(['pymol', '-qc'])  # Launch PyMOL in quiet and command-line mode
    pdb_files = glob.glob(path + "/*.pdb")
    references = []
    for i, file in enumerate(pdb_files):
        pymol.cmd.load(file, str(i))
        atom_count = pymol.cmd.count_atoms(f"chain {chain} and model {str(i)}")

        if atom_count == 0:
            print(f"Skipping {file} as it does not contain chain {chain}")
            continue
        if i == 0:
            ref_name = "reference_1"
            print(file)
            pymol.cmd.select("reference_1", f"chain {chain} and model {str(i)}")
            references.append(ref_name)
            pdb_name = os.path.basename(file).replace(".pdb", f"{ref_name}.pdb")
            pymol.cmd.save(os.path.join(output, pdb_name), f"model {str(i)}")
            continue
        pymol.cmd.select("current", f"chain {chain} and model {str(i)}")

        found_cluster = False

        for reference in references:
            rmsd = pymol.cmd.align(str(reference), "current")
            rmsd_value = rmsd[0]
            print(rmsd_value)
            if rmsd_value <= rmsd_threshold:
                pdb_name = os.path.basename(file).replace(".pdb", f"{reference}.pdb")
                pymol.cmd.save(os.path.join(output, pdb_name), f"model {str(i)}")
                found_cluster = True
                break
        if not found_cluster:
            num_cluster = len(references) + 1
            ref_name_new = f"reference_{num_cluster}"
            pymol.cmd.select(ref_name_new,  f"chain {chain} and model {str(i)}")
            references.append(ref_name_new)
            pdb_name = os.path.basename(file).replace(".pdb", f"{ref_name_new}.pdb")
            pymol.cmd.save(os.path.join(output, pdb_name), f"model {str(i)}")

    pymol.cmd.reinitialize()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--chain_to_align", type=str, help="chain identifier for chains to align", default="A", required=True)
    parser.add_argument("-T", "--rmsd_threshold", type=float, help="RMSD Threshold Difference", default=0.3)
    parser.add_argument("-i", "--input_path", type=str, help="path of the input directory", required=True)
    parser.add_argument("-o", "--output_path", type=str, help="path of the output directory", required=True)

    args = parser.parse_args()
    create_directory(args.output_path)
    align_small(args.input_path, args.output_path, args.chain_to_align, args.rmsd_threshold)
