import argparse
import pymol
import os
import glob


def align(path, output, chain, rmsd_threshold):
    pymol.finish_launching(['pymol', '-qc'])  # Launch PyMOL in quiet and command-line mode
    pdb_files = glob.glob(path + "/*.pdb")
    for i, file in enumerate(pdb_files):
        pymol.cmd.load(file, str(i))
        pymol.cmd.h_add( "all" )
        atom_count = pymol.cmd.count_atoms(f"chain {chain} and model {str(i)}")
        if atom_count == 0:
            print(f"Skipping {file} as it does not contain chain {chain}")
            continue
        if i == 0:
            pymol.cmd.select("reference", f"chain {chain} and model {str(i)}")
            continue
        pymol.cmd.select("current", f"chain {chain} and model {str(i)}")
        rmsd = pymol.cmd.cealign("reference", "current")
        rmsd_value  = rmsd["RMSD"]

        if rmsd_value <= rmsd_threshold:
            print("saved")
            pdb_name = os.path.basename(file).replace(".pdb", "_aligned.pdb")
            pymol.cmd.save(os.path.join(output, pdb_name), f"model {str(i)}")
            pymol.cmd.delete(f"model {str(i)}")
            pymol.cmd.remove("current")
        else:
            print(f"Alignment RMSD {rmsd_value} for {file} exceeds threshold {rmsd_threshold}, so it is not included in the aligment")
    pymol.cmd.reinitialize()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--chain_to_align", type=str, help="chain identifier for chains to align", default="A", required=True)
    parser.add_argument("-T", "--rmsd_threshold", type=float, help="RMSD Threshold for filtering.", default=3)
    parser.add_argument("-i", "--input_path", type=str, help="path of the input directory", required=True)
    parser.add_argument("-o", "--output_path", type=str, help="path of the output directory", required=True)

    args = parser.parse_args()

    align(args.input_path, args.output_path, args.chain_to_align, args.rmsd_threshold)
