import os
import numpy as np
from glob import glob
import argparse


def load_and_shuffle(particle_name, base_dir):
    """
    Loads all .data files inside <base_dir>/<particle_name>,
    shuffles rows consistently across all files, and returns a dict.
    """
    base_dir = os.path.abspath(base_dir)
    folder = os.path.join(base_dir, particle_name)

    if not os.path.isdir(folder):
        raise FileNotFoundError(f"Folder not found: {folder}")

    data_files = sorted(glob(os.path.join(folder, "*.data")))
    if not data_files:
        raise FileNotFoundError(f"No .data files found in {folder}")

    arrays = [np.loadtxt(f) for f in data_files]

    nrows = {arr.shape[0] for arr in arrays}
    if len(nrows) != 1:
        raise ValueError("All .data files must have the same number of rows.")

    n = arrays[0].shape[0]
    perm = np.random.permutation(n)

    shuffled = {
        os.path.basename(f): arr[perm]
        for f, arr in zip(data_files, arrays)
    }

    return shuffled


def save_shuffled(shuffled_dict, output_dir):
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    for filename, arr in shuffled_dict.items():
        outpath = os.path.join(output_dir, filename)
        np.savetxt(outpath, arr)
        print(f"Saved: {outpath}")


def main():
    parser = argparse.ArgumentParser(
        description="Shuffle .data files consistently across rows"
    )
    parser.add_argument("--particle", required=True,
                        help="Particle name (e.g. Positron)")
    parser.add_argument("--input-dir", required=True,
                        help="Absolute path to base input directory")
    parser.add_argument("--output-dir", required=True,
                        help="Absolute path to output directory")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed (e.g. Condor PROCESS)")

    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)

    shuffled = load_and_shuffle(args.particle, args.input_dir)
    save_shuffled(shuffled, args.output_dir)


if __name__ == "__main__":
    main()
