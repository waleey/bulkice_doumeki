import os
import numpy as np
from glob import glob

def load_and_shuffle(particle_name, base_dir="."):
    """
    Loads all .data files inside <base_dir>/<particle_name>,
    shuffles rows consistently across all files, and returns a dict.
    """
    folder = os.path.join(base_dir, particle_name)
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"Folder not found: {folder}")

    # Find .data files
    data_files = sorted(glob(os.path.join(folder, "*.data")))
    if not data_files:
        raise FileNotFoundError("No .data files found in folder.")

    # Load all files into a list of arrays
    arrays = [np.loadtxt(f) for f in data_files]

    # Check all files have the same number of rows
    nrows = {arr.shape[0] for arr in arrays}
    if len(nrows) != 1:
        raise ValueError("All .data files must have the same number of rows.")

    n = arrays[0].shape[0]

    # Create a single random permutation
    perm = np.random.permutation(n)

    # Apply the permutation to each array
    shuffled = {os.path.basename(f): arr[perm] for f, arr in zip(data_files, arrays)}

    return shuffled


def save_shuffled(shuffled_dict, output_dir):
    """
    Saves the shuffled arrays back to disk inside output_dir.
    """
    os.makedirs(output_dir, exist_ok=True)

    for filename, arr in shuffled_dict.items():
        outpath = os.path.join(output_dir, filename)
        np.savetxt(outpath, arr)
        print(f"Saved: {outpath}")


if __name__ == "__main__":
    particle = "Positron"  # change or pass via argparse
    shuffled = load_and_shuffle(particle)

    # Example: save to a new folder
    save_shuffled(shuffled, f"{particle}")

