import os
from tqdm import tqdm

from pmt_locator import *
import awkward as ak


class DataFormater:
    def __init__(self, 
                 sig_dir, 
                 rad_dir,
                 noi_dir,
                 out_dir,
                 events_per_file = 1000,
                 apply_filter = True,
                 data_augmentation = False,
                 generate_darknoise = False):
        
        self.sig_dir = sig_dir
        self.rad_dir = rad_dir
        self.noi_dir = noi_dir
        self.out_dir = out_dir
        self.events_per_file = events_per_file
        self.apply_filter = apply_filter
        self.data_augmentation = data_augmentation
        self.generate_darknoise = generate_darknoise

        self.num_sig = get_files_in_directory(sig_dir)
        self.num_rad = get_files_in_directory(rad_dir)
        self.num_noi = get_files_in_directory(noi_dir)

        self.num_batches = self.num_sig // self.events_per_file

        assert (self.num_sig % self.events_per_file == 0), "Warning: Number of files is not cleanly divideable by the number of events per output file."

        if not data_augmentation:
            assert (self.num_sig == self.num_rad) and (self.num_sig == self.num_noi), "Error: Number of signal files unequals number of background files."
        
        print("=======================================")
        print(f"Signal files from {self.sig_dir}")
        print(f"Radioactive decay files from {self.rad_dir}")
        print(f"Dark noise files from {self.noi_dir}")
        print(f"Output directory {self.out_dir}\n")
        print(f"Number of infiles/events: {self.num_sig}")
        print(f"Number of events per outfile: {self.events_per_file}")
        print(f"Number of batches: {self.num_batches}")
        print("=======================================")        

    def load_file(self, path):
        """Loads data from a given path.
        Args:
            path ('str'): Path to file.

        Returns:
            np.2darray: Returns 2D np.array with time in 1st column and PMT in 2nd column.
        """
        #Dark noise files are saved as structured numpy binaries with np.savez() and must be loaded with np.load()
        if "noi" in path:
            data = np.load(path)
        #For binary files from Geant4 np.fromfile is used.
        else:
            data = np.fromfile(path, dtype=[('time', 'f4'), ('pmt', 'u1')])
        return data
    
    def load_pmt_positions(self):
        pmt_locator = PMTLocator()
        self.pmt_lookup = pmt_locator.get_pmt_location_cartesian()
    
    def combine_data(self):
        """Combines signal and background (radioactive decay + dark noise) data and sorts the entries by time.
        Data is stored as an awkward array and can be further processed along the line.
        ToDo: Batch load data, operations like filtering, electronics time smear, etc are applied on batches
        """

        comb_time = []
        comb_pmt = []
        comb_pos = []

        # loop through all files in the signal directory
        for i in self.batch_idx:
            sig_path = os.path.join(self.sig_dir, f"sig_{i}.bin")
            rad_path = os.path.join(self.rad_dir, f"rad_{i}.bin")
            noi_path = os.path.join(self.noi_dir, f"noi_{i}.npz")

            sig = self.load_file(sig_path)
            rad = self.load_file(rad_path)
            noi = self.load_file(noi_path)

            # combine time
            c_time = np.concatenate((sig["time"], rad["time"], noi["time"]))
            c_pmt = np.concatenate((sig["pmt"], rad["pmt"], noi["pmt"]))

            # Get PMT locations using lookup table
            c_pos = np.array([self.pmt_lookup[pmt_id] for pmt_id in c_pmt])

            # sort arrays by arrival time
            c_time, c_pmt, c_pos = sort_by(c_time, c_pmt, c_pos)

            comb_time.append(c_time)
            comb_pmt.append(c_pmt)
            comb_pos.append(c_pos)

        self.time = ak.Array(comb_time)
        self.pmt = ak.Array(comb_pmt)
        self.pos = ak.Array(comb_pos)

        # combine time and position information (events_per_file * var * 4 * float)
        self.hit = ak.concatenate([self.time[..., np.newaxis], self.pos], axis = -1)

    def save(self, idx):
        out_path = os.path.join(self.out_dir, f"data_{idx}.pq")
        ak.to_parquet(self.hit, out_path)

    def reset(self):
        self.time = 0
        self.pmt = 0
        self.pos = 0
        self.hit = 0

    def run(self):
        self.load_pmt_positions()
        idx = np.arange(self.num_sig) # index array

        for b in tqdm(range(self.num_batches)):
            self.batch_idx = idx[b*self.events_per_file:(b+1)*self.events_per_file]
            self.combine_data()
            self.save(idx = b)
            self.reset()

def get_files_in_directory(dir):
    """Returns number of files in directory. Does not count sub-directories.

    Args:
        dir (str): Path to directory.

    Returns:
        int: Number of files in directory.
    """
    return len([name for name in os.listdir(dir) if os.path.isfile(os.path.join(dir,name))])       


def sort_by(key_array, *arrays):
    """
    Sort multiple arrays according to the sorting order of key_array.

    Args:
        key_array (np.ndarray): Array to determine sort order.
        *arrays (np.ndarray): Arrays to be sorted in the same order as key_array.

    Returns:
        Tuple: Sorted key_array and sorted versions of *arrays.
    """
    idx = key_array.argsort()
    return (key_array[idx],) + tuple(arr[idx] for arr in arrays)

""" ToDo: Filter and Smearing"""