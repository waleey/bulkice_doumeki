import os
import tarfile
from tqdm import tqdm
from scipy.stats import truncnorm, gumbel_r
from pmt_locator import *
import awkward as ak


class DataFormater:
    def __init__(self, 
                 sig_dir, 
                 rad_dir,
                 noi_dir,
                 out_dir,
                 events_per_out = 1000,
                 files_per_tar = 10_000,
                 pmt_transit_mode = None,
                 trigger_window = 20, # in ns
                 trigger_pmts = 2,
                 load_mode = "tar",
                 verbose = False):
        
        self.sig_dir = os.path.join(sig_dir, load_mode)
        self.rad_dir = os.path.join(rad_dir, load_mode)
        self.noi_dir = os.path.join(noi_dir, load_mode)
        self.out_dir = out_dir

        self.events_per_out = events_per_out
        self.files_per_tar = files_per_tar
        self.pmt_transit_mode = pmt_transit_mode
        self.trigger_window = trigger_window
        self.trigger_pmts = trigger_pmts
        self.load_mode = load_mode
        self.verbose = verbose

        self.pmt_transit_time = 43 # in ns
        self.pmt_transit_time_spread_norm = 1.854 # in ns, see M. Unland PhD thesis, mean of Gauss sigma from Tab. 6.1
        self.pmt_transit_time_spread_fisher_tippett = 1.944 # in ns, see M. Unland PhD thesis, mean of Fisher-Tippet sigma from Tab. 6.1
        
        if self.verbose:
            print("=======================================")
            print(f"Signal files from {self.sig_dir}")
            print(f"Radioactive decay files from {self.rad_dir}")
            print(f"Dark noise files from {self.noi_dir}")
            print(f"Output directory {self.out_dir}")
            print(f"Number of events per outfile: {self.events_per_out}")
            print(f"Trigger setting: dt_trig = {self.trigger_window} ns, n_trig = {self.trigger_pmts}")
            print("=======================================")        

    def pmt_transit_distribution(self):

        loc = self.pmt_transit_time

        if self.pmt_transit_mode == "truncnorm":
            scale = self.pmt_transit_time_spread_norm
            a_trunc = -10 # in ns
            b_trunc = +100 # in ns
            a, b = (a_trunc - loc) / scale, (b_trunc - loc) / scale

            self.pmt_transit_dist = truncnorm(a, b, loc, scale)

        elif self.pmt_transit_mode == "fisher_tippett":
            scale = self.pmt_transit_time_spread_fisher_tippett

            self.pmt_transit_dist = gumbel_r(loc, scale)      

        else:
            raise ValueError(f"{self.pmt_transit_mode} not a valid PMT transit time distribution.")


    def load_file(self, path, mode = "tar"):
        """Wrapper for raw file and tar file loader.

        Args:
            path ("str"): Path to raw/tar file.
            mode (str, optional): Data loading mode. Defaults to "tar".

        Returns:
            np.2darray: Output of load_raw or load_tar depending on mode.
        """
        if mode == "tar":
            return self.load_tar(path)
        else:
            return self.load_raw(path)

    def load_raw(self, path):
        """Loads data from a given path.
        Args:
            path ('str'): Path to file.

        Returns:
            np.2darray: Returns 2D np.array with time in 1st column and PMT in 2nd column.
        """
        return np.fromfile(path, dtype=[('time', 'f8'), ('pmt', 'u1')])
    
    def load_tar(self, tar_path):
        """Loads single tar file and extracts single files.

        Args:
            tar_path ('str'): Path to tar file.

        Returns:
            np.2darray: Returns 2D np.array with time in 1st column and PMT in 2nd column.
        """

        mode = tar_path.split("/")[-1].split("_")[0]

        with tarfile.open(tar_path, "r:") as tar:
            f = tar.extractfile(f"{mode}_{self.file_idx}.bin")
            if f is not None:
                binary_data = f.read()
                data = np.frombuffer(binary_data, dtype=[('time', 'f8'), ('pmt', 'u1')])
                return data
            else:
                raise ValueError("f is None") 
    
    def load_pmt_positions(self):
        pmt_locator = PMTLocator()
        self.pmt_lookup = pmt_locator.get_pmt_location_cartesian()
    
    def combine_event(self, i):
        """Combines signal and background (radioactive decay + dark noise) per event sorts the entries by time.
        """
        self.file_idx = i

        # data path
        if self.load_mode == "raw":
            self.tar_idx = None
            sig_path = os.path.join(self.sig_dir, f"sig_{self.file_idx}.bin")
            rad_path = os.path.join(self.rad_dir, f"rad_{self.file_idx}.bin")
            noi_path = os.path.join(self.noi_dir, f"noi_{self.file_idx}.bin")
        else:
            self.tar_idx = i // 10_000 # batch index
            sig_path = os.path.join(self.sig_dir, f"sig_10k_batch_{self.tar_idx}.tar")
            rad_path = os.path.join(self.rad_dir, f"rad_10k_batch_{self.tar_idx}.tar")
            noi_path = os.path.join(self.noi_dir, f"noi_10k_batch_{self.tar_idx}.tar")

        # load
        self.sig = self.load_file(sig_path, mode = self.load_mode)
        self.rad = self.load_file(rad_path, mode = self.load_mode)
        self.noi = self.load_file(noi_path, mode = self.load_mode)

        # combine
        # src: sig = 1, bkg (rad, noi) = 0
        c_time = np.concatenate((self.sig["time"], self.rad["time"], self.noi["time"])) # in ns
        c_pmt = np.concatenate((self.sig["pmt"], self.rad["pmt"], self.noi["pmt"]))
        c_src = np.concatenate((np.ones(self.sig.size), np.zeros(self.rad.size), np.zeros(self.noi["time"].size))).astype(np.uint8)

        # add transit time
        if self.pmt_transit_mode is not None:
            c_time = self.pmt_transit(c_time)

        # sort arrays by arrival time
        c_time, c_pmt, c_src = sort_by(c_time, c_pmt, c_src)

        return c_time, c_pmt, c_src

    def filter(self, time, pmt, src):
        """Applies filter function.

        Args:
            time (np.array): sorted times
            pmt (np.array): sorted pmts
            src (np.array): sorted sources

        Returns:
            tuple of np.arrays: filtered time, pmt, src array
        """
        f_time, f_pmt, f_src = [], [], []

        i = 0
        n = len(time)

        while i < n:
            window_end = time[i] + self.trigger_window
            j = i
            while j < n and time[j] < window_end:
                j += 1

            window_pmts = pmt[i:j]
            window_src = src[i:j]
            num_hit = j - i
            num_pmt = len(np.unique(window_pmts))
            unique_src = np.unique(window_src)

            if num_hit > 1 and num_pmt >= self.trigger_pmts:
                f_time.append(time[i:j])
                f_pmt.append(window_pmts)
                if len(unique_src) == 2:
                    f_src.append(1)
                else:
                    f_src_val = unique_src[0]
                    f_src.append(0 if f_src_val == 0 else 2)

            i = j  # move pointer ahead

        if f_time:
            return np.concatenate(f_time), np.concatenate(f_pmt), np.array(f_src)
        else:
            return None, None, None


    def batch(self):
        """Calls combine_event and filter function. Loops through batch size.
        """

        time, pmt, src = [], [], []
        sig_eff, bkg_eff = [], []

        # loop through all files in the signal directory
        for i in tqdm(self.batch_idx, desc="Combining and Filtering"):
            
            c_time, c_pmt, c_src = self.combine_event(i)
            f_time, f_pmt, f_src = self.filter(c_time, c_pmt, c_src)

            time.append(f_time)
            pmt.append(f_pmt)
            src.append(f_src)

        self.time = ak.Array(time)
        self.pmt = ak.Array(pmt)
        self.src = ak.Array(src)

        self.ntime = self.normalize_time()

        self.sig_eff = ak.Array(sig_eff)
        self.bkg_eff = ak.Array(bkg_eff)

        # combine time and pmt information (events_per_out * var * 2 * float)
        # position and direction are added in the data loader to keep memory consumption low
        self.hit = ak.concatenate([self.ntime[..., np.newaxis], self.pmt[..., np.newaxis]], axis = -1)

    def normalize_time(self):

        m_time = ak.mean(self.time) # mean time
        n_time = self.time - m_time # shift time relative to mean
        n_time = n_time/1E9 # rescale from ns to s

        return n_time
    
    def pmt_transit(self, time_arrival):
        """Simulates the PMT transit by adding a transit time and pulling a transit time spread from a truncated normal distribution.

        Args:
            time_arrival (np.array): Arrival time of photons on the PMT. Output of bulkice_doumeki.

        Returns:
            _type_: _description_
        """
        # add transit time and jitter, loc = self.pmt_transit_time
        transit_time_spread = self.pmt_transit_dist.rvs(size = time_arrival.size) # in ns
        time_transit = time_arrival + transit_time_spread
        # round to 1 ns
        time_transit = np.round(time_transit, 0)

        return time_transit
    
    def save(self, idx):
        """Save hit data (time and positions) to parquet.

        Args:
            idx (int): Batch index.
        """
        out_path = os.path.join(self.out_dir, f"data_{idx}.pq")
        ak.to_parquet(self.hit, out_path)

    def reset(self):
        """Explicitly reset parameters for batch mode. 
        """
        self.time = 0
        self.pmt = 0
        self.src = 0
        self.hit = 0

    def run(self, idx_low, idx_high):
        """Main function. Loads data, combines it, performs operations on the data and saves it in a parquet format. 
        """

        # initiate the pmt transite time distribution
        if self.pmt_transit_mode is not None:
            self.pmt_transit_distribution()

        #self.load_pmt_positions()
        self.idx_low = idx_low
        self.idx_high = idx_high
        self.idx = np.arange(self.idx_low, self.idx_high) # complete index array of run

        assert ((self.idx_high - self.idx_low) % self.events_per_out == 0), f"Warning: File range ({self.idx_low} - {self.idx_high}) is not cleanly divideable by the number of events per output file {self.events_per_outfile}."
        self.num_batches = (self.idx_high - self.idx_low) // self.events_per_out

        # batch mode, combine data is only applied to batch_idx
        for b in range(self.num_batches):
            self.batch_idx = self.idx[b*self.events_per_out:(b+1)*self.events_per_out] # index array of batch
            if self.verbose: print(f"Output batch {b+1} / {self.num_batches} ... [{self.batch_idx[0]}, {self.batch_idx[-1]}])")
            self.batch()
            self.save(idx = b)
            self.reset()

    def filter_test(self, time, pmt, src):
        """Applies filter function.

        Args:
            time (np.array): sorted times
            pmt (np.array): sorted pmts
            src (np.array): sorted sources

        Returns:
            tuple of np.arrays: filtered time, pmt, src array
        """
        
        f_time = [] # time filtered list
        t_max = []
        f_pmt = [] # pmt filtered list
        f_src = [] # source filtered list

        f_ns = 0 # number of signal hits after filter
        f_nb = 0 # number of background hits after filter

        f_sing = 0 # number of hits that cluster as singlets
        f_mult = 0 # number of hits that cluster as multiplets
    
        while len(time) > 0:

            # temporal coincidence
            mask_time = time - time[0] < self.trigger_window
            num_hit = np.sum(mask_time)
            # spatial coincidence, require at least two PMTs
            num_pmt = len(np.unique(pmt[mask_time]))
            # trigger source (bkg only: 0, sig+bkg: 1, sig only: 2)
            num_src = np.unique(src[mask_time])
            
            if num_hit == 1:
                f_sing += 1
            else:
                f_mult += np.sum(mask_time)

            if num_hit > 1 and num_pmt >= self.trigger_pmts:
                f_time.append(time[mask_time])
                f_pmt.append(pmt[mask_time])
                t_max.append(time[mask_time][-1]-time[mask_time][0])

                # mark trigger as bkg only (0), mixed sig, bkg (1), sig only (2)
                if len(num_src) == 2: # case sig + bkg
                    f_src.append(1)
                else:
                    if num_src == 0:
                        f_src.append(0)
                    elif num_src == 1:
                        f_src.append(2)
                
                # number of sig, bkg hits that survived filtering
                f_ns += np.sum(src[mask_time])
                f_nb += np.sum(np.logical_not(src[mask_time]))

            time = time[num_hit:]
            pmt = pmt[num_hit:]
            src = src[num_hit:]

        if len(f_time) > 0:
            return np.concatenate(f_time), np.concatenate(f_pmt), np.array(f_src), f_ns, f_nb, f_sing, f_mult, t_max
        else:
            return None, None, None, 0, 0, 0, 0, None

    def combine_test(self, i):
        """Combines signal and background (radioactive decay + dark noise) per event sorts the entries by time.
        """
        self.file_idx = i

        # data path
        if self.load_mode == "raw":
            sig_path = os.path.join(self.sig_dir, f"sig_{self.file_idx}.bin")
            rad_path = os.path.join(self.rad_dir, f"rad_{self.file_idx}.bin")
            noi_path = os.path.join(self.noi_dir, f"noi_{self.file_idx}.bin")
        else:
            self.tar_idx = i // 10_000 # batch index
            sig_path = os.path.join(self.sig_dir, f"sig_10k_batch_{self.tar_idx}.tar")
            rad_path = os.path.join(self.rad_dir, f"rad_10k_batch_{self.tar_idx}.tar")
            noi_path = os.path.join(self.noi_dir, f"noi_10k_batch_{self.tar_idx}.tar")

        # load
        self.sig = self.load_file(sig_path, mode = self.load_mode)
        self.rad = self.load_file(rad_path, mode = self.load_mode)
        self.noi = self.load_file(noi_path, mode = self.load_mode)

        # combine
        # src: sig = 1, bkg (rad, noi) = 0
        c_time = np.concatenate((self.sig["time"], self.rad["time"], self.noi["time"])) # in ns
        c_pmt = np.concatenate((self.sig["pmt"], self.rad["pmt"], self.noi["pmt"]))
        c_src = np.concatenate((np.ones(self.sig.size), np.zeros(self.rad.size), np.zeros(self.noi["time"].size))).astype(np.uint8)

        # add transit time
        if self.pmt_transit_mode is not None:
            c_time = self.pmt_transit(c_time)

        # sort arrays by arrival time
        c_time, c_pmt, c_src = sort_by(c_time, c_pmt, c_src)

        return c_time, c_pmt, c_src

    def batch_test(self):
        """Calls combine_event and filter function. Loops through batch size.
        """

        time, pmt, src = [], [], []
        sig_eff, mult_eff, bkg_eff = [], [], []
        t_max = []

        # loop through all files in the signal directory
        for i in tqdm(self.batch_idx, desc="Combining and Filtering"):
            
            c_time, c_pmt, c_src = self.combine_test(i)

            sort = np.argsort(self.sig["time"])
            _,_,n_src,_,_,n_sing,n_mult, _ = self.filter_test(self.sig["time"][sort], self.sig["pmt"][sort], np.ones(self.sig["time"].size))
            
            f_time, f_pmt, f_src, f_ns, f_nb, _,_,tmax = self.filter_test(c_time, c_pmt, c_src)

            time.append(f_time.value) # to value, ak.Array(list) cannot handle units
            pmt.append(f_pmt)
            src.append(f_src)
            
            # original, un-filtered number of sig, bkg
            ns = np.sum(c_src)
            nb = np.sum(np.logical_not(c_src))

            # sig, mult, bkg efficiency: events that pass filter / all events 
            sig_eff.append(f_ns/ns)
            if n_mult > 0 : mult_eff.append(f_ns/n_mult)
            else: mult_eff.append(0)
            bkg_eff.append(f_nb/nb)
            t_max.append(tmax)
        print(time)
        self.time = ak.Array(time)
        self.pmt = ak.Array(pmt)
        self.src = ak.Array(src)

        # To Do
        self.ntime = self.normalize_time()
        self.t_max = ak.Array(t_max)

        self.sig_eff = ak.Array(sig_eff)
        self.mult_eff = ak.Array(mult_eff)
        self.bkg_eff = ak.Array(bkg_eff)

        # combine time and pmt information (events_per_out * var * 2 * float)
        # position and direction are added in the data loader to keep memory consumption low
        self.hit = ak.concatenate([self.ntime[..., np.newaxis], self.pmt[..., np.newaxis]], axis = -1)

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

