import numpy as np
import scipy.constants as c
from scipy.interpolate import InterpolatedUnivariateSpline
import pandas as pd
import urllib
from io import StringIO
import os
import re
import json
from tqdm import tqdm

from pmt_locator import PMTLocator
pmtlocator = PMTLocator()
ploc_car = pmtlocator.get_pmt_location_cartesian()
ploc_sph = pmtlocator.get_pmt_location_spherical()


import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import matplotlib.patches as patches


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 12


class Decay:

    def __init__(self, dir, isolist, num_particles, mode = "exp", delimiter = ","):

        self.dir = dir
        self.isolist = isolist
        self.sim_num_particles = num_particles
        self.mode = mode
        self.delimiter = delimiter

        
    def get_sim_path(self):
        self.sim_path = self.dir + "/veri_" + self.iso + "_" + self.mode + ".dat"

    def get_hit_path(self):
        self.hit_path = self.dir + "/hits_" + self.iso + "_" + self.mode + ".dat"

    def get_url(self):
        refname = self.flip_isoname()
        self.url = "https://www-nds.iaea.org/relnsd/v1/data?fields=decay_rads&nuclides="+refname

    def get_data(self):
        print("Loading simulation data...")
        self.get_sim_data()
        print("Loading reference data from https://www-nds.iaea.org ...")
        self.get_ref_data()
        if self.return_hit_time: 
            print("Sorting and filtering photon hits...")
            self.get_hit_data()

    def get_sim_data(self):
        # simulatation data
        self.sim_df = pd.read_csv(self.sim_path, sep = self.delimiter)
        self.sim_alpha = self.sim_df.loc[self.sim_df["ParticleName"] == "alpha"]
        self.sim_electron = self.sim_df.loc[self.sim_df["ParticleName"] == "e-"]
        self.sim_positron = self.sim_df.loc[self.sim_df["ParticleName"] == "e+"]
        self.sim_gamma = self.sim_df.loc[self.sim_df["ParticleName"] == "gamma"]
        self.sim_photon = self.sim_df.loc[self.sim_df["ParticleName"] == "opticalphoton"]
        self.sim_particles = [self.sim_alpha, self.sim_electron, self.sim_positron, self.sim_gamma, self.sim_photon]

        # rescale energies (MeV -> keV)
        for p in self.sim_particles:
            p.iloc[:,2] *= 1000
    
    def get_ref_data(self):
        # reference data
        # cut-off values
        self.int_cutoff = 1E-3 # intensity cutoff in %
        self.hl_cutoff = 1E-8 # half-life cutoff in s

        df_ref_alpha = self.read_df_from_url("alpha") # alpha data
        if not df_ref_alpha.empty:
            self.ref_alpha = self.get_energies_intensities(df_ref_alpha, "alpha")
        else: self.ref_alpha = np.array([])
        
        df_ref_electron = self.read_df_from_url("electron") # beta- data
        if not df_ref_electron.empty:
            self.ref_electron = self.get_energies_intensities(df_ref_electron, "electron")
        else: self.ref_electron = np.array([])
        
        df_ref_positron = self.read_df_from_url("positron") # beta+ data
        if not df_ref_positron.empty:
            self.ref_positron = self.get_energies_intensities(df_ref_positron, "positron")
        else: self.ref_positron = np.array([])

        df_ref_gamma = self.read_df_from_url("gamma") # gamma-ray data
        if not df_ref_gamma.empty:
            # filter also half life to get rid of meta stable states
            self.ref_gamma = self.get_energies_intensities(df_ref_gamma, "gamma")
        else: self.ref_gamma = np.array([])
        
        df_ref_xray = self.read_df_from_url("xray") # X-ray data
        if not df_ref_xray.empty:
            self.ref_gamma = np.concatenate([self.ref_gamma, self.get_energies_intensities(df_ref_xray, "xray")])

        self.ref_particles = [self.ref_alpha, self.ref_electron, self.ref_positron, self.ref_gamma]

    def get_hit_data(self):
        # load dataframe
        self.hit_df = pd.read_csv(self.hit_path, sep = "\t")
        # name columns
        self.hit_df.columns = ["event_id", "time", "energy", "pmt_id", "hit_pos_x", "hit_pos_y", "hit_pos_z", 
                               "ver_pos_x", "ver_pos_y", "ver_pos_z", "parent_id", "survival", "angle"]
                
        # st = sort time
        self.hit_df_st, self.hit_time_st, self.hit_diff_time_st = self.get_hit_time(sorting=["time"])
        # st = sort time, QE biased
        self.hit_df_st_qe, self.hit_time_st_qe, self.hit_diff_time_st_qe = self.get_hit_time(sorting=["time"], biasing = True)
        # set = sort event time
        self.hit_df_set, self.hit_time_set, self.hit_diff_time_set = self.get_hit_time(sorting=["event_id", "time"])
        # sept = sort event pmt time
        self.hit_df_sept, self.hit_time_sept, self.hit_diff_time_sept = self.get_hit_time(sorting=["event_id", "pmt_id", "time"])
        # sept = sort pmt time
        self.hit_df_spt, self.hit_time_spt, self.hit_diff_time_spt = self.get_hit_time(sorting=["pmt_id", "time"])

    def get_hit_time(self, sorting, biasing = False):
        # filter only photons that survive QE biasing
        if biasing:
            df = self.hit_df[self.hit_df.iloc[:,-2]==1]
        else:
            df = self.hit_df
        # sorted dataframe
        df_sort = df.sort_values(by=sorting)
        # sorted hit times
        hit_time = df_sort.time

        if sorting == ["event_id", "pmt_id", "time"]:
            hit_unique_events = np.unique(df_sort.event_id) # unique events
            hit_unique_pmts = np.unique(df_sort.pmt_id) # unique events
            hit_diff_time = {
                f"PMT_{pmt_id}": np.concatenate([np.diff(hit_time[df_sort.event_id == event_id][df_sort.pmt_id == pmt_id]) for event_id in hit_unique_events])
            for pmt_id in hit_unique_pmts
            }
        elif sorting == ["pmt_id", "time"]:
            hit_unique_pmts = np.unique(df_sort.pmt_id) # unique events
            hit_diff_time = {
                f"PMT_{pmt_id}": np.diff(hit_time[df_sort.pmt_id == pmt_id]) for pmt_id in hit_unique_pmts}
        elif sorting == ["event_id", "time"]:
            hit_unique_events = np.unique(df_sort.event_id) # unique events
            # time difference between two hits of same event
            hit_diff_time = [np.diff(hit_time[df_sort.event_id == event_id]) for event_id in hit_unique_events]
            hit_diff_time = np.concatenate(hit_diff_time)
        elif sorting == ["time"]:
            # time difference between time sorted hits (no matter the event)
            hit_diff_time = np.diff(hit_time)
        return df_sort, hit_time, hit_diff_time

    def run(self, return_hit_times = False, flag_comb_hits = False):
        
        self.return_hit_time = return_hit_times
        self.flag_comb_hits = flag_comb_hits

        if self.flag_comb_hits:
            self.hit_path = self.dir + "/hits_" + self.isolist[0] + "Chain" + "_" + self.mode + ".dat"

            self.get_hit_data()
            plot_summary_combined_data(self)

        else:
            for iso in self.isolist:
                self.iso = iso
                self.get_sim_path()
                self.get_hit_path()
                self.get_url()            

                print("---------------------------\n{}\n---------------------------".format(self.iso))

                self.get_data()
                self.plot()        

    def plot(self):
        plot_decay_radiation(self)

        if self.return_hit_time:
            plot_example_trace(self)
            plot_time_diff(self)

    def flip_isoname(self):
        # Match and rearrange the pattern: letters followed by digits
        if self.iso.endswith("m"): 
            self.metastable = True
            name = self.iso[:-1]
        else: 
            self.metastable = False
            name = self.iso
        return re.sub(r'([a-zA-Z]+)(\d+)', r'\2\1', name).lower()
    
    def cutoff(self, df, particle):

        if particle == "electron" or particle == "positron":
            int_mask = df["intensity_beta"] > self.int_cutoff
        else:
            int_mask = df["intensity"] > self.int_cutoff
        
        if particle == "gamma":
            hl_mask = np.where(np.isnan(df["start_level_hl"]), 0, df["start_level_hl"]) < self.hl_cutoff
        else:
            hl_mask = np.ones_like(int_mask, dtype=bool)

        return (int_mask*hl_mask).values
    
    def get_energies(self, df, particle):

        mask = self.cutoff(df, particle)

        if particle == "electron" or particle == "positron":
            energies = df["mean_energy"][mask].values
        else:
            energies = df["energy"][mask].values
        
        return energies
    
    def get_intensities(self, df, particle):

        mask = self.cutoff(df, particle)

        if particle == "electron" or particle == "positron":
            intensities = df["intensity_beta"][mask].values
        else:
            intensities = df["intensity"][mask].values
        
        return intensities / 100 #from % to decimal
    
    def get_energies_intensities(self, df, particle):

        energies = self.get_energies(df, particle)
        intensities = self.get_intensities(df, particle)

        touple = np.array([energies, intensities]).T

        return np.unique(touple, axis = 0)
    
    def read_df_from_url(self, particle):

        # check particle type
        if particle == "alpha":
            full_url = self.url + "&rad_types=a"
        elif particle == "electron":
            full_url = self.url + "&rad_types=bm"
        elif particle == "positron":
            full_url = self.url + "&rad_types=bp"
        elif particle == "gamma":
            full_url = self.url + "&rad_types=g"
        elif particle == "xray":
            full_url = self.url + "&rad_types=x"

        req = urllib.request.Request(full_url)
        req.add_header('User-Agent', 
            'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
        res = urllib.request.urlopen(req)
        cont =  res.read().decode('utf-8').strip()
        file = StringIO(cont)

        df = pd.read_csv(file)

        if (particle == "gamma" or particle == "xray" or particle == "alpha") and df.size :
            if self.metastable: # only excited (meta-stable) states
                df = df.query("p_energy!=0")
            else:  # only ground states 
                df = df.query("p_energy==0")        

        return df   

class SigBkgPlotter:

    def __init__(self, 
                 dirname, 
                 filebase, 
                 type,
                 plotdir, 
                 filename = None,
                 simvol = None,
                 actvol = None, 
                 energy = None, 
                 chain = None,
                 time = None,
                 time_window = None, 
                 prob = None,
                 pos_den = None,
                 num_simpos = None,
                 mask = None,
                 depth = None,
                 pmtres = None):
        
        self.dirname = dirname
        self.filebase = filebase
        self.type = type
        self.plotdir = plotdir
        self.filename = filename
        self.simvol = simvol
        self.actvol = actvol
        self.energy = energy
        self.chain = chain
        self.time = time
        self.prob = prob
        self.pos_den = pos_den
        self.num_simpos = num_simpos
        self.mask = mask
        self.depth = depth   

        if time_window is not None:
            self.time_window = time_window
        else:
            self.time_window = 20

        if pmtres is not None:
            self.pmtres = pmtres
        else:
            self.pmtres = 5

        if self.type == "signal":
            assert ((self.simvol is not None) and (self.energy is not None) and (self.time is not None)), f"'type = signal' requires simvol (={self.simvol}), energy (={self.energy}) and time (={self.time}) to be defined."
            
            self.rescaling = self.num_simpos / (self.pos_den * self.actvol ** 3)
            
            self.suffix = f"_{self.simvol}m_{self.energy}MeV_{self.time}s"
            if self.filename is None: self.filename = self.filebase + self.suffix + ".dat"

        
        elif self.type == "background":
            assert ((self.simvol is not None) and (self.time is not None) and (self.prob is not None)), f"'type = background' requires simvol (={self.simvol}) and time (={self.time}) to be defined."
            
            self.rescaling = self.time

            if self.chain is not None:
                self.suffix = f"_{self.chain}_{self.simvol}m_{self.time}s"
            else: 
                self.suffix = f"_t={self.time}s_p={self.prob}%_V={self.simvol}m"
            if self.filename is None: self.filename = self.filebase + self.suffix + ".dat"

        
        elif self.type == "reference":
            assert ((self.mask is None) and (self.time is not None)), f"'type = reference' requires mask (={self.mask}) to be undefined and time (={self.time}) to be defined."
            
            self.rescaling = self.time
            
            self.suffix = f"_{self.time}s"
            if self.filename is None: self.filename = self.filebase + ".dat"

        else:
            raise ValueError(f"type = {self.type} is not supported.")

        self.path = os.path.join(self.dirname, self.filename)

        if self.mask is not None:
            self.suffix = self.suffix + "_qe"

        print(f"Loading data from file '{self.filename}' from directory '{self.dirname}'.")

    def load_data(self):

        self.df = load_data(self.path, self.type, self.mask, self.time)

    def combine_data(self, numfiles):

        df = pd.DataFrame() # empty dataframe

        for i in range(numfiles):
            load_path = os.path.splitext(self.path) + f"_{str(i)}.dat"

            if os.stat(load_path).st_size != 0: # check that file is not empty
                new_df = pd.read_csv(load_path, delimiter="\t", header = None) # new dataframe
                new_df.insert(loc = 0, column = "Run ID", value=np.ones(new_df.shape[0], dtype = int)*i)
                df = pd.concat([df, new_df], ignore_index=True)
            else:
                continue

        df.columns = ["run_id", "event_id", "time", "energy", "pmt_id", "hit_pos_x", "hit_pos_y", "hit_pos_z", 
                      "ver_pos_x", "ver_pos_y", "ver_pos_z", "parent_id", "survival", "angle"]
        
        df.to_csv(self.path, header=None, index=None, sep='\t', mode="w") # save file

        self.df = df

    def run(self):
        self.load_data()
        self.plot()

    def plot(self):

        plot_position_cathesian(self.df, mode = "hit")
        plt.savefig(self.plotdir + "vertex_position_cartesian" + self.suffix + ".png")
        plt.show()
        plt.close()
        """
        a = 1/0

        plot_trigger_duration(self.df, self.type, self.rescaling, time_window = self.time_window)
        plt.savefig(self.plotdir + "trigger_duration" + self.suffix + ".png", bbox_inches = "tight")
        plt.show()
        plt.close()

        a = 1/0

        plot_spectrum(self.df, depth = self.depth)
        plt.savefig(self.plotdir + "spectrum" + self.suffix + ".png")
        plt.show()
        plt.close()

        plot_position_cathesian(self.df, mode = "vertex")
        plt.savefig(self.plotdir + "vertex_position_cartesian" + self.suffix + ".png")
        plt.show()
        plt.close()

        plot_position_cathesian(self.df, mode = "hit")
        plt.savefig(self.plotdir + "hit_position_cartesian" + self.suffix + ".png")
        plt.show()
        plt.close()

        plot_position_polar(self.df, mode = "vertex")
        plt.savefig(self.plotdir + "vertex_position_spherical" + self.suffix + ".png")
        plt.show()
        plt.close()

        plot_position_polar(self.df, mode = "hit")
        plt.savefig(self.plotdir + "hit_position_spherical" + self.suffix + ".png")
        plt.show()
        plt.close()

        plot_trigger_rate_tw(self.df, rescaling = self.rescaling)
        plt.savefig(self.plotdir + "trigger_rate_tw" + self.suffix + ".png")
        plt.show()
        plt.close()

        if self.type != "reference":

            plot_trigger_pmt_tw(self.df, rescaling = self.rescaling, mode = "realistic", pmt_resolution = self.pmtres)
            plt.savefig(self.plotdir + "trigger_pmt_tw" + self.suffix + ".png", bbox_inches = "tight")
            plt.show()
            plt.close()
            
            plot_trigger_pmt_mode(self.df, self.rescaling, time_window = self.time_window)
            plt.savefig(self.plotdir + "trigger_pmt_mode" + self.suffix + ".png", bbox_inches = "tight")
            plt.show()
            plt.close()

        if self.type == "signal":

            if self.mask is None:
                plot_simulation_volume(self.df)
                plt.savefig(self.plotdir + "simvol" + self.suffix + ".png")
                plt.show()
                plt.close()

            plot_vertex(self.df, self.rescaling, time_window = self.time_window)
            plt.savefig(self.plotdir + "vertex" + self.suffix + ".png", bbox_inches = "tight")
            plt.show()
            plt.close()

            plot_vertex_tw(self.df, rescaling = self.rescaling)
            plt.savefig(self.plotdir + "vertex_tw" + self.suffix + ".png", bbox_inches = "tight")
            plt.show()
            plt.close()
        """
def load_data(path, type, mask, time):

        if type == "reference":
            df = pd.read_csv(path, sep="\s+", comment="#", header=None, names=["diff_time", "pmt_id"], dtype={"dt": np.float64, "pmt_id": np.int16})
            df.diff_time *= 1E9 # transform s -> ns
            df.insert(loc =2, column = "time", value = np.cumsum(df.diff_time)) # add time column that is cumulative sum of diff time
            df = df[df.time <= time * 1E9]
            return df
        
        else: # signal or background data from bulkice doumeki
            df = pd.read_csv(path, delimiter="\t", header = None)
            
            # really crappy condition, definitely needs rework
            if df.shape[1] == 6:
                df.columns = ["run_id", "event_id", "time", "energy", "pmt_id", "survival"]
                for r in np.unique(df.run_id):
                    run_mask = (df.run_id == r).values
                    df.loc[run_mask, "time"] += r*1E9
                df = df[df.time <= time * 1E9]
 
            else:
                if type == "background": # add extra column for run ID for background files
                    df.insert(loc = 0, column = "run_id", value=np.ones(df.shape[0], dtype = int))
                df.columns = ["run_id", "event_id", "time", "energy", "pmt_id", "hit_pos_x", "hit_pos_y", "hit_pos_z", 
                        "ver_pos_x", "ver_pos_y", "ver_pos_z", "parent_id", "survival", "angle"]
            
            if mask is not None:
                df = df[df.survival == 1] # keep photons that survive QE biasing
                
            return df

def get_coincidence_counts(df_in, mode , pmtcount = "pessimistic", return_vertex = True,
                           time_window = None, pmt_resolution = None, verbose = False):
    """_summary_

    Args:
        df_in (dataframe): incoming dataframe
        mode (str): Sort by "event" or "time" 
        pmtcount (str, optional): "non-cumulativ", "optimistic", "pessimistic", "realistic" counting of coincidences. 
        "non-cumulativ" is non-cumulative and only counts the number of events that fulfill the condition. 
        E.g. 5 hits on 3 PMTs would count as 1 only for PMT=3.
        "optimistic" is cumulative and counts the number of hits (one event can have several hits).
        E.g. 5 hits on 3 PMTs would count as 5 for PMT>=1, PMT>=2, and PMT>=3.
        "pessimistic" is cumulative and counts the number of unique PMTs (assuming that we cannot hold apart
        hits on the same PMT no matter the time resolution.)
        E.g. 5 hits on 3 PMTs would count as 3 for PMT>=1, PMT>=2, and PMT>=3.
        "realistic" is cumulative and counts the number of hits observing the time resolution of the sensor.
        E.g. 5 hits on 3 PMTs (3 on PMT1, 1 on PMT2, 1 on PMT3), and the first two hits on PMT1 arrive 
        within the time resolution of the PMT, would counts as 4 (1+1 for PMT1, 1 for PMT2, 1 for PMT3)
        for PMT>=1, PMT>=2, and PMT>=3. NOT IMPLEMENTED!
        Defaults to "pessimistic".
        return_vertex (bool, optional): Returns mean vertex position of triggered hits. Defaults to True.
        time_window (float, optional): Size of timewindow in ns. Defaults to None.
        pmt_resolution (float, optional): PMT resolution. Needed for pmtcount = "pessimistic". Defaults to None.
        verbose (bool, optional): Verbose for printing. Defaults to False.

    Returns:
        counts (list): Number of counts for siglet, dublets, multiplets
        vertex_radius (list): Mean vertex radius for given trigger (event, time)
        max_coin (list): Maximum coincidence.
        pmts (np.array): Histogram of coincidence PMTs.
    """

    if pmtcount == "realistic":
        assert  pmt_resolution is not None, "Realistic PMT counting scenario (pmtcount = 'realistic') requires pmt_resolution to be defined."

    df = df_in.copy(deep=True)
    
    counts = [0, 0, 0] # coincidence counts, =1, =2, >2
    pmts = np.zeros(24) # counts how many pmts were triggered

    vertex_radius = [[] for _ in range(3)] # coincidence vertex radius, =1, =2, >2, useful to constrain simulation volume
    dt_first_last = []

    if mode == "event":
        
        max_coinc = [0, 0, 0] # Run ID, Event ID, coincidence count

        unique_run = np.unique(df.run_id)

        # loop over all unique runs
        for r in tqdm(unique_run):
            mask_run = (df.run_id == r).values # run mask
            unique_event = np.unique(df[mask_run].event_id)

            # loop over all unique events
            for e in unique_event:
                mask_event = (df[mask_run].event_id == e).values # event mask for given run

                mdf = df[mask_run][mask_event] # masked dataframe

                coinc = mdf.shape[0] # number of coincidences

                unique_pmts = np.unique(mdf.pmt_id) # unique pmts
                num_unique_pmts = unique_pmts.size # number of pmts that "saw" the event

                if num_unique_pmts > 1:
                    times = mdf.time.copy().values
                    times.sort()
                    dt_first_last.append(times[-1]-times[0])
                
                if pmtcount == "non-cumulative":
                    pmts[num_unique_pmts-1] += 1
                elif pmtcount == "optimistic":
                    pmts[:num_unique_pmts] += coinc
                elif pmtcount == "pessimistic":
                    pmts[:num_unique_pmts] += num_unique_pmts
                elif pmtcount == "realistic":
                    cnt = 0 # initialize counts
                    for u in unique_pmts: # loop over all pmts that received hits
                        pmt_mask = (mdf.pmt_id == u).values # pmt mask for given event and run
                        mmdf = mdf[pmt_mask] # masked masked dataframe
                        pmt_time = mmdf.time.values
                        cnt += count_time_windows(pmt_time, pmt_resolution) # add resolution weighted counts
                    pmts[:num_unique_pmts] += cnt
                else:
                    raise ValueError(f"{pmtcount} not a valid argument. Chose between 'non-cumulative', 'optimistic', 'pessimistic' and 'realistic'.")

                if return_vertex:
                    ver_pos_x = mdf.ver_pos_x # vertex x position
                    ver_pos_y = mdf.ver_pos_y # vertex y position
                    ver_pos_z = mdf.ver_pos_z # vertex z position
                    # vertex radius
                    ver_rad = np.sqrt(ver_pos_x**2 + ver_pos_y**2 + ver_pos_z**2)

                # store data depending on whether the hit is a singlet, doublet or multiplet (>2)
                if coinc == 1: # singlet
                    counts[0] += 1
                    if return_vertex: vertex_radius[0].append(ver_rad.median())

                elif coinc == 2: # doublet
                    counts[1] += 1
                    if return_vertex: vertex_radius[1].append(ver_rad.median())

                elif coinc > 2: # more than doublet/multiplet
                    counts[2] += 1
                    if return_vertex: vertex_radius[2].append(ver_rad.median())

                if coinc > max_coinc[2]:
                    max_coinc[0] = r
                    max_coinc[1] = e
                    max_coinc[2] = coinc

        if verbose: print(f"Max Coincidence: Run ID {max_coinc[0]}, Event ID {max_coinc[1]}, Size {max_coinc[0]}")
    
    if mode == "time":
        
        perc = 0 # percent
        psteps = 100 # number of steps
        pbar = tqdm(total = psteps) # process bar
        df_size = df.shape[0] # original size of data frame

        max_coinc = 0

        df = df.sort_values(by="time")

        while df.shape[0] > 0:
            if ((df_size - df.shape[0]) > perc * df_size/psteps):
                pbar.update(1)
                perc += 1
            
            mask_time = df.time.values - df.time.values[0] < time_window
            mdf = df[mask_time] # masked dataframe
            coinc = mdf.shape[0] # number of coincidences

            unique_pmts = np.unique(mdf.pmt_id) # unique pmts
            num_unique_pmts = np.unique(mdf.pmt_id).size # number of pmts that "saw" the event

            if coinc > 1 and num_unique_pmts > 1:
                mdf.time.values.sort()
                dt_first_last.append(mdf.time.values[-1]-mdf.time.values[0])
            
            if pmtcount == "non-cumulative":
                pmts[num_unique_pmts-1] += 1
            elif pmtcount == "optimistic":
                pmts[:num_unique_pmts] += coinc
            elif pmtcount == "pessimistic":
                pmts[:num_unique_pmts] += num_unique_pmts
            elif pmtcount == "realistic":
                cnt = 0 # initialize counts
                for u in unique_pmts: # loop over all pmts that received hits
                    pmt_mask = (mdf.pmt_id == u).values # pmt mask for given event and run
                    mmdf = mdf[pmt_mask] # masked masked dataframe
                    pmt_time = mmdf.time.values
                    cnt += count_time_windows(pmt_time, pmt_resolution) # add resolution weighted counts
                pmts[:num_unique_pmts] += cnt

            if return_vertex:
                ver_pos_x = mdf.ver_pos_x # vertex x position
                ver_pos_y = mdf.ver_pos_y # vertex y position
                ver_pos_z = mdf.ver_pos_z # vertex z position
                # vertex radius
                ver_rad = np.sqrt(ver_pos_x**2 + ver_pos_y**2 + ver_pos_z**2)

            if coinc == 1: # singlet
                counts[0] += 1
                if return_vertex: vertex_radius[0].append(ver_rad.median())

            elif coinc == 2: # doublet
                counts[1] += 1
                if return_vertex: vertex_radius[1].append(ver_rad.median())

            elif coinc > 2: # more than doublet/multiplet
                counts[2] += 1
                if return_vertex: vertex_radius[2].append(ver_rad.median())

            if coinc > max_coinc:
                max_coinc = coinc
            
            # delete first len(coinc) rows
            df = df.iloc[coinc:]
        
        if verbose: print(f"Max Coincidence: Size {max_coinc}")

    if verbose: print(f"Singlet: {counts[0]}, Doublet: {counts[1]}, Multiplet: {counts[2]}")

    return counts, vertex_radius, max_coinc, pmts, dt_first_last

def get_trigger_vertex(df_in, sort, time_window = None):
    """Returns the average vertex position for event sorted or time sorted hit series. 
    In the time-sorted case (sort = "time") a time_window needs to be asigned.

    Args:
        df_in (dataframe): incoming dataframe
        sort (str): Sort by "event" or "time" 
        time_window (float, optional): Size of timewindow in ns. Defaults to None.

    Returns:
        vertex_radius (nested np.array): Mean vertex radius for given trigger (event, time)
    """

    if sort == "time":
        assert  time_window is not None, "Time-sorted mode (pmtcount = 'time') requires time_window to be defined."

    df = df_in.copy(deep=True)
    
    vertex_radius = [[] for _ in range(3)] # coincidence vertex radius, =1, =2, >2, useful to constrain simulation volume

    if sort == "event":
        
        unique_run = np.unique(df.run_id)

        # loop over all unique runs
        for r in tqdm(unique_run):
            mask_run = (df.run_id == r).values # run mask
            unique_event = np.unique(df[mask_run].event_id)

            # loop over all unique events
            for e in unique_event:
                mask_event = (df[mask_run].event_id == e).values # event mask for given run

                mdf = df[mask_run][mask_event] # masked dataframe

                coinc = mdf.shape[0] # number of coincidences
                
                ver_pos_x = mdf.ver_pos_x # vertex x position
                ver_pos_y = mdf.ver_pos_y # vertex y position
                ver_pos_z = mdf.ver_pos_z # vertex z position
                # vertex radius
                ver_rad = np.sqrt(ver_pos_x**2 + ver_pos_y**2 + ver_pos_z**2)

                # store data depending on whether the hit is a singlet, doublet or multiplet (>2)
                if coinc == 1: # singlet
                    vertex_radius[0].append(ver_rad.median())

                elif coinc == 2: # doublet
                    vertex_radius[1].append(ver_rad.median())

                elif coinc > 2: # more than doublet/multiplet
                    vertex_radius[2].append(ver_rad.median())
    
    elif sort == "time":
        
        # progress bar for while loop
        perc = 0 # percent
        psteps = 100 # number of steps
        pbar = tqdm(total = psteps) # process bar
        df_size = df.shape[0] # original size of data frame

        df = df.sort_values(by="time")

        while df.shape[0] > 0:
            if ((df_size - df.shape[0]) > perc * df_size/psteps):
                pbar.update(1)
                perc += 1
            
            mask_time = df.time.values - df.time.values[0] < time_window
            mdf = df[mask_time] # masked dataframe
            coinc = mdf.shape[0] # number of coincidences
        
            ver_pos_x = mdf.ver_pos_x # vertex x position
            ver_pos_y = mdf.ver_pos_y # vertex y position
            ver_pos_z = mdf.ver_pos_z # vertex z position
            # vertex radius
            ver_rad = np.sqrt(ver_pos_x**2 + ver_pos_y**2 + ver_pos_z**2)

            if coinc == 1: # singlet
                vertex_radius[0].append(ver_rad.median())

            elif coinc == 2: # doublet
                vertex_radius[1].append(ver_rad.median())

            elif coinc > 2: # more than doublet/multiplet
                vertex_radius[2].append(ver_rad.median())
            
            # delete first len(coinc) rows
            df = df.iloc[coinc:]

    else:
        raise ValueError(f"{sort} not a valid argument. Chose between 'event' and 'time'")

    # transform list of lists into array of arrays
    vertex = np.array([np.array(v) for v in vertex_radius], dtype = object)

    return vertex

def get_trigger_rate(df_in, sort, rescaling, time_window = None):
    """Returns rate for given coincidence. No PMT information is considered.
    In the time-sorted case (sort = "time") a time_window needs to be assigned.
    Rescaling is applied to correct rates to positron density or background rate.

    Args:
        df_in (dataframe): incoming dataframe
        sort (str): Sort by "event" or "time" 
        rescaling (float): Rescaling factor to correct positron density, background rate, etc.
        time_window (float, optional): Size of timewindow in ns. Defaults to None.

    Returns:
        rate (np.array): Rate for a given coincidence.
    """

    if sort == "time":
        assert  time_window is not None, "Time-sorted mode (pmtcount = 'time') requires time_window to be defined."

    df = df_in.copy(deep=True)
    
    rate = np.array([]) # coincidence counts, =1, =2, =3, ...

    if sort == "event":
        
        unique_run = np.unique(df.run_id)

        # loop over all unique runs
        for r in tqdm(unique_run):
            mask_run = (df.run_id == r).values # run mask
            unique_event = np.unique(df[mask_run].event_id)

            # loop over all unique events
            for e in unique_event:
                mask_event = (df[mask_run].event_id == e).values # event mask for given run

                mdf = df[mask_run][mask_event] # masked dataframe
                coinc = mdf.shape[0] # number of coincidences

                if rate.size < coinc: # extend array by coincidence - size of array and increase counter by one
                    rate = np.append(rate,np.zeros(coinc-rate.size), axis=0)
                
                rate[coinc-1] += 1
    
    elif sort == "time":

        # progress bar for while loop        
        perc = 0 # percent
        psteps = 100 # number of steps
        pbar = tqdm(total = psteps) # process bar
        df_size = df.shape[0] # original size of data frame

        df = df.sort_values(by="time")

        while df.shape[0] > 0:
            if ((df_size - df.shape[0]) > perc * df_size/psteps):
                pbar.update(1)
                perc += 1
            
            mask_time = df.time.values - df.time.values[0] < time_window
            mdf = df[mask_time] # masked dataframe
            coinc = mdf.shape[0] # number of coincidences

            if rate.size < coinc: # extend array by coincidence - size of array and increase counter by one
                rate = np.append(rate,np.zeros(coinc-rate.size), axis=0)
                
            rate[coinc-1] += 1            
            # delete first len(coinc) rows
            df = df.iloc[coinc:]

    else:
        raise ValueError(f"{sort} not a valid argument. Chose between 'event' and 'time'")

    rate = rate/rescaling
    return rate

def get_trigger_pmt(df_in, sort , rescaling = 1, mode = "realistic", time_window = None, pmt_resolution = None):
    """Returns coincidence rate for multiplicities from 1-24. The mode defines how coincidences are count.
    In the time-sorted case (sort = "time") a time_window needs to be assigned.
    In the realistic mode (mode = "realistic") a pmt_resolution needs to be assigned.

    Args:
        df_in (dataframe): incoming dataframe
        sort (str): Sort by "event" or "time" 
        mode (str, optional): "non-cumulativ", "optimistic", "pessimistic", "realistic" counting of coincidences. 
        "non-cumulativ" is non-cumulative and only counts the number of events that fulfill the condition. 
        E.g. 5 hits on 3 PMTs would count as 1 only for PMT=3.
        "optimistic" is cumulative and counts the number of hits (one event can have several hits).
        E.g. 5 hits on 3 PMTs would count as 5 for PMT>=1, PMT>=2, and PMT>=3.
        "pessimistic" is cumulative and counts the number of unique PMTs (assuming that we cannot hold apart
        hits on the same PMT no matter the time resolution.)
        E.g. 5 hits on 3 PMTs would count as 3 for PMT>=1, PMT>=2, and PMT>=3.
        "realistic" is cumulative and counts the number of hits observing the time resolution of the sensor.
        E.g. 5 hits on 3 PMTs (3 on PMT1, 1 on PMT2, 1 on PMT3), and the first two hits on PMT1 arrive 
        within the time resolution of the PMT, would counts as 4 (1+1 for PMT1, 1 for PMT2, 1 for PMT3)
        for PMT>=1, PMT>=2, and PMT>=3. NOT IMPLEMENTED!
        Defaults to "realistic".
        rescaling (float): Rescaling factor to correct positron density, background rate, etc.
        time_window (float, optional): Size of timewindow in ns. Defaults to None.
        pmt_resolution (float, optional): PMT resolution. Needed for pmtcount = "pessimistic". Defaults to None.

    Returns:
        pmt (list): Rate for a given multiplicity depending on the mode.
    """

    if sort == "time":
        assert  time_window is not None, "Time-sorted mode (pmtcount = 'time') requires time_window to be defined."


    if mode == "realistic":
        assert  pmt_resolution is not None, "Realistic PMT counting scenario (pmtcount = 'realistic') requires pmt_resolution to be defined."

    df = df_in.copy(deep=True)
    
    pmt = np.zeros(24) # counts how many pmts were triggered

    if sort == "event":
        
        unique_run = np.unique(df.run_id)

        # loop over all unique runs
        for r in tqdm(unique_run):
            mask_run = (df.run_id == r).values # run mask
            unique_event = np.unique(df[mask_run].event_id)

            # loop over all unique events
            for e in unique_event:
                mask_event = (df[mask_run].event_id == e).values # event mask for given run

                mdf = df[mask_run][mask_event] # masked dataframe
                coinc = mdf.shape[0] # number of coincidences

                unique_pmts = np.unique(mdf.pmt_id) # unique pmts
                num_unique_pmts = unique_pmts.size # number of pmts that "saw" the event
                
                if mode == "non-cumulative":
                    pmt[num_unique_pmts-1] += 1
                elif mode == "optimistic":
                    pmt[:num_unique_pmts] += coinc
                elif mode == "pessimistic":
                    pmt[:num_unique_pmts] += num_unique_pmts
                elif mode == "realistic":
                    cnt = 0 # initialize counts
                    for u in unique_pmts: # loop over all pmts that received hits
                        pmt_mask = (mdf.pmt_id == u).values # pmt mask for given event and run
                        mmdf = mdf[pmt_mask] # masked masked dataframe
                        pmt_time = mmdf.time.values
                        cnt += count_time_windows(pmt_time, pmt_resolution) # add resolution weighted counts
                    pmt[:num_unique_pmts] += cnt
                else:
                    raise ValueError(f"{mode} not a valid argument. Chose between 'non-cumulative', 'optimistic', 'pessimistic' and 'realistic'.")
    
    elif sort == "time":

        # progress bar for while loop        
        perc = 0 # percent
        psteps = 100 # number of steps
        pbar = tqdm(total = psteps) # process bar
        df_size = df.shape[0] # original size of data frame

        df = df.sort_values(by="time")

        while df.shape[0] > 0:
            if ((df_size - df.shape[0]) > perc * df_size/psteps):
                pbar.update(1)
                perc += 1
            
            mask_time = df.time.values - df.time.values[0] < time_window
            mdf = df[mask_time] # masked dataframe
            coinc = mdf.shape[0] # number of coincidences

            unique_pmts = np.unique(mdf.pmt_id) # unique pmts
            num_unique_pmts = unique_pmts.size # number of pmts that "saw" the event
            
            if mode == "non-cumulative":
                pmt[num_unique_pmts-1] += 1
            elif mode == "optimistic":
                pmt[:num_unique_pmts] += coinc
            elif mode == "pessimistic":
                pmt[:num_unique_pmts] += num_unique_pmts
            elif mode == "realistic":
                cnt = 0 # initialize counts
                for u in unique_pmts: # loop over all pmts that received hits
                    pmt_mask = (mdf.pmt_id == u).values # pmt mask for given event and run
                    mmdf = mdf[pmt_mask] # masked masked dataframe
                    pmt_time = mmdf.time.values
                    cnt += count_time_windows(pmt_time, pmt_resolution) # add resolution weighted counts
                pmt[:num_unique_pmts] += cnt
            else:
                raise ValueError(f"{mode} not a valid argument. Chose between 'non-cumulative', 'optimistic', 'pessimistic' and 'realistic'.")
            
            # delete first len(coinc) rows
            df = df.iloc[coinc:]

    else:
        raise ValueError(f"{sort} not a valid argument. Chose between 'event' and 'time'")
    
    pmt = pmt/rescaling
    return pmt


def get_trigger_duration(df_in, sort, time_window = None):
    """Returns time difference between first and last hit for triggers that involve more than 1 PMT.
    In the time-sorted case (sort = "time") a time_window needs to be assigned.

    Args:
        df_in (dataframe): incoming dataframe
        sort (str): Sort by "event" or "time" 
        time_window (float, optional): Size of timewindow in ns. Defaults to None.

    Returns:
        duration (list): Rate for a given multiplicity depending on the mode.
    """

    if sort == "time":
        assert  time_window is not None, "Time-sorted mode (pmtcount = 'time') requires time_window to be defined."

    df = df_in.copy(deep=True)
    
    duration_strict = [] # list of time difference between first and last entry, condition: at least two PMTs triggered
    duration_soft = [] # list of time difference between first and last entry, condition: at least two hits

    if sort == "event":
        
        unique_run = np.unique(df.run_id)

        # loop over all unique runs
        for r in tqdm(unique_run):
            mask_run = (df.run_id == r).values # run mask
            unique_event = np.unique(df[mask_run].event_id)

            # loop over all unique events
            for e in unique_event:
                mask_event = (df[mask_run].event_id == e).values # event mask for given run

                mdf = df[mask_run][mask_event] # masked dataframe
                coinc = mdf.shape[0] # number of coincidences
                num_unique_pmts = np.unique(mdf.pmt_id).size # number of pmts that "saw" the event

                if num_unique_pmts > 1: # condition: at least two PMT are triggered
                    times = mdf.time.copy().values # make copy to avoid sorting the dataframe
                    times.sort()               
                    duration_strict.append(times[-1]-times[0])
                
                if coinc > 1: # condition: at least two hits
                    times = mdf.time.copy().values # make copy to avoid sorting the dataframe
                    times.sort()               
                    duration_soft.append(times[-1]-times[0])
    
    elif sort == "time":

        # progress bar for while loop        
        perc = 0 # percent
        psteps = 100 # number of steps
        pbar = tqdm(total = psteps) # process bar
        df_size = df.shape[0] # original size of data frame

        df = df.sort_values(by="time")

        while df.shape[0] > 0:
            if ((df_size - df.shape[0]) > perc * df_size/psteps):
                pbar.update(1)
                perc += 1
            
            mask_time = df.time.values - df.time.values[0] < time_window
            mdf = df[mask_time] # masked dataframe
            coinc = mdf.shape[0] # number of coincidences
            num_unique_pmts = np.unique(mdf.pmt_id).size # number of pmts that "saw" the event

            if num_unique_pmts > 1: # condition: at least two PMT are triggered
                times = mdf.time.copy().values # make copy to avoid sorting the dataframe
                times.sort()               
                duration_strict.append(times[-1]-times[0])
            
            if coinc > 1: # condition: at least two hits
                times = mdf.time.copy().values # make copy to avoid sorting the dataframe
                times.sort()               
                duration_soft.append(times[-1]-times[0])       

            # delete first len(coinc) rows
            df = df.iloc[coinc:]

    else:
        raise ValueError(f"{sort} not a valid argument. Chose between 'event' and 'time'")

    return np.array(duration_strict), np.array(duration_soft)


def count_time_windows(time_series, resolution):
    count = 1  # Start with the first window
    window_start = time_series[0]

    for time in time_series:
        if time - window_start > resolution:
            count += 1  # Open a new window
            window_start = time  # Reset window start
    return count

# return maximum of nested 2D list
def get_max_2D_nested_list(list):
    lmax = 0
    for l in list:
        if len(l) == 0: # avoid zero size sublists
            continue
        else:
            tmax = np.max(l) # temporary maximum
            if (tmax > lmax):
                lmax = tmax
    return lmax

def cherenkov_noramlized(wvl):
    # returns normalized Cherenkov spectrum for an array of wavelength
    return 1/wvl**2 * (wvl[-1]*wvl[0]/(wvl[-1]-wvl[0]))

def absorption_normalized(wvl, depth, file = None):

    # load depth dependent ice properties table
    if file is None:
        file = "../../data/Materials/IceCubeICE.dat"
        
    with open(file, "r") as f:
        data = json.load(f)

    depth_list = np.array(data["jDepth_spice"])
    abs_400nm_list = 1/np.array(data["ja400inv_spice"])
    abs_400nm = abs_400nm_list[depth_list == depth]
    #abs_400nm_interpolate = InterpolatedUnivariateSpline(depth_list, abs_400nm_list, k = 3, ext = 3)
    
    kappa = 1.084106802940
    A = 6954.090332031250
    B = 6617.754394531250

    # see 2013 SPICE paper
    abs_dust = get_absorption_dust(wvl, abs_400nm, kappa)
    dtau = get_temperature(depth)-get_temperature(1730)
    abs_tot = abs_dust + A * np.exp(-B/wvl) * (1 + 0.01 * dtau)

    # absorption length = 1/ absorption coefficient
    abs_len_tot = 1/abs_tot

    # normalize
    abs_len_tot = abs_len_tot / ((wvl[1]-wvl[0])*abs_len_tot.sum())

    return abs_len_tot

def get_absorption_dust(wvl, abs_400nm, kappa):
    return abs_400nm * (wvl / 400) ** (-kappa)

def get_temperature(depth):
    return 221.5 - 0.00045319 * depth + 5.822E-6 * depth**2


######################################################

####   #       ###   #####  #####  ###  #   #   ###
#   #  #      #   #    #      #     #   ##  #  #   #
####   #      #   #    #      #     #   # # #  #
#      #      #   #    #      #     #   #  ##  #  ##
#      #####   ###     #      #    ###  #   #   ###

######################################################

def plot_spectrum(df, qefile = None, depthfile = None, depth = None):

    if qefile is None:
        qe_path = "./../../mdom/InputFile/TA0001_HamamatsuQE.data"

    wvl = (c.Planck * c.speed_of_light) / (df.energy * c.electron_volt) * 1E9 # wvl = (hc/E) * (1E9 for nm)
    qe = pd.read_csv(qe_path, sep=" ", header = None)
    qe.columns = ["wvl", "pdf"]
    qe.pdf = qe.pdf/((qe.wvl[1]-qe.wvl[0])*qe.pdf.sum())

    ch_pdf = cherenkov_noramlized(qe.wvl.values)

    if depth is not None:
        abs_pdf = absorption_normalized(qe.wvl.values, depth = depth, file = depthfile)
        comb_pdf = (qe.pdf * ch_pdf * abs_pdf)
    else:
        comb_pdf = (qe.pdf * ch_pdf)

    comb_pdf = comb_pdf/((qe.wvl[1]-qe.wvl[0])*comb_pdf.sum())

    fig, ax = plt.subplots(1,1)
    ax.hist(wvl, bins = qe.wvl, histtype = "step", density = True, lw = 2, label = "Simulation")
    ax.step(qe.wvl, qe.pdf, lw = 2, ls = "--", color = "grey", label = "QE Hamamatsu")
    ax.step(qe.wvl, ch_pdf, lw = 2, ls = "--", color = "black", label = "Cherenkov")
    if depth is not None: ax.step(qe.wvl, abs_pdf, lw = 2, ls = "--", color = "midnightblue", label = "Absorption")
    ax.step(qe.wvl, comb_pdf, lw = 2, ls = "--", label = "Combined")

    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Density")
    ax.legend()
    plt.tight_layout()


def plot_decay_radiation(self):

    fig, ax = plt.subplots(2,3, figsize = (15,10))
    ax = ax.flatten()
    texts = [r"$\alpha$", r"$e^-$", r"$e^+$", r"$\gamma$", r"$\lambda$"]

    # plot spectrum for alpha, beta, gamma and optical photons
    for i in range(6):
        if i < 5: 
            p = self.sim_particles[i]

            energy = p["Energy"]
            if i > 3: # for optical photons plot wavelength [nm]
                energy *= 1E3 # keV -> eV
                energy = (c.Planck * c.speed_of_light) / (energy * c.electron_volt) # wvl = (hc/E)
                energy *= 1E9 # m -> nm

            y, bins = np.histogram(energy, bins = 100)
            y = y / self.sim_num_particles
            x = (bins[1:]+bins[:-1])/2 
            ax[i].step(x,y, where="mid")

            if (i < 4) and (self.ref_particles[i].size): 
                if np.sum(y) == 0:
                    ysim_min = 1
                else: 
                    ysim_min = 0.5* y[np.where(y == 0, False, True)].min()
                ymin = np.minimum(0.5*np.min(self.ref_particles[i][:,1]), ysim_min)
                ax[i].vlines(x = self.ref_particles[i][:,0], ymin = 0, ymax = self.ref_particles[i][:,1], color = "red", lw = 0.5, ls = "--")
                ax[i].set_ylim(ymin,1)

            ax[i].text(x = 0.5, y = 0.8, s = texts[i], color = "k", transform = ax[i].transAxes, ha = "center", va = "center")
            ax[i].set_xlabel("Energy [keV]")
            if i == 4: ax[i].set_xlabel("Wavelength [nm]")
            ax[i].set_ylabel("Density")
            ax[i].set_yscale("log")
            ax[i].grid(True)
        else:
            p_name, p_count = np.unique(self.sim_df["ParticleName"], return_counts=True)

            ax[i].barh(p_name, p_count, color = "grey")
            ax[i].set_xscale("log")
            ax[i].set_xlim(0.5,)
            ax[i].set_yticklabels([])
            for y, label in zip(range(len(p_name)), p_name):
                if len(p_name) < 15: size = 10
                elif len(p_name) < 30: size = 6
                else: size = 4
                ax[i].text(10, y, label, va="center", ha="left", color="black", size = size)

    plt.show()
    plt.close()


def plot_time_diff(self):

    fig, ax = plt.subplots(1,1)
    ax.hist(self.hit_diff_time_st, density=True, bins = 100, histtype="step", label = "all events", lw = 2)
    ax.hist(self.hit_diff_time_set, density=True, bins = 100, histtype="step", label = "per event", lw = 2)
    for key in self.hit_diff_time_sept.keys():
        ax.hist(self.hit_diff_time_sept[key], density=True, bins = 100, histtype="step", alpha = 0.1, color = "grey", lw = 2)
        ax.hist(self.hit_diff_time_spt[key], density=True, bins = 100, histtype="step", alpha = 0.1, color = "red", lw = 2)

    ax.set_xscale("log")
    ax.set_xlim(1E0, 1E9)
    ax.set_yscale("log")
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Normalized Counts")
    ax.legend()

    plt.tight_layout()
    plt.show()
    plt.close()


def plot_example_trace(self, num = 5):
    num = 5
    colors = plt.cm.plasma(np.linspace(0,1,num))

    times = []
    for i in range(num):
        times.append(self.hit_df_set.time[self.hit_df_set.event_id == i])

    # Main plot, draw times for num events               
    fig, ax = plt.subplots(1,1, figsize = (10,2))
    for i in range(num):
        ax.vlines(times[i], ymin=0, ymax=1, alpha = 0.5, color = colors[i])

    ax.set_xlim(0,1E9)
    ax.set_ylim(0,1)
    ax.set_yticks([])
    ax.set_xlabel("Time [ns]")


    # Add Zoom-in on a specific event
    focus_event = 3  # Choose which event to zoom into
    focus_times = times[focus_event]

    x_min, x_max = 0.99999*focus_times.min(), 1.00001*focus_times.max()

    # Draw a box around the zoomed-in sequence in ax
    rect = patches.Rectangle(
        (x_min, 0),  # Bottom-left corner (x_min, ymin)
        x_max - x_min,  # Width of the event region
        1,  # Height (full y-range)
        linewidth=1, edgecolor="black", facecolor="none", zorder=10
    )
    ax.add_patch(rect)  # Add the rectangle to the main plot

    # Create inset axes
    width=ax.bbox.transformed(fig.gca().transAxes).width
    height=ax.bbox.transformed(fig.gca().transAxes).height
    axins = inset_axes(ax, width="100%", height="40%", bbox_to_anchor=(0, 1.05, 1, 0.4), bbox_transform=ax.transAxes, borderpad=0) # Adjust size & location


    # Plot only the focused event in the inset
    axins.vlines(focus_times, ymin=0, ymax=1, alpha=0.8, color=colors[focus_event])

    # Set zoomed-in x range
    axins.set_xlim(x_min, x_max)
    axins.set_ylim(0,1)
    axins.set_xticks([])
    axins.set_yticks([])

    # Add a box & connecting lines
    mark_inset(ax, axins, loc1=4, loc2=3, fc="none", ec="black", lw=1)  # loc1 & loc2 define diagonal connection points


    plt.tight_layout()
    plt.show()
    plt.close()

def plot_summary_combined_data(self):
    
    fig, ax = plt.subplots(2,1)

    # Plots hit times of all time-sorted optical photons
    ax[0].vlines(self.hit_time_st / 1E9, ymin = 0, ymax = 1, lw = 0.1, alpha = 0.2)
    ax[0].axvline(1, color = "k", lw = 2, ls = "--")
    ax[0].set_ylim(0,1)
    ax[0].set_xlabel("Time [s]")
    ax[0].set_yticks([])

    # Plots histogram of time difference for time-sorted optical photons
    # for all photons and QE biased
    ax[1].hist(self.hit_diff_time_st / 1E6, bins = 100, histtype="step", label = "all hits", lw = 2)
    ax[1].hist(self.hit_diff_time_st_qe / 1E6, bins = 100, histtype="step", label = "QE biased", lw = 2, ls = "--")

    ax[1].set_yscale("log")
    ax[1].set_xlabel("Time Difference [ms]")
    ax[1].set_ylabel("Counts")
    ax[1].legend()

    plt.tight_layout()
    plt.show()
    plt.close()

def plot_position_cathesian(df, mode = "vertex"):

    if mode == "vertex":
        x = df.ver_pos_x
        y = df.ver_pos_y
        z = df.ver_pos_z

    elif mode == "hit":
        x = df.hit_pos_x
        y = df.hit_pos_y
        z = df.hit_pos_z

    fig = plt.figure(figsize = (16,4))

    ax1 = fig.add_subplot(1,4,1, projection='3d')
    ax1.scatter(x,y,z, alpha = 0.05)
    if mode == "hit": ax1.scatter(ploc_car[:,0],ploc_car[:,1],ploc_car[:,2], color = "black")
    ax1.scatter(0,0,0, color = "red")
    ax1.view_init(elev=30., azim=45)
    ax1.set_xlabel("x [m]")
    ax1.set_ylabel("y [m]")
    ax1.set_zlabel("z [m]")
    ax1.set_box_aspect(None, zoom=0.8)

    ax2 = fig.add_subplot(1,4,2)
    ax2.scatter(x,y, alpha = 0.1)
    if mode == "hit": ax2.scatter(ploc_car[:,0],ploc_car[:,1], color = "black")
    ax2.scatter(0,0, color = "red", zorder = 10)
    ax2.set_xlabel("x [m]")
    ax2.set_ylabel("y [m]")

    ax3 = fig.add_subplot(1,4,3)
    ax3.scatter(x,z, alpha = 0.1)
    if mode == "hit": ax3.scatter(ploc_car[:,0],ploc_car[:,2], color = "black")
    ax3.scatter(0,0, color = "red", zorder = 10)
    ax3.set_xlabel("x [m]")
    ax3.set_ylabel("z [m]")

    ax4 = fig.add_subplot(1,4,4)
    ax4.scatter(y,z, alpha = 0.1)
    if mode == "hit": ax4.scatter(ploc_car[:,1],ploc_car[:,2], color = "black")
    ax4.scatter(0,0, color = "red", zorder = 10)
    ax4.set_xlabel("y [m]")
    ax4.set_ylabel("z [m]")
    plt.tight_layout()


def plot_position_polar(df, mode = "vertex", bins = 20):

    if mode == "vertex":
        x = df.ver_pos_x
        y = df.ver_pos_y
        z = df.ver_pos_z

    elif mode == "hit":
        x = df.hit_pos_x
        y = df.hit_pos_y
        z = df.hit_pos_z

    # convert cartesian to spherical coordinates
    r = np.sqrt(x**2 + y**2 + z**2) # radial distance
    phi = np.arctan2(y, x) # azimuth
    theta = np.arccos(z/r) # zenith

    # plot
    fig = plt.figure(figsize = (10,4))

    ax1 = fig.add_subplot(1,2,1)
    ax1.hist(r, bins = bins, histtype = "step")
    ax1.set_xlabel("r [m]")
    ax1.set_ylabel("Counts")

    ax2 = fig.add_subplot(1,2,2, polar = True)
    sc = ax2.scatter(phi, np.rad2deg(theta), c=r, alpha = 0.5)
    if mode == "hit": ax2.scatter(ploc_sph[:,1],np.rad2deg(ploc_sph[:,2]), color = "black")
    ax2.set_ylim(0,180)
    ax2.set_yticks([0,45,90,135,180])

    cax = fig.add_axes([1,0.1,0.03,0.9])
    cbar = fig.colorbar(sc, cax=cax, label="r [m]")#, ticks=[0,10,20,30])
    cbar.solids.set(alpha = 1)

    plt.tight_layout()

def plot_simulation_volume(df, bins = 20):

    # QE biased photons
    mask = df["survival"] == 1
    # vertex radius
    ver_radius = np.sqrt(df.ver_pos_x**2 + df.ver_pos_y**2 + df.ver_pos_z**2)

    # histogramming
    xrange = (np.min(ver_radius), np.max(ver_radius))
    y, edges = np.histogram(ver_radius, bins=bins, range=xrange)
    y_masked, edges = np.histogram(ver_radius[mask], bins=bins, range=xrange)

    x = (edges[1:]+edges[:-1])/2

    # cumulative distribution
    cum_y = np.cumsum(y)/np.sum(y)
    cum_y_masked = np.cumsum(y_masked)/np.sum(y_masked)

    # threshold value
    cum_thresh = 0.9
    idx = np.argmin(np.where(cum_y > cum_thresh, cum_y, np.inf))
    idx_masked = np.argmin(np.where(cum_y_masked > cum_thresh, cum_y_masked, np.inf))

    # plot
    fig, ax = plt.subplots(1,2, figsize = (10,4))

    ax[0].step(x, y, where = "mid", lw = 2, color = "C0", label="All")
    ax[0].errorbar(x, y, yerr = np.sqrt(y), ls = "", lw = 2, color = "C0", capsize = 5)
    ax[0].step(x, y_masked, where = "mid", lw = 2, ls = "--", color = "C1", label="QE biased")
    ax[0].errorbar(x, y_masked, yerr = np.sqrt(y_masked), ls = "", lw = 2, color = "C1", capsize = 5)
    ax[0].set_xlabel(r"Vertex distance [m]")
    ax[0].set_ylabel("Counts")
    ax[0].set_xlim(0,)
    #ax[0].set_ylim(0,)
    ax[0].grid()
    ax[0].set_yscale("log")
    ax[0].legend()

    # Add a second x-axis
    ax0_top = ax[0].twiny()  # Twin x-axis on the bottom subplot
    ax0_top.set_xlabel("Simulation box size [m]")  # Label for second x-axis

    # Adjust the second x-axis ticks to be 1/sqrt(3) of the first x-axis
    factor = 1 / np.sqrt(3)
    x_min, x_max = ax[0].get_xlim()
    trans_min, trans_max = x_min * factor, x_max * factor

    # Set whole number ticks on the upper axis
    upper_ticks = np.linspace(np.round(trans_min,0), np.round(trans_max,0) + 1, 5, endpoint=True)  # Whole numbers
    ax0_top.set_xticks(upper_ticks)


    ax[1].step(x, cum_y, lw = 2, label="All")
    ax[1].step(x, cum_y_masked, lw = 2, ls = "--", label="QE biased")
    ax[1].axhline(cum_thresh, color="grey", ls = "--", lw = 2)

    if idx != idx_masked:
        ax[1].axvline(x[idx], color="C0", ls = "--", lw = 2)
        ax[1].axvline(x[idx_masked], color="C1", ls = "--", lw = 2)
    else:
        ax[1].axvline(x[idx], color="grey", ls = "--", lw = 2)

    ax[1].set_xlabel(r"Vertex distance [m]")
    ax[1].set_ylabel("Cumulative Counts")
    ax[1].set_xlim(0,)
    ax[1].grid()

    ax1_top = ax[1].twiny()  # Twin x-axis on the bottom subplot
    ax1_top.set_xlabel("Simulation box size [m]")  # Label for second x-axis
    ax1_top.set_xticks(upper_ticks)

    plt.tight_layout()


def plot_trigger_rate_tw(df, rescaling):

    fig, ax = plt.subplots(1,1, sharex="row", sharey="col", figsize = (5,5))

    tws = np.array([5, 10, 20, 50]) # time windows in ns

    rate0 = get_trigger_rate(df, sort = "event", rescaling = rescaling)
    rate0 = np.array([rate0[0], rate0[1], np.sum(rate0[2:])]) # multiplet is cumulative for multiplicity > 2

    for i, tw in enumerate(tws):

        r = get_trigger_rate(df, sort = "time", time_window=tw, rescaling = rescaling)
        rate = np.array([r[0], r[1], np.sum(r[2:])]) # multiplet is cumulative for multiplicity > 2

        if i == 0:
            ax.scatter(tw, rate[0], label = "singlet", marker = "o", s = 20, color = "C0")
            ax.scatter(tw, rate[1], label = "doublet", marker = "x", s = 20, color = "C1")
            ax.scatter(tw, rate[2], label = "multiplet", marker = "^", s = 20, color = "C2")
        else:
            ax.scatter(tw, rate[0], marker = "o", s = 20, color = "C0")
            ax.scatter(tw, rate[1], marker = "x", s = 20, color = "C1")
            ax.scatter(tw, rate[2], marker = "^", s = 20, color = "C2")

    ax.axhline(rate0[0], color = "C0", lw = 2, alpha = 0.5)
    ax.axhline(rate0[1], color = "C1", lw = 2, alpha = 0.5)
    ax.axhline(rate0[2], color = "C2", lw = 2, alpha = 0.5)

    ax.set_xlabel("Time Window [ns]")
    ax.set_ylabel("Rate [Hz]")
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.legend()
    ax.grid()

def plot_vertex(df, rescaling, bins = 50, time_window = 20):

    verrad0 = get_trigger_vertex(df, sort = "event")
    verrad = get_trigger_vertex(df, sort = "time", time_window=time_window)
    labels = ["singlet", "doublet", "multiplet"]
    colors = ["C0", "C1", "C2"]

    # find the maximum vertex distance    
    xmax = get_max_2D_nested_list(verrad0)
    xrange = (0, xmax)

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    ax = ax.ravel()

    for i in range(3):
        y0, edge = np.histogram(verrad0[i], bins = bins, range = xrange)
        y, _ = np.histogram(verrad[i], bins = bins, range = xrange)
        y0, y = y0/rescaling, y/rescaling
        x = (edge[1:]+edge[:-1])/2

        ax[0].step(x, y0, lw = 2, color = colors[i], label = labels[i])
        ax[1].step(x, y, lw = 2, color = colors[i])

    ax[0].text(x = 0.5, y = 0.9, s = "event-sorted", color = "k", transform = ax[0].transAxes, ha = "center", va = "center")
    ax[1].text(x = 0.5, y = 0.9, s = r"time-sorted, $\Delta t$ = {:.0f} ns".format(time_window), color = "k", transform = ax[1].transAxes, ha = "center", va = "center")
    
    for i in range(2):
        ax[i].set_xlabel("Vertex distance [m]")
        ax[i].set_ylabel("Rate [Hz]")
        ax[i].set_yscale("log")
        ax[i].grid()

    fig.legend(ncol = 3, bbox_to_anchor=(0.75, 1.1))
    plt.tight_layout()

def plot_vertex_tw(df, rescaling, bins = 50):

    verrad0 = get_trigger_vertex(df, sort = "event")
    colors = ["C0", "C1", "C2"]
    labels = ["singlet", "doublet", "multiplet"]
    tws = np.array([10, 50, 1000, 5000]) # time windows in ns

    # find the maximum vertex distance    
    xmax = get_max_2D_nested_list(verrad0)
    xrange = (0, xmax)

    fig, ax = plt.subplots(2, 2, figsize=(8, 6), sharex="row", sharey="col")
    ax = ax.ravel()

    for i, tw in enumerate(tws):

        verrad = get_trigger_vertex(df, sort = "time", time_window=tw)
        for j in range(3):

            y0, edge = np.histogram(verrad0[j], bins = bins, range = xrange)
            y, _ = np.histogram(verrad[j], bins = bins, range = xrange)
            y0, y = y0/rescaling, y/rescaling
            x = (edge[1:]+edge[:-1])/2


            if i==0: 
                ax[i].step(x, y, lw = 2, color = colors[j], label = labels[j])
        
            else:
                ax[i].step(x, y0, lw = 2, alpha = 0.5, color = colors[j])
                ax[i].step(x, y, lw = 2, color = colors[j])

        ax[i].text(x = 0.5, y = 0.9, s = r"$\Delta t$ = {:.0f} ns".format(tw), color = "k", transform = ax[i].transAxes, ha = "center", va = "center")
    
        ax[i].set_yscale("log")
        ax[i].grid()    
    
    fig.supxlabel("Vertex distance [m]")
    fig.supylabel("Events")
    fig.legend(ncol = 3, bbox_to_anchor=(0.75, 1.05))
    plt.tight_layout()

def plot_trigger_pmt_tw(df, rescaling, mode = "realistic", pmt_resolution = 5):

    pmt0 = get_trigger_pmt(df, sort = "event", rescaling = rescaling, mode = mode, pmt_resolution = pmt_resolution)

    tws = np.array([5, 10, 20, 50]) # time windows in ns
    labels = [r"$\Delta t$=5 ns",r"$\Delta t$=10 ns",r"$\Delta t$=20 ns",r"$\Delta t$=50 ns"] # labels
    colors = plt.cm.viridis(np.arange(len(tws))/(len(tws)-1))

    fig, ax = plt.subplots(2, 1, figsize=(5, 5), gridspec_kw={'height_ratios': [3, 1]})

    ax[0].step(np.arange(24), pmt0, where = "mid", color = "k", label = "Event", alpha = 1, lw = 2)
    
    diffmax = 0
    for i, tw in enumerate(tws):

        pmt = get_trigger_pmt(df, sort = "time", rescaling = rescaling, mode = mode, time_window=tw, pmt_resolution = pmt_resolution)
        ax[0].step(np.arange(24), pmt, where = "mid", color = colors[i], label = labels[i].format(tw), 
                zorder = len(tws)-i, alpha = 1, lw = 2)
        
        if i > 0:
            diff = pmt-pmt0
            diffmax = np.maximum(diffmax, np.max(np.abs(diff)))
        
            ax[1].step(np.arange(24), diff, where = "mid", color = colors[i], 
                    zorder = len(tws)-i, alpha = 1, lw = 2)

    for i in range(2):
        ax[i].set_xlim(0,)
        ax[i].set_xlabel("PMT coincidence")
        ax[i].set_yscale("symlog")
        ax[i].grid()

    ylim = -10 ** np.ceil(np.log10(diffmax))
    ax[0].set_ylim(0,)
    ax[0].set_ylabel("Rate [Hz]")
    ax[1].set_yticks([-ylim,0,ylim])
    ax[1].set_ylabel(r"Residual")


    fig.legend(ncol = 3, bbox_to_anchor=(1, 1.15))
    plt.tight_layout()

def plot_trigger_pmt_mode(df, rescaling, time_window = 20, pmt_resolution = 5):

    modes = ["non-cumulative", "optimistic", "pessimistic", "realistic"]

    tws = np.array([5, 10, 20, 50]) # time windows in ns
    colors = plt.cm.viridis(np.arange(len(tws))/(len(tws)-1))
    lss = ["-.", "--", ":", "-"]

    fig, ax = plt.subplots(1, 2, figsize = (10,4))

    for i, md in enumerate(modes):
        # PMT distributions
        pmt0 = get_trigger_pmt(df, sort = "event", rescaling = rescaling, mode = md, pmt_resolution = pmt_resolution) # event-sorted
        pmt = get_trigger_pmt(df, sort = "time", rescaling = rescaling, mode = md, time_window=time_window, pmt_resolution = pmt_resolution) # time-sorted

        ax[0].step(np.arange(24), pmt0, where = "mid", color = colors[i], zorder = len(tws)-i, label = modes[i], alpha = 1, lw = 2, ls = lss[i])
        ax[1].step(np.arange(24), pmt, where = "mid", color = colors[i], zorder = len(tws)-i, alpha = 1, lw = 2, ls = lss[i])

    ax[0].text(x = 0.5, y = 0.9, s = "event-sorted", color = "k", transform = ax[0].transAxes, ha = "center", va = "center")
    ax[1].text(x = 0.5, y = 0.9, s = r"time-sorted, $\Delta t$ = {:.0f} ns".format(time_window), color = "k", transform = ax[1].transAxes, ha = "center", va = "center")

    for i in range(2):
        ax[i].set_ylim(0,)
        ax[i].set_xlabel("PMT multiplicity")
        ax[i].set_ylabel("Rate [Hz]")
        ax[i].set_yscale("symlog")
        ax[i].grid()

    fig.legend(ncol = 2, bbox_to_anchor=(0.6, 1.15))
    plt.tight_layout()

def plot_trigger_duration(df, type, rescaling, time_window = 20, bins = 100):

    if type == "reference":
        duration_time_strict, duration_time_soft = get_trigger_duration(df, sort = "time", time_window=time_window)

        dt_min = np.min(np.concatenate([duration_time_strict, duration_time_strict]))
        dt_max = np.max(np.concatenate([duration_time_strict, duration_time_soft]))

        dt_range = [None, np.logspace(np.log10(dt_min), np.log10(dt_max), num=(bins+1), endpoint=True)]
        dt = [None, None, duration_time_soft, duration_time_strict]

        i0 = 1

    else:

        duration_event_strict, duration_event_soft = get_trigger_duration(df, sort = "event")
        duration_time_strict, duration_time_soft = get_trigger_duration(df, sort = "time", time_window=time_window)

        q_event_90, q_event_99 = np.percentile(duration_event_strict, [90, 99])

        if type == "signal":
            dt_event_min = np.min(np.concatenate([duration_event_strict,  duration_event_soft, duration_time_strict, duration_time_strict]))
            dt_event_max = np.max(np.concatenate([duration_event_strict,  duration_event_soft, duration_time_strict, duration_time_soft]))

            dt_time_min = dt_event_min
            dt_time_max = dt_event_max

        else:
            dt_event_min = np.min(np.concatenate([duration_event_strict,  duration_event_soft]))
            dt_event_max = np.max(np.concatenate([duration_event_strict,  duration_event_soft]))

            dt_time_min = np.min(np.concatenate([duration_time_strict, duration_time_strict]))
            dt_time_max = np.max(np.concatenate([duration_time_strict, duration_time_soft]))

        dt_event_range = np.logspace(np.log10(dt_event_min), np.log10(dt_event_max), num=(bins+1), endpoint=True)
        dt_time_range = np.logspace(np.log10(dt_time_min), np.log10(dt_time_max), num=(bins+1), endpoint=True)

        dt_range = [dt_event_range, dt_time_range]
        dt = [duration_event_soft, duration_event_strict, duration_time_soft, duration_time_strict]

        i0 = 0

    fig, ax = plt.subplots(1,2, figsize = (10,4))

    for i in range(i0,len(ax)):
        y_soft, edge_soft = np.histogram(dt[2*i], bins = dt_range[i])
        y_strict, edge_strict = np.histogram(dt[2*i+1], bins = dt_range[i])

        y_soft = y_soft/rescaling
        y_strict = y_strict/rescaling

        x_soft = (edge_soft[1:]+edge_soft[:-1])/2
        x_strict = (edge_strict[1:]+edge_strict[:-1])/2

        ax[i].step(x_soft, y_soft, label = r"$N_{\text{hit}} > 1$" if i == 1 else None)
        ax[i].step(x_strict, y_strict, label = r"$N_{\text{hit}} > 1,\ N_{\text{PMT}} > 1$" if i == 1 else None)

        ax[i].set_xlabel("Time [ns]")
        ax[i].set_ylabel("Rate [Hz]")

        ax[i].set_xscale("log")
        ax[i].set_yscale("log")
        ax[i].grid()

    if type != "reference":
        ax[0].axvline(q_event_90, color = "black", ls = "--", label = "90%")
        ax[0].axvline(q_event_99, color = "black", ls = ":", label = "99%")
    
    ax[0].text(x = 0.5, y = 0.9, s = "event-sorted", color = "k", transform = ax[0].transAxes, ha = "center", va = "center")
    ax[1].text(x = 0.5, y = 0.9, s = r"time-sorted, $\Delta t$ = {:.0f}".format(time_window), color = "k", transform = ax[1].transAxes, ha = "center", va = "center")

    fig.legend(ncol = 2, bbox_to_anchor=(0.65, 1.1))


def plot_trigger_rate_energy(params):

    # read parameters
    dirname, energy_range, simvol, time, mask = params

    markers = ["o", "s", "^"]
    colors = ["C0", "C1", "C2"]
    labels = ["singlet", "doublet", "multiplet"]

    fig, ax = plt.subplots(1,1,figsize = (6,4))

    rates = np.zeros((3,len(energy_range)))

    for e, energy in enumerate(energy_range): # loop over energy

        filename = f"simvol_{simvol}m_{energy}MeV_{time}s.dat"

        print(f"Loading data from file '{filename}' from directory '{dirname}'.")

        path = os.path.join(dirname, filename)
        df = load_data(path = path, type = "signal", mask = mask) # load data

        rate = get_trigger_rate(df, sort = "event")
        rate =np.array([rate[0], rate[1], np.sum(rate[2:])]) # multiplet is cumulative for multiplicity > 2

        
        rates[:,e] =  rate

    for k in range(3): # loop over singlet, doublet, multiplet
        ax.plot(energy_range, rate[k], color = colors[k], lw = 2, 
                marker = markers[k], linestyle = "-", ms = 7, label = labels[k]) # event-based rate
        
    ax.set_xlabel("Energy [MeV]")
    ax.set_xticks([10,15,20])
    ax.set_ylabel("Rate [Hz]")
    ax.set_yscale("log")
    ax.grid()

    fig.legend(ncol = 3, bbox_to_anchor=(0.85, 1.1))

    fig.tight_layout()

def plot_gun_volume(params):

    # read parameters
    dirname, energy, simvol, time, mask = params
    
    # load file
    filename = f"simvol_{simvol}m_{energy}MeV_{time}s.dat"
    path = os.path.join(dirname, filename)
    df = load_data(path = path, type = "signal", mask = mask)

    # compute coincidence counts for event mode
    vertex_radius = get_trigger_vertex(df, sort = "event")

    radius_range = np.array([1,2,5,10,20,50,100,200,500]) # gun radius range
    volume_range = 4/3 * np.pi * radius_range**3 # gun volume

    density_positrons = 0.42969 # positron density per m**3
    num_positrons = density_positrons * volume_range # number of positrons in gun volume
    num_photons = np.zeros((3,radius_range.size)) # empty array

    for i, radius in enumerate(radius_range): # loop over radii
        for j in range(3): # loop over singlet, doublet, multiplet
            num_photons[j,i] = np.sum(vertex_radius[j] < radius)


    eff_volume = volume_range * num_photons/num_positrons # effective volume = gun volume * num_ph / num_positrons

    # plot
    markers = ["o", "s", "^"]
    colors = ["C0", "C1", "C2"]

    labels = ["singlet", "doublet", "multiplet"]

    fig, ax = plt.subplots(1,1)

    for i in range(3):
        ax.plot(volume_range, eff_volume[i], label = labels[i], marker = markers[i])

    ax.set_xlabel(r"Gun Volume [m$^3$]")
    ax.set_ylabel(r"Effective Volume [m$^3$]")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, axis = "y", which="major", linestyle = "-", alpha = 1)
    ax.legend()

    ax2 = ax.twiny()  # Twin x-axis on the bottom subplot
    ax2.set_xlabel("Gun Radius [m]")  # Label for second x-axis
    ax2.set_xticks(np.log10(volume_range),radius_range) # Correctly set the axis ticks
    ax2.set_xlim(np.log10(ax.get_xlim())) # Match the axis limits
    ax2.grid(True, which="major", axis = "x", linestyle = "-", alpha = 1)

    plt.tight_layout()
