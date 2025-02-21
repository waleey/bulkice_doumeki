import numpy as np
import pandas as pd
import urllib
from io import StringIO
import os
import re

from plthelper import *

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
        self.hit_df.columns = ["Event ID", "Time [ns]", "Energy [eV]", "PMT ID", "HIT POS X [m]", "HIT POS Y [m]", "HIT POS Z [m]", "VER POS X [m]", "VER POS Y [m]", "VER POS Z [m]", "?", "?", "?"]
        
        # st = sort time
        self.hit_df_st, self.hit_time_st, self.hit_diff_time_st = self.get_hit_time(sorting=["Time [ns]"])
        # st = sort time, QE biased
        self.hit_df_st_qe, self.hit_time_st_qe, self.hit_diff_time_st_qe = self.get_hit_time(sorting=["Time [ns]"], biasing = True)
        # set = sort event time
        self.hit_df_set, self.hit_time_set, self.hit_diff_time_set = self.get_hit_time(sorting=["Event ID", "Time [ns]"])
        # sept = sort event pmt time
        self.hit_df_sept, self.hit_time_sept, self.hit_diff_time_sept = self.get_hit_time(sorting=["Event ID", "PMT ID", "Time [ns]"])
        # sept = sort pmt time
        self.hit_df_spt, self.hit_time_spt, self.hit_diff_time_spt = self.get_hit_time(sorting=["PMT ID", "Time [ns]"])

    def get_hit_time(self, sorting, biasing = False):
        # filter only photons that survive QE biasing
        if biasing:
            df = self.hit_df[self.hit_df.iloc[:,-2]==1]
        else:
            df = self.hit_df
        # sorted dataframe
        df_sort = df.sort_values(by=sorting)
        # sorted hit times
        hit_time = df_sort["Time [ns]"]

        if sorting == ["Event ID", "PMT ID", "Time [ns]"]:
            hit_unique_events = np.unique(df_sort["Event ID"]) # unique events
            hit_unique_pmts = np.unique(df_sort["PMT ID"]) # unique events
            hit_diff_time = {
                f"PMT_{pmt_id}": np.concatenate([np.diff(hit_time[df_sort["Event ID"] == event_id][df_sort["PMT ID"] == pmt_id]) for event_id in hit_unique_events])
            for pmt_id in hit_unique_pmts
            }
        elif sorting == ["PMT ID", "Time [ns]"]:
            hit_unique_pmts = np.unique(df_sort["PMT ID"]) # unique events
            hit_diff_time = {
                f"PMT_{pmt_id}": np.diff(hit_time[df_sort["PMT ID"] == pmt_id]) for pmt_id in hit_unique_pmts}
        elif sorting == ["Event ID", "Time [ns]"]:
            hit_unique_events = np.unique(df_sort["Event ID"]) # unique events
            # time difference between two hits of same event
            hit_diff_time = [np.diff(hit_time[df_sort["Event ID"] == event_id]) for event_id in hit_unique_events]
            hit_diff_time = np.concatenate(hit_diff_time)
        elif sorting == ["Time [ns]"]:
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