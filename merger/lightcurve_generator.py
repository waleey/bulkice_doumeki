import numpy as np
import astropy.units as u
from tqdm import tqdm

def get_mean_square_energy(mean_energy, alpha):
    return (2+alpha)/(1+alpha) * mean_energy **2

class GenerateLightcurve:
    
    def __init__(self, parameter_range, format, flavor_same, time_const):
        self.time_start = parameter_range["time_start"].to(u.s).value
        self.time_end = parameter_range["time_end"].to(u.s).value
        self.time_steps = parameter_range["time_steps"]
        self.lumi_low = parameter_range["lumi_low"].to(u.erg/u.s).value
        self.lumi_high = parameter_range["lumi_high"].to(u.erg/u.s).value
        self.mean_energy_low = parameter_range["mean_energy_low"].to(u.MeV).value
        self.mean_energy_high = parameter_range["mean_energy_high"].to(u.MeV).value
        self.alpha_low = parameter_range["alpha_low"]
        self.alpha_high = parameter_range["alpha_high"]

        self.format = format
        self.flavor_same = flavor_same
        self.time_const = time_const

        self.time = np.linspace(start = self.time_start, stop = self.time_end, num = self.time_steps)

    def sample(self, num_files):

        print("---- SAMPLE LIGHTCURVE PARAMETERS FROM RANGE ----")
        print("Luminosity: {:.2e} - {:.2e}".format(self.lumi_low, self.lumi_high))
        print("Mean Energy: {:.2e} - {:.2e}".format(self.mean_energy_low, self.mean_energy_high))
        print("Alpha: {:.2e} - {:.2e}".format(self.alpha_low, self.alpha_high))
        
        self.num_files = num_files
        # if flavors differ, generate random entries for each one, i.e. three time more
        if self.flavor_same is True:
            self.samples = self.num_files
        else:
            self.samples = 3 * self.num_files

        # sample luminosity, mean energy and alpha from uniform distribution
        self.luminosity_sample = np.random.uniform(low = self.lumi_low, high = self.lumi_high, size = self.samples)
        self.mean_energy_sample = np.random.uniform(low = self.mean_energy_low, high = self.mean_energy_high, size = self.samples)
        self.alpha_sample = np.random.uniform(low = self.alpha_low, high = self.alpha_high, size = self.samples)

        if self.format == "mean_squared_energy":
            self.mean_squared_energy_sample = get_mean_square_energy(self.mean_energy_sample, self.alpha_sample)
        
        # for time constant evolution copy sampled data by len(time)
        if self.time_const:
            self.luminosity_data = np.tile(self.luminosity_sample, (self.time_steps, 1))
            self.mean_energy_data = np.tile(self.mean_energy_sample, (self.time_steps, 1))

            if self.format == "mean_squared_energy":
                self.mean_squared_energy_data = np.tile(self.mean_squared_energy_sample, (self.time_steps, 1))        
            else: 
                self.alpha_data = np.tile(self.alpha_sample, (self.time_steps, 1))
        else:
            raise ValueError("{} not supported yet. Only time constant evolution supported.".format(self.time_const))

    def write(self, filebase):
        
        print("---- WRITE DATA TO FILE ----")
        print(f"Directory: {filebase}")

        for i in tqdm(range(self.num_files)):

            filename = filebase + f"_{i}.txt"

            if self.flavor_same:
                flavor_data = np.array([self.mean_energy_data[:,i], self.mean_squared_energy_data[:,i], self.luminosity_data[:,i]])
                data = np.concatenate([self.time[np.newaxis,:], flavor_data, flavor_data, flavor_data], axis = 0)

            else:
                flavor_data = np.array([[self.mean_energy_data[:,3*i+0], self.mean_squared_energy_data[:,3*i+0], self.luminosity_data[:,3*i+0]],
                                        [self.mean_energy_data[:,3*i+1], self.mean_squared_energy_data[:,3*i+1], self.luminosity_data[:,3*i+1]],
                                        [self.mean_energy_data[:,3*i+2], self.mean_squared_energy_data[:,3*i+2], self.luminosity_data[:,3*i+2]]])
                self.flavor_data = flavor_data
                data = np.concatenate([self.time[np.newaxis,:], flavor_data[0], flavor_data[1], flavor_data[2]], axis = 0)
            
            np.savetxt(filename, data.T, delimiter=", ", fmt='%g')

    def run(self, num_files, filebase):

        self.sample(num_files)
        self.write(filebase)

        return 0
    
#Define variables
time_start = 0 * u.s
time_end = 1 * u.s
time_steps = 11

flavour_same = True # flavour has same luminosity, mean energy and mean squared energy
format = "mean_squared_energy" # "alpha", "mean_squared_energy"

lumi_low = 1E51 * u.erg/u.s 
lumi_high = 1E53 * u.erg/u.s

mean_energy_low = 10 * u.MeV
mean_energy_high = 15 * u.MeV

alpha_low = 2
alpha_high = 4

parameter_range = {"time_start": time_start,
                   "time_end": time_end,
                   "time_steps": time_steps,
                   "lumi_low": lumi_low,
                   "lumi_high": lumi_high,
                   "mean_energy_low": mean_energy_low,
                   "mean_energy_high": mean_energy_high,
                   "alpha_low": alpha_low,
                   "alpha_high": alpha_high}

num_files = 10
filebase = "/home/jakob/software/doumeki/bulkice_doumeki/analysis/files/input_sntools/gamma/gamma"

genlc = GenerateLightcurve(parameter_range=parameter_range, format = format, flavor_same = flavour_same, time_const = True)
genlc.run(num_files = num_files, filebase = filebase)    

