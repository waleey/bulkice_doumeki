import os
import argparse
import time
import datetime
from tqdm import tqdm

import numpy as np
from scipy.stats import poisson

def parseCommandLine():
    parser = argparse.ArgumentParser()

    parser.add_argument('time_low', type=float, help='Start time t0')    
    parser.add_argument('time_high', type=float, help='End time t1')    
    parser.add_argument('mean_rate', type=float, help='Mean rate per PMT')    
    parser.add_argument('num_pmts', type=int, help='Number of PMTs')
    parser.add_argument('num_files', type=int, help='Number of files')
    parser.add_argument('dir', type=str, help='Path to directory')    

    args = parser.parse_args()

    return args.time_low, args.time_high, args.mean_rate, args.num_pmts, args.num_files, args.dir

class DarkNoiseGenerator:
    def __init__(self, time_low, time_high, mean_rate, num_pmts):
        self.time_low = time_low
        self.time_high = time_high
        self.mean_rate = mean_rate
        self.num_pmts = num_pmts

        assert time_high > time_low, "Error: time_high <= time_low!"

    def generate(self):
        self.noise_per_pmt = poisson.rvs(mu = self.mean_rate, size = self.num_pmts)

        self.time = np.random.uniform(low = self.time_low, high = self.time_high, size = np.sum(self.noise_per_pmt)).astype(dtype = np.float64)
        self.pmt = np.repeat(np.arange(self.num_pmts), self.noise_per_pmt).astype(dtype = np.uint8)

        self.data = np.array([self.time, self.pmt])

    def reset(self):
        self.noise_per_pmt = 0
        self.time = 0
        self.pmt = 0
        self.data = 0

    def save(self, idx):
        path = os.path.join(self.noi_dir, f"noi_{str(idx)}")
        np.savez(path, time = self.time, pmt = self.pmt)
        #self.data.tofile(path, format())

    def run(self, noi_dir, num_files):
        self.noi_dir = noi_dir
        for i in tqdm(range(num_files)):
            self.generate()
            self.save(idx = i)
            self.reset()

if __name__ == "__main__":

    time0 = time.time()

    time_low, time_high, mean_rate, num_pmts, num_files, dir = parseCommandLine()

    generator = DarkNoiseGenerator(time_low=time_low, time_high=time_high, mean_rate=mean_rate, num_pmts=num_pmts)
    generator.run(noi_dir=dir, num_files=num_files)

    runtime = time.time()-time0
    print("Run time: {}".format(str(datetime.timedelta(seconds=runtime))))