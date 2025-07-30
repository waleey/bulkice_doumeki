import argparse
from data_formater import *

def parseCommandLine():
    parser = argparse.ArgumentParser()

    parser.add_argument('idx_low', type=int, help='Lower file index')    
    parser.add_argument('idx_high', type=int, help='Higher file index')       

    args = parser.parse_args()

    return args.idx_low, args.idx_high

if __name__ == "__main__":

    # path to directories
    sig_dir = "/home/jakob/software/doumeki/bulkice_doumeki/files/sig"
    rad_dir = "/home/jakob/software/doumeki/bulkice_doumeki/files/rad"
    noi_dir = "/home/jakob/software/doumeki/bulkice_doumeki/files/noi"
    out_dir = "/home/jakob/software/doumeki/bulkice_doumeki/files/comb"

    events_per_out = 10_000
    files_per_tar = 10_000

    trigger_window = 20 # in ns
    trigger_pmts = 2
    pmt_transit_mode = "fisher_tippett"

    load_mode = "tar"
    verbose = True

    idx_low, idx_high = parseCommandLine()

    # Call DataFormater class
    dformater = DataFormater(sig_dir, rad_dir, noi_dir, out_dir, 
                            events_per_out, files_per_tar, 
                            pmt_transit_mode, 
                            trigger_window, trigger_pmts, 
                            load_mode, verbose)

    # Run combine and filtering on subset of data
    dformater.run(idx_low, idx_high)