import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c

def plot_hits_event(events):

    unique, counts = np.unique(events, return_counts=True)

    fig, ax = plt.subplots(1,1, figsize = (10,4))
    ax.step(np.arange(len(unique)), counts)
    ax.set_xlabel("Triggered event [#]")
    ax.set_ylabel("Number of hits from the \n same event per sensor [#]")
    ax.grid(axis = "y")
    plt.gca().set_ylim(bottom=0)
    plt.tight_layout()

def plot_spectrum(energy):

    wvl = (c.Planck * c.speed_of_light) / (energy * c.electron_volt) * 1E9 # wvl = (hc/E) * (1E9 for nm)

    fig, ax = plt.subplots(1,1)
    ax.hist(wvl, bins = 10, histtype = "step", density = True)
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Density")
    plt.tight_layout()

def plot_lightcurve(time):

    y, bins = np.histogram(time, bins = 25)
    x = (bins[:-1] + bins[1:])/2
    t = x / 1E6 # in ms

    fig, ax = plt.subplots(1,1)

    ax.step(t,y, where = "mid")
    ax.set_xlabel("Time [ms]")
    ax.errorbar(t, y, yerr=(np.sqrt(y)), marker = "o", capsize = 5, ls = "none", color = "C0", alpha = 0.5)
    ax.set_ylabel("Counts/bin/sensor")
    plt.tight_layout()

def plot_hits_pmt(pmts):
    y, bins = np.histogram(pmts, bins = 24, range = (-0.5,23.5))
    x = (bins[:-1]+bins[1:])/2

    fig, ax = plt.subplots(1,1)
    ax.step(x, y, where = "mid")
    ax.errorbar(x, y, yerr=(np.sqrt(y)), marker = "o", capsize = 5, ls = "none", color = "C0", alpha = 0.5)
    ax.axvline(x = 11, color = "k", ls = "--")
    ax.text(x = 0.5, y = 0.5, s = "Equator", color = "k", rotation = 90, transform = ax.transAxes, ha = "center", va = "center")
    ax.text(x = 0.25, y = 0.9, s = "Upper Hemisphere", transform = ax.transAxes, ha = "center", va = "center")
    ax.text(x = 0.75, y = 0.9, s = "Lower Hemisphere", transform = ax.transAxes, ha = "center", va = "center")
    ax.set_xlim(-0.5,23.5)
    ax.set_xlabel("PMT ID")
    ax.set_ylabel("Counts")
    plt.tight_layout()

def plot_time_event(time_event):
    
    fig, ax = plt.subplots(1,1)

    for i in range(time_event.size):
        tevent = time_event[i]
        ax.step(np.arange(1,tevent.size+1), tevent, alpha = 0.75, where = "mid")

    ax.set_xlabel("Number of hits from the same event per sensor [#]")
    ax.set_ylabel("Time difference [ns]")
    plt.tight_layout()

def plot_coincidence_rate(coincidence_rate):

    fig, ax = plt.subplots(1,1)
    ax.step(np.arange(1,25), coincidence_rate, where = "mid")
    ax.set_xlim(1,25)
    ax.set_yscale("log")
    ax.set_yticks([1E3,1E2,1E1,1,1E-1,1E-2,1E-3])
    ax.set_xlabel("Number of coincidence PMTs")
    ax.set_ylabel("Time integrated rate [Hz]")
    ax.grid()
    plt.tight_layout()

def hist_time_diff(time_diff_pmt):

    ran = np.linspace(0, np.array(np.max(time_diff_pmt)).max(), 10)

    fig, ax = plt.subplots(1,1)

    for i in range(24):
        time_diff = time_diff_pmt[i]
        if time_diff == []:
            continue
        else:
            ax.hist(time_diff, bins = ran, histtype="step", label=">{:.0f}".format(i))

    ax.set_xlabel("Time difference [ns]")
    ax.set_ylabel("Counts")
    ax.legend(ncol = 2)
    plt.tight_layout() 