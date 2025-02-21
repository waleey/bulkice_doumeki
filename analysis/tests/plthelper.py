import numpy as np
import scipy.constants as c

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import matplotlib.patches as patches

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 12

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
        times.append(self.hit_df_set["Time [ns]"][self.hit_df_set["Event ID"] == i])

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