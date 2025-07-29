# Bulk-Ice DOUMEKI

This is a GEANT4 based simulation of different optical modules currently in use and to be deployed in future in IceCube Neutrino Observatory in South Pole. Currently available optical modules for simulation are MDOM, LOM16, LOM18, PDOM, D-Egg, WOM (currently under development). Except for WOM, the optical module simulation was initially written by M. Unland and C. Lozano, and later they were modified by me. The simulation also contains detailed depth dependent ice properties under antarctic ice sheet and it can simulation 40*40*40 cubic meter of ice at different temperature. Currently, it can simulate positron and electron flux from CCSN neutrinos, background radioactivity inside pressure vessel of MDOM, LOMs, and D-Egg, and photon wave with different Zenith Angle. It also accepts SNEWPY neutrino flux models in a python program called "merger" and use sntools to simulate the positron and electron flux from ibd and enees interactions. 


## Quick Start Guide

-Make sure you have Geant4 installed in your local machine.  

### Prerequisite for compiling bulkice\_doumeki: 

- **Set up env.sh file (must do before compilation):**
  - Go to "/path\_to\_bulkice\_doumeki/mdom/"
  - Open env.sh file in writing mode.
  - Set G4BUILD variable to your Geant4 installation directory
  - Save and exit.
  - Type ". env.sh" on your terminal to set up the environment
- **Set up build directory and compile:**
  - Go to "/path\_to\_bulkice\_doumeki/mdom/"
  - Type "mkdir build"
  - Type "cd build"
  - Type "cmake **..**"
  - Type "make -jN" or "make"

This will generate a binary executable under the "/path\_to\_bulkice\_doumeki/mdom/build/" directory.

- **Set Up the Output Folder**
  - Go to "/path\_to\_bulkice\_doumeki/mdom/build/"
  - Type "mkdir output". (You can name your output folder however you want).
  - Changing the output folder to a different directory will require changing the code, and will be explained later along with how to change the output fields.
  - Note: by default, Output will be dumped in a folder in the build directory.

## Simulating Events

Events can be simulated in two different ways. One is using the Visualization driver OpenGL, where you can see the detector geometry and the particle interaction, but it has limited capabilities in terms of visualizing a large number of particles. Another way is running it in batch mode, where you can simulate a realistic number of events easily.

## Visualization

Visualization is currently under development, so limited functionality is available. You can only visualize a plane wave of optical photons at an user-defined angle and their interaction with the WOM module. In order to do that:

- Make sure you are in the build directory and the program is compiled and the executable is generated.
- An example run could be **"./bulkice\_doumeki wom vis 45"**. Here, the simulated photon wave will make a 45 degree angle with the detector surface.

## Batch Mode

### Simulating IBD/ENEES/ALL/Radioactive\_Background Noise
    - Type **"./bulkice\_doumeki [om model] [interaction channel] [depth index] [output folder] [run id]"**
    - Available OM Models: [dom, mdom, lom16, lom18, pmt, degg, wom]
    - Available interaction channels: [ibd, enees, all, radioactivity]
    - Example run: " **./bulkice\_doumeki mdom ibd 88 output 0**"

### Simulating Photon Waves for Measuring Angular Sensitivity

#### Simulating photon waves at a single angle

  - Type "**./bulkice\_doumeki [om model] opticalphoton [depth index] [output folder] [run id] [distance from the detector center (m)] [angle (degree)]"**
  - Available OM Models: [dom, mdom, lom16, lom18, pmt, degg, wom]
  - Example run: " **./bulkice\_doumeki wom opticalphoton 88 output 0 2 45**".

#### Simulating photon waves within a range of angles

  - Type "**./bulkice\_doumeki [om model] opticalphoton [depth index] [output folder] [run id] [distance from the detector center (m)] [start angle (degree)] [final angle (degree)] [increment (degree)]"**
  - Available OM Models: [dom, mdom, lom16, lom18, pmt, degg, wom]
  - Example run: " **./bulkice\_doumeki wom opticalphoton 88 output 0 2 0 180 0.5**".

## Explanation of Input Parameters

**om model:**

Available Optical Module models in the simulation. For now, they are MDOM, WOM, DOM, DEGG, LOM16, and LOM18

**interaction channel:**

Possible interaction the input particles might have in the volume. If we are injecting positron flux, the possible interaction is **ibd** (Inverse Beta Decay). If we are injecting an electron flux, the possible interaction would be **enees** (Electron-Neutrino Electron Elastic Scattering). If someone wants both interactions to happen simultaneously, the channel would be **all**.

One could be interested in studying the angular sensitivity of the detector to optical photons. Therefore, they have to use the **opticalphoton** interaction channel. More on this later.

**depth index:**

The ice property varies with depth, and each depth is denoted by a depth index in the range [0, 108]. For example, 2.2Km depth has an index of 88. The depth indexes for each depth are given in this table:

**DEPTH INDEX TABLE**

**output folder**

By default, the output folder would be in the build directory of the simulation. If the user wants to dump the output data to a separate directory, they can provide the full path to the directory on the command line. An example with a full path would be: " **./bulkice\_doumeki mdom ibd 88 /home/waly/dump/ 0".**

**run id:**

Each run can be assigned an unique run id by the user. It is mainly to keep track of the files while running multiple runs in batch mode. One can set it to whatever integer number they like.

**distance from detector center:**

This option is available for simulating plane wave of photons only. It let's the user define a distance (in meters) from the center of the detector from where the photon wave will be generated.

**angle, start angle finish angle, increment:**

This is for simulating a plane wave of photons where the wave vector makes a specific angle with the detector surface.

If one wants to simulate a plane wave only at one angle, they might specify only the angle in degree.

If one simulates a range of angles with a specific step size, they might specify that start angle, final angle and the increment (step size). All are in degrees.

# Merger
Connects sntools to bulkice_doumeki. 
usage: <code> python3 merge.py [-h] [-t START_TIME] [-T END_TIME] progenitorModel inputFormat outfileS distance omModel simType depthIndex outputFolderG runID</code>

Arguments:
<ul>
  <li>progenitorModel:      name of the progenitor file to be used in sntools</li>
  <li>inputFormat:          format of the progenitor file to be used in sntools</li>
  <li>outfileS:             name of the output file for sntools</li>
  <li>distance:              distance to the progenitor from earth in kpc</li>
  <li>omModel:               Optical Module model to be used: (dom, mdom, lom18, lom16, pmt)</li>
  <li>simType:               Simulation type: (ibd, enees, all, radioactivity)</li>
  <li>depthIndex:            Simulation depth index: (0, 1, ....., 108)</li>
  <li>outputFolderG:         Output folder for bulkice_doumeki</li>
  <li>runID:                 Run ID for each simulation run in bulkice_doumeki</li></ul>

Options:<ul>
  <li>-h, --help:            show this help message and exit</li>
  <li>-t START_TIME, --start_time START_TIME:
                        Simulation start time passed to sntools [ms]</li>
  <li>-T END_TIME, --end_time END_TIME:
                        Simulation start time passed to sntools [ms]</li></ul>


**Advanced Features:**

**\<To be available soon\>**


*NB: * I am still in the process of building and debugging it. Especially, running the simulation in parallel with sntools to have new positron and electron flux each time is still under development and will be updated soon. Nevertheless, please feel free to share any suggestion, criticism, or thoughts. My email: wkarim@u.rochester.edu slack: Waly M Z Karim  

## A Note about Input Particle Files

When running the `opticalphoton`, `ibd`, and other input types, the particle energies and positions are read from files in the folder
```
mdom/InputFile
```
## A Note About Output Files

Output files are .dat files consisting of columns of numbers separated by tabs. The columns in order are: the runID, the time of the hit in nanoseconds, the energy of the photon, the PMT ID, the x,y, and z coordinates of the photon in meters, the x,y,x coordinates of the vertex in meters, the positron ID, and a 1 or 0 representing whether or not the hit passed the quantum efficiency of the PMT. For `opticalphoton` simulations, there is an additional column giving the angle. 

