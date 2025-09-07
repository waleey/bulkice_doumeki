from Event import Event
from NUANCEReader import NUANCEReader
from WritePrimaries import WritePrimaries
from g4tools import G4tools
from stools import stools
import argparse
import sys
import yaml
from pathlib import Path
import numpy as np
"""
This is the main function that calls sntools and geant4.
All variables in merge() function that is related to
calling sntools has upper case 'S' at the end and for
Geant4 it's G.
"""
def parseCommandLine():
    parser = argparse.ArgumentParser()

    parser.add_argument('--progenitorModel', help='name of the progenitor file to be used in sntools')
    parser.add_argument('--inputFormat', help='format of the progenitor file to be used in sntools')
    parser.add_argument('--outfileS', help='name of the output file for sntools')
    parser.add_argument('--distance', help='distance to the progenitor from earth in kpc')
    parser.add_argument('--omModel', help='Optical Module model to be used: [dom, mdom, lom18, lom16, pmt]')
    parser.add_argument('--simType', help='Simulation type: [ibd, enees, all, radioactivity]')
    parser.add_argument('--depthIndex', help='Simulation depth index: [0, 1, ....., 108]')
    parser.add_argument('--outputFolderG', help = 'Output folder for bulkice_doumeki')
    parser.add_argument('-t', '--start_time', dest='start_time', default='0',
                         help='Simulation start time passed to sntools [ms]')
    parser.add_argument('-T', '--end_time', dest='end_time', default='None',
                         help='Simulation end time passed to sntools [ms]')
    parser.add_argument('--runID', help = 'Run ID for each simulation run in bulkice_doumeki', default = 0)
    parser.add_argument('--transformation', help = 'add neutrino flavor transformation ', default = 'NoTransformation', required = False)
    parser.add_argument('--ndoms', '-n', help = 'number of optical modules to simulate', default = 1, required = True, type = int)
    args = parser.parse_args()

    return args
#    return args.progenitorModel, args.outfileS, args.distance, args.omModel, args.simType, args.depthIndex, args.outputFolderG, args.runID

def merge():
    args = parseCommandLine()
#    progenitorModelS, outfileS, distanceS, omModelG, simTypeG, depthIndex, outputFolderG, runIDG = parseCommandLine()
    basefolderS = '/Users/walu/icecube/sntools/fluxes/' #you need to change this path
    basefolderG = '/Users/walu/icecube/bulkice_doumeki/mdom/build/'  #goes back to the build folder!
   
    #initializing modules to call sntools and bulkice_doumeki
    stool=stools(args.progenitorModel,args.inputFormat,args.distance,args.outfileS,basefolderS,args.start_time,args.end_time, args.transformation)
    useStool = True 
    """
    if(args.simType == 'ibd'):
        stool.setChannel('ibd') #IBD events specified
        stool.callSntools()

    elif(args.simType == 'enees'): #ENEES event specified
        stool.setChannel('es')
        stool.callSntools()

    elif(args.simType == 'all'): #no channel specified. Will run every possible neutrino interaction in sntools
        stool.callSntools()

    elif(args.simType == 'radioactivity'): #no need to call sntools as radioactivity is studied separately.
        useStool = False
    
    else:
        print(f"invalid simulation type {args.simType} requested. Aborting...")
        sys.exit(0)
    """

    """
    Write the output of sntools to input files of bulkice_doumeki.
    If you need to save the output somewhere else, change baseFolderW.
    """
    baseFolderW = '/Users/walu/icecube/bulkice_doumeki/mdom/build/merger/inputfile_dump/' #goes back to to the InputFile dir.
    
    if(useStool):
        with NUANCEReader(args.outfileS) as reader:
            events = reader.get_events()
        writer = WritePrimaries(events, baseFolderW, args.ndoms)

        if(args.simType == 'ibd'):
            #print("writes ibd")
            writer.writePositron()
            writer.writeNeutron()

        elif(args.simType == 'enees'):
            writer.writeElectron()

        else:
            writer.writePositron()
            writer.writeNeutron()
            writer.writeElectron()

    """
    runnign bulkice_doumeki on iteration
    for each unique dom id
    """
    del reader, events, writer
    #loading all positron data
    if(args.simType == 'ibd'):
        energy_positron = np.loadtxt(baseFolderW + 'Positron/pos20002nkibd_energy.data')
        x_positron = np.loadtxt(baseFolderW + 'Positron/pos20002nkibd_x.data')
        y_positron = np.loadtxt(baseFolderW + 'Positron/pos20002nkibd_y.data')
        z_positron = np.loadtxt(baseFolderW + 'Positron/pos20002nkibd_z.data')
        dirX_positron = np.loadtxt(baseFolderW + 'Positron/pos20002nkibd_ax.data')
        dirY_positron = np.loadtxt(baseFolderW + 'Positron/pos20002nkibd_ay.data')
        dirZ_positron = np.loadtxt(baseFolderW + 'Positron/pos20002nkibd_az.data')
        inTime_positron = np.loadtxt(baseFolderW + 'Positron/pos20002nkibd_time.data')
        dom_ids_positron = np.loadtxt(baseFolderW + 'Positron/pos20002nkibd_domid.data', dtype = int)
       
    #loading all electron data
    if(args.simType == 'enees'):
        energy_electron = np.loadtxt(baseFolderW + 'Electron/e20002nkibd_energy.data')
        x_electron = np.loadtxt(baseFolderW + 'Electron/e20002nkibd_x.data')
        y_electron = np.loadtxt(baseFolderW + 'Electron/e20002nkibd_y.data')
        z_electron = np.loadtxt(baseFolderW + 'Electron/e20002nkibd_z.data')
        dirX_electron = np.loadtxt(baseFolderW + 'Electron/e20002nkibd_ax.data')
        dirY_electron = np.loadtxt(baseFolderW + 'Electron/e20002nkibd_ay.data')
        dirZ_electron = np.loadtxt(baseFolderW + 'Electron/e20002nkibd_az.data')
        inTime_electron = np.loadtxt(baseFolderW + 'Electron/e20002nkibd_time.data')
        dom_ids_electron = np.loadtxt(baseFolderW + 'Electron/e20002nkibd_domid.data', dtype = int)

    #iterating over unique dom ids
    ndoms = args.ndoms
    for domid in range(ndoms):
        #writing the primaries in InputFile
        if(args.simType == 'ibd'):
            indices = np.where(dom_ids_positron == domid)[0]
            energy = energy_positron[indices]
            x = x_positron[indices]
            y = y_positron[indices]
            z = z_positron[indices]
            dirX = dirX_positron[indices]
            dirY = dirY_positron[indices]
            dirZ = dirZ_positron[indices]
            inTime = inTime_positron[indices]
            true_dom_ids = dom_ids_positron[indices]
            temp_write(energy, x, y, z, dirX, dirY, dirZ, inTime, true_dom_ids, particle_type = 'Positron')

        elif(args.simType == 'enees'):
            indices = np.where(dom_ids_electron == domid)[0]
            energy = energy_electron[indices]
            x = x_electron[indices]
            y = y_electron[indices]
            z = z_electron[indices]
            dirX = dirX_electron[indices]
            dirY = dirY_electron[indices]
            dirZ = dirZ_electron[indices]
            inTime = inTime_electron[indices]
            true_dom_ids = dom_ids_electron[indices]
            temp_write(energy, x, y, z, dirX, dirY, dirZ, inTime, true_dom_ids, particle_type = 'Electron')

        elif(args.simType == 'all'):
            indices = np.where(dom_ids_positron == domid)[0]
            energy = energy_positron[indices]
            x = x_positron[indices]
            y = y_positron[indices]
            z = z_positron[indices]
            dirX = dirX_positron[indices]
            dirY = dirY_positron[indices]
            dirZ = dirZ_positron[indices]
            inTime = inTime_positron[indices]
            true_dom_ids = dom_ids_positron[indices]
            temp_write(energy, x, y, z, dirX, dirY, dirZ, inTime, true_dom_ids, particle_type = 'Positron')

            indices = np.where(dom_ids_electron == domid)[0]
            energy = energy_electron[indices]
            x = x_electron[indices]
            y = y_electron[indices]
            z = z_electron[indices]
            dirX = dirX_electron[indices]
            dirY = dirY_electron[indices]
            dirZ = dirZ_electron[indices]
            inTime = inTime_electron[indices]
            true_dom_ids = dom_ids_electron[indices]
            temp_write(energy, x, y, z, dirX, dirY, dirZ, inTime, true_dom_ids, particle_type = 'Electron')

          
    
        bulkice=G4tools(args.omModel,args.simType,args.depthIndex,args.outputFolderG, domid ,basefolderG)
        bulkice.callG4()


def temp_write(energy, x, y, z, dirX, dirY, dirZ, inTime, dom_ids, particle_type = 'Positron'):
    folder = '/Users/walu/icecube/bulkice_doumeki/mdom/InputFile/' #goes back to to the InputFile dir.
    dtypes = [['energy', energy], ['x', x], ['y', y], ['z', z], ['ax', dirX], ['ay', dirY], ['az', dirZ], ['time', inTime], ['domid', dom_ids]]
    if(particle_type == 'Positron'):
        pref = 'pos20002nkibd_'
    elif(particle_type == 'Neutron'):
        pref = 'neu20002nkibd_'
    elif(particle_type == 'Electron'):
        pref = 'e20002nkibd_'
    for dtype in dtypes:
        file = open(folder + particle_type + '/' + pref + dtype[0] + '.data', 'w')
        for data in dtype[1]:
            file.write(str(data) + "\n")
        file.close()

    del dtypes


if __name__ == "__main__":
    merge()