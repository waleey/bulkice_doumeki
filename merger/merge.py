from Event import Event
from NUANCEReader import NUANCEReader
from WritePrimaries import WritePrimaries
from g4tools import G4tools
from stools import stools
import argparse
import sys

import time
import datetime

"""
This is the main function that calls sntools and geant4.
All variables in merge() function that is related to
calling sntools has upper case 'S' at the end and for
Geant4 it's G.
"""
def parseCommandLine():
    parser = argparse.ArgumentParser()

    #parser.add_argument('progenitorModel', help='name of the progenitor file to be used in sntools')
    parser.add_argument('nameS', help='name of the lightcurve file to be used in sntools')    
    #parser.add_argument('outfileS', help='name of the output file for sntools')
    parser.add_argument('distance', help='distance to the progenitor from earth in kpc')
    parser.add_argument('omModel', help='Optical Module model to be used: [dom, mdom, lom18, lom16, pmt]')
    parser.add_argument('simType', help='Simulation type: [ibd, enees, all, radioactivity]')
    parser.add_argument('depthIndex', help='Simulation depth index: [0, 1, ....., 108]')
    #parser.add_argument('outputFolderG', help = 'Output folder for bulkice_doumeki')
    parser.add_argument('runID', help = 'Run ID for each simulation run in bulkice_doumeki')
    args = parser.parse_args()

    #return args.progenitorModel, args.outfileS, args.distance, args.omModel, args.simType, args.depthIndex, args.outputFolderG, args.runID
    return args.nameS, args.distance, args.omModel, args.simType, args.depthIndex, args.runID
    

def merge():
    #progenitorModelS, outfileS, distanceS, omModelG, simTypeG, depthIndex, outputFolderG, runIDG = parseCommandLine()
    nameS, distanceS, omModelG, simTypeG, depthIndex, runIDG = parseCommandLine()

    doumekiFolder = '/home/jakob/software/doumeki/bulkice_doumeki/' # change here!!!

    basefolderG = doumekiFolder + 'mdom/build/' #goes back to the build folder!

    if simTypeG != "radioactivity":
        outTypeG = "signal/"
    else:
        outTypeG = "background/"

    infileS = doumekiFolder + 'analysis/files/input_sntools/gamma/' + nameS + '_' + runIDG + '.txt'
    outfileS = doumekiFolder + 'analysis/files/output_sntools/gamma/' + nameS + '_' + runIDG + '.kin'
    outputFolderG = doumekiFolder + 'analysis/files/output_geant4/' + outTypeG
   
    #initializing modules to call sntools and bulkice_doumeki
    #stool = stools(progenitorModelS, distanceS, infileS, outfileS, basefolderS)
    stool = stools(distanceS, infileS, outfileS)
    bulkice = G4tools(omModelG, simTypeG, depthIndex, outputFolderG, runIDG, basefolderG)
    useStool = True
     
    if(simTypeG == 'ibd'):
        stool.setChannel('ibd') #IBD events specified
        stool.callSntools()

    elif(simTypeG == 'enees'): #ENEES event specified
        stool.setChannel('es')
        stool.callSntools()

    elif(simTypeG == 'all'): #no channel specified. Will run every possible neutrino interaction in sntools
        stool.callSntools()

    elif(simTypeG == 'radioactivity'): #no need to call sntools as radioactivity is studied separately.
        useStool = False
    
    else:
        print(f"invalid simulation type {simTypeG} requested. Aborting...")
        sys.exit(0)
    

    """
    Write the output of sntools to input files of bulkice_doumeki.
    If you need to save the output somewhere else, change baseFolderW.
    """
    baseFolderW = doumekiFolder + 'analysis/files/input_geant4/gamma/' #goes back to to the InputFile dir.
    nameW = "gamma_" + runIDG

    if(useStool):
        with NUANCEReader(outfileS) as reader:
            events = reader.get_events()

        writer = WritePrimaries(events, baseFolderW, nameW)

        if(simTypeG == 'ibd'):
            #print("writes ibd")
            writer.writePositron()
            writer.writeNeutron()

        elif(simTypeG == 'enees'):
            writer.writeElectron()

        else:
            writer.writePositron()
            writer.writeNeutron()
            writer.writeElectron()
    
    #Calling Geant4 here
    bulkice.callG4()

if __name__ == "__main__":
    time0 = time.time()
    merge()
    runtime = time.time()-time0
    print("Run time: {}".format(str(datetime.timedelta(seconds=runtime))))



