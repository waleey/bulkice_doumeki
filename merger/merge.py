from Event import Event
from NUANCEReader import NUANCEReader
from WritePrimaries import WritePrimaries
from g4tools import G4tools
from stools import stools
import argparse
import sys
import yaml
from pathlib import Path
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
    bulkice=G4tools(args.omModel,args.simType,args.depthIndex,args.outputFolderG,args.runID,basefolderG)
    useStool = True 
     
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
    Write the output of sntools to input files of bulkice_doumeki.
    If you need to save the output somewhere else, change baseFolderW.
    """
    baseFolderW = '/Users/walu/icecube/bulkice_doumeki/mdom/InputFile/' #goes back to to the InputFile dir.
    
    if(useStool):
        with NUANCEReader(args.outfileS) as reader:
            events = reader.get_events()
        writer = WritePrimaries(events, baseFolderW)

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
    
    #Calling Geant4 here
    #bulkice.callG4()

if __name__ == "__main__":
    merge()

