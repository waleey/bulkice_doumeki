from Event import Event
from NUANCEReader import NUANCEReader
from WritePrimaries import WritePrimaries
from g4tools import G4tools
from stools import stools
import argparse
import sys
"""
This is the main function that calls sntools and geant4.
All variables in merge() function that is related to
calling sntools has upper case 'S' at the end and for
Geant4 it's G.
"""
def parseCommandLine():
    parser = argparse.ArgumentParser()

    parser.add_argument('progenitorModel', help='name of the progenitor file to be used in sntools')
    parser.add_argument('outfileS', help='name of the output file for sntools')
    parser.add_argument('distance', help='distance to the progenitor from earth in kpc')
    parser.add_argument('omModel', help='Optical Module model to be used: [dom, mdom, lom18, lom16, pmt]')
    parser.add_argument('simType', help='Simulation type: [ibd, enees, all, radioactivity]')
    parser.add_argument('depthIndex', help='Simulation depth index: [0, 1, ....., 108]')
    parser.add_argument('outputFolderG', help = 'Output folder for bulkice_doumeki')
    parser.add_argument('runID', help = 'Run ID for each simulation run in bulkice_doumeki')
    args = parser.parse_args()

    return args.progenitorModel, args.outfileS, args.distance, args.omModel, args.simType, args.depthIndex, args.outputFolderG, args.runID
    

def merge():
    progenitorModelS, outfileS, distanceS, omModelG, simTypeG, depthIndex, outputFolderG, runIDG = parseCommandLine()
    basefolderS = '/home/waly/snewpy/SNEWPY_models/Nakazato_2013/' #you need to change this path
    basefolderG = '/home/waly/bulkice_doumeki/mdom/build/' #you need to change this path
   
    #initializing modules to call sntools and bulkice_doumeki
    stool = stools(progenitorModelS, distanceS, outfileS, basefolderS)
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
    baseFolderW = '/home/waly/bulkice_doumeki/mdom/InputFile/' #will change this later
    
    if(useStool):
        with NUANCEReader(outfileS) as reader:
            events = reader.get_events()

        writer = WritePrimaries(events, baseFolderW)

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
    merge()



