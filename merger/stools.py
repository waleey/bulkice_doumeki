import subprocess
import os

class stools:

    def __init__(self, distance, inFile, outFile):
        #self.progenitorModel = progenitorModel
        self.distance = str(distance)
        self.inFile = inFile #baseFolder + self.progenitorModel
        self.outFile = outFile
        self.channel = None
        #self.baseFolder = baseFolder

    def callSntools(self):
        print('Called sntools....')
        print(f'+++ Reading neutrino data from: {os.path.dirname(self.inFile)} +++')
        try:
            if self.channel is not None:
                subprocess.check_call(['sntools', self.inFile, '--format', 'gamma' '--channel', self.channel, '--detector', 'IceCube', '--distance', self.distance, '--output', self.outFile])

            else:
                subprocess.check_call(['sntools', self.inFile, '--format', 'gamma', '--detector', 'IceCube', '--distance', self.distance, '--output', self.outFile])
        except subprocess.CalledProcessError as e:
            print(f"error running package {e}")

    def getChannel(self):
        return self.channel
    def setChannel(self, channel):
        self.channel = channel
