import subprocess

class stools:

    def __init__(self, progenitorModel, distance, outFile, baseFolder):
        self.progenitorModel = progenitorModel
        self.distance = str(distance)
        self.outFile = outFile
        self.channel = None
        self.baseFolder = baseFolder
        self.inFile = baseFolder + self.progenitorModel

    def callSntools(self):
        print('Called sntools....')
        try:
            if self.channel is not None:
                subprocess.check_call(['sntools', self.inFile, '--format', 'gamma', '--channel', self.channel, '--detector', 'IceCube', '--distance', self.distance, '--output', self.outFile,'--maxworkers','4','--verbose','--endtime','1'])

            else:
                subprocess.check_call(['sntools', self.inFile, '--format', 'gamma', '--detector', 'IceCube', '--distance', self.distance, '--output', self.outFile,'--maxworkers','4','--verbose']) #,'--endtime','10'])
        except subprocess.CalledProcessError as e:
            print(f"error running package {e}")

    def getChannel(self):
        return self.channel
    def setChannel(self, channel):
        self.channel = channel
