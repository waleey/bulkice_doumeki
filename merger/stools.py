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
                subprocess.check_call(['sntools', self.inFile, '--format', 'SNEWPY-Nakazato_2013', '--channel', self.channel, '--detector', 'IceCube', '--distance', self.distance, '--output', self.outFile])

            else:
                subprocess.check_call(['sntools', self.inFile, '--format', 'SNEWPY-Nakazato_2013', '--detector', 'IceCube', '--distance', self.distance, '--output', self.outFile])
        except subprocess.CalledProcessError as e:
            print(f"error running package {e}")

    def getChannel(self):
        return self.channel
    def setChannel(self, channel):
        self.channel = channel
