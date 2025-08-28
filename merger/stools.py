import subprocess

class stools:

    def __init__(self, progenitorModel,inputFileFormat, distance, outFile, baseFolder,start_time,end_time, transformation = None):
        self.progenitorModel = progenitorModel
        self.format=inputFileFormat
        self.distance = str(distance)
        self.outFile = outFile
        self.channel = None
        self.baseFolder = baseFolder
        self.inFile = baseFolder + self.progenitorModel
        self.start=start_time
        self.end=end_time
        self.transformation = transformation
    def callSntools(self):
        print('Called sntools....')
        try:
            if self.channel is not None:
                if self.end != 'None':
                    subprocess.check_call(['sntools', self.inFile, '--format', self.format, '--channel', self.channel, '--detector', 'IceCube', '--distance', self.distance, '--output', self.outFile,'--maxworkers','4','--verbose','--starttime',self.start,'--endtime',self.end, '--transformation', self.transformation])
                else:
                    subprocess.check_call(['sntools', self.inFile, '--format', self.format, '--channel', self.channel, '--detector', 'IceCube', '--distance', self.distance, '--output', self.outFile,'--maxworkers','4','--verbose','--starttime',self.start, '--transformation', self.transformation])
            else:
                if self.end != 'None':
                    subprocess.check_call(['sntools', self.inFile, '--format', self.format, '--detector', 'IceCube', '--distance', self.distance, '--output', self.outFile,'--maxworkers','4','--verbose','--starttime',self.start,'--endtime',self.end, '--transformation', self.transformation])
                else:
                    subprocess.check_call(['sntools', self.inFile, '--format', self.format, '--detector', 'IceCube', '--distance', self.distance, '--output', self.outFile,'--maxworkers','4','--verbose','--starttime',self.start, '--transformation', self.transformation])        
        except subprocess.CalledProcessError as e:
            print(f"error running package {e}")

    def getChannel(self):
        return self.channel
    def setChannel(self, channel):
        self.channel = channel
