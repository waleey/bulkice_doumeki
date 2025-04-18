import subprocess

class G4tools:

    def __init__(self, omModel, interaction, depthIndex, outputFolder, runID, baseFolder):
        self.omModel = omModel
        self.interaction = interaction
        self.baseFolder = baseFolder
        self.depthIndex = depthIndex
        self.outputFolder = outputFolder
        self.runID = runID
        self.executable = "./bulkice_doumeki"

    def callG4(self):
        print('called bulkice_doumeki')
        try:
            #subprocess.check_call(['. env.sh'])
            subprocess.check_call([self.executable, self.omModel, self.interaction, self.depthIndex, self.outputFolder, self.runID],
                                  cwd = self.baseFolder)

        except subprocess.CalledProcessError as e:
            print(f"error running package {e}")


    def GetChannel(self):
        return self.interaction
