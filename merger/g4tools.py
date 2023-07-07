import subprocess

class G4tools:

    def __init__(self, omModel, interaction, runMacro, baseFolder):
        self.omModel = omModel
        self.interaction = interaction
        self.baseFolder = baseFolder
        self.runMacro = self.baseFolder + runMacro
        self.executable = baseFolder + "./bulkice_doumeki"

    def callG4(self):
        print('called bulkice_doumeki')
        try:
            #subprocess.check_call(['. env.sh'])
            subprocess.check_call([self.executable, self.omModel, self.interaction, self.runMacro])

        except subprocess.CalledProcessError as e:
            print(f"error running package {e}")


    def GetChannel(self):
        return self.interaction
