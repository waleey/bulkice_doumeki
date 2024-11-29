from Event import Event
import os

class WritePrimaries:

    def __init__(self, events, baseFolder, name):
        self.events = events
        self.baseFolder = baseFolder
        self.name = name
        self.energy = []
        self.x = []
        self.y = []
        self.z = []
        self.dirX = []
        self.dirY = []
        self.dirZ = []
        self.inTime = []
        #self.dtypes = [['energy', self.energy], ['x', self.x], ['y', self.y], ['z', self.z], ['ax', self.dirX], ['ay', self.dirY], ['az', self.dirZ], ['time', self.inTime]]


    def writePositron(self):
        pdgID = -11
        reactionMode = -1001001
        #self.Reset()
        self.setupArrays(reactionMode, pdgID)
        self.dtypes = [['energy', self.energy], ['x', self.x], ['y', self.y], ['z', self.z], ['ax', self.dirX], ['ay', self.dirY], ['az', self.dirZ], ['time', self.inTime]]

        # Ensure the Positron directory exists
        positron_dir = os.path.join(self.baseFolder, 'Positron')
        if not os.path.exists(positron_dir):
            os.makedirs(positron_dir)

        for dtype in self.dtypes:
            file = open(self.baseFolder + 'Positron/pos_' + self.name + '_' + dtype[0] + '.data', 'w')
            for data in dtype[1]:
                file.write(str(data) + "\n")
            file.close()
        print(f"+++ Writing positron data to {self.baseFolder} +++")

    def writeNeutron(self):
        pdgID = 2112
        reactionMode = -1001001
        #self.Reset()
        self.setupArrays(reactionMode, pdgID)
        self.dtypes = [['energy', self.energy], ['x', self.x], ['y', self.y], ['z', self.z], ['ax', self.dirX], ['ay', self.dirY], ['az', self.dirZ], ['time', self.inTime]]

        # Ensure the Neutron directory exists
        neutron_dir = os.path.join(self.baseFolder, 'Neutron')
        if not os.path.exists(neutron_dir):
            os.makedirs(neutron_dir)

        for dtype in self.dtypes:
            file = open(self.baseFolder + 'Neutron/neu_' + self.name + '_' + dtype[0] + '.data', 'w')
            for data in dtype[1]:
                file.write(str(data) + "\n")
            file.close()
        print(f"+++ Writing neutron data to {self.baseFolder} +++")
    
    def writeElectron(self):
        pdgID = 11
        reactionMode = 98
        #self.Reset()
        self.setupArrays(reactionMode, pdgID)
        self.dtypes = [['energy', self.energy], ['x', self.x], ['y', self.y], ['z', self.z], ['ax', self.dirX], ['ay', self.dirY], ['az', self.dirZ], ['time', self.inTime]]

        # Ensure the Electron directory exists
        electron_dir = os.path.join(self.baseFolder, 'Electron')
        if not os.path.exists(electron_dir):
            os.makedirs(electron_dir)

        for dtype in self.dtypes:
            file = open(self.baseFolder + 'Electron/ele_' + self.name + '_' + dtype[0] + '.data', 'w')
            for data in dtype[1]:
                file.write(str(data) + "\n")
            file.close()
        print(f"+++ Writing electron data to {self.baseFolder} +++")
    
    def setupArrays(self, reactionMode, pdgID):
        self.Reset()

        for event in self.events:
            if event.reaction_mode == reactionMode:
                for track in event.tracks:
                    if(track['id'] == pdgID and track['state'] == 0):
                        self.energy.append(track['E'])
                        self.dirX.append(track['x'])
                        self.dirY.append(track['y'])
                        self.dirZ.append(track['z'])

                self.x.append((event.vertex[0] - 2000)/100)
                self.y.append((event.vertex[1] - 2000)/100)
                self.z.append((event.vertex[2] - 2000)/100)
                self.inTime.append(event.vertex[3])

    def Reset(self):
        self.energy = []
        self.x = []
        self.y = []
        self.z = []
        self.dirX = []
        self.dirY = []
        self.dirZ = []
        self.inTime = []


