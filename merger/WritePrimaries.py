from Event import Event


class WritePrimaries:

    def __init__(self, events, baseFolder):
        self.events = events
        self.baseFolder = baseFolder
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

        for dtype in self.dtypes:
            file = open(self.baseFolder + 'Positron/pos20002nkibd_' + dtype[0] + '.data', 'w')
            for data in dtype[1]:
                file.write(str(data) + "\n")
            file.close()
        print(f"positron data has been written to {self.baseFolder}")
    def writeNeutron(self):
        pdgID = 2112
        reactionMode = -1001001
        #self.Reset()
        self.setupArrays(reactionMode, pdgID)
        self.dtypes = [['energy', self.energy], ['x', self.x], ['y', self.y], ['z', self.z], ['ax', self.dirX], ['ay', self.dirY], ['az', self.dirZ], ['time', self.inTime]]

        for dtype in self.dtypes:
            file = open(self.baseFolder + 'Neutron/neu20002nkibd_' + dtype[0] + '.data', 'w')
            for data in dtype[1]:
                file.write(str(data) + "\n")
            file.close()
        print(f"neutron data has been written to {self.baseFolder}")
    
    def writeElectron(self):
        pdgID = 11
        reactionMode = 98
        #self.Reset()
        self.setupArrays(reactionMode, pdgID)
        self.dtypes = [['energy', self.energy], ['x', self.x], ['y', self.y], ['z', self.z], ['ax', self.dirX], ['ay', self.dirY], ['az', self.dirZ], ['time', self.inTime]]

        for dtype in self.dtypes:
            file = open(self.baseFolder + 'Electron/e20002nkibd_' + dtype[0] + '.data', 'w')
            for data in dtype[1]:
                file.write(str(data) + "\n")
            file.close()
        print(f"electron data has been written to {self.baseFolder}")
    
    def setupArrays(self, reactionMode, pdgID):
        self.Reset()

        for event in self.events:
            if event.reaction_mode == reactionMode:
                for track in event.tracks[:10]:
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


