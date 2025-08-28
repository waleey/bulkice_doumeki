from Event import Event
import numpy as np

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
        """
        re-weighting positron energy
        distribution according to
        effective volume in C. Lozano's
        Thesis
        """
        self.energy, self.x, self.y, self.z, self.dirX, self.dirY, self.dirZ, self.inTime = self.energy_reweighting(self.energy, self.inTime)
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

        """
        re-weighting positron energy
        distribution according to
        effective volume in C. Lozano's
        Thesis
        """

        self.energy, self.x, self.y, self.z, self.dirX, self.dirY, self.dirZ, self.inTime = self.energy_reweighting(self.energy, self.inTime)

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

    def energy_reweighting(self, energy_dist, time_dist, gen_vol_side = 40, energy_bin_width = 0.5):
        gen_vol = gen_vol_side ** 3
        energy_bins = np.arange(np.min(energy_dist), np.max(energy_dist), energy_bin_width)
        N, edges = np.histogram(energy_dist, bins=energy_bins)
        centers = 0.5 * (edges[1:] + edges[:-1])
        widths = edges[1:] - edges[:-1] #assumes uniform energy binning
        dNdE = N / (gen_vol * widths) #m^-3 MeV^-1

        veff_lozano = 127.9 * centers 
        #weighting the dNdE by the effective volume
        dNdE_veff = dNdE * veff_lozano #m^-3 MeV^-1

        N_center = dNdE_veff * 1.5e4 * widths
        radius_eff = (3 * veff_lozano / (4 * np.pi))**(1/3)  # m

        # Ensure N_center is integer
        N_center = np.rint(N_center).astype(int)

        # Side length of cubic volume (±20m around origin)
        L = gen_vol_side
        halfL = L / 2.0

        # Storage for all particles (across bins)
        all_x, all_y, all_z = [], [], []

        for N, r_eff in zip(N_center, radius_eff):
            x_vals, y_vals, z_vals = [], [], []
            
            while len(x_vals) < N:
                # Propose uniform random points in cube
                x = np.random.uniform(-halfL, halfL, size=N)
                y = np.random.uniform(-halfL, halfL, size=N)
                z = np.random.uniform(-halfL, halfL, size=N)
                
                # Radial distance
                r = np.sqrt(x**2 + y**2 + z**2)
                
                # Accept only those within radius_eff
                mask = r <= r_eff
                
                x_vals.extend(x[mask].tolist())
                y_vals.extend(y[mask].tolist())
                z_vals.extend(z[mask].tolist())
            
            # Trim to exactly N and extend the global lists
            all_x.extend(x_vals[:N])
            all_y.extend(y_vals[:N])
            all_z.extend(z_vals[:N])

        # Convert to single 1D numpy arrays
        all_x = np.array(all_x)
        all_y = np.array(all_y)
        all_z = np.array(all_z)

        all_E = []

        for N, mu, sigma in zip(N_center, centers, widths):
            # Sample N energies from Gaussian for this bin
            E_vals = np.random.normal(loc=mu, scale=sigma, size=N)
            all_E.extend(E_vals.tolist())

        # Convert to numpy array
        all_E = np.array(all_E)

        # Resample arrival times directly from empirical distribution
        all_time = np.random.choice(time_dist, size=len(all_E), replace=True)

        # Number of particles
        N_total = len(all_E)

        # Sample random directions isotropically
        phi = np.random.uniform(0, 2*np.pi, N_total)         # azimuth
        costheta = np.random.uniform(-1, 1, N_total)         # cos(theta)
        theta = np.arccos(costheta)

        # Convert to Cartesian components
        ax = np.sin(theta) * np.cos(phi)
        ay = np.sin(theta) * np.sin(phi)
        az = np.cos(theta)

        return all_E, all_x, all_y, all_z, ax, ay, az, all_time


