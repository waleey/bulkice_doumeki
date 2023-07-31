import numpy as np

"""
Each event in sntools is considered as an object here.
events reaction_mode, vertex, tracks, and info are managed by Event class.
"""
track_dtype = [('id', int), ('E', float), ('x', float), ('y', float), ('z', float), ('state', int)] 
info_dtype = [('not-used0', int), ('not-used1', float), ('flux', float)]

class Event(object):    
    def __init__(self, reaction_mode, vertex, tracks, info, idx=0):
        """NUANCE Format event Object based on docs found at
        http://neutrino.phy.duke.edu/nuance-format/
        
        Parameters:
        -----------
        
        reaction_mode : int
            The particle interaction this event represents
            
        vertex : np.ndarray of float
            Array representing position and time of this event 
            idx  :  Quantity
              0  :  x-coord. [cm]
              1  :  y-coord. [cm]
              2  :  z-coord. [cm]
              3  :  time [ms] (Can differ, set by simulation time)
            
        tracks : tuple or list of np.ndarray
            List of tracks produced during this event (size may be variable)
            Each track is an array with elements representing the following:
            idx  :  Quantity
              0  :  PDG Id. Code (Can differ based on target)
              1  :  Total particle energy [MeV]
              2  :  x-direction cosine.
              4  :  y-direction cosine.
              5  :  z-direction cosine.
              6  :  State of particle (-1: initial, -2: intermediary, 0: final)
        
        info : np.ndarray
            Information on Event for reweighting to different fluxes
            idx  : Quantity
              0  : Defunct
              1  : Defunct
              2  : Neutrino flux in units matching simulation
        """
        self.reaction_mode = reaction_mode
        self.vertex = vertex
        self.tracks = tracks
        self.info = info
        self.idx = idx
        
    def __repr__(self):
        """String represenation of an Event
        """
        repr_str = f"{get_reaction_mode(self.reaction_mode)} Interaction Event: "
        
        # Show initial particles (pcl)
        pcls = []
        for track in self.tracks[self.tracks['state'] == -1]:
            pcls.append(str(Particle.from_pdgid(track['id'])))
        repr_str += ", ".join(pcls)
        repr_str += " --> "
        
        # Show final particles (pcl)
        pcls = []
        for track in self.tracks[self.tracks['state'] == 0]:
            pcls.append(str(Particle.from_pdgid(track['id'])))
        repr_str += ", ".join(pcls)
        return repr_str
 
    @classmethod
    def from_text(cls, lines):
        """Create event object from lines of text in NUANCE format
        """
        reaction_mode = int(lines[0].split()[1])
        vertex = np.array([float(var) for var in lines[1].split()[1:]])

        tracks = np.empty(0, dtype=track_dtype)
        info = np.empty(0, dtype=info_dtype)
        for line in lines[2:]:
            if 'track' in line:
                track = np.array(tuple(float(var) for var in line.split()[1:]), dtype=track_dtype)
                tracks = np.append(tracks, track)
            elif 'info' in line:
                info_line = np.array(tuple(float(var) for var in line.split()[1:]), dtype=info_dtype)
                info = np.append(info, info_line)
        return cls(reaction_mode, vertex, tracks, info)
   
def get_reaction_mode(mode, short=True):
    """Retrieve human-readable interaction type from NUANCE reaction mode
    
    Parameters
    ----------
    mode : int
        Index of interaction mode from the table provided at
        http://neutrino.phy.duke.edu/nuance-format/
        
    short : bool
        Switch to print short (if True) interaction name or
        full (if false) interaction name
    """
    interaction = (None, None)
    if mode == 1:
        interaction = ("CCQE", "Charged-current quasi-elastic")
    elif mode == 2:
        interaction = ("NCQE", "Neutral-current quasi-elastic")
    elif mode in range(3,17):
        interaction = ("SPrP", "Single-pion resonant production")
    elif mode in range(17,91):
        interaction = ("SPnrP", "Single-pion non-resonant production")
    elif mode == 91:
        interaction = ("CCDI", "Charged-current deep-inelastic")
    elif mode == 92:
        interaction = ("NCDI", "Neutral-current deep-inelastic")
    elif mode in range(93, 95):
        interaction = ('N/a', 'Not Used')
    elif mode == 95:
        interaction = ("CSQE", "Cabibbo-suppressed quasi-elastic")
    elif mode == 96:
        interaction = ("NCc/dPP", "Neutral current coherent/diffractive pion production")
    elif mode == 97:
        interaction = ("CCc/dPP", "Charged current coherent/diffractive pion production")
    elif mode == 98:
        interaction = ("ES", "Elastic scattering from electrons")
    elif mode == 99:
        interaction = ("IMD", "Inverse muon decay (electron target)")
    return interaction[0] if short else interaction[1]
