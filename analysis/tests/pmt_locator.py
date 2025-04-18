import numpy as np
import json
from astropy import units as u


class PMTLocator:

    def __init__(self,
                 num_pmts = 24,
                 units = False,
                 mdom_json_path = None,
                 pmt_json_path = None):
        
        self.num_pmts = num_pmts
        self.units = units
        self.mdom_json_path = mdom_json_path
        self.pmt_json_path = pmt_json_path

        # set default path
        if self.mdom_json_path is None:
            self.mdom_json_path = "../../data/SegmentedModules/om_mDOM.dat" 
        
        if self.pmt_json_path is None:
            self.pmt_json_path = "../../data/PMTs/pmt_Hamamatsu_R15458.dat"

    def load_parameters(self):
        """Loads json files and defines geometry variables.
        """

        # open and load json files
        mdom_file = open(self.mdom_json_path)
        mdom_data = json.load(mdom_file)

        pmt_file = open(self.pmt_json_path)
        pmt_data = json.load(pmt_file)
            
        # number of PMTs
        self.num_pol_pmts = mdom_data['jNrPolarPMTs']
        self.num_eq_pmts = mdom_data['jNrEqPMTs']
        self.num_tot_pmts = (self.num_pol_pmts + self.num_eq_pmts) * 2

        # geometrical measures
        self.glass_outer_radius = mdom_data['jGlassOutRad']['jValue'] * u.mm
        self.glass_thickness = mdom_data['jGlassThick']['jValue'] * u.mm
        self.glass_inner_radius = self.glass_outer_radius - self.glass_thickness
        self.gel_thinkness_front_pmt = mdom_data['jGelThicknessFrontPMT']['jValue'] * u.mm
        self.pmt_offset = pmt_data['jOuterShape']['jOutRad']['jValue'] * u.mm + pmt_data['jOuterShape']['jSpherePos_y']['jValue'] * u.mm - pmt_data['jOuterShape']['jEllipsePos_y']['jValue'] * u.mm
        self.eq_pmt_r_offset = mdom_data['jEqPMTrOffset']['jValue'] * u.mm
        self.theta_pol =  mdom_data['jThetaPolar']['jValue'] * u.deg
        self.theta_eq =  mdom_data['jThetaEquatorial']['jValue'] * u.deg

    def select_pmt(self, pmt_id):
        """Get spherical coordinates of mDOM PMTs given the PMT ID used internally within doumeki. 
        See OMSimMDOM.cc and OMSimPMTConstruction.cc for more information.

        Args:
            pmt_id (int): bulkice_doumeki PMT ID

        Returns:
            r (np.float)
            phi (np.float)
            theta (np.float)
        """
        r = self.glass_inner_radius - self.gel_thinkness_front_pmt - self.pmt_offset

        # 4 upper polar PMTs
        if ((pmt_id >= 0) and (pmt_id <= self.num_pol_pmts - 1)):
            theta = self.theta_pol
            phi = pmt_id * 360. * u.deg / self.num_pol_pmts

        # 8 upper equitorial PMTs
        if ((pmt_id >= self.num_pol_pmts) and (pmt_id <= self.num_pol_pmts + self.num_eq_pmts - 1)):
            theta = self.theta_eq
            phi = (pmt_id - self.num_pol_pmts + 0.5) * 360. * u.deg / self.num_eq_pmts
            r += self.eq_pmt_r_offset

        # 8 lower equitorial PMTs
        if ((pmt_id >= self.num_pol_pmts + self.num_eq_pmts) and (pmt_id <= self.num_pol_pmts + 2 * self.num_eq_pmts - 1)):
            theta = 180. * u.deg - self.theta_eq
            phi = (pmt_id - (self.num_pol_pmts + self.num_eq_pmts) + 0.5) * 360. * u.deg / self.num_eq_pmts
            r += self.eq_pmt_r_offset

        # 8 lower polar PMTs
        if ((pmt_id >= self.num_pol_pmts + 2 * self.num_eq_pmts) and (pmt_id <= self.num_tot_pmts - 1)):
            theta = 180. * u.deg - self.theta_pol
            phi = (pmt_id - (self.num_pol_pmts + 2 * self.num_eq_pmts)) * 360. * u.deg / self.num_pol_pmts
            
        r = r.to(u.m) # convert mm to m
        # return unitless values
        if not self.units:
            r = r.value
            phi = np.deg2rad(phi.value)
            theta = np.deg2rad(theta.value)

        return r, phi, theta
    
    def get_pmt_location_spherical(self):
        """Returns PMT locations in spherical coordinates (r, phi, theta).

        Returns:
            np.array(): PMT locations in spherical coordinates.
        """

        self.load_parameters()
        pmt_location_spherical = np.zeros((self.num_pmts, 3))
        for i in range(self.num_pmts):
            r, phi, theta = self.select_pmt(i)
            pmt_location_spherical[i] = np.array([r,phi,theta])

        return pmt_location_spherical
    
    def get_pmt_location_cartesian(self):
        """Returns PMT locations in cartesian coordinates (x, y, z).

        Returns:
            np.array(): PMT locations in cartesian coordinates.
        """

        self.load_parameters()
        pmt_location_cartesian = np.zeros((self.num_pmts, 3))
        for i in range(self.num_pmts):
            r, phi, theta = self.select_pmt(i)
            x, y, z = spherical2cartesian(r, phi, theta)
            pmt_location_cartesian[i] = np.array([x,y,z])

        return pmt_location_cartesian
    
def spherical2cartesian(r, phi, theta):
    """Coordinate transformation from spherical (r, phi, theta) to cartesian (x, y, z) frame.

    Args:
        r (np.float): Radius
        phi (np.float): Azimuth
        theta (np.float): Zenith

    Returns:
        x: x
        y: y
        z: z
    """

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z