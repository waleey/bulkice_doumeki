#ifndef OMSimLOM16_h
#define OMSimLOM16_h 1
#include "abcDetectorComponent.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "OMSimInputData.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"

extern G4double gRefCone_angle;


class LOM16 : public abcDetectorComponent
{
public:
    LOM16(OMSimInputData* pData, G4bool pPlaceHarness = false);
    void Construction();
    G4double mCylinderAngle;
    G4double mGlassOutRad;
    G4double mCylHigh;
    G4String mDataKey = "om_LOM16";

  

private:
    OMSimPMTConstruction* mPMTManager;

    //functions
    void GetSharedData();
    G4UnionSolid* PressureVessel(const G4double pOutRad, G4String pSuffix);

    //Lom specific functions
    void PlaceDummySupportStructure(G4LogicalVolume* lInnerVolumeLogical);
    void PlaceCADSupportStructure(G4LogicalVolume* lInnerVolumeLogical);
    G4LogicalVolume* mSupportStructureLogical;
    G4LogicalVolume* CreateEquatorBand(const G4double pOutRad);
    
    //for gelpad and PMT creation
    void PlacePMTsAndGelpads(G4UnionSolid* lGelSolid,G4LogicalVolume* lGelLogical);
    void SetPMTAndGelpadPositions();
    void CreateGelpadLogicalVolumes(G4UnionSolid* lGelSolid);
    void PlacePMTs(G4LogicalVolume* lInnerVolumeLogical);
    void PlaceGelpads(G4LogicalVolume* lInnerVolumeLogical);

    //selection variables
    G4bool mPlaceHarness = true;
    G4bool mHarnessUnion = true; //it should be true for the first module that you build, and then false

    //vectors for positions and rotations
    std::vector<G4ThreeVector> mPMTPositions;
    std::vector<G4ThreeVector> mGelpadPositions;
    std::vector<G4double> mPMT_theta;
    std::vector<G4double> mPMT_phi;


    

    //Shared data from jSON file
    G4double mGlassThick;
    G4double mGelThicknessFrontPMT;
    G4double mGelThickness;
    G4double mEqPMTrOffset;
    G4double mEqPMTzOffset;
    G4double mThetaPolar;
    G4double mThetaEquatorial;
    G4int mNrPolarPMTs;
    G4int mNrEqPMTs;
    G4double mGlassInRad;
    G4int mTotalNrPMTs;
    G4double mEqTiltAngle;
    G4double mPolEqPMTPhiPhase;
    G4double mPolPadOpeningAngle;
    G4double mEqPadOpeningAngle;

    //from PMTConstruction class (not readable directly...needs to be changed)
    G4double mTotalLenght;
    G4double mOutRad;
    G4double mSpherePos_y; 
    G4double mEllipsePos_y; 

    //from PMT manager
    G4double mPMToffset;
    G4double mMaxPMTRadius;   

    
    //helper variables
    std::stringstream converter;
    std::stringstream converter2;
    G4Transform3D lTransformers;

    //for gelpads
    G4double jGelPadDZ;
    G4RotationMatrix* lRot = new G4RotationMatrix();

    G4double lPMT_theta, lPMT_phi, lPMT_x, lPMT_y, PMT_z, GelPad_x, GelPad_y, GelPad_z; //PMT
    G4double lPMT_theta_cone, lPMT_phi_cone, PMT_rho, GelPad_rho; //cones

    //logical of gelpads
    std::vector<G4LogicalVolume*> mGelPad_logical;


};

#endif
//
