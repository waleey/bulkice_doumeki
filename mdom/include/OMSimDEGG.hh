// #ifndef DEGG_Test_h
// #define DEGG_Test_h 1
#include "abcDetectorComponent.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "OMSimInputData.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
class dEGGHarness;
class dEGG : public abcDetectorComponent
{
public:
    
    dEGG(OMSimInputData* pData,G4bool pPlaceHarness=true);
    G4double mOut_sphere_r_max;
    void Construction();
private:
    OMSimPMTConstruction *mPMTManager;
    dEGGHarness* mHarness;
    
    

    G4int mOutSegments1;
    G4double mOutSphereRadiusMax;
    G4double mOutSphereDtheta;
    G4double mOutTransformZ;
    G4double mOutTorusRadius1;
    G4double mOutTorusRadius2;
    G4double mOutCenterOfTorusRadius1;
    G4int mOutSegments2;
    G4double mOutCenterOfTorusRadius2;
    G4double mOutCenterOfTorusZ2;
    G4double mOutTorusZmin2;
    G4double mOutTorusZmax2;
    G4double mOutTorusZ0;
    G4double mOutTorusTransformZ;

    G4int mInnSegments1;
    G4double mInnSphereRmax;
    G4double mInnSphereDtheta;
    G4double mInnTransformZ;
    G4double mInnTorusRadius1;
    G4double mInnTorusRadius2;
    G4double mInnCenterOfTorusRadius1;
    G4int mInnSegments2;
    G4double mInnCenterOfTorusRadius2;
    G4double mInnCenterOfTorusZ2;
    G4double mInnTorusZmin2;
    G4double mInnTorusZmax2;
    G4double mInnTorusZ0;
    G4double mInnTorusTransformZ;
    
    G4double mPmtDistance;
    G4double mGelHeight;
    G4double mMainBoardRmin;
    G4double mMainBoardRmax;
    G4double mMainBoardDz;
    G4double mHVBoardRmin;
    G4double mHVBoardRmax;
    G4double mHVBoardDz;
    G4double mMainBoardPosition;
    G4double mHVBoardPosition;
    G4double mGelOffset;
    G4String mDataKey = "om_DEGG";
    
    void GetSharedData();
    void PlaceGel();
    void PlacePMT();
    void PlaceCADSupportStructure(G4LogicalVolume* lInnerVolumeLogical);
    
    G4LogicalVolume* lgelsolid;
    G4LogicalVolume* lgelsolid1;
    G4LogicalVolume* lInnerVolumeLogical;
    
    G4VSolid* CreateEggSolid(G4int segments_1,
                            G4double r_max,
                            G4double d_theta,
                            G4double sphere_transform_z,
                            G4double torus_r_1,
                            G4double centeroftorus_r_1,
                            G4int segments_2,
                            G4double torus_r_2,
                            G4double centeroftorus_r_2,
                            G4double centeroftorus_z_2,
                            G4double torus2_zmin,
                            G4double torus2_zmax,
                            G4double torus_z0_2,
                            G4double torus1_transform_z);
    G4VSolid* Egg_Inner(G4int segments_1,
                            G4double r_max,
                            G4double d_theta,
                            G4double sphere_transform_z,
                            G4double torus_r_1,
                            G4double centeroftorus_r_1,
                            G4int segments_2,
                            G4double torus_r_2,
                            G4double centeroftorus_r_2,
                            G4double centeroftorus_z_2,
                            G4double torus2_zmin,
                            G4double torus2_zmax,
                            G4double torus_z0_2,
                            G4double torus1_transform_z);
    G4VSolid*Egg_Innertube(G4int segments_2,
                            G4double torus_r_2,
                            G4double centeroftorus_r_2,
                            G4double centeroftorus_z_2,
                            G4double torus2_zmin,
                            G4double torus2_zmax,
                            G4double torus_z0_2);

};