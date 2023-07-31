#ifndef OMSimMDOM_h
#define OMSimMDOM_h 1
#include "abcDetectorComponent.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "OMSimInputData.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"

extern G4double gRefCone_angle;

class mDOMHarness;

class mDOM : public abcDetectorComponent
{
public:
    // K.H. added default constructor for delivered class
    mDOM() {}
    mDOM(OMSimInputData* pData, G4bool pPlaceHarness = true);
    ~mDOM() {}
    // K.H. virtual added
    virtual void Construction();
    G4double mCylinderAngle;
    G4double mGlassOutRad;
    G4double mCylHigh;
    G4String mDataKey = "om_mDOM";
    G4int mNrTotalLED;
    std::vector<G4Transform3D> mLEDTransformers; //coordinates from center of the module
    std::vector<std::vector<G4double>> mLED_AngFromSphere; //stores rho (mm),theta (deg),phi (deg) of each LED from the center of its corresponding spherical part. Useful to run the particles.
    G4double getGlassOutRad() {return mGlassOutRad; }
    G4double getGlassInRad() {return mGlassInRad; }
    inline G4UnionSolid* GetOuterSolid() { return lGlassSolid; }
    inline G4UnionSolid* GetInnerSolid() { return lGelSolid; }

protected:
    OMSimPMTConstruction* mPMTManager;
    mDOMHarness* mHarness;
    void GetSharedData();
    G4SubtractionSolid* EquatorialReflector(G4VSolid* pSupportStructure, G4Cons* pReflCone, G4double pAngle, G4String pSuffix);
    void SetPMTPositions();
    G4UnionSolid* PressureVessel(const G4double pOutRad, G4String pSuffix);
    G4SubtractionSolid* SubstractHarnessPlug(G4VSolid* pSolid);
    std::tuple<G4SubtractionSolid*, G4UnionSolid*> SupportStructure();
    std::tuple<G4SubtractionSolid*, G4UnionSolid*, G4UnionSolid*, G4Tubs*>  LedFlashers(G4VSolid* lSupStructureSolid);
    void SetLEDPositions();

    // K.H. separated code to generate logicals from Construction function
    void GenerateLogicals();

    G4bool mPlaceHarness = true;
    G4bool mHarnessUnion = true; //it should be true for the first module that you build, and then false
    std::vector<G4ThreeVector> mPMTPositions;
    std::vector<G4RotationMatrix> mPMTRotations;
    std::vector<G4RotationMatrix> mPMTRotPhi;
    std::vector<G4ThreeVector> mReflectorPositions;

    G4LogicalVolume *lGlassLogical;
    G4LogicalVolume *lSupStructureLogical;
    G4LogicalVolume *lGelLogical;
    G4LogicalVolume *lRefConePolarLogical;
    G4LogicalVolume *lRefconeEqUpCutLogical;
    G4LogicalVolume *lRefconeEqLoCutLogical;
    G4LogicalVolume *lLEDHoleAirLogical;
    G4LogicalVolume *lLEDGlassTopLogical;
    G4LogicalVolume *lLEDLogical;

    //Shared data from jSON file
    G4double mGlassThick;
    G4double mGelThicknessFrontPMT;
    G4double mGelThickness;
    G4double mEqPMTrOffset;
    G4double mEqPMTzOffset;
    G4double mRefConeHalfZ;
    G4double mRefConeSheetThickness;
    G4double mThetaPolar;
    G4double mThetaEquatorial;
    G4int mNrPolarPMTs;
    G4int mNrEqPMTs;
    G4double mGlassInRad;
    G4double mRefConeAngle;
    G4int mTotalNrPMTs;
    G4double mPMToffset;
    G4double mRefConeIdealInRad;
    G4double mSupStructureRad;

private:
    G4UnionSolid* lGlassSolid;
    G4UnionSolid* lGelSolid;

};

#endif
//
