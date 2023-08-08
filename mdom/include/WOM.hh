#ifndef WOM_HH_INCLUDED
#define WOM_HH_INCLUDED
/**
*Very rough draft of possible WOM Simulation
*Waly M Z Karim
*8/8/2023
**/
#include "G4VSolid.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4MultiUnion.hh"
#include "G4UnionSolid.hh"
#include "G4Tubs.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Sphere.hh"
#include "G4EllipticalCone.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

class WOM
{
public:
    WOM(G4LogicalVolume* );
    ~WOM();

    void PlaceIt();

   // inline G4MultiUnion* GetOuterSolid() { return fGlassSolid; }
    //inline G4MultiUnion* GetInnerSolid() { return fGelSolid; }

private:
    void Construction();
    void GetSharedData();
    void GenerateLogicals();
    G4MultiUnion* PressureVessel(const G4String&, G4double, G4double);
    G4VSolid* WOMTube(const G4String&, G4double);
    // PMT Construction
    G4VSolid* PMTConstruction();

    void SubtractHarnessPlug();
    void SetPMTPosition();
    void SelLEDPosition();
    void ConstructMaterial();

    G4LogicalVolume* fLogicMother;
    //cylindrical tube
    G4double fGlassTubeOuterRad;
    G4double fGlassTubeInnerRad;
    G4double fGlassTubeHalfLength;


    //Spherical Cap
    G4double fGlassCapOuterRad;
    G4double fGlassCapInnerRad;

    //Vessel glass
    G4MultiUnion* fGlassSolid;
    //G4UnionSolid* fGlassSolid;
    G4MultiUnion* fGelSolid;
    G4LogicalVolume* fGlassLogical;
    G4LogicalVolume* fGelLogical;

    //WOM tube
    G4double fWomTubeHalfLength;
    G4double fWomTubeOuterRad;
    G4double fWomTubeInnerRad;

    G4VSolid* fWomTubeSolid;
    G4VSolid* fWomInsideSolid;

    G4LogicalVolume* fWomTubeLogical;
    G4LogicalVolume* fWomTubeInsideLogical;

    //PMTs
    G4double fPMTGlobalZ;

    G4LogicalVolume* fPMTBodyLogical1;
    G4LogicalVolume* fPMTBodyLogical2;
    G4LogicalVolume* fPMTCathodeLogical1;
    G4LogicalVolume* fPMTCathodeLogical2;

    //PMT Solids
    G4VSolid* fPMTSolid1;
    G4VSolid* fPMTSolid2;
    G4VSolid* fPMTHolder;
    G4VSolid* fPMTCone;
    G4VSolid* fPMTTube;
    G4VSolid* fPMTBoard;

    //PMT Holder Parameters
    G4double fPMTHolderRad;
    G4double fPMTHolderZ;
    G4double fPMTHolderHalfLength;
    //PMT Cone Parameters
    G4double fPMTConeHeight;
    G4double fPMTConeCut;
    G4double fPMTConeRatio;
    G4double fPMTConeZ;
    //PMT Tube Parameters
    G4double fPMTTubeRad;
    G4double fPMTTubeZ;
    G4double fPMTTubeHalfLength;
    //PMT Board parameters
    G4double fPMTBoardRad;
    G4double fPMTBoardZ;
    G4double fPMTBoardHalfLength;

    G4double fPMTOffset;

    //material construction
    G4Material* air; //general material for now. Will change later.



};


#endif // WOM_HH_INCLUDED
