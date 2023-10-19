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

#include "OMSimInputData.hh"

class WOM
{
public:
    WOM(G4LogicalVolume*, OMSimInputData* );
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
    G4VSolid* WOMPaint(const G4String&);
    // PMT Construction
    G4VSolid* PMTConstruction();
    G4VSolid* PMTCathodeConstruction();

    void SubtractHarnessPlug();
    void SetPMTPosition();
    void SelLEDPosition();
    void ConstructMaterial();

    OMSimInputData* fData;

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

    //WOM paint
    G4double fWOMPaintOuterRad;
    G4double fWOMPaintInnerRad;
    G4double fWOMPaintThickness;
    G4VSolid* fWOMPaintSolid;

    G4LogicalVolume* fWOMPaintLogical;
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

    //PMT Photo Cathode
    G4VSolid* fPMTCathode;
    G4double fPMTCathodeRad;
    G4double fPMTCathodeHalfLength;
    G4double fPMTCathodeZ;

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
    G4Material* glassMaterial;
    G4Material* fillerMaterial;
    G4Material* paintMaterial;
    G4Material* tubeMaterial;
    G4Material* tubeInsideMaterial;
    G4Material* pmtAbsorber;
    G4Material* pmtPhotoCathode;

    G4Material* air; //keeping this as dummy for PMT solids



};


#endif // WOM_HH_INCLUDED
