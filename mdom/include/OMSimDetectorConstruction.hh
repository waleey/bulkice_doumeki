#ifndef OMSimDetectorConstruction_h
#define OMSimDetectorConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
//#include "G4Orb.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Box.hh"

#include "OMSimInputData.hh"
#include "OMSimPMTConstruction.hh"
#include "OMSimPMTConstruction.hh"
#include "OMSimMDOM.hh"
#include "OMSimPDOM.hh"
#include "OMSimLOM16.hh"
#include "OMSimLOM18.hh"
#include "OMSimDEGG.hh"
#include "OMSimRadioactivityData.hh"
#include "WOM.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include <cmath>
class OMSimDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    OMSimDetectorConstruction(G4int, G4double);
    OMSimDetectorConstruction();
    ~OMSimDetectorConstruction();
    G4VPhysicalVolume *Construct();

   /* inline void SetOMModel(G4int model) { fDOM = model; }
    inline void SetWorldSize(G4double worldSize) { fworldSize = worldSize; }*/

    G4ThreeVector DrawFromVolume();


private:
    //G4Orb *mWorldSolid;
    G4Box *mWorldSolid;
    G4LogicalVolume *mWorldLogical;
    G4VPhysicalVolume *mWorldPhysical;
    void ConstructWorld();
    void ConstructWorldMat();
    OMSimInputData *mData;
    OMSimRadioactivityData* radData;

    OMSimPMTConstruction* mPMTManager;
    mDOM* fMDOM;
    pDOM* fPDOM;
    LOM16* fLOM16;
    LOM18* fLOM18;
    dEGG* fDEGG;
    WOM* fWOM;

    //vars for drawing primaries from within a volume
    G4VSolid* fOuterSolid;
    G4VSolid* fInnerSolid;
    G4double fGlassWeight;
    G4double glassInRad;
    G4double glassOutRad;

    //world material
    G4int z;
    G4int nelements;
    G4double a;
    G4double density;
    G4Element* H;
    G4Element* O;
    G4Material* ice;
    G4int fDOM;
    G4double fworldSize;

};

#endif
//
