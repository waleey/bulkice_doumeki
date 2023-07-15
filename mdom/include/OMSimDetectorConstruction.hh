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
#include "OMSimInputData.hh"
#include "OMSimPMTConstruction.hh"
#include "OMSimRadioactivityData.hh"
#include "OMSimMDOM.hh"
#include "OMSimPDOM.hh"
#include "OMSimLOM16.hh"
#include "OMSimLOM18.hh"
#include "OMSimDEGG.hh"
#include "OMSimParticleSetup.hh"

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

    inline void SetOMModel(G4int model) { fDOM = model; }
    inline void SetWorldSize(G4double worldSize) { fworldSize = worldSize; }

    void DrawFromVolume();
    inline std::vector<std::vector<G4ThreeVector>>& GetPositions() { return isotopePos; }


private:
    //G4Orb *mWorldSolid;
    G4Box *mWorldSolid;
    G4LogicalVolume *mWorldLogical;
    G4VPhysicalVolume *mWorldPhysical;
    void ConstructWorld();
    void ConstructWorldMat();
    OMSimInputData *mData;

    OMSimPMTConstruction* mPMTManager;
    mDOM* fMDOM;
    pDOM* fPDOM;
    LOM16* fLOM16;
    LOM18* fLOM18;
    dEGG* fDEGG;

    //vars for drawing primaries from within a volume
    G4VSolid* fSolid;
    G4double fGlassWeight;
    G4double glassInRad;
    G4double glassOutRad;

    std::vector<G4ThreeVector> K40_pos;
    std::vector<G4ThreeVector> U238_pos;
    std::vector<G4ThreeVector> U235_pos;
    std::vector<G4ThreeVector> Th232_pos;

    std::vector<std::vector<G4ThreeVector>> isotopePos {K40_pos, U238_pos, U235_pos, Th232_pos};

    void GenerateInVolume(G4int);

    enum {K40, U238, U235, Th232};
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

    //ice properties
    G4int numEntries;
    std::vector<G4double> energyice;
    std::vector<G4double> rindexice;
    std::vector<G4double> absorptionice;
    std::vector<G4double> mieice;

    G4MaterialPropertiesTable *mptice;

};

#endif
//
