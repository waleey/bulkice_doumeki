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

class OMSimDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    OMSimDetectorConstruction(G4int, G4double);
    ~OMSimDetectorConstruction();
    G4VPhysicalVolume *Construct();

private:
    //G4Orb *mWorldSolid;
    G4Box *mWorldSolid;
    G4LogicalVolume *mWorldLogical;
    G4VPhysicalVolume *mWorldPhysical;
    void ConstructWorld();
    void ConstructWorldMat();
    OMSimInputData *mData;
    OMSimPMTConstruction* mPMTManager;

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
