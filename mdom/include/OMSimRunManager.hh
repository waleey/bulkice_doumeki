#ifndef OMSIMRUNMANAGER_HH_INCLUDED
#define OMSIMRUNMANAGER_HH_INCLUDED

#include "G4RunManager.hh"

#include "G4VisExecutive.hh"
#include "G4RayTracer.hh"
#include <vector>

#include "OMSimDetectorConstruction.hh"
#include "OMSimPhysicsList.hh"
#include "OMSimPrimaryGeneratorAction.hh"
#include "OMSimRunAction.hh"
#include "OMSimEventAction.hh"
#include "OMSimTrackingAction.hh"
#include "OMSimSteppingAction.hh"
#include "OMSimSteppingVerbose.hh"
#include "OMSimRadioactivityData.hh"

class OMSimRunManager
{
public:
    OMSimRunManager();
    OMSimRunManager(G4int, G4double, G4String&);
    ~OMSimRunManager();

    void Initialize();
    void BeamOn();
    void OpenFile();
    void CloseFile();

    /*inline void SetZenithAngle(G4double angle) { fZenithAngle = angle; }
    inline void SetStartAngle(G4double angle) { fStartAngle = angle; }
    inline void SetFinalAngle(G4double angle)   { fFinalAngle = angle; }
    inline void SetAngleIncrement(G4double angle) { fAngleIncrement = angle; }
    inline void SetDistance(G4double distance) { fDistance = distance; }*/
    

private:
    G4int fpmtModel;
    G4double fworldSize;
    G4String fInteraction;
//  mandatory classes
    G4RunManager* fRunManager;
    OMSimDetectorConstruction* fDetectorConstruction;
    OMSimPrimaryGeneratorAction* fPrimaryGenerator;
    G4VisManager* fVisManager;
    OMSimPhysicsList* fPhysicsList;
    OMSimRunAction* fRunAction;
    OMSimEventAction* fEventAction;
    OMSimTrackingAction* fTrackingAction;
    OMSimSteppingAction* fSteppingAction;
    OMSimRadioactivityData* fRadData;

//  Select Action Type
    G4int fActionType;

//  Generator functions
    void GeneratePositron();
    void GenerateNeutron();
    void GenerateElectron();
    void GenerateK40();
    void GenerateU238();
    void GenerateU235();
    void GenerateTh232();
    void GeneratePhoton();
    void GenerateGamma();
    void GenerateToVisualize();

    /**
    *Activity is given in Bq/Kg unit
    **/
    std::vector<G4double> K40Activity {0, 61, 61, 61, 61, 2.548, 61};
    std::vector<G4double> U238Activity {0, 4.61, 4.61, 4.61, 4.61, 0.2146, 4.61};
    std::vector<G4double> U235Activity {0, 0.60, 0.60, 0.60, 0.60, 0, 0.60};
    std::vector<G4double> Th232Activity {0, 1.28, 1.28, 1.28, 1.28, 0.1373, 1.28};
    std::vector<G4double> glassWeight {10., 13.0, 10.0, 10.0, 10.0, 13.0, 13.0}; //for now, weight of Degg is assumed to be equal to MDOM. Will be changed later.

    enum {Positron, Neutron, Electron, K40, U238, U235, Th232, Photon, Visualization, wave, beam, Gamma}; //make sure this order is the same in primary generator class enum
    enum {pmt, mdom, dom, lom16, lom18, degg, wom};

    G4double startingtime;
    G4double finishtime;

    //optical photon simulation variables
  /*  G4double fZenithAngle;
    G4double fStartAngle;
    G4double fFinalAngle;
    G4double fAngleIncrement;
    G4double fDistance;*/
};

#endif // OMSIMRUNMANAGER_HH_INCLUDED
