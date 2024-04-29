#ifndef OMSimPrimaryGeneratorAction_h
#define OMSimPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
//#include "G4ParticleTable.hh"

#include "OMSimPositronAction.hh"
#include "OMSimNeutronAction.hh"
#include "OMSimElectronAction.hh"
#include "OMSimK40Action.hh"
#include "OMSimU238Action.hh"
#include "OMSimU235Action.hh"
#include "OMSimTh232Action.hh"
#include "OMSimPhotonAction.hh"


class PrimaryGeneratorMessenger;

class OMSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	OMSimPrimaryGeneratorAction();
	OMSimPrimaryGeneratorAction(G4String&);
	~OMSimPrimaryGeneratorAction();


public:
	void GeneratePrimaries(G4Event* anEvent);
	bool ParticleExist();
	void LoadData();
	void SetPosition(G4ThreeVector&);
	inline void SetActionType(G4int actionType) { fActionType = actionType; }
	void SetAngle(G4double angle);
	void SetupMuons(G4int);

    //Setters for messenger class
	inline void SetParticleName(G4String& name) {fParticleName = name; }
    inline void SetParticleDistance(G4double distance) {fDistance = distance; }
    inline void SetParticleAngle(G4double angle) {fAngle = angle; }
    inline void SetParticleEnergy(G4double energy) {fEnergy = energy; }
    inline void SetBeamPosition(G4double z) {fZPosition = z; }
    inline void SetStartEnergy(G4double energy) {fStartEnergy = energy; }
    inline void SetFinalEnergy(G4double energy) {fFinalEnergy = energy; }

private:
    void GenerateToVisualize();
    void GenerateWave();
    void GenerateBeam();

	OMSimPositronAction* fPositronAction;
	OMSimNeutronAction* fNeutronAction;
	OMSimElectronAction* fElectronAction;
	OMSimK40Action* fK40Action;
	OMSimU238Action* fU238Action;
	OMSimU235Action* fU235Action;
	OMSimTh232Action* fTh232Action;
	OMSimPhotonAction* fPhotonAction;
	G4ParticleGun *fParticleGun;
	G4int fActionType;

	G4String fParticleName;
	G4double fDistance;
	G4double fAngle;
	G4double fEnergy;
	G4double fStartEnergy;
	G4double fFinalEnergy;
	G4double fZPosition;
	G4String fInteraction;
    G4double fPhotonAngle;

    PrimaryGeneratorMessenger* fGeneratorMessenger;

	enum {Positron, Neutron, Electron, K40, U238, U235, Th232, Photon, Visualization, wave, beam, muon};
};



#endif
