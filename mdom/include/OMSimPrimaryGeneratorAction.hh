#ifndef OMSimPrimaryGeneratorAction_h
#define OMSimPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
//#include "G4ParticleTable.hh"

#include "OMSimPositronAction.hh"
#include "OMSimNeutronAction.hh"
#include "OMSimElectronAction.hh"
#include "OMSimK40Action.hh"
#include "OMSimU238Action.hh"
#include "OMSimU235Action.hh"
#include "OMSimTh232Action.hh"


class OMSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	OMSimPrimaryGeneratorAction();
	~OMSimPrimaryGeneratorAction();


public:
	void GeneratePrimaries(G4Event* anEvent);
	bool ParticleExist();
	void LoadData();
	void SetPosition(G4ThreeVector&);
	inline void SetActionType(G4int actionType) { fActionType = actionType; }

private:

	OMSimPositronAction* fPositronAction;
	OMSimNeutronAction* fNeutronAction;
	OMSimElectronAction* fElectronAction;
	OMSimK40Action* fK40Action;
	OMSimU238Action* fU238Action;
	OMSimU235Action* fU235Action;
	OMSimTh232Action* fTh232Action;
	G4ParticleGun *fParticleGun;
	G4int fActionType;

	enum {Positron, Neutron, Electron, K40, U238, U235, Th232};
};



#endif
