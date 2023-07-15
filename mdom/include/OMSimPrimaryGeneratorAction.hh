#ifndef OMSimPrimaryGeneratorAction_h
#define OMSimPrimaryGeneratorAction_h 1

//#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "OMSimParticleSetup.hh"

#include <vector>
#include <string>
//#include "globals.hh"

/*class G4GeneralParticleSource;
class G4Event;*/


class OMSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	OMSimPrimaryGeneratorAction(G4String&, G4int);
	OMSimPrimaryGeneratorAction(G4String&, G4int, OMSimParticleSetup*);
	~OMSimPrimaryGeneratorAction();


public:
	void GeneratePrimaries(G4Event* anEvent);

private:
	//G4GeneralParticleSource* particleSource;

	//creating particle gun and make it read from sntools output files

	//void SetUpEnergyAndPosition();

	G4ParticleGun *fParticleGun;
	std::string finteraction_channel;
	G4int fomModel;
	OMSimParticleSetup* fParticleSetup;
};



#endif
