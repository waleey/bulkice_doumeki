#include "OMSimPrimaryGeneratorAction.hh"


#include "G4Event.hh"
#include "G4ParticleTypes.hh"

#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdlib.h>

extern G4double gworldsize;


OMSimPrimaryGeneratorAction::OMSimPrimaryGeneratorAction(G4String& interaction_channel, G4int omModel) : finteraction_channel(interaction_channel), fomModel(omModel)
{

	fParticleGun = new G4ParticleGun(1);
	G4cout << ":::::::::::::::Particle Gun Created:::::::::::" << G4endl;
}

OMSimPrimaryGeneratorAction::~OMSimPrimaryGeneratorAction()
{
	delete fParticleGun;
}

void OMSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	OMSimParticleSetup* fParticleSetup = new OMSimParticleSetup(fParticleGun, anEvent, fomModel);

   try
    {
        if(finteraction_channel == "ibd")
        {
	    std::cout << "************Only particles from IBD is considered************* " << std::endl;
	    fParticleSetup -> GeneratePositron();
            fParticleSetup -> GenerateNeutron();
        }
        else if(finteraction_channel == "enees")
        {
            std::cout << "************Only particles from ENEES is considered***********" << std::endl;
            fParticleSetup -> GenerateElectron();
            fParticleSetup -> GenerateNeutron();
        }
        else if(finteraction_channel == "all")
        {
            std::cout << "**********Particles from all interactions are considered********" << std::endl;
            fParticleSetup -> GeneratePositron();
            fParticleSetup -> GenerateElectron();
            fParticleSetup -> GenerateNeutron();
        }
        /**
        *This is the block that shoots the radioactive element from inside the glass materials
        **/

        else if(finteraction_channel == "radioactivity")
        {
            std::cout << "********Only particles from RadioactiveDecay is being simulated**********" << std::endl;
            fParticleSetup -> GenerateK40();
            fParticleSetup -> GenerateTh238();
            fParticleSetup -> GenerateU238();
            fParticleSetup -> GenerateU235();
        }
        else
        {
            throw finteraction_channel;
        }
    }
	catch(...)
	{
        std::cout << "Invalid Interaction Channel! Aborting..." << std::endl;
	}

	/*fParticleSetup -> GenerateK40();
    fParticleSetup -> GenerateTh238();
    fParticleSetup -> GenerateU238();
    fParticleSetup -> GenerateU235();*/




}
