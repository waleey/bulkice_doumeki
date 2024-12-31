#include "OMSimDecayChainAction.hh"


OMSimDecayChainAction::OMSimDecayChainAction()   
{
    fRadData = new OMSimRadioactivityData();
}

OMSimDecayChainAction::OMSimDecayChainAction(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun)
{
    fRadData = new OMSimRadioactivityData();
}

OMSimDecayChainAction::~OMSimDecayChainAction()
{
    delete fRadData;
}


void OMSimDecayChainAction::GenerateIsotope() { 

    // Find and configure the ion
    fisotope = G4IonTable::GetIonTable()->GetIon(fZ, fA, fexcitationEnergy * keV);
    if (!fisotope) {
        std::cerr << "Error: Could not find ion Z=" << fZ << ", A=" << fA << ", E=" << fexcitationEnergy *1000 << " keV" << std::endl;
        return;
    }
}

void OMSimDecayChainAction::GeneratePrimaries(G4Event* anEvent) {
    
    G4double energy = 0. * keV;
    G4double charge = 0. * eplus;
    // random direction
    G4ThreeVector orientation = fRadData -> SetupOrientation();
    // random time
    G4double initialTime = fRadData -> GetInitialTime() * s;

    try
    {
        fParticleGun -> SetParticlePosition(fPosition);
    }
    catch(...)
    {
        std::cerr << "Error in particle Position!!" << std::endl
        << "OMSimDecayChainAction::GeneratePrimaries()" << std::endl;
    }
    
    // Set up the particle gun
    fParticleGun -> SetParticleDefinition(fisotope);
    fParticleGun -> SetParticleCharge(charge);
    fParticleGun -> SetParticleTime(initialTime);
    fParticleGun -> SetParticleMomentumDirection(orientation);
    fParticleGun -> SetParticleEnergy(energy);
    fParticleGun -> GeneratePrimaryVertex(anEvent);
}

