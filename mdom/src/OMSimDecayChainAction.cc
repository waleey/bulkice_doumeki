#include "OMSimDecayChainAction.hh"
#include "CLHEP/Random/RandFlat.h"

extern G4bool gVerbose;
extern G4double gSimulationTime;

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
    fisotope = G4IonTable::GetIonTable()->GetIon(fZ, fA, fexcitationEnergy * keV, ftotalAngularMomentum);
    if (!fisotope) {
        std::cerr << "Error: Could not find ion Z=" << fZ << ", A=" << fA << ", E=" << fexcitationEnergy / keV << " keV, J=" << ftotalAngularMomentum << std::endl;
        return;
    }
}

void OMSimDecayChainAction::GeneratePrimaries(G4Event* anEvent) {
    
    G4double energy = 0. * keV;
    G4double charge = 0. * eplus;
    // random direction
    G4ThreeVector orientation = fRadData -> SetupOrientation();
    G4double timeWindow = gSimulationTime;
    G4double meanLifeTime = fisotope -> GetPDGLifeTime();

    G4double initialTime = CLHEP::RandFlat::shoot(ftimeLow, ftimeHigh);
    
    if (gVerbose){
    std::cout << "+++ (DECAY):"
              << " ||| Mean Life Time [ns] = " << meanLifeTime
              << " ||| Time Window [ns] = " << timeWindow
              << " ||| Time Low [ns] = " << ftimeLow
              << " ||| Time High [ns] = " << ftimeHigh
              << " ||| Initial Time [ns] : " << initialTime << std::endl;
    }
    // now that initial time is set, set life time to zero
    fisotope -> SetPDGLifeTime(0 * ns);
    
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

