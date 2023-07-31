#include "OMSimTh232Action.hh"

OMSimTh232Action::OMSimTh232Action(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun)
{
    fRadData = new OMSimRadioactivityData();
}
OMSimTh232Action::~OMSimTh232Action()
{
    delete fRadData;
}
void OMSimTh232Action::GeneratePrimaries(G4Event* anEvent)
{
    G4int Z = 90;
    G4int A = 232;
    G4double energy = 0. * keV;
    G4double charge = 0. * eplus;
    G4ParticleDefinition* particle = G4IonTable::GetIonTable() -> GetIon(Z, A, energy);
    particle -> SetPDGLifeTime(0. * ns);
    G4ThreeVector orientation = fRadData -> SetupOrientation();
    G4double initialTime = fRadData -> GetInitialTime() * s;

    try
    {
        fParticleGun -> SetParticlePosition(fPosition);
    }
    catch(...)
    {
        std::cerr << "Error in particle Position!!" << std::endl
        << "OMSimU235Action::GeneratePrimaries()" << std::endl;
    }
    fParticleGun -> SetParticleMomentumDirection(orientation);
    fParticleGun -> SetParticleEnergy(energy);
    fParticleGun -> SetParticleCharge(charge);
    fParticleGun -> SetParticleDefinition(particle);
    fParticleGun -> SetParticleTime(initialTime);
    fParticleGun -> GeneratePrimaryVertex(anEvent);
}
