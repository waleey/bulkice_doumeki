#include "OMSimU235Action.hh"

OMSimU235Action::OMSimU235Action(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun)
{
    fRadData = new OMSimRadioactivityData();
}
OMSimU235Action::~OMSimU235Action()
{
    delete fRadData;
}
void OMSimU235Action::GeneratePrimaries(G4Event* anEvent)
{
    G4int Z = 92;
    G4int A = 235;
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
