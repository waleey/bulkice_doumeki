#include "OMSimPhotonAction.hh"
#include <cmath>

extern G4double gDistance;

OMSimPhotonAction::OMSimPhotonAction(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun),
        fZenith(0)
{
    fRadData = new OMSimRadioactivityData();
}
OMSimPhotonAction::~OMSimPhotonAction()
{
    delete fRadData;
}
void OMSimPhotonAction::GeneratePrimaries(G4Event* anEvent)
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable -> FindParticle("opticalphoton");

    G4ThreeVector position = GeneratePosition();
    G4ThreeVector direction = GenerateDirection(fZenith);
    G4double energy = (fRadData -> RandomGen(1.775, 6.212)) * eV;

    fParticleGun -> SetParticlePosition(position);
    fParticleGun -> SetParticleMomentumDirection(direction);
    fParticleGun -> SetParticleEnergy(energy);
    fParticleGun -> SetParticleDefinition(particle);
    fParticleGun -> GeneratePrimaryVertex(anEvent);

}

G4ThreeVector OMSimPhotonAction::GeneratePosition()
{
    G4double x = (fRadData -> RandomGen(0., gDistance));
    G4double y = sqrt(pow(gDistance, 2) - pow(x, 2));
    G4double z = (fRadData -> RandomGen(-600, 600));
    x = x * m;
    y = y * m;
    z = z * mm;
    std::cout << "distance: " << gDistance << std::endl <<
    "x: " << x / m<< std::endl
    << "y: " << y / m << std::endl
    << "z: " << z / m << std::endl;

    return G4ThreeVector(x * m, y * m, z * mm);
}

G4ThreeVector OMSimPhotonAction::GenerateDirection(G4double angle)
{
    G4double uz = cos(angle * deg);
    G4double alpha = (fRadData -> RandomGen(0, 2 * CLHEP::pi));
    G4double ux = (- sin(angle * deg) * cos(alpha));
    G4double uy = (- sin(angle * deg) * sin(alpha)) ;

    std::cout << "angle: " << angle << std::endl
    << "alpha: " << alpha / deg << std::endl
    << "ux: " << ux << std::endl
    << "uy: " << uy  << std::endl
    << "uz: " << uz  << std::endl;


    return G4ThreeVector(ux, uy, uz);
}
