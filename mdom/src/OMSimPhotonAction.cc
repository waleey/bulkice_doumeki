#include "OMSimPhotonAction.hh"
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <stdlib.h>

extern G4double gDistance;
extern G4double gAngle;

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
    G4ThreeVector direction = GenerateDirection(fZenith);
    fParticleGun -> SetParticleDefinition(particle);
    std::vector<G4double> energy;

/*
    G4String fileName = "../InputFile/energy_temp_photon_Ch.dat";
    std::ifstream file(fileName, std::ios::binary);
    G4double temp;
    if(!file.is_open())
    {
        G4cout << "Failed to open your photon file." << G4endl;
    }
    else
    {
        while(file.read((char*)&temp, sizeof(temp)))
        {
            energy.push_back(temp * eV);
        }
        file.close();
    }

    G4int nPhotons = energy.size();
*/

    G4double initialEnergy = 1 * eV;
    G4double finalEnergy = 8 * eV;
    G4int nPhotons = 10000000;

    std::cout << "Generating " << nPhotons << " between " << initialEnergy / eV << " eV and " << finalEnergy / eV << " eV." << std::endl;
    for(int i = 0; i < nPhotons; i++)
    {
        G4ThreeVector position = GeneratePosition();
        //G4double energyP = 3.55 * eV;
        G4double energyP = fRadData -> RandomGen(initialEnergy, finalEnergy);
        fParticleGun -> SetParticlePosition(position);
        fParticleGun -> SetParticleMomentumDirection(direction);
        fParticleGun -> SetParticleEnergy(energyP);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
    }





   /* std::cout <<
    "x: " << position.x() / m<< std::endl
    << "y: " << position.y() / m << std::endl
    << "z: " << position.z() / m << std::endl;*/



    //G4double energy = (fRadData -> RandomGen(1.7, 6.2)) * eV;

    //std::cout << "Energy: " << energy / eV<< std::endl;




    /*OMSimRadioactivityData* radData = new OMSimRadioactivityData();

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName = "opticalphoton";
    G4ParticleDefinition* particle = particleTable -> FindParticle(particleName);

    G4double theta = radData -> RandomGen(0, 2 * CLHEP::pi);
    G4double x = 0.5 * sin(theta);
    G4double y = 0.5 * cos(theta);
    G4double z = (radData -> RandomGen(-500, 500));
    x = x * m;
    y = y * m;
    z = z * mm;
std::cout <<
    "x: " << x / m<< std::endl
    << "y: " << y / m << std::endl
    << "z: " << z / m << std::endl;


    G4ThreeVector position(x, y, z);
    G4double angle = 90;
    G4double uz = cos(angle * deg);
    G4double alpha = (radData -> RandomGen(0, 2 * CLHEP::pi));
    G4double ux = (- sin(angle * deg) * cos(alpha * rad));
    G4double uy = (- sin(angle * deg) * sin(alpha * rad)) ;
    G4ThreeVector direction(ux, uy, uz);

    std::cout << "angle: " << angle << std::endl
    << "alpha: " << alpha / deg << std::endl
    << "ux: " << ux << std::endl
    << "uy: " << uy  << std::endl
    << "uz: " << uz  << std::endl;

    fParticleGun -> SetParticleDefinition(particle);
    fParticleGun -> SetParticleEnergy((radData -> RandomGen(1.7, 6.2) * eV));
    fParticleGun -> SetParticlePosition(position);
    fParticleGun -> SetParticleMomentumDirection(direction);

    delete radData;*/



}

G4ThreeVector OMSimPhotonAction::GeneratePosition()
{
   /* G4double theta = fRadData -> RandomGen(0, 2 * CLHEP::pi);
    G4double x = gDistance * sin(theta);
    G4double y = gDistance * cos(theta);
    G4double z = (fRadData -> RandomGen(-3000, 3000));
    x = x * m;
    y = y * m;
    z = z * mm;*/
    /*std::cout << "distance: " << gDistance << std::endl <<
    "x: " << x / m<< std::endl
    << "y: " << y / m << std::endl
    << "z: " << z / m << std::endl;*/

    G4double r = fRadData -> RandomGen(-1, 1);
    G4double y = - gDistance * sin(gAngle * deg) + r * cos(gAngle * deg);
    G4double z = - gDistance * cos(gAngle * deg) - r * sin(gAngle * deg);
    G4double x = 0;

    x = x * m;
    y = y * m;
    z = z * m;

    //std::cout << "+++ (PHOTON-ACTION): Position x=" << x << " m, y=" << y << " m, z=" << z << " m." << std::endl;

    return G4ThreeVector(x, y, z);
}

G4ThreeVector OMSimPhotonAction::GenerateDirection(G4double angle)
{
    /*G4double uz = cos(angle * deg);
    G4double alpha = (fRadData -> RandomGen(0, 2 * CLHEP::pi));
    G4double ux = (- sin(angle * deg) * cos(alpha * rad));
    G4double uy = (- sin(angle * deg) * sin(alpha * rad)) ;*/

   /* std::cout << "angle: " << angle << std::endl
    << "alpha: " << alpha / deg << std::endl
    << "ux: " << ux << std::endl
    << "uy: " << uy  << std::endl
    << "uz: " << uz  << std::endl;*/

    G4double uz = cos(angle * deg);
    G4double uy = sin(angle * deg);
    G4double ux = 0;

    return G4ThreeVector(ux, uy, uz);
}
