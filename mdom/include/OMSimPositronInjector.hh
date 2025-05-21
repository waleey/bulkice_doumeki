#ifndef OMSIMPOSITRONINJECTOR_HH_INCLUDED
#define OMSIMPOSITRONINJECTOR_HH_INCLUDED

#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "OMSimRadioactivityData.hh"

#include <fstream>
#include <vector>
#include <iostream>

class OMSimPositronInjector
{
public:
    OMSimPositronInjector(G4ParticleGun*);
    ~OMSimPositronInjector();

    void GeneratePrimaries(G4Event*);
    inline bool PositronExist() { return fParticleExist; }
    void LoadData();


private:

    G4double neutronMass = 939.565 * MeV;
    G4double protonMass = 938.272 * MeV;
    G4double electronMass = 0.511 * MeV;
    G4double IBDThresholdEnergy = (neutronMass-protonMass) + electronMass;
    G4double a0 = -0.1;

    G4ParticleGun* fParticleGun;
    bool fParticleExist;
    G4int fParticleNum;
    G4int fIdx;
    G4int fPerc;

    OMSimRadioactivityData* fRadData;

    std::vector<G4double> energy;
    std::vector<G4double> fX;
    std::vector<G4double> fY;
    std::vector<G4double> fZ;
    std::vector<G4double> alpha_X;
    std::vector<G4double> alpha_Y;
    std::vector<G4double> alpha_Z;
    std::vector<G4double> inTime;
    std::vector<std::vector<G4double>> fPositronData {energy, fX, fY, fZ, alpha_X, alpha_Y, alpha_Z, inTime};
    std::vector<std::string>  dtypes {"energy", "x", "y", "z", "ax", "ay", "az", "time"};
    enum {ENERGY, X, Y, Z, AX, AY, AZ, TIME};

    G4double SampleEnergy(G4double, G4double);
    G4double SampleZenith(G4double, G4double, G4String, G4int);
    G4ThreeVector rotateFrameFromeTo(const G4ThreeVector&, const G4ThreeVector&, const G4ThreeVector&);

};


#endif // OMSIMPOSITRONINJECTOR_HH_INCLUDED
