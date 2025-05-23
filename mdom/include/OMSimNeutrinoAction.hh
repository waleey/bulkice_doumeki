#ifndef OMSIMNEUTRINOACTION_HH_INCLUDED
#define OMSIMNEUTRINOACTION_HH_INCLUDED

#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "OMSimRadioactivityData.hh"

#include <fstream>
#include <vector>
#include <iostream>

class OMSimNeutrinoAction
{
public:
    OMSimNeutrinoAction(G4ParticleGun*);
    ~OMSimNeutrinoAction();

    void GeneratePrimaries(G4Event*);
    void LoadData();
    inline bool NeutrinoExist() { return fParticleExist; }

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

    G4double SampleEnergy(G4double, G4double);
    G4double SampleZenith(G4double, G4double, G4String, G4int);
    G4ThreeVector rotateFrameFromeTo(const G4ThreeVector&, const G4ThreeVector&, const G4ThreeVector&);

};


#endif // OMSIMNEUTRINOACTION_HH_INCLUDED