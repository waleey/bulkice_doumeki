#ifndef OMSIMPHOTONACTION_HH_INCLUDED
#define OMSIMPHOTONACTION_HH_INCLUDED

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "OMSimRadioactivityData.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"

class OMSimPhotonAction
{
public:
    OMSimPhotonAction(G4ParticleGun*);
    ~OMSimPhotonAction();

    void GeneratePrimaries(G4Event* );
    inline void SetAngle(G4double angle) { fZenith = angle; }

private:
    G4ThreeVector GeneratePosition();
    G4ThreeVector GenerateDirection(G4double zenith);

    G4ParticleGun* fParticleGun;
    OMSimRadioactivityData* fRadData;

    G4double fZenith;
};


#endif // OMSIMPHOTONACTION_HH_INCLUDED
