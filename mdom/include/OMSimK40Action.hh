#ifndef OMSIMK40ACTION_HH_INCLUDED
#define OMSIMK40ACTION_HH_INCLUDED

#include "G4ThreeVector.hh"
#include "OMSimRadioactivityData.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"

class OMSimK40Action
{
public:
    OMSimK40Action(G4ParticleGun*);
    ~OMSimK40Action();

    void GeneratePrimaries(G4Event*);
    inline void SetPosition(G4ThreeVector& position) { fPosition = position; }

private:
    G4ThreeVector fPosition;
    G4ParticleGun* fParticleGun;
    OMSimRadioactivityData* fRadData;
};

#endif // OMSIMK40ACTION_HH_INCLUDED
