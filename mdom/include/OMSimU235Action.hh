#ifndef OMSIMU235ACTION_HH_INCLUDED
#define OMSIMU235ACTION_HH_INCLUDED

#include "G4ThreeVector.hh"
#include "OMSimRadioactivityData.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"

class OMSimU235Action
{
public:
    OMSimU235Action(G4ParticleGun*);
    ~OMSimU235Action();

    void GeneratePrimaries(G4Event*);
    inline void SetPosition(G4ThreeVector& position) { fPosition = position; }

private:
    G4ThreeVector fPosition;
    G4ParticleGun* fParticleGun;
    OMSimRadioactivityData* fRadData;
};


#endif // OMSIMU235ACTION_HH_INCLUDED
