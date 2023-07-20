#ifndef OMSIMU238ACTION_HH_INCLUDED
#define OMSIMU238ACTION_HH_INCLUDED

#include "G4ThreeVector.hh"
#include "OMSimRadioactivityData.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"

class OMSimU238Action
{
public:
    OMSimU238Action(G4ParticleGun*);
    ~OMSimU238Action();

    void GeneratePrimaries(G4Event*);
    inline void SetPosition(G4ThreeVector& position) { fPosition = position; }

private:
    G4ThreeVector fPosition;
    G4ParticleGun* fParticleGun;
    OMSimRadioactivityData* fRadData;
};


#endif // OMSIMU238ACTION_HH_INCLUDED
