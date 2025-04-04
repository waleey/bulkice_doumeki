#ifndef OMSIMDECAYCHAINACTION_HH_INCLUDED
#define OMSIMDECAYCHAINACTION_HH_INCLUDED

#include "G4ThreeVector.hh"
#include "OMSimRadioactivityData.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"

class OMSimDecayChainAction
{
public:
    OMSimDecayChainAction();
    OMSimDecayChainAction(G4ParticleGun*);
    ~OMSimDecayChainAction();

    void GenerateIsotope();
    void GeneratePrimaries(G4Event*);
    static std::vector<std::tuple<std::string, G4int, G4int, G4double, G4int, G4double>> GetDecayChain();

    inline void SetPosition(G4ThreeVector& position) { fPosition = position; }
    inline void SetZ(G4int z) { fZ = z; }
    inline void SetA(G4int a) { fA = a; }
    inline void SetTimeLow(G4double timeLow) {ftimeLow = timeLow; }
    inline void SetTimeHigh(G4double timeHigh) {ftimeHigh = timeHigh; }
    inline void SetExcitationEnergy(G4double excitationEnergy) { fexcitationEnergy = excitationEnergy; }
    inline void SetTotalAngularMomentum(G4int totalAngularMomentum) { ftotalAngularMomentum = totalAngularMomentum; }
    inline void SetPDGLifeTime(G4double meanLifeTime) {fisotope -> SetPDGLifeTime(meanLifeTime); }
    inline G4double GetPDGLifeTime() {return fisotope -> GetPDGLifeTime(); }

private:
    G4ParticleDefinition* fisotope;
    G4ThreeVector fPosition;
    G4int fZ;
    G4int fA;
    G4double ftimeLow;
    G4double ftimeHigh;
    G4double fexcitationEnergy;
    G4int ftotalAngularMomentum;
    G4ParticleGun* fParticleGun;
    OMSimRadioactivityData* fRadData;
};


#endif // OMSIMDECAYCHAINACTION_HH_INCLUDED
