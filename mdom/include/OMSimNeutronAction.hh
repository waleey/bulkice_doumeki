#ifndef OMSIMNEUTRONACTION_HH_INCLUDED
#define OMSIMNEUTRONACTION_HH_INCLUDED

#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>
#include <vector>
#include <iostream>

class OMSimNeutronAction
{
public:
    OMSimNeutronAction(G4ParticleGun*);
    ~OMSimNeutronAction();

    void GeneratePrimaries(G4Event*);
    inline bool NeutronExist() { return fParticleExist; }
    void LoadData();

private:
    G4ParticleGun* fParticleGun;
    bool fParticleExist;
    G4int fParticleNum;
    G4int fIdx;

    std::vector<G4double> energy;
    std::vector<G4double> fX;
    std::vector<G4double> fY;
    std::vector<G4double> fZ;
    std::vector<G4double> alpha_X;
    std::vector<G4double> alpha_Y;
    std::vector<G4double> alpha_Z;
    std::vector<G4double> inTime;
    std::vector<std::vector<G4double>> fNeutronData {energy, fX, fY, fZ, alpha_X, alpha_Y, alpha_Z, inTime};
    std::vector<std::string>  dtypes {"energy", "x", "y", "z", "ax", "ay", "az", "time"};
    enum {ENERGY, X, Y, Z, AX, AY, AZ, TIME};
};



#endif // OMSIMNEUTRONACTION_HH_INCLUDED
