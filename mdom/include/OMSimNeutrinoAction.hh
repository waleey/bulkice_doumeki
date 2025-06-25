#ifndef OMSIMNEUTRINOACTION_HH_INCLUDED
#define OMSIMNEUTRINOACTION_HH_INCLUDED

#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "OMSimRadioactivityData.hh"

#include <fstream>
#include <vector>
#include <iostream>
#include <numeric>
#include "spline.h"

extern G4double gSimulationTime;
extern G4double gworldsize;
extern G4double gNeutrinoFlux;
extern G4double gNeutrinoMeanEnergy;
extern G4double gNeutrinoEnergyPinch;

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
    G4double iceDensity = 0.92 * g/(cm*cm*cm);
    G4double waterMolarMass = 18 * g / mole;
    G4double avogadroNumber = 6.022E23 / mole;

    G4double IBDThresholdEnergy = (std::pow(neutronMass + electronMass, 2) - std::pow(protonMass, 2))/(2 * protonMass);
    G4double a0 = -0.1;

    G4ParticleGun* fParticleGun;
    bool fParticleExist;
    G4int fParticleNum;
    G4int fIdx;
    G4int fPerc;

    // copy global variables
    G4double fSimulationTime = gSimulationTime;
    G4double fworldsize = gworldsize * m;
    G4double fNeutrinoFlux = gNeutrinoFlux;
    G4double fNeutrinoMeanEnergy = gNeutrinoMeanEnergy;
    G4double fNeutrinoEnergyPinch = gNeutrinoEnergyPinch;

    tk::spline fInverseCDF;

    OMSimRadioactivityData* fRadData;

    G4double SampleNeutrinoEnergy();
    G4double SamplePositronZenith(G4double, G4double, G4String, G4int);

    void InverseCDFNeutrinoSpectrumCrossSection();
    std::vector<double> CDFNeutrinoSpectrumCrossSection(const std::vector<double>&);
    std::pair<std::vector<double>, std::vector<double>> NormalizedPDFNeutrinoSpectrumCrossSection(G4double, G4double, G4int);

    std::pair<G4double, G4double> NumberNeutrinoInteractions();
    G4double IntegralNeutrinoSpectrumCrossSection(G4double, G4double, G4int);
    G4double NeutrinoBlackBodySpectrum(G4double, G4double, G4double);
    
    G4double IBDCrossSection(G4double, G4int);
    G4double IBDCrossSectionZeroth(G4double);
    G4double IBDCrossSectionFirst(G4double);

    G4double IBDPositronEnergyZeroth(G4double);
    G4double IBDPositronMomentumZeroth(G4double);
    G4double IBDPositronVeZeroth(G4double);
    G4double IBDPositronEnergyFirst(G4double, G4double, G4double);

    G4double TargetNumber();
    G4double SimulationVolume();

    G4ThreeVector rotateFrameFromeTo(const G4ThreeVector&, const G4ThreeVector&, const G4ThreeVector&);

};


#endif // OMSIMNEUTRINOACTION_HH_INCLUDED