#ifndef OMSIMRADIOACTIVITYDATA_HH_INCLUDED
#define OMSIMRADIOACTIVITYDATA_HH_INCLUDED

#include "G4Types.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <iostream>
#include <boost/math/distributions/poisson.hpp>
#include <random>
#include <cmath>

class OMSimRadioactivityData
{
public:

    OMSimRadioactivityData();
    ~OMSimRadioactivityData();

    G4ThreeVector SetupPosition();
    G4ThreeVector SetupOrientation();

    void SetTimeWindow(G4double timeWindow);
    void SetActivity(G4double activity) { factivity = activity; }
    static void SetGlassRad(G4double, G4double);

    static G4double fglassOutRad;
    static G4double fglassInRad;
    static G4double ftimeWindow;
    static G4int fomModel;
    G4int GetNumDecay();
    G4double GetTimeWindow();
    G4double GetInitialTime();



private:

    G4double finitialTime;
    G4int fnumDecay;
    G4double factivity;
    G4ThreeVector fdecayPosition;
    G4ThreeVector fdecayOrientation;

    //Optical Module Specific
    G4double fglassWeight;
    G4double fglassRadius;
    G4double fglassTheta;
    G4double fglassPhi;

    void GeneratePosition();
    void GenerateOrientation();
    void GenerateFlatInitialTime();
    void GenerateExpInitialTime();
    G4double RandomGen(G4double, G4double);
    void InversePoisson();
};


#endif // OMSIMRADIOACTIVITYDATA_HH_INCLUDED
