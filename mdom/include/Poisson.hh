#ifndef POISSON_HH_INCLUDED
#define POISSON_HH_INCLUDED

#include <iostream>
#include <boost/math/distributions/poisson.hpp>
#include <random>
#include "G4Types.hh"

class Poisson
{
public:
    Poisson(G4double timeWindow, G4double activity) : mtimeWindow(timeWindow), mactivity(activity), mnumDecay(0), mstepSize(0), mminLim(0), mmaxLim(0) {}
    ~Poisson() {}

    G4int GetNumDecay();
    void SetTimeWindow(G4double timeWindow) { mtimeWindow = timeWindow ; }
    void SetActivity(G4double activity) { mactivity = activity ; }


private:

    G4double mtimeWindow;
    G4double mactivity;
    G4double mnumDecay;
    G4double mstepSize;
    G4double mminLim;
    G4double mmaxLim;

    G4double RandomGen();
    void InversePoisson();

};


#endif // POISSON_HH_INCLUDED
