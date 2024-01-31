#ifndef PRIMARYGENERATORMESSENGER_HH_INCLUDED
#define PRIMARYGENERATORMESSENGER_HH_INCLUDED

#include "G4SystemOfUnits.hh"
#include "G4UImessenger.hh"
/*
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "OMSimPrimaryGeneratorAction.hh"
*/

class OMSimPrimaryGeneratorAction;

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;


class PrimaryGeneratorMessenger : public G4UImessenger
{

    public:
        PrimaryGeneratorMessenger(OMSimPrimaryGeneratorAction*);
        ~PrimaryGeneratorMessenger() override;

        void SetNewValue(G4UIcommand*, G4String) override;

    private:
        OMSimPrimaryGeneratorAction* fPrimaGen;

        G4UIdirectory* fDirectory;

        G4UIcmdWithAString* fActionType;
        G4UIcmdWithAString* fParticleType;
        G4UIcmdWithADouble* fDistance;
        G4UIcmdWithADoubleAndUnit* fEnergy; //This produces a monoenergetic beam
        G4UIcmdWithADoubleAndUnit* fStartEnergy; //this samples energy from a uniform distribution
        G4UIcmdWithADoubleAndUnit* fFinalEnergy; //this samples energy from a uniform distribution
        G4UIcmdWithADouble* fZPosition;
        G4UIcmdWithADouble* fAngle; //incident angle of the beam.



};
#endif // PRIMARYGENERATORYMESSENGER_HH_INCLUDED
