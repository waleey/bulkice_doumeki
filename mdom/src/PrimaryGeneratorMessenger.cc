#include "PrimaryGeneratorMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "OMSimPrimaryGeneratorAction.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(OMSimPrimaryGeneratorAction* primaGen)
                          : fPrimaGen(primaGen)
{
    fDirectory = new G4UIdirectory("/particleGen/");
    fDirectory -> SetGuidance("Particle generator action setup commands");

    //setting up the action type.
    fActionType = new G4UIcmdWithAString("/particleGen/shape", this);
    fActionType -> SetGuidance("Setting up the gun shape");
    fActionType -> SetGuidance("[usage] /particleGen/shape wave for particle wave");
    fActionType -> SetGuidance("[usage] /particleGen/shape beam for particle beam");
    fActionType -> SetParameterName("beamShape", false);


    //setting up the particle type.
    fParticleType = new G4UIcmdWithAString("/particleGen/particle", this);
    fParticleType -> SetGuidance("Setting up the particle type.");
    fParticleType -> SetGuidance("[usage] /particleGen/particle opticalphoton for photons");
    fParticleType -> SetGuidance("[usage] /particleGen/particle e+ for positrons");
    fParticleType -> SetGuidance("[usage] /particleGen/particle e- for electrons");
    fParticleType -> SetParameterName("ParticleType", true);
    fParticleType -> SetGuidance("If no particle type is chosen, optical photons would be used.");
    fParticleType -> SetDefaultValue("opticalphoton");

    //setting up the angle for plane waves
    fAngle = new G4UIcmdWithADouble("/particleGen/angle/", this);
    fAngle -> SetGuidance("Angle of incidence in degrees..");
    fAngle -> SetGuidance("Angle must be within 0 to 180 degrees.");
    fAngle -> SetGuidance("[usage] /particleGen/angle [your angle]");
    fAngle -> SetParameterName("incAngle", true);
    fAngle -> SetRange("incAngle >= 0 && incAngle <= 180");
    fAngle -> SetGuidance("If no angle is provided, 45 degree will be used.");
    fAngle -> SetDefaultValue(45.);

    //setting up the distance from the center
    fDistance = new G4UIcmdWithADouble("/particleGen/distance/", this);
    fDistance -> SetGuidance("Distance of the source from the center.");
    fDistance -> SetGuidance("Distance must be within 0 to 5 in meter. (user can change the limit from the code later)");
    fDistance -> SetGuidance("[usage] /particleGen/distance [distance]");
    fDistance -> SetParameterName("Distance", true);
    fDistance -> SetRange("Distance >= 0 && Distance <= 5");
    fDistance -> SetDefaultValue(2);

    //setting up the energy range
    fStartEnergy = new G4UIcmdWithADoubleAndUnit("/particleGen/startEnergy/", this);
    fStartEnergy -> SetGuidance("lower limit in energy range in eV.");
    fStartEnergy -> SetGuidance("[usage] /particleGen/startEnergy [energy] [unit]");
    fStartEnergy -> SetParameterName("StartEnergy", true);
    fStartEnergy -> SetGuidance("If no energy is provided, 0 eV will be assumed.");
    fStartEnergy -> SetDefaultValue(0);
    fStartEnergy -> SetDefaultUnit("eV");

    //setting up the energy range
    fFinalEnergy = new G4UIcmdWithADoubleAndUnit("/particleGen/finalEnergy/", this);
    fFinalEnergy -> SetGuidance("upper limit in energy range in eV.");
    fFinalEnergy -> SetGuidance("[usage] /particleGen/finalEnergy [energy] [unit]");
    fFinalEnergy -> SetParameterName("FinalEnergy", true);
    fFinalEnergy -> SetGuidance("If no energy value is provided, 0 eV will be assumed");
    fFinalEnergy -> SetDefaultValue(0);
    fFinalEnergy -> SetDefaultUnit("eV");

    //setting up energy for mono-energetic beam
    fEnergy = new G4UIcmdWithADoubleAndUnit("/particleGen/fEnergy", this);
    fEnergy -> SetGuidance("Setting up energy for monoenergetic beams.");
    fEnergy -> SetGuidance("[usage] /particleGen/fEnergy [energy] [unit]");
    fEnergy -> SetParameterName("Energy", true);
    fEnergy -> SetGuidance("If no energy value is provided, 3.5 eV will be assumed");
    fEnergy -> SetDefaultValue(3.5);
    fEnergy -> SetDefaultUnit("eV");

    //z position of the pencil beam
    fZPosition = new G4UIcmdWithADouble("/particleGen/z/", this);
    fZPosition -> SetGuidance("z coordiante of the pencil beam in cm");
    fZPosition -> SetGuidance("[usage] /particleGen/z [coordinate]");
    fZPosition -> SetParameterName("Z", true);
    fZPosition -> SetRange("Z >= -550 && Z <= 550");
    fZPosition -> SetGuidance("If no value is provided, 0 cm will be assumed.");
    fZPosition -> SetDefaultValue(0);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{

}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    /**
    Needs implementation!!!!
    **/
    if(command == fActionType)
    {
        if(newValue == "wave")
        {
            fPrimaGen -> SetActionType(8); //need to check action number
        }
        else if(newValue == "beam")
        {
            fPrimaGen -> SetActionType(9); //need to check action number.
        }
        else
        {
            G4cerr << "Invalid action type provided!! Aborting..." << G4endl;
            exit(0);
        }
    }
    if (command == fParticleType)
    {
        fPrimaGen -> SetParticleName(newValue);
    }
    if(command == fAngle)
    {
        fPrimaGen -> SetParticleAngle(fAngle -> GetNewDoubleValue(newValue));
    }
    if(command == fDistance)
    {
        fPrimaGen -> SetParticleDistance(fDistance -> GetNewDoubleValue(newValue));
    }
    if(command == fStartEnergy)
    {
        fPrimaGen -> SetStartEnergy(fStartEnergy -> GetNewDoubleValue(newValue));
    }
    if(command == fFinalEnergy)
    {
        fPrimaGen -> SetFinalEnergy(fFinalEnergy -> GetNewDoubleValue(newValue));
    }
    if(command == fEnergy)
    {
        fPrimaGen -> SetParticleEnergy(fEnergy -> GetNewDoubleValue(newValue));
    }
    if(command == fZPosition)
    {
        fPrimaGen -> SetBeamPosition(fZPosition -> GetNewDoubleValue(newValue));
    }
    else
    {
        G4cerr << "Invalid photon action command given. Aborting..." << G4endl;
        exit(0);
    }
}
