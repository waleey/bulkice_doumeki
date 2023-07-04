#ifndef OMSimPrimaryGeneratorAction_h
#define OMSimPrimaryGeneratorAction_h 1

//#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "OMSimParticleSetup.hh"

#include <vector>
#include <string>
//#include "globals.hh"

/*class G4GeneralParticleSource;
class G4Event;*/


class OMSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	OMSimPrimaryGeneratorAction(G4String&, G4int);
	~OMSimPrimaryGeneratorAction();


public:
	void GeneratePrimaries(G4Event* anEvent);

private:
	//G4GeneralParticleSource* particleSource;

	//creating particle gun and make it read from sntools output files

	//void SetUpEnergyAndPosition();

	G4ParticleGun *fParticleGun;
	std::string finteraction_channel;
	G4int fomModel;
    /*G4int numParticles;

    std::vector<G4double> energy;
    std::vector<G4double> fX;
    std::vector<G4double> fY;
    std::vector<G4double> fZ;
    std::vector<G4double> alpha_X;
    std::vector<G4double> alpha_Y;
    std::vector<G4double> alpha_Z;
    std::vector<G4double> inTime;
    std::vector<std::vector<G4double>> data {energy, fX, fY, fZ, alpha_X, alpha_Y, alpha_Z, inTime};
    std::string filePath = "/home/waly/bulkice_doumeki/mdom/InputFile/20002nkibd_"; //CHANGE THE FILE PATH
    std::vector<std::string>  dtypes {"energy", "x", "y", "z", "ax", "ay", "az", "time"};
    enum {ENERGY, X, Y, Z, AX, AY, AZ, TIME};
    */
};



#endif
