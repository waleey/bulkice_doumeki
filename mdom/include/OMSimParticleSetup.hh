#ifndef OMSIMPARTICLESETUP_HH_
#define OMSIMPARTICLESETUP_HH_

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"

#include "Poisson.hh"


#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <cmath>

class OMSimParticleSetup
{
    public:
    OMSimParticleSetup(G4ParticleGun* ParticleGun, G4Event* anEvent, G4int omModel);
    ~OMSimParticleSetup() {}
    void GeneratePositron();
    void GenerateNeutron();
    void GenerateElectron();

    /**
    *The decay chains are simulated independently
    *the photon produced due to Cerenkov effects and scintillation are mixed
    *Each daughter produced from Radioactive process will have a global time of 0
    *This is equivalent to making the lifetime of the radioactive particles 0
    **/
    void GenerateK40();
    void GenerateTh238();
    void GenerateU238();
    void GenerateU235();
    void SetGlassWeight(G4double glassWeight) { fglassweight = glassWeight; }

    static G4double glassInRad;
    static G4double glassOutRad;
    static void SetGlassRad(G4double, G4double);

    /**
    **/

    private:
    void SetupPositron();
    void SetupNeutron();
    void SetupElectron();

    void SetupRadioactiveParticles();


    G4int fomModel;
    /*OMSimInputData* mData;
    G4String mDataKey;*/
    G4int pos_count;
    G4int neu_count;
    G4int e_count;
    G4ParticleGun* fParticleGun;
    G4Event* fEvent;
    G4ParticleTable* fParticleTable = G4ParticleTable::GetParticleTable();
    void DataReader(std::string filepath, std::vector<std::vector<G4double>>& data);
    double RandomGen(double minLim, double maxLim);

    /**
    *Setting up dummy parameters for radioactive particles
    *for now, U235, U238, Th238, and K40 decays are simulated
    *it is assumed that they are emitted from a spherical surface
    *no information about their decay time, activity, and quantity are accruately known to me

    the activity is set to a dummy value for now.
    Will be updates later following M.Unland's thesis.
    Also, the timewindo is dummy for now.
    **/
    int radioactiveParticleNum = 500; //just a dummy number. More realistic approach will be used later.
    double radioactiveParticleEnergy;
    double radioactiveParticleCharge;
    G4ThreeVector radioactiveParticlePosition;
    G4ThreeVector radioactiveParticleOrientation;
    G4ParticleDefinition* radioactiveParticle;
    G4double k40activity;
    G4double u238activity;
    G4double u235activity;
    G4double th238activity;
    G4double timeWindow;
    G4double fglassweight;
    Poisson* poisson;
    /**
    **/
    std::vector<G4double> energy;
    std::vector<G4double> fX;
    std::vector<G4double> fY;
    std::vector<G4double> fZ;
    std::vector<G4double> alpha_X;
    std::vector<G4double> alpha_Y;
    std::vector<G4double> alpha_Z;
    std::vector<G4double> inTime;
    std::vector<std::vector<G4double>> pos_data {energy, fX, fY, fZ, alpha_X, alpha_Y, alpha_Z, inTime};
    std::vector<std::vector<G4double>> neu_data {energy, fX, fY, fZ, alpha_X, alpha_Y, alpha_Z, inTime};
    std::vector<std::vector<G4double>> e_data{energy, fX, fY, fZ, alpha_X, alpha_Y, alpha_Z, inTime};
    std::string pos_filePath = "/home/waly/bulkice_doumeki/mdom/InputFile/Positron/pos20002nkibd_"; //CHANGE THE FILE PATH
    std::string neu_filePath = "/home/waly/bulkice_doumeki/mdom/InputFile/Neutron/neu20002nkibd_"; //CHANGE THE FILE PATH
    std::string e_filePath = "/home/waly/bulkice_doumeki/mdom/InputFile/Electron/e20002nkibd_"; //CHANGE THE FILE PATH
    std::vector<std::string>  dtypes {"energy", "x", "y", "z", "ax", "ay", "az", "time"};
    enum {ENERGY, X, Y, Z, AX, AY, AZ, TIME};
};

#endif // OMSIMPARTICLESETUP_HH_
