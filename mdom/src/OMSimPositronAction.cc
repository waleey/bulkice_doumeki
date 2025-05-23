#include "OMSimPositronAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Positron.hh"
#include <cmath>

extern G4double gworldsize;
extern G4bool gVerbose;
extern G4int gRunID;
extern G4double gPositronDensity;

OMSimPositronAction::OMSimPositronAction(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun),
        fParticleExist(true),
        fParticleNum(0),
        fIdx(0),
        fPerc(0)
{
    fRadData = new OMSimRadioactivityData();
}
OMSimPositronAction::~OMSimPositronAction()
{

}
void OMSimPositronAction::GeneratePrimaries(G4Event* anEvent)
{
    
    try
    {
        std::cout << "Positron x: " << fPositronData.at(X).at(fIdx) << std::endl
    << "Positron y: " << fPositronData.at(Y).at(fIdx) << std::endl
    << "Positron z: " << fPositronData.at(Z).at(fIdx) << std::endl;
    }
    catch(...)
    {
        std::cerr << "Vector out of range. " << std::endl
        << "OMSimPositronAction::GeneratePrimaries()" << std::endl;
        exit(0);
    }

    fParticleGun -> SetParticlePosition(G4ThreeVector(fPositronData.at(X).at(fIdx) * m, fPositronData.at(Y).at(fIdx) * m, fPositronData.at(Z).at(fIdx) * m));
    fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(fPositronData.at(AX).at(fIdx), fPositronData.at(AY).at(fIdx), fPositronData.at(AZ).at(fIdx)));
    fParticleGun -> SetParticleEnergy(fPositronData.at(ENERGY).at(fIdx) * MeV);
    fParticleGun -> SetParticleTime(fPositronData.at(TIME).at(fIdx) * ms);
    fParticleGun -> SetParticleDefinition(G4Positron::PositronDefinition());
    fParticleGun -> GeneratePrimaryVertex(anEvent); 
    
    // prints progress for every percentage point
    if (fIdx > (fParticleNum * fPerc/100)){
        std::cerr << fPerc << " percent done" << std::endl;
        fPerc++;
    }

    fIdx++;
    if(fIdx == fParticleNum)
    {
        std::cout << "Total " << fIdx << " positrons are generated!" << std::endl;
        fParticleExist = false;
    }
}
void OMSimPositronAction::LoadData()
{
    using namespace std;
    G4String filePath = "../../analysis/files/input_geant4/gamma/Positron/pos_gamma_" + std::to_string(gRunID) + "_"; //will change soon
    G4double temp;
    G4String fileName;

    for(unsigned int i = 0; i < fPositronData.size(); i++)
    {
        temp = 0;
        fileName = filePath + dtypes.at(i) + ".data";

        ifstream file(fileName);

        if(!file.is_open())
        {
            std::cout << "Failed to open" + fileName << std::endl;

        }

        while(file >> temp)
        {
            fPositronData.at(i).push_back(temp);
        }

        file.close();
    }
    fParticleNum = fPositronData.at(0).size();
    std::cout << "Number of positrons: " << fParticleNum << std::endl;
}
