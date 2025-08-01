#include "OMSimPositronAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Positron.hh"
extern G4String inputFolder;

OMSimPositronAction::OMSimPositronAction(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun),
        fParticleExist(true),
        fParticleNum(0),
        fIdx(0)
{

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
    try
    {
        fParticleGun -> SetParticlePosition(G4ThreeVector(fPositronData.at(X).at(fIdx) * m, fPositronData.at(Y).at(fIdx) * m, fPositronData.at(Z).at(fIdx) * m));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(fPositronData.at(AX).at(fIdx), fPositronData.at(AY).at(fIdx), fPositronData.at(AZ).at(fIdx)));
        fParticleGun -> SetParticleEnergy(fPositronData.at(ENERGY).at(fIdx) * MeV);
        fParticleGun -> SetParticleTime(fPositronData.at(TIME).at(fIdx) * ms);
        fParticleGun -> SetParticleDefinition(G4Positron::PositronDefinition());
        fParticleGun -> GeneratePrimaryVertex(anEvent);
    }
    catch(...)
    {
        std::cerr << "Vector out of range. " << std::endl
        << "OMSimPositronAction::GeneratePrimaries()" << std::endl;
        exit(0);
    }
    fIdx++;
    if(fIdx == fParticleNum)
    {
        std::cout << "Total " << fIdx << " poistrons are generated!" << std::endl;
        fParticleExist = false;
    }
}
void OMSimPositronAction::LoadData()
{
    using namespace std;
    G4String filePath = "../"+inputFolder+"/Positron/pos20002nkibd_"; //will change soon.
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

}
