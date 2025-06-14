#include "OMSimGammaAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"

OMSimGammaAction::OMSimGammaAction(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun),
        fParticleExist(true),
        fParticleNum(0),
        fIdx(0)
{

}
OMSimGammaAction::~OMSimGammaAction()
{

}
void OMSimGammaAction::GeneratePrimaries(G4Event* anEvent)
{
    try
    {
        std::cout << "Gamma x: " << fGammaData.at(X).at(fIdx) << std::endl
        << "Gamma y: " << fGammaData.at(Y).at(fIdx) << std::endl
        << "Gamma z: " << fGammaData.at(Z).at(fIdx) << std::endl;
    }
    catch(...)
    {
        std::cerr << "Vector out of range. " << std::endl
        << "OMSimGammaAction::GeneratePrimaries()" << std::endl;
        exit(0);
    }
    try
    {
        fParticleGun -> SetParticlePosition(G4ThreeVector(fGammaData.at(X).at(fIdx) * m, fGammaData.at(Y).at(fIdx) * m, fGammaData.at(Z).at(fIdx) * m));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(fGammaData.at(AX).at(fIdx), fGammaData.at(AY).at(fIdx), fGammaData.at(AZ).at(fIdx)));
        fParticleGun -> SetParticleEnergy(fGammaData.at(ENERGY).at(fIdx) * MeV);
        fParticleGun -> SetParticleTime(fGammaData.at(TIME).at(fIdx) * ms);
        fParticleGun -> SetParticleDefinition(G4Gamma::GammaDefinition());
        fParticleGun -> GeneratePrimaryVertex(anEvent);
    }
    catch(...)
    {
        std::cerr << "Vector out of range. " << std::endl
        << "OMSimGammaAction::GeneratePrimaries()" << std::endl;
        exit(0);
    }
    fIdx++;
    if(fIdx == fParticleNum)
    {
        std::cout << "Total " << fIdx << " gammas are generated!" << std::endl;
        fParticleExist = false;
    }
}
void OMSimGammaAction::LoadData()
{
    using namespace std;
    G4String filePath = "../InputFile/Gamma/gamma_"; //will change soon.
    G4double temp;
    G4String fileName;

    for(unsigned int i = 0; i < fGammaData.size(); i++)
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
            fGammaData.at(i).push_back(temp);
        }

        file.close();
    }
    fParticleNum = fGammaData.at(0).size();

}