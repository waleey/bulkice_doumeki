#include "OMSimElectronAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"

extern G4int gRunID;

OMSimElectronAction::OMSimElectronAction(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun),
        fParticleExist(true),
        fParticleNum(0),
        fIdx(0)
{

}
OMSimElectronAction::~OMSimElectronAction()
{

}
void OMSimElectronAction::GeneratePrimaries(G4Event* anEvent)
{
    try
    {
        std::cout << "Electron x: " << fElectronData.at(X).at(fIdx) << std::endl
    << "Electron y: " << fElectronData.at(Y).at(fIdx) << std::endl
    << "Electron z: " << fElectronData.at(Z).at(fIdx) << std::endl;
    }
    catch(...)
    {
        std::cerr << "Vector out of range. " << std::endl
        << "OMSimElectronAction::GeneratePrimaries()" << std::endl;
        exit(0);
    }

    fParticleGun -> SetParticlePosition(G4ThreeVector(fElectronData.at(X).at(fIdx) * m, fElectronData.at(Y).at(fIdx) * m, fElectronData.at(Z).at(fIdx) * m));
    fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(fElectronData.at(AX).at(fIdx), fElectronData.at(AY).at(fIdx), fElectronData.at(AZ).at(fIdx)));
    fParticleGun -> SetParticleEnergy(fElectronData.at(ENERGY).at(fIdx) * MeV);
    fParticleGun -> SetParticleTime(fElectronData.at(TIME).at(fIdx) * ms);
    fParticleGun -> SetParticleDefinition(G4Electron::ElectronDefinition());
    fParticleGun -> GeneratePrimaryVertex(anEvent);

  /*  std::cout << "Electron x: " << fElectronData.at(X).at(fIdx) /m << std::endl
    << "Electron y: " << fElectronData.at(Y).at(fIdx) / m << std::endl
    << "Electron z: " << fElectronData.at(Z).at(fIdx) / m << std::endl;*/

    fIdx++;
    if(fIdx == fParticleNum)
    {
        std::cout << "Total " << fIdx << " electrons are generated!" << std::endl;
        fParticleExist = false;
    }
}
void OMSimElectronAction::LoadData()
{
    using namespace std;
    G4String filePath = "../../analysis/files/input_geant4/gamma/Electron/ele_gamma_" + std::to_string(gRunID) + "_";
    G4double temp;
    G4String fileName;

    for(unsigned int i = 0; i < fElectronData.size(); i++)
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
            fElectronData.at(i).push_back(temp);
        }

        file.close();
    }

    fParticleNum = fElectronData.at(0).size();

}
