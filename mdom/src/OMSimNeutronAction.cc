#include "OMSimNeutronAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Neutron.hh"

extern G4int gRunID;

OMSimNeutronAction::OMSimNeutronAction(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun),
        fParticleExist(true),
        fParticleNum(0),
        fIdx(0)
{

}
OMSimNeutronAction::~OMSimNeutronAction()
{

}
void OMSimNeutronAction::GeneratePrimaries(G4Event* anEvent)
{
    try
    {
        std::cout << "Neutron x: " << fNeutronData.at(X).at(fIdx) << std::endl
    << "Neutron y: " << fNeutronData.at(Y).at(fIdx) << std::endl
    << "Neutron z: " << fNeutronData.at(Z).at(fIdx) << std::endl;
    }
    catch(...)
    {
        std::cerr << "Vector out of range. " << std::endl
        << "OMSimNeutronAction::GeneratePrimaries()" << std::endl;
        exit(0);
    }

    fParticleGun -> SetParticlePosition(G4ThreeVector(fNeutronData.at(X).at(fIdx) * m, fNeutronData.at(Y).at(fIdx) * m, fNeutronData.at(Z).at(fIdx) * m));
    fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(fNeutronData.at(AX).at(fIdx), fNeutronData.at(AY).at(fIdx), fNeutronData.at(AZ).at(fIdx)));
    fParticleGun -> SetParticleEnergy(fNeutronData.at(ENERGY).at(fIdx) * MeV);
    fParticleGun -> SetParticleTime(fNeutronData.at(TIME).at(fIdx) * ms);
    fParticleGun -> SetParticleDefinition(G4Neutron::NeutronDefinition());
    fParticleGun -> GeneratePrimaryVertex(anEvent);

   /* std::cout << "Neutron x: " << fNeutronData.at(X).at(fIdx) /m << std::endl
    << "Neutron y: " << fNeutronData.at(Y).at(fIdx) / m << std::endl
    << "Neutron z: " << fNeutronData.at(Z).at(fIdx) / m << std::endl; */

    fIdx++;

    if(fIdx == fParticleNum)
    {
        std::cout << "Total " << fIdx << " neutrons are generated!" << std::endl;
        fParticleExist = false;
    }
}
void OMSimNeutronAction::LoadData()
{
    using namespace std;
    G4String filePath = "../../analysis/files/input_geant4/gamma/Neutron/neu_gamma_" + std::to_string(gRunID) + "_"; //will change soon
    G4double temp;
    G4String fileName;

    for(unsigned int i = 0; i < fNeutronData.size(); i++)
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
            fNeutronData.at(i).push_back(temp);
        }

        file.close();
    }
    fParticleNum = fNeutronData.at(0).size();
}
