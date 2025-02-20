#include "OMSimPositronAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Positron.hh"

extern G4double gworldsize;
extern G4bool gVerbose;
OMSimPositronAction::OMSimPositronAction(G4ParticleGun* particleGun)
    :   fParticleGun(particleGun),
        fParticleExist(true),
        fParticleNum(0),
        fIdx(0)
{
    fRadData = new OMSimRadioactivityData();
}
OMSimPositronAction::~OMSimPositronAction()
{

}
void OMSimPositronAction::GeneratePrimaries(G4Event* anEvent)
{
    /*
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
    */
    
    G4double positionRandom[3] = {fRadData -> RandomGen(-gworldsize, gworldsize),fRadData -> RandomGen(-gworldsize, gworldsize),fRadData -> RandomGen(-gworldsize, gworldsize)};
    G4ThreeVector orientationRandom = fRadData -> SetupOrientation();
    G4double timeRandom = fRadData -> RandomGen(0, 1);

    fParticleGun -> SetParticlePosition(G4ThreeVector(positionRandom[0] * m, positionRandom[1] * m, positionRandom[2] * m));
    fParticleGun -> SetParticleMomentumDirection(orientationRandom);
    fParticleGun -> SetParticleEnergy(10 * MeV);
    fParticleGun -> SetParticleTime(timeRandom * s);
    fParticleGun -> SetParticleDefinition(G4Positron::PositronDefinition());
    fParticleGun -> GeneratePrimaryVertex(anEvent);

    if (fIdx % 1000 == 0){
        std::cerr << (fIdx / 1000) << " percent done" << std::endl;
    }

    try
    {
        if (true)//(gVerbose)
        {
            std::cout << "Positron x: " << positionRandom[0] << std::endl
                    << "Positron y: " << positionRandom[1] << std::endl
                    << "Positron z: " << positionRandom[2] << std::endl;
        }
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
        std::cout << "Total " << fIdx << " positrons are generated!" << std::endl;
        fParticleExist = false;
    }
}
void OMSimPositronAction::LoadData()
{
    using namespace std;
    G4String filePath = "../InputFile/Positron/pos20002nkibd_"; //will change soon.
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
    fParticleNum = 12500;
}
