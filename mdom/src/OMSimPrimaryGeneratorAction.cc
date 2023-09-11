#include "OMSimPrimaryGeneratorAction.hh"
#include "OMSimDetectorConstruction.hh"
#include "OMSimRadioactivityData.hh"

//#include "G4ParticleTypes.hh"

OMSimPrimaryGeneratorAction::OMSimPrimaryGeneratorAction()
    :   fPositronAction(0),
        fNeutronAction(0),
        fElectronAction(0),
        fK40Action(0),
        fU238Action(0),
        fU235Action(0),
        fTh232Action(0),
        fPhotonAction(0),
        fParticleGun(0),
        fActionType(0),
        fPhotonAngle(0)
{

	fParticleGun = new G4ParticleGun(1);
	G4cout << ":::::::::::::::Particle Gun Created:::::::::::" << G4endl;

	fPositronAction = new OMSimPositronAction(fParticleGun);
	fNeutronAction = new OMSimNeutronAction(fParticleGun);
	fElectronAction = new OMSimElectronAction(fParticleGun);
	fK40Action = new OMSimK40Action(fParticleGun);
	fU238Action = new OMSimU238Action(fParticleGun);
	fU235Action = new OMSimU235Action(fParticleGun);
	fTh232Action = new OMSimTh232Action(fParticleGun);
	fPhotonAction = new OMSimPhotonAction(fParticleGun);
}
OMSimPrimaryGeneratorAction::OMSimPrimaryGeneratorAction(G4String& interaction)
    :   fPositronAction(0),
        fNeutronAction(0),
        fElectronAction(0),
        fK40Action(0),
        fU238Action(0),
        fU235Action(0),
        fTh232Action(0),
        fPhotonAction(0),
        fParticleGun(0),
        fActionType(0),
        fInteraction(interaction)
{

	fParticleGun = new G4ParticleGun(1);
	G4cout << ":::::::::::::::Particle Gun Created:::::::::::" << G4endl;
	fPositronAction = new OMSimPositronAction(fParticleGun);
	fNeutronAction = new OMSimNeutronAction(fParticleGun);
	fElectronAction = new OMSimElectronAction(fParticleGun);
	fK40Action = new OMSimK40Action(fParticleGun);
	fU238Action = new OMSimU238Action(fParticleGun);
	fU235Action = new OMSimU235Action(fParticleGun);
	fTh232Action = new OMSimTh232Action(fParticleGun);
	fPhotonAction = new OMSimPhotonAction(fParticleGun);
}

OMSimPrimaryGeneratorAction::~OMSimPrimaryGeneratorAction()
{
	delete fPositronAction;
	delete fNeutronAction;
	delete fElectronAction;
	delete fK40Action;
	delete fU238Action;
	delete fU235Action;
	delete fTh232Action;
	delete fPhotonAction;
	delete fParticleGun;
}

void OMSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	/*if(fInteraction == "vis")
	{
        fActionType = Visualization;
	}*/
	switch(fActionType)
	{
        case Positron:
            std::cout << "Generating Positrons!" << std::endl;
            fPositronAction -> GeneratePrimaries(anEvent);
            break;
        case Neutron:
            std::cout << "Generating Neutrons!" << std::endl;
            fNeutronAction -> GeneratePrimaries(anEvent);
            break;
        case Electron:
            std::cout << "Generating Electrons!" << std::endl;
            fElectronAction -> GeneratePrimaries(anEvent);
            break;
        case K40:
            std::cout << "Generating K40!" << std::endl;
            fK40Action -> GeneratePrimaries(anEvent);
            break;
        case U238:
            std::cout << "Generating U238!" << std::endl;
            fU238Action -> GeneratePrimaries(anEvent);
            break;
        case U235:
            std::cout << "Generating U235!" << std::endl;
            fU235Action -> GeneratePrimaries(anEvent);
            break;
        case Th232:
            std::cout << "Generating Th232!" << std::endl;
            fTh232Action -> GeneratePrimaries(anEvent);
            break;
       /* case Visualization:
            //std::cout << "Visualization called!" << std::endl;
            GenerateToVisualize();
            fParticleGun -> GeneratePrimaryVertex(anEvent);
            break;*/
        case Photon:
            //std::cout << "Generating Photon wave!" << std::endl;
            fPhotonAction -> GeneratePrimaries(anEvent);
            break;
        default:
            std::cerr << "Invalid action type in OMSimPrimaryGeneratorAction::GeneratePrimaries() " << std::endl
            << "Aborting...." << std::endl;
            exit(0);
	}
}
bool OMSimPrimaryGeneratorAction::ParticleExist()
{
    switch(fActionType)
    {
        case Positron:
            return fPositronAction -> PositronExist();
        case Neutron:
            return fNeutronAction -> NeutronExist();
        case Electron:
            return fElectronAction -> ElectronExist();
        default:
            std::cerr << "Invalid action type in OMSimPrimaryGeneratorAction::ParticleExist() " << std::endl
            << "Aborting...." << std::endl;
            exit(0);
    }
}
void OMSimPrimaryGeneratorAction::LoadData()
{
    switch(fActionType)
    {
        case Positron:
            fPositronAction -> LoadData();
            break;
        case Neutron:
            fNeutronAction -> LoadData();
            break;
        case Electron:
            fElectronAction -> LoadData();
            break;
        default:
            std::cerr << "Invalid action type in OMSimPrimaryGeneratorAction::ParticleExist() " << std::endl
            << "Aborting...." << std::endl;
            exit(0);
    }
}
void OMSimPrimaryGeneratorAction::SetPosition(G4ThreeVector& position)
{
    switch(fActionType)
    {
        case K40:
            fK40Action -> SetPosition(position);
            break;
        case U238:
            fU238Action -> SetPosition(position);
            break;
        case U235:
            fU235Action -> SetPosition(position);
            break;
        case Th232:
            fTh232Action -> SetPosition(position);
            break;
        default:
            std::cerr << "Invalid action type in OMSimPrimaryGeneratorAction::SetPosition() " << std::endl
            << "Aborting...." << std::endl;
            exit(0);
    }
}
void OMSimPrimaryGeneratorAction::GenerateToVisualize()
{

    OMSimRadioactivityData* radData = new OMSimRadioactivityData();

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName = "opticalphoton";
    G4ParticleDefinition* particle = particleTable -> FindParticle(particleName);

    G4ThreeVector position(0., -100 * mm, 350.0 * mm);
    G4ThreeVector direction(radData -> SetupOrientation());

    fParticleGun -> SetParticleDefinition(particle);
    fParticleGun -> SetParticleEnergy(5 * eV);
    fParticleGun -> SetParticlePosition(position);
    fParticleGun -> SetParticleMomentumDirection(direction);

    delete radData;

}
void OMSimPrimaryGeneratorAction::SetAngle(G4double angle)
{
    fPhotonAction -> SetAngle(angle);
}
