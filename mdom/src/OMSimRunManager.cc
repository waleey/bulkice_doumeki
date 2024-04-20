#include "OMSimRunManager.hh"
#include "OMSimAnalysisManager.hh"

#define G4VIS_USE 1

extern OMSimAnalysisManager gAnalysisManager;
extern G4String	ghitsfilename;
extern G4String	gHittype;
extern G4int gNumCherenkov;
extern G4int gNumScint;
extern G4int gEvent;
extern G4double gDistance;
extern G4double gZenithAngle;
extern G4double gStartAngle;
extern G4double gFinalAngle;
extern G4double gAngleIncrement;
extern G4bool gMultipleAngle;
extern G4int gPhotonNotAbsorbed;
extern G4bool gVis;
G4double gAngle = 0;

//G4int gIdx = 0;
G4double OMSimRadioactivityData::ftimeWindow = 60.0; //for now just running for 1 sec.
OMSimRunManager::OMSimRunManager()
{

    fRadData = new OMSimRadioactivityData();

}

OMSimRunManager::OMSimRunManager(G4int pmtModel, G4double worldSize, G4String& interactionChannel)
    :   fpmtModel(pmtModel),
        fworldSize(worldSize),
        fInteraction(interactionChannel),
        fRunManager(0),
        fDetectorConstruction(0),
        fPrimaryGenerator(0),
        fVisManager(0),
        fPhysicsList(0),
        fRunAction(0),
        fEventAction(0),
        fTrackingAction(0),
        fSteppingAction(0),
        fRadData(0),
        fActionType(0)
{
    fRunManager = new G4RunManager();
    fDetectorConstruction = new OMSimDetectorConstruction(fpmtModel, fworldSize);
    fRunManager -> SetUserInitialization(fDetectorConstruction);
    fPhysicsList = new OMSimPhysicsList();
    fRunManager -> SetUserInitialization(fPhysicsList);
    if(fInteraction == "vis")
    {
        fPrimaryGenerator = new OMSimPrimaryGeneratorAction(fInteraction);
    }
    else
    {
        fPrimaryGenerator = new OMSimPrimaryGeneratorAction();
    }
    fRunManager -> SetUserAction(fPrimaryGenerator);
    fRunAction = new OMSimRunAction();
    fRunManager -> SetUserAction(fRunAction);
    fEventAction = new OMSimEventAction();
     fRunManager -> SetUserAction(fEventAction);
    fTrackingAction = new OMSimTrackingAction();
     fRunManager -> SetUserAction(fTrackingAction);
    fSteppingAction = new OMSimSteppingAction();
    fRunManager -> SetUserAction(fSteppingAction);
    fRadData = new OMSimRadioactivityData();
    #ifdef G4VIS_USE
    fVisManager = new G4VisExecutive();
    #endif

    fRadData = new OMSimRadioactivityData();
}

OMSimRunManager::~OMSimRunManager()
{
    std::cout << "::::::::::::::this is the end of a run:::::::::::::"<< std::endl;

    if(gVis)
    {
    	delete fVisManager;
    }

    delete fRunManager;
    delete fRadData;
}

void OMSimRunManager::Initialize()
{
    #ifdef G4VIS_USE
    fVisManager -> RegisterGraphicsSystem(new G4RayTracer);
    fVisManager -> Initialize();
    std::cerr << " ------------------------------- " << std::endl
         << " ---- VisManager created! ---- " << std::endl
         << " ------------------------------- " << std::endl;
    std::cerr << std::endl;
    #endif

    std::cerr << "about to initialize runManager" << std::endl;
    fRunManager -> Initialize();
    std::cerr << "initialize runManager succeed" << std::endl;

}

void OMSimRunManager::BeamOn()
{
    if(fInteraction == "ibd")
    {
        //GenerateNeutron();
        GeneratePositron();

    }
    else if(fInteraction == "enees")
    {
        GenerateElectron();
    }
    else if(fInteraction == "all")
    {
        GeneratePositron();
        //GenerateNeutron();
        GenerateElectron();
    }
    else if(fInteraction == "radioactivity")
    {
        switch(fpmtModel)
        {
            case mdom:
                GenerateK40();
                GenerateU238();
                GenerateU235();
                GenerateTh232();
                break;
            case degg:
                GenerateK40();
                GenerateU238();
                GenerateTh232();
                break;
            default:
                std::cerr << "Radioactivity is not defined for this OM model yet. Aborting...." << std::endl;
                exit(0);
        }
    }
    /*else if(fInteraction == "vis")
    {
        GenerateToVisualize();
    }*/
    else if(fInteraction == "opticalphoton")
    {
        /**
        optical photon wave simulation
        **/
        GeneratePhoton();
    }
    else
    {
        std::cerr << "Invalid interaction channel. Aborting..." << std::endl;
        exit(0);
    }
}

void OMSimRunManager::GeneratePositron()
{
    fPrimaryGenerator -> SetActionType(Positron);
    fPrimaryGenerator -> LoadData();
   while(fPrimaryGenerator -> ParticleExist())
    {
        gEvent++;
        fRunManager -> BeamOn(1);
    }
}

void OMSimRunManager::GenerateNeutron()
{
    fPrimaryGenerator -> SetActionType(Neutron);
    fPrimaryGenerator -> LoadData();
    while(fPrimaryGenerator -> ParticleExist())
    {
        gEvent++;
        fRunManager -> BeamOn(1);
    }
}

void OMSimRunManager::GenerateElectron()
{
    fPrimaryGenerator -> SetActionType(Electron);
    fPrimaryGenerator -> LoadData();
    while(fPrimaryGenerator -> ParticleExist())
    {
        gEvent++;
        fRunManager -> BeamOn(1);
    }
}

void OMSimRunManager::GenerateK40()
{
    fPrimaryGenerator -> SetActionType(K40); //This member function needs to be implemented
    G4double activity = K40Activity.at(fpmtModel);
    G4double timeWindow = OMSimRadioactivityData::ftimeWindow / 60;
    G4double meanRate = activity * timeWindow * glassWeight.at(fpmtModel); //multiplied by mdom vessel weight of 13. Will be made flexible later.
    G4int numParticle = fRadData -> GetNumDecay(meanRate);
    G4ThreeVector point;
    for(int i = 0; i < numParticle; i++)
    {
        gEvent++;
        point = fDetectorConstruction -> DrawFromVolume();
        fPrimaryGenerator -> SetPosition(point);
        fRunManager -> BeamOn(1); //each primary will contain a separate event.
    }

}

void OMSimRunManager::GenerateU238()
{
    fPrimaryGenerator -> SetActionType(U238); //This member function needs to be implemented
    G4double activity = U238Activity.at(fpmtModel);
    G4double timeWindow = OMSimRadioactivityData::ftimeWindow / 60;
    G4double meanRate = activity * timeWindow * glassWeight.at(fpmtModel); //multiplied by mdom vessel weight of 13. Will be made flexible later.
    G4int numParticle = fRadData -> GetNumDecay(meanRate);
    G4ThreeVector point;
    for(int i = 0; i < numParticle; i++)
    {
        gEvent++;
        point = fDetectorConstruction -> DrawFromVolume();
        fPrimaryGenerator -> SetPosition(point); //not implemented yet.
        fRunManager -> BeamOn(1); //each primary will contain a separate event.
    }

}

void OMSimRunManager::GenerateU235()
{
    fPrimaryGenerator -> SetActionType(U235); //This member function needs to be implemented
    G4double activity = U235Activity.at(fpmtModel);
    G4double timeWindow = OMSimRadioactivityData::ftimeWindow / 60;
    G4double meanRate = activity * timeWindow * glassWeight.at(fpmtModel); //multiplied by mdom vessel weight of 13. Will be made flexible later.
    G4int numParticle = fRadData -> GetNumDecay(meanRate);
    G4ThreeVector point;
    for(int i = 0; i < numParticle; i++)
    {
        gEvent++;
        point = fDetectorConstruction -> DrawFromVolume();
        fPrimaryGenerator -> SetPosition(point); //not implemented yet.
        fRunManager -> BeamOn(1); //each primary will contain a separate event.
    }
}

void OMSimRunManager::GenerateTh232()
{
    fPrimaryGenerator -> SetActionType(Th232); //This member function needs to be implemented
    G4double activity = Th232Activity.at(fpmtModel);
    G4double timeWindow = OMSimRadioactivityData::ftimeWindow / 60;
    G4double meanRate = activity * timeWindow * glassWeight.at(fpmtModel); //multiplied by mdom vessel weight of 13. Will be made flexible later.
    G4int numParticle = fRadData -> GetNumDecay(meanRate);
    G4ThreeVector point;
    for(int i = 0; i < numParticle; i++)
    {
        gEvent++;
        point = fDetectorConstruction -> DrawFromVolume();
        fPrimaryGenerator -> SetPosition(point); //not implemented yet.
        fRunManager -> BeamOn(1); //each primary will contain a separate event.
    }
}
/*void OMSimRunManager::GenerateToVisualize()
{
    fPrimaryGenerator -> SetActionType(Visualization);

}*/
void OMSimRunManager::GeneratePhoton()
{
    fPrimaryGenerator -> SetActionType(Photon);
    for(G4double angle = gStartAngle; angle <= gFinalAngle; angle += gAngleIncrement)
    {
        gAngle = angle;
        fPrimaryGenerator -> SetAngle(angle);
        fRunManager -> BeamOn(1);
    }
}
void OMSimRunManager::OpenFile()
{
    G4cout << ":::::::::This is the beginning of Run Action::::::::" << G4endl;
    startingtime = clock() / CLOCKS_PER_SEC;
    G4String filename = ghitsfilename + ".dat";
	gAnalysisManager.datafile.open(filename, std::ios::out /*| std::ios::app*/);
	gNumCherenkov = 0;
	gNumScint = 0;
}

void OMSimRunManager::CloseFile()
{
    G4cout << "::::::::::::This is the end of Run Action:::::::::::" << G4endl;

	std::cout << "Cherenkov Produced: " << gNumCherenkov << std::endl
	<< "Scintillation Photons Produced: " << gNumScint << std::endl
	<< "Photon Reached WLS but Not Absorbed: " << gPhotonNotAbsorbed << std::endl;


// 	Close output data file
    gAnalysisManager.datafile.close();
    gAnalysisManager.Reset();
    double finishtime = clock() / CLOCKS_PER_SEC;
//std::cout << "Total Positron Count in World Volume:: " << gPosCount << std::endl;
    std::cout << "Total Event Produced: " << gEvent << std::endl;
    G4cout << "Computation time: " << finishtime - startingtime << " seconds." << G4endl;
}
