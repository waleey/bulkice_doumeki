#include "OMSimRunManager.hh"
#include "OMSimAnalysisManager.hh"
#include "OMSimUtils.hh"
#include "OMSimPMTQE.hh"


#define G4VIS_USE 1

extern OMSimAnalysisManager gAnalysisManager;
extern G4String	ghitsfilename;
extern std::ofstream gBinaryHitFile;
//extern std::fstream gRadioDecayFile;
//extern std::fstream gPositronData;
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
extern G4String gHittype;
extern G4double gRadioactiveLeakContainment;
extern G4double gSimulationTime;
extern G4bool gVerbose;
G4double gAngle = 0;

//G4int gIdx = 0;
G4double OMSimRadioactivityData::ftimeWindow = gSimulationTime / s;
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
    OMSimPMTQE* qe = new OMSimPMTQE(); // create QE table instance
    qe -> ReadQeTable(); // load QE table once
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
    fSteppingAction = new OMSimSteppingAction(qe);
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
    if (fInteraction == "neu")
    {
        GenerateNeutrino();
    }
    else if(fInteraction == "ibd")
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
                //GenerateK40();
                //GenerateU238();
                //GenerateU235();
                //GenerateTh232();
                
                //GenerateDecayChain("K40");
                //GenerateDecayChain("U238");
                //GenerateDecayChain("U235");
                GenerateDecayChain("Th232");
                break;
            case lom18:
                GenerateK40();
                GenerateU238();
                GenerateU235();
                GenerateTh232();
                break;
            case lom16:
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
            case wom:
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

void OMSimRunManager::GenerateNeutrino()
{
    fPrimaryGenerator -> SetActionType(Neutrino);
    fPrimaryGenerator -> LoadData();
    while(fPrimaryGenerator -> ParticleExist())
    {
        gEvent++;
        fRunManager -> BeamOn(1);
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
    G4double timeWindow = OMSimRadioactivityData::ftimeWindow;
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
    G4double timeWindow = OMSimRadioactivityData::ftimeWindow;
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
    G4double timeWindow = OMSimRadioactivityData::ftimeWindow;
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
    G4double timeWindow = OMSimRadioactivityData::ftimeWindow;
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

void OMSimRunManager::GenerateDecayChain(G4String chainName)
{
    fPrimaryGenerator -> SetActionType(DecayChain);
    std::string filename = "../InputFile/";
    G4double activity;

    // Select decay chain
    if (chainName == "U238"){
        filename += "U238DecayChain.txt";
        activity = U238Activity.at(fpmtModel) * glassWeight.at(fpmtModel) * 1/s;
    }
    else if (chainName == "U235"){
        filename += "U235DecayChain.txt";
        activity = U235Activity.at(fpmtModel) * glassWeight.at(fpmtModel) * 1/s;
    }
    else if (chainName == "Th232"){
        filename += "Th232DecayChain.txt";
        activity = Th232Activity.at(fpmtModel) * glassWeight.at(fpmtModel) * 1/s;
    }
    else if (chainName == "K40"){
        filename += "K40DecayChain.txt";
        activity = K40Activity.at(fpmtModel) * glassWeight.at(fpmtModel) * 1/s;
    }
    else
    {
        std::cerr << "Invalid decay chain name. Aborting..." << std::endl;
        exit(0);
    }
    
    // Time window
    G4double timeWindow = gSimulationTime;

    G4double prob1 = 1-gRadioactiveLeakContainment; // probability for isotope to decay in the next timeWindow duration
    G4double prob2 = gRadioactiveLeakContainment; // probability that isotope decayed in leak time window
    G4double lifeTimeThresh = - timeWindow / std::log(1-prob1);
    G4double leakTimeThresh = timeWindow * std::log(1-prob2)/std::log(1-prob1);
    G4double leakTimeWindow;

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "Simulating Full " << chainName << " Decay Chain" << std::endl;
    std::cout << chainName << " Activity [Bq/Hz]: " << activity * s << std::endl;
    std::cout << "Time window [s]: " << timeWindow / s<< std::endl;
    std::cout << "p1 [%]: " << prob1 * 100 << ", p2 [%]: " << prob2 * 100 << std::endl;    
    std::cout << "Threshold Mean Life [s]: " << lifeTimeThresh / s << std::endl;
    std::cout << "Threshold Leak Time [s]: " << leakTimeThresh / s << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

    G4ThreeVector point;
    std::vector<std::string> skipIsotopes; // Set to store isotope names

    // load decay chain configuration file
    std::vector<ElementData> decayChain = loadDataFromFile(filename);
    // loop over the decay chain
    for (size_t i = 0; i < decayChain.size(); i++)
    {
        const auto& [isotopeName, Z, A, excitationEnergy, totalAngularMomentum, branchingRatio] = decayChain[i]; 
        
        if (std::find(skipIsotopes.begin(), skipIsotopes.end(), isotopeName) != skipIsotopes.end())
        {
            std::cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n|" 
                  << "   0 decays " 
                  << isotopeName 
                  << "\n| --> Skipping because " << isotopeName << " is included in skip list."
                  << "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" 
                  << std::endl;
            continue;
        }        

        G4double meanRate;
        G4double timeLow;
        G4double timeHigh;

        G4double maxLifeTime = 0;
        std::vector<std::string> skipGroup; // sub-group of isotope names to avoid
        size_t j = i + 1;

        while (j < decayChain.size())
        {
            const auto& [nextIsotopeName, nextZ, nextA, nextExcitationEnergy, nextTotalAngularMomentum, nextBranchingRatio] = decayChain[j]; 
            

            fPrimaryGenerator -> SetZ(nextZ);
            fPrimaryGenerator -> SetA(nextA);
            fPrimaryGenerator -> SetExcitationEnergy(nextExcitationEnergy);
            fPrimaryGenerator -> SetTotalAngularMomentum(nextTotalAngularMomentum);
            fPrimaryGenerator -> GenerateIsotope();
            G4double nextLifeTime = fPrimaryGenerator -> GetPDGLifeTime(); // safe original PDG life time
            if (nextIsotopeName == "Pa234m") // hard coded life time for Pa234m, can be loaded from data?
            {
                nextLifeTime = 100.325 * s;
                //std::cout << "Pa234m lifetime " << nextLifeTime / ns << " ns!" << std::endl;
            }
            if (nextLifeTime < lifeTimeThresh)
            {
                skipGroup.push_back(nextIsotopeName);
                maxLifeTime = std::max(maxLifeTime, nextLifeTime);
                j++;
            }
            else
            {
                //skipGroup.push_back(nextIsotopeName); // last daughter to fulfill criteria has semi-stable daugher lifetime
                break; // stop if we reach a "semi-stable" isotope
            }
        }
        
        if (skipGroup.size() != 0)
        {
            leakTimeWindow = -maxLifeTime * std::log(1-prob2); // duration in for which 99% of daughters decay, take max of leak chain
            meanRate = activity * (timeWindow+leakTimeWindow); // needs to move inside the isotope loop
            timeLow = -leakTimeWindow; // lower time window
            timeHigh = timeWindow; // higher time window
            skipIsotopes.insert(skipIsotopes.end(), skipGroup.begin(), skipGroup.end()); // when leak time is simulated, avoid injecting daughters to prevent double counting
            
            // print all isotopes that will be skipped
            std::cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
                        << "| Simulate leakage in time window [" << timeLow / s << ", " << timeHigh / s << "] s."
                        << std::endl;
            for (const auto& iso : skipGroup)
            {
                std::cout << "| Adding " << iso << " to skip list." << std::endl;
            }
            std::cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
                        << std::endl;
        }
        else
        {
            meanRate = activity * timeWindow; // needs to move inside the isotope loop
            timeLow = 0; // lower time window
            timeHigh = timeWindow; // higher time window
        }
        
        // configure primary generator/isotope
        fPrimaryGenerator -> SetZ(Z);
        fPrimaryGenerator -> SetA(A);
        fPrimaryGenerator -> SetExcitationEnergy(excitationEnergy);
        fPrimaryGenerator -> SetTotalAngularMomentum(totalAngularMomentum);
        fPrimaryGenerator -> GenerateIsotope();
        fPrimaryGenerator -> SetTimeLow(timeLow); // set ftimeLow in fradData, used inside OMSimDecayChainAction class
        fPrimaryGenerator -> SetTimeHigh(timeHigh); // set ftimeHigh in fradData

        G4double const originalLifeTime = fPrimaryGenerator -> GetPDGLifeTime(); // safe original PDG life time
        //fPrimaryGenerator -> SetPDGLifeTime(0 * ns);

        // number of decays pulled from Poisson distribution
        G4int numDecays = CLHEP::RandPoisson::shoot(meanRate);
        if (branchingRatio < 1.0) //use rejection sampling for branching ratios != 1
        {
            G4int counter = 0;
            for (G4int i = 0; i < numDecays; ++i)
            {
                G4double rand = CLHEP::RandFlat::shoot(0.0,1.0);
                if (branchingRatio > rand)
                {
                    ++counter;
                }
            }
            numDecays = counter;
        }
        std::cout << "----------------------------------------------------------------------\n|   " 
                  << numDecays 
                  << " decays ("
                  << meanRate
                  << " exp. decays) " 
                  << isotopeName 
                  << " (Z=" 
                  << Z 
                  << ", A=" 
                  << A 
                  << ", E=" 
                  << excitationEnergy
                  << " keV, J=" 
                  << totalAngularMomentum * 1/2 
                  << ", tau="
                  << originalLifeTime / ns
                  << " ns)\n----------------------------------------------------------------------" 
                  << std::endl;
        for (G4int i = 0; i < numDecays; ++i)
        {
            // Generate a random point for the decay
            point = fDetectorConstruction -> DrawFromVolume();
            // Set particle position    
            fPrimaryGenerator -> SetPosition(point);
            // Process one event with RunManager
            fRunManager -> BeamOn(1);
            // restore original PDG life time, inside loop, because life-time is needed to compute initial decay time
            fPrimaryGenerator -> SetPDGLifeTime(originalLifeTime);
        }
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

    if (gHittype != "binary")
    {   G4String filename = ghitsfilename + ".dat";
        std::cout << "filename :" << filename << std::endl;
        gAnalysisManager.datafile.open(filename, std::ios::out /*| std::ios::app*/);
        gNumCherenkov = 0;
        gNumScint = 0;
    }
    else
    {
        G4String binaryName = ghitsfilename + ".bin";
        std::cout << "filename :" << binaryName << std::endl;
        gBinaryHitFile.open(binaryName, std::ios::binary);
    }
    //gRadioDecayFile.open("radiodecay_output.csv", std::ios::out /*| std::ios::app*/);
    //gRadioDecayFile << "TrackID,ParticleName,Energy\n";

    //gPositronData.open("positron_data.csv", std::ios::out);
    //gPositronData << "E_nu,E_p0,E_p1,th_uni,th_0,th_1,phi,posx,posy,poz\n";
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
    
    gBinaryHitFile.close();

    //gRadioDecayFile.close();
    
    //gPositronData.close();
}
