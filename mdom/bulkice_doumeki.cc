#include <iostream>
#include <sstream>
#include <string>
#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4Cerenkov.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4DecayPhysics.hh"

#define G4VIS_USE 1
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4RayTracer.hh"
#endif

#include "OMSimDetectorConstruction.hh"
#include "OMSimPhysicsList.hh"
#include "OMSimPrimaryGeneratorAction.hh"
#include "OMSimRunAction.hh"
#include "OMSimEventAction.hh"
#include "OMSimTrackingAction.hh"
#include "OMSimSteppingAction.hh"
#include "OMSimSteppingVerbose.hh"
#include "OMSimAnalysisManager.hh"
//#include "OMSimPMTQE.hh"

//setting up the external variables
G4int           gGlass = 1;
G4int           gGel = 1;
G4double        gRefCone_angle = 51;
G4int           gConeMat = 1;
G4int           gHolderColor = 1;
//G4int           gDOM = 0; // 0 : single PMT
G4int           gDOM = 1; // 1 : mdom
//G4int           gDOM = 2; // 2 : pdom
G4int           gPMT = 3; // i don't know what is it
G4bool          gPlaceHarness = true;
G4int           gHarness = 1;
G4int           gRopeNumber = 1;
G4double        gworldsize = 0.5*m;

G4bool          gCADImport = true;
G4String        gHittype = "individual"; // seems like individual records each hit per pmt
G4bool          gVisual = true; // may be visualization on?
G4int           gEnvironment = 2; // 1 is ICeCUbe ice without any ice property, 2 is with property, 0 is air.
G4String        ghitsfilename = "/mnt/c/Users/Waly/bulkice_doumeki/hit_";
//G4String        ghitsfilename = "hit.dat";
G4int           gcounter = 0;
G4int           gPosCount = 0;
G4String        gQEFile = "/home/waly/bulkice_doumeki/mdom/InputFile/TA0001_HamamatsuQE.data";
enum {PMT, MDOM, DOM, LOM16, LOM18};

OMSimAnalysisManager gAnalysisManager;

void ParseCommandLine(int argc, char** argv, G4String& macroname, G4int& PMT_model, G4double& worldsize, G4String& interaction_channel)
{
    if(argc == 4 || argc == 2)
    {
        std::string model = argv[1];

        if(model == "dom")
        {
            PMT_model = 2;
        }
        else if(model == "mdom")
        {
            PMT_model = 1;
        }
        else if(model == "lom16")
        {
            PMT_model = 3;
        }
        else if(model == "lom18")
        {
            PMT_model = 4;
        }
        else if(model == "pmt")
        {
            PMT_model = 0;
        }
        else
        {
            std::cerr << "ERROR! Invalid OM Model. Aborting..." << std::endl;
            std::cout << "Usage: " << std::endl;
            std::cout << "Available OM Models: [dom, mdom, lom16, lom18, pmt]" << std::endl;
            exit(0);
        }

        if(argc == 4)
        {
            macroname = argv[3];
            interaction_channel = argv[2];
            worldsize = 20;
            ghitsfilename += model + "_" + interaction_channel + ".dat";
        }
        else
        {
            worldsize = 0.25;
            ghitsfilename += model + ".dat";
        }

    }
    else
    {
        std::cerr << "Invalid number of argument given. Aborting...." << std::endl;
        std::cout << "Usage: " << std::endl;
        std::cout << "./bulkice_doumeki" << " " << "[OM Model]" << " " << "[interaction channel]" << " " <<  "[macro name] (for batch mode operation)" << std::endl;
        std::cout << "./bulkice_doumeki" << " " << "[OM Model]" << " (for visualization)" << std::endl;
        std::cout << "Available OM Models: [dom, mdom, lom16, lom18, pmt]" << std::endl;
        std::cout << "Available interaction channels: [ibd, enees, all, radioactivity]" << std::endl;
        exit(0);
    }
}

std::vector<G4String> explode(G4String s, char d) {
        std::vector<G4String> o;
        int i,j;
        i = s.find_first_of("#");
        if (i == 0) return o;
        while (s.size() > 0) {
                i = s.find_first_of(d);
                j = s.find_last_of(d);
                o.push_back(s.substr(0, i));
                if (i == j) {
                        o.push_back(s.substr(j+1));
                        break;
                }
                s.erase(0,i+1);
        }
        return o;
    }

std::vector<G4String> explode(char* cs, char d) {
        std::vector<G4String> o;
        G4String s = cs;
        return explode(s,d);
}

int main(int argc, char** argv)
{



    G4double world_size(0);
    G4String macroname;
    G4int PMT_model(0);
    G4String interaction_channel;
    ParseCommandLine(argc, argv, macroname, PMT_model, world_size, interaction_channel);


    G4RunManager* runmanager = new G4RunManager();

    //Generating Physics List
    G4VModularPhysicsList* physicsList = new FTFP_BERT;
    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
    /***
    Will add multithreading and update the optical physics code later!!!
    **/
    physicsList->RegisterPhysics(opticalPhysics);
    physicsList -> RegisterPhysics(new G4RadioactiveDecayPhysics);
    physicsList -> RegisterPhysics(new G4DecayPhysics);

    runmanager -> SetUserInitialization(new OMSimDetectorConstruction(PMT_model, world_size));
    runmanager -> SetUserInitialization(physicsList);
    //runmanager -> SetUserInitialization(new OMSimPhysicsList);
    runmanager -> SetUserAction(new OMSimPrimaryGeneratorAction(interaction_channel, PMT_model));
    runmanager -> SetUserAction(new OMSimEventAction);
    runmanager -> SetUserAction(new OMSimRunAction);
    runmanager -> SetUserAction(new OMSimSteppingAction);
    runmanager -> SetUserAction(new OMSimTrackingAction);

    #ifdef G4VIS_USE
  // initialize visualization package
    G4VisManager* vismanager = new G4VisExecutive();
    vismanager -> RegisterGraphicsSystem(new G4RayTracer);
  //G4VisManager* visManager= new J4VisManager;
    vismanager -> Initialize();
    std::cerr << " ------------------------------- " << std::endl
         << " ---- VisManager created! ---- " << std::endl
         << " ------------------------------- " << std::endl;
    std::cerr << std::endl;
    #endif

    std::cerr << "about to initialize runManager" << std::endl;
    runmanager-> Initialize();
    std::cerr << "initialize runManager succeed" << std::endl;

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    G4UIExecutive* ui = 0;

    if ( argc == 4 ) {
    // batch mode
        std::cerr << ":::::::::::::::::::Batch Mode Called:::::::::::::::::" << std::endl;
        G4String command = "/control/execute " + macroname;
        std::cout << "command " << command << std::endl;
        UImanager->ApplyCommand(command);
        }
    else {
    // interactive mode
        std::cerr << "interactive mode called" << std::endl;
        ui = new G4UIExecutive(argc, argv);
        UImanager->ApplyCommand("/control/execute ray_trace.mac");
        ui-> SessionStart();
        delete ui;
    }

  //-----------------------
  // terminating...
  //-----------------------

    #ifdef G4VIS_USE
    delete vismanager;
    #endif

    delete runmanager;
    std::cout << "::::::::::::::this is the end:::::::::::::"<< std::endl;
    return 0;

}
