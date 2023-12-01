#include <iostream>
#include <sstream>
#include <string>
#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4SystemOfUnits.hh"

#define G4VIS_USE 1
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4RayTracer.hh"
#endif

#include "OMSimAnalysisManager.hh"
#include "OMSimRunManager.hh"

#include <ctime>
#include <string>


//setting up the external variables
G4int           gGlass = 0;
G4int           gGel = 0;
G4double        gRefCone_angle = 51;
G4int           gConeMat = 1;
G4int           gHolderColor = 1;
G4int           gDOM = 1;
G4int           gPMT = 0;
G4bool          gPlaceHarness = true;
G4int           gHarness = 1;
G4int           gRopeNumber = 1;
G4double        gworldsize = 0.5 * m;
G4double        gElectronFactor = 9.5;
G4bool          gCADImport = true;
G4String        gHittype = "individual"; // seems like individual records each hit per pmt
G4bool          gVisual = true; // may be visualization on?
G4int           gEnvironment = 2; // 1 is ICeCUbe ice without any ice property, 2 is with property, 0 is air.
G4String        ghitsfilename = "/mnt/c/Users/Waly/bulkice_doumeki/"; //change location of your output file
//G4double        gSimulatedTime = 1.0;
G4int           gcounter = 0;
G4int           gPosCount = 0;
G4String        gQEFile = "../InputFile/TA0001_HamamatsuQE.data";
G4int           gNumCherenkov = 0;
G4int           gNumScint = 0;
G4int           gEvent = 0;
G4int           gDepthIndex = 88; //default value
G4int           gRunID = 0; //Run ID for each run
G4double        gPaintThickness = 0;
/**
variables for optical photon waves
**/
G4bool          gVis = false;
G4double        gDistance = 1 * m; //default distance for photon from symmetry axis
G4double        gZenithAngle = 0; //default zenith angle for photon wave
G4double        gStartAngle = 0;
G4double        gFinalAngle = 180;
G4double        gAngleIncrement = 10;
G4bool          gMultipleAngle = false;
G4bool          gWriteZenithAngle = false;
G4bool          gPhotonSim = false;

enum {PMT, MDOM, DOM, LOM16, LOM18, DEGG, WOM};

OMSimAnalysisManager gAnalysisManager;

void help() //needs change
{
    std::cout << "Usage: " << std::endl;

    std::cout << "./bulkice_doumeki" << " " << "[OM Model]" << " " << "[interaction channel]" << " " <<  "[depth index]" << " " << "[output folder] " << "[run id]"
    << " (for simulating supernovae neutrino flux)" << std::endl;

    std::cout << " " << std::endl;
    //std::cout << "./bulkice_doumeki" << " " << "[OM Model]" << "vis" << " (for visualization)" << std::endl;
    std::cout << "./bulkice_doumeki" << " " << "[OM Model]" << " " << "opticalphoton" << " " << "[depth index]" << " " << "[output folder] " << "[run id] "<< "[distance (m)] "
    << "[zenith angle (degree)]" << " " << " (for photon wave with a single zenith angle)" << std::endl;

    std::cout << " " << std::endl;

    std::cout << "./bulkice_doumeki " << "[OM Model] " << "opticalphoton " << "[depth index] " << "[output folder] " << "[run id] " << "[distance (m] "
    << "[start angle (degree)] " << "[final angle (degree)] " << "[increment (degree)] "
    << "(for photon wave with multiple zenith angle)" << std::endl;

    std::cout << " " << std::endl;

    std::cout << "Available OM Models: [dom, mdom, lom16, lom18, pmt, degg, wom]" << std::endl;

    std::cout << "Available interaction channels: [ibd, enees, all, radioactivity]" << std::endl;

    std::cout << "Available depth index: [0, 1, ......, 108]" << std::endl;

    std::cout << "To see associated depth with each depth index, go to build -> data -> Materials -> IceCubeIce.dat -> jDepth_spice" << std::endl;
}
void ParseCommandLine(int argc, char** argv, G4int& PMT_model, G4double& worldsize, G4String& interaction_channel)
{
    if(argc == 6 || argc ==  9/*change this t0 9 when done with paint thickness test */ || argc == 11 /*change this to 10 when done with paint thickness test*/ || argc == 2)
    {
        G4String model = argv[1];


        if(model == "help")
        {
            help();
            exit(0);
        }
        interaction_channel = argv[2];
        gDepthIndex = atoi(argv[3]);

        if(gDepthIndex < 0 || gDepthIndex > 108)
        {
            std::cerr << "Invalid depth index." << std::endl;
            help();
            exit(0);
        }

        G4String outputFolder = argv[4];
        gRunID = atoi(argv[5]);
        worldsize = 20 * m;
        ghitsfilename += outputFolder + "/";

        if(model == "dom")
        {
            std::cout << "*****DOM simulation selected*****" << std::endl;
            PMT_model = 2;
            gPMT = 5;
        }
        else if(model == "mdom")
        {
            std::cout << "*****MDOM simulation selected*****" << std::endl;
            PMT_model = 1;
            gPMT = 0;
        }
        else if(model == "lom16")
        {
            std::cout << "*****LOM16 simulation selected*****" << std::endl;
            PMT_model = 3;
            gPMT = 3;
        }
        else if(model == "lom18")
        {
            std::cout << "*****LOM18 simulation selected*****" << std::endl;
            PMT_model = 4;
            gPMT = 3;
        }
        else if(model == "pmt")
        {
            std::cout << "*****Single PMT simulation selected*****" << std::endl;
            PMT_model = 0;
        }
        else if(model == "degg")
        {
            std::cout << "*****DEGG simulation selected*****" << std::endl;
            PMT_model = 5;
            gPMT = 4;
        }
        else if(model == "wom")
        {
            std::cout << "*****WOM simulation selected*****" << std::endl;
            PMT_model = 6;
            gPMT = 5;
        }
        else
        {
            std::cout << "Invalid OM Model selected!" << std::endl;
            help();
            exit(0);
        }

        if(interaction_channel == "opticalphoton")
        {
            gDistance = atof(argv[6]);
            gPhotonSim = true;

            if(argc == 9)
            {
                gZenithAngle = atof(argv[7]);
                gStartAngle = gZenithAngle;
                gFinalAngle = gZenithAngle;
                gAngleIncrement = 10; //some dummy number, doesn't affect the output.
                gMultipleAngle = false;
                gPaintThickness = atof(argv[8]); // remove this when done with paint thickness experiment.

                ghitsfilename += model + "_" + interaction_channel + "_" + std::to_string(gZenithAngle) + "_"+ std::to_string(gRunID);
            }
            else
            {
                gMultipleAngle = true;
                gStartAngle = atof(argv[7]);
                gFinalAngle = atof(argv[8]);
                gAngleIncrement = atof(argv[9]);
                gPaintThickness = atof(argv[10]); //will remove later


                ghitsfilename += model + "_" + interaction_channel + "_" + std::to_string(gStartAngle) + "_"
                + std::to_string(gFinalAngle) + "_"+ std::to_string(gAngleIncrement) + "_"+ std::to_string(gRunID);
            }
        }
        else
        {
            /**
            neutrino flux simulations here
            **/
            ghitsfilename += model + "_" + interaction_channel + "_" + std::to_string(gRunID);
        }
    }
    else if(argc == 4)
    {
        PMT_model = 6;
        gPMT = 5;
        gVis = true;
        interaction_channel = argv[2];
        worldsize = 20;
        gZenithAngle = atof(argv[3]);
    }
    else
    {
        std::cerr << "Invalid number of argument given. Aborting...." << std::endl;
        help();
        exit(0);
    }
}
  /*  if(argc == 6 || argc == 3)
    {
        std::string model = argv[1];
        interaction_channel = argv[2];

        if(model == "dom")
        {
            std::cout << "*****DOM simulation selected*****" << std::endl;
            PMT_model = 2;
            gPMT = 5;
        }
        else if(model == "mdom")
        {
            std::cout << "*****MDOM simulation selected*****" << std::endl;
            PMT_model = 1;
            gPMT = 0;
        }
        else if(model == "lom16")
        {
            std::cout << "*****LOM16 simulation selected*****" << std::endl;
            PMT_model = 3;
            gPMT = 3;
        }
        else if(model == "lom18")
        {
            std::cout << "*****LOM18 simulation selected*****" << std::endl;
            PMT_model = 4;
            gPMT = 3;
        }
        else if(model == "pmt")
        {
            std::cout << "*****Single PMT simulation selected*****" << std::endl;
            PMT_model = 0;
        }
        else if(model == "degg")
        {
            std::cout << "*****DEGG simulation selected*****" << std::endl;
            PMT_model = 5;
            gPMT = 4;
        }
        else if(model == "wom")
        {
            std::cout << "*****WOM simulation selected*****" << std::endl;
            PMT_model = 6;
            gPMT = 5;
        }
        else if(model == "help")
        {
            help();
            exit(0);
        }
        else
        {
            std::cout << "Invalid OM Model selected!" << std::endl;
            help();
            exit(0);
        }

        if(argc == 6)
        {
            //macroname = argv[3];
            G4String outputFolder = argv[4];
            gDepthIndex = atoi(argv[3]);
            if(gDepthIndex < 0 || gDepthIndex > 108)
            {
                std::cerr << "Invalid depth index." << std::endl;
                help();
                exit(0);
            }
            gRunID = atoi(argv[5]);
            worldsize = 20;
            ghitsfilename += outputFolder+ "/" + model + "_" + interaction_channel + "_" + std::to_string(gRunID);
        }
        else
        {
            //worldsize = 0.25;
            worldsize = 2.5;
            ghitsfilename += model;
        }

    }
    else
    {
        std::cerr << "Invalid number of argument given. Aborting...." << std::endl;
        help();
        exit(0);
    }
}*/

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

std::vector<double> readColumnDouble (G4String fn, int col) {
    std::vector<double>	values;
    unsigned int c;
    double	a;
    c = col;
    std::ifstream	infile;
    std::vector<G4String> n;
    char l[256];
    G4String l2;
    infile.open(fn);
    while (infile.good() && !infile.eof()) {
        infile.getline(l,255);
        l2 = l;
        n = explode(l2,'\t');
        if (n.size()>=c) {
            a = atof(n.at(c-1));
            values.push_back(a);
        }
    }
    infile.close();

    return values;//values enthält den c. Wert aus fn (aus jeder Spalte,welche  nach 255 zeichen oder durch \n beendet wird?)
}

int main(int argc, char** argv)
{
    /*std::time_t currentTime = std::time(nullptr);
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%H%M%S", std::localtime(&currentTime));
    std::string time(buffer);*/




    G4double world_size(0);
    G4String macroname;
    G4int PMT_model(0);
    G4String interaction_channel;
    ParseCommandLine(argc, argv, PMT_model, world_size, interaction_channel);
    //ghitsfilename +=  "_" + time;

    OMSimRunManager* runManager = new OMSimRunManager(PMT_model, world_size, interaction_channel);
    runManager -> Initialize();

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    G4UIExecutive* ui = 0;

    if (!gVis) // needs change
    {
    // batch mode
       /* std::cerr << ":::::::::::::::::::Batch Mode Called:::::::::::::::::" << std::endl;
        G4String command = "/control/execute " + macroname;
        std::cout << "command " << command << std::endl;
        UImanager->ApplyCommand(command);*/
        runManager -> OpenFile();
        runManager -> BeamOn();
        runManager -> CloseFile();
        }
    else {
    // interactive mode
        std::cerr << "interactive mode called" << std::endl;
        ui = new G4UIExecutive(argc, argv);
        //UImanager->ApplyCommand("/control/execute ray_trace.mac");
        UImanager -> ApplyCommand("/control/execute vis.mac");
        ui-> SessionStart();
        delete ui;
    }

  //-----------------------
  // terminating...
  //-----------------------

   /* #ifdef G4VIS_USE
    delete vismanager;
    #endif

    delete runmanager;*/
    delete runManager;
    std::cout << "::::::::::::::this is the end:::::::::::::"<< std::endl;
    return 0;

}
