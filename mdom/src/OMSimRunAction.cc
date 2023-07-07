#include "OMSimRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include <ctime>
#include <sys/time.h>

#include "OMSimAnalysisManager.hh"
#include <time.h>
#include <sys/time.h>
extern G4String	ghitsfilename;
extern G4String	gHittype;
extern OMSimAnalysisManager gAnalysisManager;
extern G4int gcounter;
extern G4int gPosCount;


OMSimRunAction::OMSimRunAction(){}
OMSimRunAction::~OMSimRunAction(){}

void OMSimRunAction::BeginOfRunAction(const G4Run* Run)
{
    G4cout << ":::::::::This is the beginning of Run Action::::::::" << G4endl;
    gcounter = 0;
    gPosCount = 0;
    startingtime = clock() / CLOCKS_PER_SEC;
    G4int RunId = Run -> GetRunID();
    G4String filename = ghitsfilename.c_str() + std::to_string(RunId) + ".dat";
	gAnalysisManager.datafile.open(filename, std::ios::out /*| std::ios::app*/);
	gAnalysisManager.SetRunId(RunId);

}

void OMSimRunAction::EndOfRunAction(const G4Run*)
{
	//G4cout << "***********Total Number of Secondary Produced: " << gcounter << " **************" << G4endl;
	G4cout << "::::::::::::This is the end of Run Action:::::::::::" << G4endl;
	if (gHittype == "individual") {
		gAnalysisManager.Write(); // for K40 analysis
	}
	if (gHittype == "collective") {
		gAnalysisManager.WriteAccept(); // mainly for acceptance
	}

// 	Close output data file
gAnalysisManager.datafile.close();
gAnalysisManager.Reset();
double finishtime=clock() / CLOCKS_PER_SEC;
//std::cout << "Total Positron Count in World Volume:: " << gPosCount << std::endl;
G4cout << "Computation time: " << finishtime-startingtime << " seconds." << G4endl;
}

