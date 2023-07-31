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
extern G4int gNumCherenkov;
extern G4int gNumScint;


OMSimRunAction::OMSimRunAction(){}
OMSimRunAction::~OMSimRunAction(){}

void OMSimRunAction::BeginOfRunAction(const G4Run* Run)
{
    G4int RunId = Run -> GetRunID();
	gAnalysisManager.SetRunId(RunId);

}

void OMSimRunAction::EndOfRunAction(const G4Run*)
{

	if (gHittype == "individual") {
		gAnalysisManager.Write(); // for K40 analysis
	}
	if (gHittype == "collective") {
		gAnalysisManager.WriteAccept(); // mainly for acceptance
	}
//  setting all arrays to zero
    gAnalysisManager.Reset();
}

