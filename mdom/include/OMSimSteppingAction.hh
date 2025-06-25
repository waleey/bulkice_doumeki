#ifndef OMSimSteppingAction_h
#define OMSimSteppingAction_h 1
#include "G4Types.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4UserSteppingAction.hh"
#include "OMSimPMTQE.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern G4String gQEFile;

class OMSimSteppingAction : public G4UserSteppingAction
{
  public:
    OMSimSteppingAction(OMSimPMTQE* pmtQe);
   ~OMSimSteppingAction(){};

    void UserSteppingAction(const G4Step*);
    bool QEcheck(G4double lambda);

  private:
    OMSimPMTQE* pmt_qe;

    void WOMCheck(const G4Step*);
    long vesselID = 0;
    long pmtBodyID = 0;
    long processID = 0;
    long tempID;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
