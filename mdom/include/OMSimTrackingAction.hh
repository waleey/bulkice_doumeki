#ifndef OMSimTrackingAction_h
#define OMSimTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4OpticalPhoton.hh"
#include "OMSimSteppingAction.hh"
class OMSimTrackingAction : public G4UserTrackingAction
{
	public:
		OMSimTrackingAction();
		~OMSimTrackingAction();

		void PreUserTrackingAction(const G4Track*);
		void PostUserTrackingAction(const G4Track*);

	private:
        G4int counter;
        G4TrackingManager* ftrackingManager;
        OMSimPMTQE* pmt_qe = new OMSimPMTQE;
};

#endif
