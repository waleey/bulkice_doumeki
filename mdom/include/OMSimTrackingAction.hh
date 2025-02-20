#ifndef OMSimTrackingAction_h
#define OMSimTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4OpticalPhoton.hh"

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
		std::map<const G4ParticleDefinition*, G4double> fOriginalLifeTimes;
		std::map<const G4ParticleDefinition*, G4double> fOriginalGlobalTimes;
};

#endif
