#include "OMSimEventAction.hh"
#include "OMSimTrackingAction.hh"
#include "OMSimRunAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"


OMSimTrackingAction::OMSimTrackingAction()
:G4UserTrackingAction(), counter(0)
{
   // ftrackingManager = new G4TrackingManager();
}

OMSimTrackingAction::~OMSimTrackingAction()
{
}

void OMSimTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

}

void OMSimTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
    /*G4TrackVector* secTracks = ftrackingManager -> GimmeSecondaries();
    size_t n_secondaries = (*secTracks).size();
    G4cout << "SIZE " << n_secondaries << G4endl;
   if(secTracks)
    {
        for(size_t i = 0; i < n_secondaries; i++)
        {
            if((*secTracks)[i] -> GetDefinition() == G4OpticalPhoton::Definition())
            {
                counter++;
            }
        }
    }

    G4cout << "***********Total Number of Secondary Produced: " << counter << " **************" << G4endl;*/
}
