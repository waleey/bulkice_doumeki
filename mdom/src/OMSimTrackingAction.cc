#include "OMSimEventAction.hh"
#include "OMSimTrackingAction.hh"
#include "OMSimRunAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "OMSimRadioactivityData.hh"

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
    /**
    *temporal correlation between parent and daughter nucleus is preserved
    *if the mean life time of daughter nucleus is within the time window
    *if not, the decay time of the daughter is set to a random time t_rnd
    * drawn from a flat distribution between [o, t_window]
    *and the daughter is forced to decay immediately with initial time t_rnd
    *This is necessary to maintain the secular equilibrium
    *among the nucleus in the decay chain of each isotope.
    **/
    G4double timeWindow = OMSimRadioactivityData::ftimeWindow * s;
    OMSimRadioactivityData* radData = new OMSimRadioactivityData();
    if(aTrack -> GetCreatorProcess())
    {
           if(!(aTrack -> GetParticleDefinition() -> GetPDGStable()))
            {
                if(aTrack -> GetParticleDefinition() -> GetPDGLifeTime() >= timeWindow)
                {
                    aTrack -> GetDefinition() -> SetPDGLifeTime((radData -> GetInitialTime()) * s);
                    /*std::cout << aTrack -> GetParticleDefinition() -> GetPDGLifeTime() << std::endl;
                    std::cout << aTrack -> GetGlobalTime() << std::endl;
                    std::cout << aTrack -> GetParticleDefinition() -> GetParticleName() << std::endl;*/
                }
            }
    }
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
