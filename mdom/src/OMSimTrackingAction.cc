#include "OMSimEventAction.hh"
#include "OMSimTrackingAction.hh"
#include "OMSimRunAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "OMSimRadioactivityData.hh"

extern std::fstream gRadioDecayFile;

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
    
    if (aTrack -> GetParticleDefinition() -> GetParticleName() == "gamma") {
        std::cout << "+++ (PRE) Gamma Energy [keV]: " << aTrack -> GetKineticEnergy() * 1000 << std::endl;
    }
    
    if (gRadioDecayFile.is_open()){
                G4int trackID = aTrack -> GetTrackID();
                G4String particleName = aTrack -> GetParticleDefinition() -> GetParticleName();
                G4double kineticEnergy = aTrack -> GetKineticEnergy();
                gRadioDecayFile << trackID << "," << particleName << "," << kineticEnergy / CLHEP::MeV << "\n";
            }

    if (aTrack -> GetParticleDefinition() -> GetPDGCharge() > 2.) // print every isotope name
    {
        std::cout << "-------------------\n" << aTrack -> GetParticleDefinition()->GetParticleName() << "\n-------------------" << std::endl; 
    }
    
    if(aTrack -> GetCreatorProcess())
    {
        // Remove the const qualifier to make the track modifiable
        G4Track* mTrack = const_cast<G4Track*>(aTrack);
        // check if particle is instable
        if(!(mTrack -> GetParticleDefinition() -> GetPDGStable()))
        {  
            //G4double excitationEnergy = ((const G4Ions*)(mTrack -> GetParticleDefinition())) -> GetExcitationEnergy(); // in MeV
            G4double meanLife = mTrack -> GetParticleDefinition() -> GetPDGLifeTime(); 
            std::cout << "+++ Mean Life Time [ns]: " << meanLife << std::endl;
            //std::cout << "+++ Excitation Energy [keV]: " << excitationEnergy * 1000 << std::endl;
            if (meanLife * 10 > timeWindow) // || excitationEnergy == 0)) // long lived states and ground states are killed
            {
                mTrack -> SetTrackStatus(fStopAndKill);
                std::cout << ".... ---- |||| STOP AND KILL |||| ---- ...." << std::endl;
            }
        }
    }
    delete radData;
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
