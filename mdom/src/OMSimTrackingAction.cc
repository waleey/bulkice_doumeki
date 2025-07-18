#include "OMSimEventAction.hh"
#include "OMSimTrackingAction.hh"
#include "OMSimRunAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "OMSimRadioactivityData.hh"
#include "OMSimSteppingAction.hh"
extern G4bool QEFilter;
extern G4String gQEFile;


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
                if(aTrack -> GetParticleDefinition() -> GetPDGLifeTime() > OMSimRadioactivityData::ftimeWindow * s)
                //if(aTrack -> GetParticleDefinition() -> GetPDGLifeTime() > 60 * s)
                {
                    aTrack -> GetDefinition() -> SetPDGLifeTime((radData -> GetInitialTime()) * s);
                    /*std::cout << aTrack -> GetParticleDefinition() -> GetPDGLifeTime() << std::endl;
                    std::cout << aTrack -> GetGlobalTime() << std::endl;
                    std::cout << aTrack -> GetParticleDefinition() -> GetParticleName() << std::endl;*/
                }
            }
    }
    delete radData;
       
    if(QEFilter){
       if (aTrack -> GetDefinition()->GetParticleName()=="opticalphoton" ){
	  if(aTrack -> GetCreatorProcess() ->GetProcessName() == "Cerenkov")
                       {//std::cout<<"Optical Photon Produced by Cherenkov Process"<<std::endl;
			//std::cout<<"by Cherenkov Process"<<std::endl;
			G4int survived;
			G4double Ekin;
			G4double lambda;
			G4double hc = 1240*nm;
			Ekin = aTrack->GetKineticEnergy() ;
                	lambda = (hc/Ekin) * nm;
			pmt_qe -> ReadQeTable();
                	//std::cout << "+++++Wavelength: " << lambda / nm<< std::endl;
                	double qe = (pmt_qe -> GetQe(lambda)) / 100;
	        	double random = CLHEP::RandFlat::shoot(0.0, 1.0);
                	//std::cout << "++++++++++QE : " << qe << "++++" << std::endl;
                	survived = (random < (qe)) ? 1 : 0;
			//std::cout<<"-------Survived?"<<survived<<"----------"<<std::endl; 
                	if (survived == 0){
                    		//std::cout << "Killed Photon" << std::endl;
		    		const_cast<G4Track*>(aTrack)->SetTrackStatus(fStopAndKill);
				}
				}}}

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
