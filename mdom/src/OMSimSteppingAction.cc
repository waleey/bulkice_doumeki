#include "OMSimSteppingAction.hh"

#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4StepPoint.hh"

#include "OMSimAnalysisManager.hh"

extern OMSimAnalysisManager gAnalysisManager;
extern G4String	gHittype;
extern G4int gcounter;
extern G4String gQEFile;
extern G4int gPosCount;


OMSimSteppingAction::OMSimSteppingAction()
{

}


void OMSimSteppingAction::UserSteppingAction(const G4Step* aStep)
{
    G4Track* aTrack = aStep->GetTrack();

    //kill particles that are stuck... e.g. doing a loop in the pressure vessel
    if ( aTrack-> GetCurrentStepNumber() > 100000) {
        G4cout << "Particle stuck   " <<  aTrack->GetDefinition()->GetParticleName()  << " " << 1239.84193/(aTrack->GetKineticEnergy()/eV)<< G4endl;
        if ( aTrack->GetTrackStatus() != fStopAndKill ) {
            aTrack->SetTrackStatus(fStopAndKill);
        }
    }

    //	Check if optical photon is about to hit a photocathode, if so, destroy it and save the hit
    if ( aTrack->GetDefinition()->GetParticleName() == "opticalphoton" ) {

        if ( aTrack->GetTrackStatus() != fStopAndKill ) {

            if( aStep -> GetPostStepPoint() -> GetPhysicalVolume() -> GetName() == "Photocathode_pv_OMSIM") {

                G4double Ekin;
                G4double t1, t2;
                G4Track* aTrack = aStep -> GetTrack();
                G4ThreeVector vertex_pos;
                vertex_pos = aTrack -> GetVertexPosition();
                G4double hc = 1240 * nm;
                G4double lambda;
                /**
                *Will be removed soon
                **/
                std::string creator = aTrack -> GetCreatorProcess() -> GetProcessName();
                Ekin = aTrack->GetKineticEnergy() ;
                lambda = (hc/Ekin) * nm;
                pmt_qe -> ReadQeTable();
                std::cout << "+++++Wavelength: " << lambda / nm<< std::endl;
                double qe = (pmt_qe -> GetQe(lambda)) / 100;
                double random = CLHEP::RandFlat::shoot(0.0, 1.0);
                std::cout << "++++++++++QE : " << qe << "++++" << std::endl;
                G4int survived = (random < (qe)) ? 1 : 0;
                std::vector<G4String> n;
                extern std::vector<G4String> explode (G4String s, char d);
                G4ThreeVector deltapos;


                n = explode(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName(),'_');
                gAnalysisManager.stats_PMT_hit.push_back(atoi(n.at(1)));

                if (gHittype == "individual") {
                    deltapos = aTrack->GetVertexPosition() - aTrack->GetPosition();
                    t1 = aTrack->GetGlobalTime() /ns;
                    t2 = aTrack->GetLocalTime();
                    Ekin = aTrack->GetKineticEnergy();
                    gAnalysisManager.stats_photon_direction.push_back(aTrack->GetMomentumDirection());
                    gAnalysisManager.stats_photon_position.push_back(aTrack->GetPosition());
                    gAnalysisManager.stats_event_id.push_back(gAnalysisManager.current_event_id);
                    gAnalysisManager.stats_photon_flight_time.push_back(t2);
                    gAnalysisManager.stats_photon_track_length.push_back(aTrack->GetTrackLength()/m);
                    gAnalysisManager.stats_hit_time.push_back(t1);
                    gAnalysisManager.stats_photon_energy.push_back(Ekin/eV);
                    gAnalysisManager.stats_event_distance.push_back(deltapos.mag()/m);
                    gAnalysisManager.stats_vertex_position.push_back(vertex_pos);
                    gAnalysisManager.stats_positron_id.push_back(aTrack -> GetParentID());
                    gAnalysisManager.stats_survived_qe.push_back(survived);
                    //Will be removed soon
                    gAnalysisManager.stats_creator.push_back(creator);
                }

                aTrack->SetTrackStatus(fStopAndKill);
               //} 	// kills counted photon to prevent scattering and double-counting
            }
        }

    }

    {

    }
}
