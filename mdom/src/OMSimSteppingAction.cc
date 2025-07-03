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
extern G4int gNumCherenkov;
extern G4int gNumScint;
extern G4bool gWOM;
extern G4bool gVerbose;
extern G4bool gPhotonQEBiasing;
extern G4double gSimulationTime;

G4int gVesselCount = 0;
G4int gPMTBodyCount = 0;
G4int gWLSCount = 0;
G4int gWLSAndHitCount = 0;


OMSimSteppingAction::OMSimSteppingAction(OMSimPMTQE* pmtQe)
{
    pmt_qe = pmtQe;
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

       //std::cout << "Photons are in: " << aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() << std::endl;
       //WOMCheck(aStep);
       std::string creator;

        if ( aTrack->GetTrackStatus() != fStopAndKill ) {

            if((aStep -> GetPostStepPoint() -> GetMaterial() -> GetName() == "RiAbs_Photocathode") && (aTrack->GetGlobalTime() >= 0) && (aTrack->GetGlobalTime() <= gSimulationTime)) {

                G4double Ekin = aTrack->GetKineticEnergy();
                G4double t1 = aTrack->GetGlobalTime();
                G4double t2 = aTrack->GetLocalTime();
                G4ThreeVector vertex_pos = aTrack -> GetVertexPosition();
		        G4double hc = 1239.84198433 * eV * nm; // hc constant (unit eV*nm)
                G4double quantumEfficiencyMax = 0.25871250; //

                G4double lambda = (hc/Ekin);
                G4double qe = (pmt_qe -> GetQe(lambda)) / 100;
                G4double random = CLHEP::RandFlat::shoot(0.0, 1.0);
                G4int survived;
                if (gPhotonQEBiasing) { survived = (random < (qe/quantumEfficiencyMax)) ? 1 : 0; }
                else { survived = (random < (qe)) ? 1 : 0; }

                if (!survived && gPhotonQEBiasing && gVerbose){ std::cout << "(STEP): Photon killed because of biasing!\n" << random << ">" << qe/quantumEfficiencyMax << std::endl; }
                if (!survived && gVerbose){ std::cout << "(STEP): Photon killed because of biasing!\n" << random << ">" << qe << std::endl; }

                if (survived && gPhotonQEBiasing && gVerbose){ std::cout << "(STEP): Photon survived because of biasing!\n" << random << "<" << qe/quantumEfficiencyMax << std::endl; }
                if (survived && gVerbose){ std::cout << "(STEP): Photon survived because of biasing!\n" << random << "<" << qe << std::endl; }

                std::vector<G4String> n;
                extern std::vector<G4String> explode (G4String s, char d);
                G4ThreeVector deltapos = aTrack->GetVertexPosition() - aTrack->GetPosition();

                if(!gWOM)
                {
                    n = explode(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName(),'_');
                    gAnalysisManager.stats_PMT_hit.push_back(atoi(n.at(1)));
                }
                else
                {
                    G4String pmtName = aStep -> GetPostStepPoint() -> GetPhysicalVolume() -> GetName();
                    if(pmtName == "0_physical")
                    {
                        gAnalysisManager.stats_PMT_hit.push_back(0);
                    }
                    else if(pmtName == "1_physical")
                    {
                        gAnalysisManager.stats_PMT_hit.push_back(1);
                    }
                    else
                    {
                        gAnalysisManager.stats_PMT_hit.push_back(100); //dummy number for invalid PMT
                    }
                }

                if (gHittype == "individual") {
                    if (gVerbose){
                        std::cout << "+++ (STEP) Optical Photon reached photocathode!" << std::endl;
                    }
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
                   //gAnalysisManager.stats_creator.push_back(creator);
                }
                
                else if (gHittype == "binary")
                {
                    if (survived)
                    {
                        gAnalysisManager.stats_PMT.push_back(static_cast<uint8_t>(atoi(n.at(1))));
                        gAnalysisManager.stats_TIME.push_back(aTrack->GetGlobalTime() /ns);
                    }
                }
                aTrack->SetTrackStatus(fStopAndKill);
                //} 	// kills counted photon to prevent scattering and double-counting

            }
        }

    }

}

void OMSimSteppingAction::WOMCheck(const G4Step* aStep)
{
    tempID = aStep -> GetTrack() -> GetTrackID();
    if(aStep -> GetTrack() -> GetTrackStatus() != fStopAndKill)
    {
        if(aStep -> GetPostStepPoint() -> GetPhysicalVolume() -> GetName() == "glassPhysical")
        {
            if(tempID != vesselID)
            {
                gVesselCount++;
                vesselID = tempID;
            }
        }

        if(aStep -> GetTrack() -> GetCreatorProcess())
        {
            if(aStep -> GetTrack() -> GetCreatorProcess() -> GetProcessName() == "WLSBoundary")
            {
                if(tempID != processID)
                {
                    gWLSCount++;
                    processID = tempID;
                }
            }
        }

    }
    if(aStep -> GetPostStepPoint() -> GetMaterial() -> GetName() == "NoOptic_Absorber")
    {
        if(tempID != pmtBodyID)
        {
            gPMTBodyCount++;
            pmtBodyID = tempID;
        }
    }

}
