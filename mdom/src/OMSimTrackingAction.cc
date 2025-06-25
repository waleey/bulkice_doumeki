#include "OMSimEventAction.hh"
#include "OMSimTrackingAction.hh"
#include "OMSimRunAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Ions.hh"

extern std::fstream gRadioDecayFile;
extern G4bool gVerbose;
extern G4bool gPhotonQEBiasing;
extern G4double gSimulationTime;
extern G4String gHittype;

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

    
    G4double timeWindow = gSimulationTime;
    
    // save track ID, particle ID and particle energy of all particles
    if (gRadioDecayFile.is_open()){
                G4int trackID = aTrack -> GetTrackID();
                G4String particleName = aTrack -> GetParticleDefinition() -> GetParticleName();
                G4double kineticEnergy = aTrack -> GetKineticEnergy();
                gRadioDecayFile << trackID << "," << particleName << "," << kineticEnergy / CLHEP::MeV << "\n";
            }


    if (gPhotonQEBiasing)
    {
        if (aTrack -> GetParticleDefinition() -> GetParticleName() == "opticalphoton")// and gVerbose)
        {

            G4Track* mTrack = const_cast<G4Track*>(aTrack); // make track modifiable

            G4double photonKineticEnergy = aTrack->GetKineticEnergy(); // kinetic energy
            G4double photonWavelength = (fhc/photonKineticEnergy); // wavelength

            if (photonWavelength < fWavelengthAcceptanceLow || photonWavelength > fWavelengthAcceptanceHigh)
            {
                if (gVerbose) { std::cout << "+++ (PRE-TRACK) Photon killed because it's outside the wavelength acceptance range!\nWavelength [nm]: " << photonWavelength / nm << ", [" << fWavelengthAcceptanceLow / nm << "," << fWavelengthAcceptanceHigh / nm << "]" << std::endl; }
                mTrack -> SetTrackStatus(fStopAndKill);
            }

            else
            {
                G4double random = CLHEP::RandFlat::shoot(0.0, 1.0); // random value between 0 and 1
                G4bool survival = (random < (fQuantumEfficiencyMax)) ? true : false; // rejection sampling

                if (!survival)
                {
                    if (gVerbose) { std::cout << "+++ (PRE-TRACK) Photon killed because of biasing!\nWavelength [nm]: " << photonWavelength / nm << ", " << random << ">" << fQuantumEfficiencyMax << std::endl; }
                    mTrack -> SetTrackStatus(fStopAndKill);
                }
                else
                {
                    if (gVerbose) { std::cout << "+++ (PRE-TRACK) Photon survived!\nWavelength [nm]: " << photonWavelength / nm << ", " << random << "<" << fQuantumEfficiencyMax << std::endl; }
                }
            }
            /*
            pmt_qe -> ReadQeTable();
            G4double qe = (pmt_qe -> GetQe(lambda)) / 100; // pdf value
            G4double random = CLHEP::RandFlat::shoot(0.0, 1.0); // random value between 0 and 1
            G4bool survived = (random < (qe)) ? true : false; // rejection sampling
            G4ThreeVector position = aTrack -> GetPosition();
            if (gVerbose)
            {
                std::cout << "+++ (PRE-TRACK) Optical photon with wavelength [nm]: " << lambda / nm << ", QE: " << qe << ", Random: " << random << ", Survived: " << survived << std::endl;
                std::cout << "                Position: x=" << position[0] / m << "m, y=" << position[1] / m << "m, z=" << position[2] / m << std::endl;
            }

            if (!survived) // kill photon
            {
                G4Track* mTrack = const_cast<G4Track*>(aTrack);
                mTrack -> SetTrackStatus(fStopAndKill);
                if (gVerbose)
                {
                    std::cout << "+++ (PRE-TRACK) QE Biasing killed optical photon!" << std::endl;
                }
            }
            else 
            {
                std::cout << "+++ (INFO) Optical photons survived! What happens next?" << std::endl;
                std::cout << "+++ (INFO) gHittype=" << gHittype << std::endl;
            }
            */
        }
    }

    if (aTrack -> GetParticleDefinition() -> GetParticleName() == "gamma" and gVerbose) // prints the energy for each gamma
    {
        std::cout << "+++ (PRE-TRACK) Gamma Energy [keV]: " << aTrack -> GetKineticEnergy() * 1000 << std::endl;
    }

    if (aTrack -> GetParticleDefinition() -> GetPDGCharge() > 2.  and gVerbose) // print every isotope name
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
            G4double excitationEnergy = ((const G4Ions*)(mTrack -> GetParticleDefinition())) -> GetExcitationEnergy(); // in MeV
            G4double kineticEnergy = mTrack -> GetKineticEnergy(); // in MeV
            // temporarily set lifetime to zero to force immediate decay
            G4double originalLifeTime = mTrack -> GetDefinition() -> GetPDGLifeTime();

            G4double globalTime = mTrack -> GetGlobalTime(); // global time of track
            G4double remainTime = timeWindow - globalTime; // remaining time
            G4double u = CLHEP::RandFlat::shoot(0.0,1.0); // random number between 0 and 1
            G4double decayTime = -log(1-u) * originalLifeTime; // sampled decay time

            if (gVerbose){
                std::cout << "+++ (PRE-TRACK) Global Time [ns]: " << globalTime 
                        << " ||| Time Window [ns]: " << timeWindow
                        << " ||| Remaining Time [ns]: " << remainTime
                        << std::endl;
                std::cout << "+++ (PRE-TRACK) Mean Life Time [ns]: " << originalLifeTime 
                        << " ||| Decay Time [ns]: " << decayTime
                        << std::endl;
                std::cout << "+++ (PRE-TRACK) Excitation Energy [keV]: " << excitationEnergy * 1000 
                        << " ||| Kinetic Energy [keV]: " << kineticEnergy * 1000
                        << std::endl;
            }
             
            if (decayTime != 0.0)
            {
                if (decayTime > remainTime)
                {
                    // if decay outside simulation window, kill immediately
                    mTrack -> SetTrackStatus(fStopAndKill);
                    if (gVerbose){
                        std::cout << ".... ---- |||| STOP AND KILL |||| ---- ...." << std::endl;
                    }
                }   
                else
                {
                    // update global time with the decay time
                    G4double newTime = globalTime + decayTime;
                    mTrack -> SetGlobalTime(newTime);

                    // Store the global time for restoration in PostUserTrackingAction
                    fOriginalGlobalTimes[mTrack->GetDefinition()] = newTime;

                    // Store the original lifetime for restoration in PostUserTrackingAction
                    fOriginalLifeTimes[mTrack->GetDefinition()] = originalLifeTime;

                    // force immediate decay
                    G4ParticleDefinition* mParticleDefinition = const_cast<G4ParticleDefinition*>(aTrack->GetParticleDefinition());
                    mParticleDefinition -> SetPDGLifeTime(0 * ns);
                }
            }
        }
    }    
}

void OMSimTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
    
    // Check if the particle definition had its lifetime temporarily changed
    auto it_lt = fOriginalLifeTimes.find(aTrack->GetParticleDefinition());
    if (it_lt != fOriginalLifeTimes.end()) {

        // remove const qualifier from GetParticleDefinition
        G4ParticleDefinition* mParticleDefinition = const_cast<G4ParticleDefinition*>(aTrack->GetParticleDefinition());

        // Restore the original lifetime
        mParticleDefinition -> SetPDGLifeTime(it_lt->second);

        // Remove the entry from the map
        fOriginalLifeTimes.erase(it_lt);

        if (gVerbose){                    
            std::cout << "+++ (POST-TRACK) Restored lifetime for "
                    << aTrack->GetParticleDefinition()->GetParticleName()
                    << " to " << it_lt->second << " ns." << std::endl;
        }
    }

    // Check if the particle definition had its global time changed
    auto it_gt = fOriginalGlobalTimes.find(aTrack->GetParticleDefinition());
    if (it_gt != fOriginalGlobalTimes.end()) {

        // remove const qualifier from GetParticleDefinition
        G4Track* mTrack = const_cast<G4Track*>(aTrack);

        // Restore the original lifetime
        mTrack -> SetGlobalTime(it_gt->second);

        // Remove the entry from the map
        fOriginalGlobalTimes.erase(it_gt);

        if (gVerbose){
            std::cout << "+++ (POST-TRACK) Updated global time for parent "
                    << aTrack->GetParticleDefinition()->GetParticleName()
                    << " to " << it_gt->second << " ns." << std::endl;
        }
        // Update global time for all secondaries
        const auto* secondaries = aTrack->GetStep()->GetSecondaryInCurrentStep();
        for (size_t i = 0; i < secondaries->size(); ++i)
        {
            G4Track* secondaryTrack = const_cast<G4Track*>((*secondaries)[i]);
            secondaryTrack->SetGlobalTime(it_gt->second);

            if (gVerbose){
                std::cout << "+++ (POST) Updated global time for daughter "
                        << secondaryTrack->GetDefinition()->GetParticleName()
                        << " to " << it_gt->second << " ns." << std::endl;
            }
        }
    }
    
}
