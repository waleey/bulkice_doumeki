#include "OMSimPhysicsList.hh"
#include "OMSimScintillation.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4StepLimiter.hh"
#include "G4SystemOfUnits.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4PhysicsListHelper.hh"

#include "G4Cerenkov.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4eplusAnnihilation.hh"
#include "G4Scintillation.hh"
#include "G4OpWLS.hh"
#include "G4OpWLS3.hh"
#include "G4OpWLS2.hh"

#include "G4LivermoreIonisationModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4BraggIonModel.hh"
#include "G4BraggModel.hh"
#include "G4hIonisation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4ionIonisation.hh"
#include "G4hMultipleScattering.hh"

#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"
#include "WLSBoundaryProcess.hh"

#include "G4StepLimiterPhysics.hh"


OMSimPhysicsList::OMSimPhysicsList():  G4VUserPhysicsList()
{
	defaultCutValue = 0.1*mm;
	SetVerboseLevel(0);
	radioactiveList = new G4RadioactiveDecayPhysics();
}

OMSimPhysicsList::~OMSimPhysicsList()
{
}

void OMSimPhysicsList::ConstructParticle()
{
	G4Gamma::GammaDefinition();
	G4OpticalPhoton::OpticalPhotonDefinition();
    G4GenericIon::GenericIonDefinition();
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4Proton::ProtonDefinition();
    G4Neutron::NeutronDefinition();
    G4MuonMinus::MuonMinusDefinition();
    G4MuonPlus::MuonPlusDefinition();
	G4NeutrinoE::NeutrinoEDefinition();
	G4AntiNeutrinoE::AntiNeutrinoEDefinition();
	G4NeutrinoMu::NeutrinoMuDefinition();
	G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
	G4Alpha::AlphaDefinition();
}

void OMSimPhysicsList::ConstructProcess()
{

	G4PhysicsListHelper* plh = G4PhysicsListHelper::GetPhysicsListHelper();
	AddTransportation();

	/**
	* Temporarily adding step limiter phyiscs
	**/

    G4StepLimiterPhysics* stepLimiterProcess = new G4StepLimiterPhysics();
    stepLimiterProcess -> ConstructProcess();

//  The Radioactive Decay Process
	radioactiveList -> ConstructProcess();

//	The Cherenkov process
	G4Cerenkov* theCerenkovProcess = new G4Cerenkov("Cerenkov");
    theCerenkovProcess->SetTrackSecondariesFirst(false);
    theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    theCerenkovProcess->SetMaxNumPhotonsPerStep(10);

//The Scintillation Process
    /*G4Scintillation* theScintillationProcess = new G4Scintillation("Scintillation");
    theScintillationProcess -> SetTrackSecondariesFirst(false);
    theScintillationProcess -> SetScintillationByParticleType(true);*/

    OMSimScintillation* theScintillationProcess = new OMSimScintillation("Scintillation");
    theScintillationProcess -> SetTrackSecondariesFirst(false);
    /**
    *Keeping the scintillation turned off
    *for now.
    **/

//	The Livermore models
	G4eIonisation* theIonizationModel = new G4eIonisation();
	theIonizationModel->SetEmModel(new G4LivermoreIonisationModel());

//  The Ion Ionization model
    G4ionIonisation* theionionisation = new G4ionIonisation();
    theionionisation -> SetEmModel(new G4BraggIonModel());

    plh -> RegisterProcess(theionionisation, G4Alpha::AlphaDefinition());
    plh -> RegisterProcess(theionionisation, G4GenericIon::GenericIon());

//  The Hadron Ionisation Process
    G4hIonisation* theHadronIon = new G4hIonisation();
    theHadronIon -> SetEmModel(new G4BraggModel());

    plh -> RegisterProcess(theHadronIon, G4Proton::ProtonDefinition());
    plh -> RegisterProcess(theHadronIon, G4Neutron::NeutronDefinition());

//  The Multiple Scattering Process
    G4hMultipleScattering* theMultScat = new G4hMultipleScattering();

    plh -> RegisterProcess(theMultScat, G4Alpha::AlphaDefinition());
    plh -> RegisterProcess(theMultScat, G4GenericIon::GenericIon());

//  PhotoElectric Process
	G4PhotoElectricEffect* thePhotoElectricEffectModel = new G4PhotoElectricEffect();
    thePhotoElectricEffectModel->SetEmModel(new G4LivermorePhotoElectricModel());

//  The Neutron Capture Process
    G4NeutronCaptureProcess* theNeutronCaptureProcess = new G4NeutronCaptureProcess();

//  Atomic Deexcitation Process
    G4UAtomicDeexcitation* deexc = new G4UAtomicDeexcitation();

    deexc -> SetFluo(true);
    deexc -> SetAuger(true);
    deexc -> SetPIXE(true);
    G4LossTableManager::Instance() -> SetAtomDeexcitation(deexc);

//  Gamma Conversion Process
	G4GammaConversion* theGammaConversionModel = new G4GammaConversion();
    theGammaConversionModel->SetEmModel(new G4LivermoreGammaConversionModel());

//	Now assign processes to generated particles
    auto theParticleIterator=GetParticleIterator();
	theParticleIterator->reset();
	while( (*theParticleIterator)() ){

		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		G4double particleMass = particle->GetPDGMass();
		G4String particleType = particle->GetParticleType();
		G4double particleCharge = particle->GetPDGCharge();
		//std::cout << "PARTICLE NAMES ARE: " << particleName << std::endl;

		if (particleName == "opticalphoton") {
			pmanager->AddDiscreteProcess(new G4OpAbsorption());
			pmanager->AddDiscreteProcess(new G4OpBoundaryProcess());
			pmanager->AddDiscreteProcess(new G4OpRayleigh);
			pmanager->AddDiscreteProcess(new G4OpMieHG);
			//pmanager->AddDiscreteProcess(new G4OpWLS);
			pmanager->AddDiscreteProcess(new WLSBoundaryProcess);
			//pmanager->AddDiscreteProcess(new G4OpWLS2);
			//pmanager->AddDiscreteProcess(new G4OpWLS3);
		}
		else if (particleName == "gamma") {
			pmanager->AddDiscreteProcess(theGammaConversionModel);
    		pmanager->AddDiscreteProcess(new G4ComptonScattering());
    		pmanager->AddDiscreteProcess(thePhotoElectricEffectModel);
		}
      	else if (particleName == "e-") {
	      	pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
	      	pmanager->AddProcess(theIonizationModel,-1,2,2);
			pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
      	}
		else if (particleName == "e+") {
			pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
			pmanager->AddProcess(new G4eIonisation,-1,2,2);
			pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
			pmanager->AddProcess(new G4eplusAnnihilation, 0,-1, 4);
		}
		if(particleName == "neutron")
		{
            pmanager -> AddProcess(theNeutronCaptureProcess);
		}
		if (theCerenkovProcess->IsApplicable(*particle)) {
  			pmanager->AddProcess(theCerenkovProcess);
			pmanager->SetProcessOrdering(theCerenkovProcess, idxPostStep);
        }
    /**
    *Scintillation is only set for electron and alpha
    *Because scintillation yield for other particles are not known.
    **/
        if(theScintillationProcess -> IsApplicable(*particle))
        {
            /*if(particle -> GetParticleName() == "e-" || particle -> GetParticleName() == "alpha" || particle -> GetParticleName() == "e+")
            {
                std::cout << "SCINTILLLLLLLLLLLLATION" << std::endl;
                pmanager -> AddDiscreteProcess(theScintillationProcess);
                pmanager -> SetProcessOrderingToLast(theScintillationProcess, idxAtRest);
                pmanager -> SetProcessOrderingToLast(theScintillationProcess, idxPostStep);
            }*/
            //pmanager -> AddProcess(theScintillationProcess);
            //pmanager -> SetProcessOrdering(theScintillationProcess, idxPostStep);
            //pmanager -> SetProcessOrdering(theScintillationProcess, idxAtRest);
            /**
            *Commented out the Scint Process
            **/
        }

	}

}

void OMSimPhysicsList::SetCuts()
{
	SetCutsWithDefault();
	if (verboseLevel>0) DumpCutValuesTable();
}


