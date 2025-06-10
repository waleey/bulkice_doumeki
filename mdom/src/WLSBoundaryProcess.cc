////////////////////////////////////////////////////////////////////////
// Optical Photon Wavelength Shifting as a surface boundary process
////////////////////////////////////////////////////////////////////////
//
// File:        WLSBoundaryProcess.cc
// Description: Discrete Process -- Wavelength Shifting of Optical Photons on surface
// Version:     1.0
// Created:     11/28/23
// Author:      Waly M Z Karim
//              (Adaptation of G4OpWLS and G4OpBoundaryProcess)
//
////////////////////////////////////////////////////////////////////////


#include "WLSBoundaryProcess.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpProcessSubType.hh"
#include "G4Poisson.hh"
#include "G4OpticalParameters.hh"
#include "G4WLSTimeGeneratorProfileDelta.hh"
#include "G4WLSTimeGeneratorProfileExponential.hh"
#include "G4GeometryTolerance.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"

#include "bits/stdc++.h"
#include <cstdlib>

extern G4int gPhotonNotAbsorbed;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
WLSBoundaryProcess::WLSBoundaryProcess(const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type),
    fInnerMaterialName("Quartz"),
    fOuterMaterialName("Filler"),
    fPaintThickness(25 * micrometer)
{
  WLSTimeGeneratorProfile = nullptr;
  Initialise();
  SetProcessSubType(fOpBoundary); //might be changed to fOpWLS
  theIntegralTable = nullptr;

  if(verboseLevel > 0)
    G4cout << GetProcessName() << " is created " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
WLSBoundaryProcess::~WLSBoundaryProcess()
{
  if(theIntegralTable)
  {
    theIntegralTable->clearAndDestroy();
    delete theIntegralTable;
  }
  delete WLSTimeGeneratorProfile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WLSBoundaryProcess::PreparePhysicsTable(const G4ParticleDefinition&) { Initialise(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WLSBoundaryProcess::Initialise()
{
  G4OpticalParameters* params = G4OpticalParameters::Instance();
  SetVerboseLevel(params->GetWLSVerboseLevel());
  UseTimeProfile(params->GetWLSTimeProfile());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* WLSBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
                                         const G4Step& aStep)
/**
* PostStepDoIT of WLSBoundaryProcess
* No change of the track should be made unless
* The photon is on the Geometric Boundary
* inner and outer materials are appropriate
* The Photon Was Absorbed
* If the photon got absorbed but not re-emitted
* it should be killed.
**/
{
  G4Material* fInnerMaterial(0);
  G4Material* fOuterMaterial(0);

  fStatus = Undefined;
  aParticleChange.Initialize(aTrack);
  const G4Step* pStep = &aStep;

  if(pStep -> GetPostStepPoint() -> GetStepStatus() == fGeomBoundary)
  {
    fInnerMaterial = pStep -> GetPreStepPoint() -> GetMaterial();
    fOuterMaterial = pStep -> GetPostStepPoint() -> GetMaterial();

  }
  else
  {
    fStatus = NotAtBoundary;
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  if(!fInnerMaterial || !fOuterMaterial)
  {
    std::cout << "Material Not Defined. Aborting.." << std::endl;
    exit(0);
  }

  G4String innerMaterial = fInnerMaterial -> GetName();
  G4String outerMaterial = fOuterMaterial -> GetName();
  if(innerMaterial == "Quartz" || innerMaterial == "Filler")
  {
    if(outerMaterial == "Quartz" || outerMaterial == "Filler")
    {

    }
  }
  else
  {
    //G4cout << "Photon not on WLS Boundary.." << G4endl;
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  const G4Material* WLSMaterial;

  if(fInnerMaterial -> GetName() == "Quartz")
  {
    WLSMaterial = fInnerMaterial;
  }
  else
  {
    WLSMaterial = fOuterMaterial;
  }

  G4double paintThickness = GetPaintThickness();
  G4MaterialPropertiesTable* MPT = aTrack.GetMaterial() -> GetMaterialPropertiesTable();
  if(!MPT)
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  G4MaterialPropertyVector* attVector = MPT -> GetProperty(kWLSABSLENGTH);
  if(!attVector)
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  G4double thePhotonEnergy = aTrack.GetDynamicParticle()->GetTotalEnergy();
  G4double attLength = abs(attVector -> Value(thePhotonEnergy)) /1000;

  //Weird issue. Geant4 thinks the unit is mm where it should be micrometer.
  //Also, correcting for the negative absorption length

  G4double probWLS = 1 - exp(-paintThickness / attLength);
  G4double randomNum = G4UniformRand();

  if(randomNum > probWLS)
  {
    //std::cout << "Photon Energy: " << thePhotonEnergy << std::endl;
    if(aTrack.GetCreatorProcess())
    {
        std::cout << "Creator Process of the Photons: " << aTrack.GetCreatorProcess() -> GetProcessName() << std::endl;
    }
    gPhotonNotAbsorbed++;
    std::cout << "absorption lengths: " << attLength / micrometer << std::endl;
    std::cout << "random Num: " << randomNum << std::endl;
    std::cout << "prob: " << probWLS << std::endl;
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  /**
  *By this point, the photon should get absorbed.
  **/
  std::vector<G4Track*> proposedSecondaries;

  aParticleChange.ProposeTrackStatus(fStopAndKill);

  if(verboseLevel > 1)
  {
    G4cout << "\n** G4OpWLS: Photon absorbed! **" << G4endl;
  }

  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

  if(!MPT->GetProperty(kWLSCOMPONENT))
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  G4int NumPhotons = 1;
  if(MPT->ConstPropertyExists(kWLSMEANNUMBERPHOTONS))
  {
    G4double MeanNumberOfPhotons = MPT->GetConstProperty(kWLSMEANNUMBERPHOTONS);
    NumPhotons                   = G4int(G4Poisson(MeanNumberOfPhotons));
    if(NumPhotons <= 0)
    {
      // return unchanged particle and no secondaries
      aParticleChange.SetNumberOfSecondaries(0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
  }

  // Retrieve the WLS Integral for this material
  // new G4PhysicsFreeVector allocated to hold CII's
  G4double primaryEnergy = aTrack.GetDynamicParticle()->GetKineticEnergy();
  G4double WLSTime       = 0.;
  G4PhysicsFreeVector* WLSIntegral = nullptr;

  WLSTime     = MPT->GetConstProperty(kWLSTIMECONSTANT);
  WLSIntegral = (G4PhysicsFreeVector*) ((*theIntegralTable)(
    aTrack.GetMaterial()->GetIndex()));

  // Max WLS Integral
  G4double CIImax       = WLSIntegral->GetMaxValue();
  G4int NumberOfPhotons = NumPhotons;

  for(G4int i = 0; i < NumPhotons; ++i)
  {
    G4double sampledEnergy;
    // Make sure the energy of the secondary is less than that of the primary
    for(G4int j = 1; j <= 100; ++j)
    {
      // Determine photon energy
      G4double CIIvalue = G4UniformRand() * CIImax;
      sampledEnergy     = WLSIntegral->GetEnergy(CIIvalue);
      if(sampledEnergy <= primaryEnergy)
        break;
    }
    // If no such energy can be sampled, return one less secondary, or none
    if(sampledEnergy > primaryEnergy)
    {
      if(verboseLevel > 1)
      {
        G4cout << " *** G4OpWLS: One less WLS photon will be returned ***"
               << G4endl;
      }
      NumberOfPhotons--;
      if(NumberOfPhotons == 0)
      {
        if(verboseLevel > 1)
        {
          G4cout
            << " *** G4OpWLS: No WLS photon can be sampled for this primary ***"
            << G4endl;
        }
        // return unchanged particle and no secondaries
        aParticleChange.SetNumberOfSecondaries(0);
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
      continue;
    }
    else if(verboseLevel > 1)
    {
      G4cout << "G4OpWLS: Created photon with energy: " << sampledEnergy
             << G4endl;
    }
    // Generate random photon direction
    G4double cost = 1. - 2. * G4UniformRand();
    G4double sint = std::sqrt((1. - cost) * (1. + cost));
    G4double phi  = twopi * G4UniformRand();
    G4double sinp = std::sin(phi);
    G4double cosp = std::cos(phi);
    G4ParticleMomentum photonMomentum(sint * cosp, sint * sinp, cost);

    G4ThreeVector photonPolarization(cost * cosp, cost * sinp, -sint);
    G4ThreeVector perp = photonMomentum.cross(photonPolarization);

    phi                = twopi * G4UniformRand();
    sinp               = std::sin(phi);
    cosp               = std::cos(phi);
    photonPolarization = (cosp * photonPolarization + sinp * perp).unit();

    // Generate a new photon:
    auto sec_dp =
      new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), photonMomentum);
    sec_dp->SetPolarization(photonPolarization);
    sec_dp->SetKineticEnergy(sampledEnergy);

    G4double secTime = pPostStepPoint->GetGlobalTime() +
                       WLSTimeGeneratorProfile->GenerateTime(WLSTime);
    G4ThreeVector secPos = pPostStepPoint->GetPosition();
    G4Track* secTrack    = new G4Track(sec_dp, secTime, secPos);

    secTrack->SetTouchableHandle(aTrack.GetTouchableHandle());
    secTrack->SetParentID(aTrack.GetTrackID());

    proposedSecondaries.push_back(secTrack);
  }

  aParticleChange.SetNumberOfSecondaries((G4int)proposedSecondaries.size());
  for(auto sec : proposedSecondaries)
  {
    aParticleChange.AddSecondary(sec);
  }
  if(verboseLevel > 1)
  {
    G4cout << "\n Exiting from G4OpWLS::DoIt -- NumberOfSecondaries = "
           << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WLSBoundaryProcess::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if(theIntegralTable)
  {
    theIntegralTable->clearAndDestroy();
    delete theIntegralTable;
    theIntegralTable = nullptr;
  }

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  std::size_t numOfMaterials           = G4Material::GetNumberOfMaterials();
  theIntegralTable                     = new G4PhysicsTable(numOfMaterials);

  // loop for materials
  for(std::size_t i = 0; i < numOfMaterials; ++i)
  {
    auto physVector = new G4PhysicsFreeVector();

    // Retrieve vector of WLS wavelength intensity for
    // the material from the material's optical properties table.
    G4MaterialPropertiesTable* MPT =
      (*materialTable)[i]->GetMaterialPropertiesTable();
    if(MPT)
    {
      G4MaterialPropertyVector* wlsVector = MPT->GetProperty(kWLSCOMPONENT);
      if(wlsVector)
      {
        // Retrieve the first intensity point in vector
        // of (photon energy, intensity) pairs
        G4double currentIN = (*wlsVector)[0];
        if(currentIN >= 0.0)
        {
          // Create first (photon energy)
          G4double currentPM  = wlsVector->Energy(0);
          G4double currentCII = 0.0;
          physVector->InsertValues(currentPM, currentCII);

          // Set previous values to current ones prior to loop
          G4double prevPM  = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN  = currentIN;

          // loop over all (photon energy, intensity)
          // pairs stored for this material
          for(std::size_t j = 1; j < wlsVector->GetVectorLength(); ++j)
          {
            currentPM = wlsVector->Energy(j);
            currentIN = (*wlsVector)[j];
            currentCII =
              prevCII + 0.5 * (currentPM - prevPM) * (prevIN + currentIN);

            physVector->InsertValues(currentPM, currentCII);

            prevPM  = currentPM;
            prevCII = currentCII;
            prevIN  = currentIN;
          }
        }
      }
    }
    theIntegralTable->insertAt(i, physVector);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double WLSBoundaryProcess::GetMeanFreePath(const G4Track& aTrack, G4double,
                                  G4ForceCondition* condition)
/**
* a dummy attLength is returned
* if the photon is on the correct volume boundary
* condition will be forced to implement this process
**/
{
  G4double attLength       = DBL_MAX;

  const G4Step* aStep = aTrack.GetStep();
  const G4String outerMat = aStep -> GetPreStepPoint() -> GetMaterial() -> GetName();
  const G4String innerMat = aStep -> GetPostStepPoint() -> GetMaterial() -> GetName();

  if(outerMat == "Quartz" || outerMat == "Filler")
  {
    if(innerMat == "Quartz" || innerMat == "Filler")
    {
        *condition = Forced;
    }
  }
  //*condition = Forced;
  return attLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WLSBoundaryProcess::UseTimeProfile(const G4String name)
{
  if(WLSTimeGeneratorProfile)
  {
    delete WLSTimeGeneratorProfile;
    WLSTimeGeneratorProfile = nullptr;
  }
  if(name == "delta")
  {
    WLSTimeGeneratorProfile = new G4WLSTimeGeneratorProfileDelta("delta");
  }
  else if(name == "exponential")
  {
    WLSTimeGeneratorProfile =
      new G4WLSTimeGeneratorProfileExponential("exponential");
  }
  else
  {
    G4Exception("G4OpWLS::UseTimeProfile", "em0202", FatalException,
                "generator does not exist");
  }
  G4OpticalParameters::Instance()->SetWLSTimeProfile(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WLSBoundaryProcess::SetVerboseLevel(G4int verbose)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetWLSVerboseLevel(verboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String& WLSBoundaryProcess::GetOuterMaterialName()
{
    return fOuterMaterialName;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String& WLSBoundaryProcess::GetInnerMaterialName()
{
    return fInnerMaterialName;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double WLSBoundaryProcess::GetPaintThickness()
{
    return fPaintThickness;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WLSBoundaryProcess::SetOuterMaterialName(G4String& name)
{
    fOuterMaterialName = name;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WLSBoundaryProcess::SetInnerMaterialName(G4String& name)
{
    fInnerMaterialName = name;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WLSBoundaryProcess::SetPaintThickness(G4double thickness)
{
    fPaintThickness = thickness;
}
