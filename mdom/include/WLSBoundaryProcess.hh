#ifndef WLSBOUNDARYPROCESS_HH_INCLUDED
#define WLSBOUNDARYPROCESS_HH_INCLUDED

////////////////////////////////////////////////////////////////////////
// Optical Photon Wavelength Shifting as a surface boundary process
////////////////////////////////////////////////////////////////////////
//
// File:        WLSBoundaryProcess.hh
// Description: Discrete Process -- Wavelength Shifting of Optical Photons on surface
// Version:     1.0
// Created:     11/28/23
// Author:      Waly M Z Karim
//              (Adaptation of G4OpWLS and G4OpBoundaryProcess)
//
////////////////////////////////////////////////////////////////////////

#include "G4VDiscreteProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4RandomTools.hh"
#include "G4OpBoundaryProcess.hh"

class G4VWLSTimeGeneratorProfile;


class WLSBoundaryProcess : public G4VDiscreteProcess
{
 public:
  explicit WLSBoundaryProcess(const G4String& processName = "WLSBoundary",
                   G4ProcessType type          = fOptical);
  virtual ~WLSBoundaryProcess();

  virtual G4bool IsApplicable(
    const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable' only for an optical photon.

  virtual void BuildPhysicsTable(
    const G4ParticleDefinition& aParticleType) override;
  // Build the WLS integral table at the right time

  virtual G4double GetMeanFreePath(const G4Track& aTrack, G4double,
                                   G4ForceCondition*) override;
  // Returns the absorption length for WLS absorption of optical
  // photons in media with a specified attenuation length.

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep) override;
  // This is the method implementing WLS for optical photons.

  virtual G4PhysicsTable* GetIntegralTable() const;
  // Returns the address of the WLS integral table.

  virtual void DumpPhysicsTable() const;
  // Prints the WLS integral table.

  virtual void UseTimeProfile(const G4String name);
  // Selects the time profile generator

  virtual void PreparePhysicsTable(const G4ParticleDefinition&) override;
  virtual void Initialise();

  void SetVerboseLevel(G4int);

  G4String& GetOuterMaterialName();
  G4String& GetInnerMaterialName();
  G4double GetPaintThickness();

  void SetOuterMaterialName(G4String&);
  void SetInnerMaterialName(G4String&);
  void SetPaintThickness(G4double);

 protected:
  G4VWLSTimeGeneratorProfile* WLSTimeGeneratorProfile;
  G4PhysicsTable* theIntegralTable;

 private:
  WLSBoundaryProcess(const WLSBoundaryProcess& right) = delete;
  WLSBoundaryProcess& operator=(const WLSBoundaryProcess& right) = delete;

  std::size_t idx_wls = 0;

  G4String fInnerMaterialName;
  G4String fOuterMaterialName;

  G4double fPaintThickness;

  G4OpBoundaryProcessStatus fStatus;

};

////////////////////
// Inline methods
////////////////////

inline G4bool WLSBoundaryProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

inline G4PhysicsTable* WLSBoundaryProcess::GetIntegralTable() const
{
  return theIntegralTable;
}

inline void WLSBoundaryProcess::DumpPhysicsTable() const
{
  std::size_t PhysicsTableSize = theIntegralTable->entries();
  G4PhysicsFreeVector* v;

  for(std::size_t i = 0; i < PhysicsTableSize; ++i)
  {
    v = (G4PhysicsFreeVector*) (*theIntegralTable)[i];
    v->DumpValues();
  }
}


#endif // WLSBOUNDARYPROCESS_HH_INCLUDED
