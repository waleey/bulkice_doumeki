//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 coFDumpFllaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
//
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4Scintillation.hh
// Description:	Discrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//                         of energy deposited by particle type
//                         Thanks to Zach Hartwig (Department of Nuclear
//                         Science and Engineeering - MIT)
//              2005-07-28 add G4ProcessType to constructor
//              2002-11-21 change to user G4Poisson for small MeanNumPotons
//              2002-11-07 allow for fast and slow scintillation
//              2002-11-05 make use of constant material properties
//              2002-05-16 changed to inherit from VRestDiscreteProcess
//              2002-05-09 changed IsApplicable method
//              1999-10-29 add method and class descriptors
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef OMSIMSCINTILLATION_HH
#define OMSIMSCINTILLATION_HH

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4UImessenger.hh"
#include "G4EmSaturation.hh"

// Class Description:
// RestDiscrete Process - Generation of Scintillation Photons.
// Class inherits publicly from G4VRestDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class OMSimScintillation : public G4VRestDiscreteProcess
{

public:

	////////////////////////////////
	// Constructors and Destructor
	////////////////////////////////

	OMSimScintillation(const G4String& processName = "Scintillation",
                                 G4ProcessType type = fElectromagnetic);
	~OMSimScintillation();

private:

        OMSimScintillation(const OMSimScintillation &right);

        //////////////
        // Operators
        //////////////

        OMSimScintillation& operator=(const OMSimScintillation &right);

public:

        ////////////
        // Methods
        ////////////

        // G4Scintillation Process has both PostStepDoIt (for energy
        // deposition of particles in flight) and AtRestDoIt (for energy
        // given to the medium by particles at rest)

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable', for any particle type except
        // for an 'opticalphoton' and for short-lived particles

	G4double GetMeanFreePath(const G4Track& aTrack,
				       G4double ,
                                       G4ForceCondition* );
        // Returns infinity; i. e. the process does not limit the step,
        // but sets the 'StronglyForced' condition for the DoIt to be
        // invoked at every step.

        G4double GetMeanLifeTime(const G4Track& aTrack,
                                 G4ForceCondition* );
        // Returns infinity; i. e. the process does not limit the time,
        // but sets the 'StronglyForced' condition for the DoIt to be
        // invoked at every step.

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
			                const G4Step&  aStep);
        G4VParticleChange* AtRestDoIt (const G4Track& aTrack,
                                       const G4Step& aStep);

        // These are the methods implementing the scintillation process.

	void SetTrackSecondariesFirst(const G4bool state);
        // If set, the primary particle tracking is interrupted and any
        // produced scintillation photons are tracked next. When all
        // have been tracked, the tracking of the primary resumes.

        void SetFiniteRiseTime(const G4bool state);
        // If set, the G4Scintillation process expects the user to have
        // set the constant material property FAST/SLOWSCINTILLATIONRISETIME.

        G4bool GetTrackSecondariesFirst() const;
        // Returns the boolean flag for tracking secondaries first.

        G4bool GetFiniteRiseTime() const;
        // Returns the boolean flag for a finite scintillation rise time.

        void SetScintillationYieldFactor(const G4double yieldfactor);
        // Called to set the scintillation photon yield factor, needed when
        // the yield is different for different types of particles. This
        // scales the yield obtained from the G4MaterialPropertiesTable.

        G4double GetScintillationYieldFactor() const;
        // Returns the photon yield factor.

        void SetScintillationExcitationRatio(const G4double excitationratio);
        // Called to set the scintillation exciation ratio, needed when
        // the scintillation level excitation is different for different
        // types of particles. This overwrites the YieldRatio obtained
        // from the G4MaterialPropertiesTable.

        G4double GetScintillationExcitationRatio() const;
        // Returns the scintillation level excitation ratio.

        G4PhysicsTable* GetFastIntegralTable() const;
        // Returns the address of the fast scintillation integral table.

        G4PhysicsTable* GetSlowIntegralTable() const;
        // Returns the address of the slow scintillation integral table.

	G4PhysicsTable* GetfirstIntegralTable() const;

	G4PhysicsTable* GetsecondIntegralTable() const;

	G4PhysicsTable* GetthirdIntegralTable() const;

        void AddSaturation(G4EmSaturation* sat) { emSaturation = sat; }
        // Adds Birks Saturation to the process.

        void RemoveSaturation() { emSaturation = NULL; }
        // Removes the Birks Saturation from the process.

        G4EmSaturation* GetSaturation() const { return emSaturation; }
        // Returns the Birks Saturation.

        void SetScintillationByParticleType(const G4bool );
        // Called by the user to set the scintillation yield as a function
        // of energy deposited by particle type

        G4bool GetScintillationByParticleType() const
        { return scintillationByParticleType; }
        // Return the boolean that determines the method of scintillation
        // production

        void DumpPhysicsTable() const;
        // Prints the fast and slow scintillation integral tables.

protected:

        void BuildThePhysicsTable();
        // It builds either the fast or slow scintillation integral table;
        // or both.

        ///////////////////////
        // Class Data Members
        ///////////////////////

        G4PhysicsTable* theSlowIntegralTable;
        G4PhysicsTable* theFastIntegralTable;
        G4PhysicsTable* firstIntegralTable;
	G4PhysicsTable* secondIntegralTable;
	G4PhysicsTable* thirdIntegralTable;
        G4bool fTrackSecondariesFirst;
        G4bool fFiniteRiseTime;

        G4double YieldFactor;
	G4double FirstAmplitude;
        G4double SecondAmplitude;
        G4double ExcitationRatio;

        G4bool scintillationByParticleType;

private:

        G4double single_exp(G4double t, G4double tau2);
        G4double bi_exp(G4double t, G4double tau1, G4double tau2);

        // emission time distribution when there is a finite rise time
        G4double sample_time(G4double tau1, G4double tau2);

        G4EmSaturation* emSaturation;

};

////////////////////
// Inline methods
////////////////////

inline
G4bool OMSimScintillation::IsApplicable(const G4ParticleDefinition& aParticleType)
{
       if (aParticleType.GetParticleName() == "opticalphoton") return false;
       if (aParticleType.IsShortLived()) return false;

       return true;
}

inline
void OMSimScintillation::SetTrackSecondariesFirst(const G4bool state)
{
	fTrackSecondariesFirst = state;
}

inline
void OMSimScintillation::SetFiniteRiseTime(const G4bool state)
{
        fFiniteRiseTime = state;
}

inline
G4bool OMSimScintillation::GetTrackSecondariesFirst() const
{
        return fTrackSecondariesFirst;
}

inline
G4bool OMSimScintillation::GetFiniteRiseTime() const
{
        return fFiniteRiseTime;
}

inline
void OMSimScintillation::SetScintillationYieldFactor(const G4double yieldfactor)
{
        YieldFactor = yieldfactor;
}

inline
G4double OMSimScintillation::GetScintillationYieldFactor() const
{
        return YieldFactor;
}

inline
void OMSimScintillation::SetScintillationExcitationRatio(const G4double excitationratio)
{
        ExcitationRatio = excitationratio;
}

inline
G4double OMSimScintillation::GetScintillationExcitationRatio() const
{
        return ExcitationRatio;
}

inline
G4PhysicsTable* OMSimScintillation::GetSlowIntegralTable() const
{
        return theSlowIntegralTable;
}

inline
G4PhysicsTable* OMSimScintillation::GetFastIntegralTable() const
{
        return theFastIntegralTable;
}

inline
G4PhysicsTable* OMSimScintillation::GetfirstIntegralTable() const
{
        return firstIntegralTable;
}

inline
G4PhysicsTable* OMSimScintillation::GetsecondIntegralTable() const
{
        return secondIntegralTable;
}

inline
G4PhysicsTable* OMSimScintillation::GetthirdIntegralTable() const
{
        return thirdIntegralTable;
}

inline
void OMSimScintillation::DumpPhysicsTable() const
{
        if (theFastIntegralTable) {
           G4int PhysicsTableSize = theFastIntegralTable->entries();
           G4PhysicsOrderedFreeVector *v;

           for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
           {
        	v = (G4PhysicsOrderedFreeVector*)(*theFastIntegralTable)[i];
        	v->DumpValues();
           }
         }
/*
        if (theSlowIntegralTable) {
           G4int PhysicsTableSize = theSlowIntegralTable->entries();
           G4PhysicsOrderedFreeVector *v;

           for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
           {
                v = (G4PhysicsOrderedFreeVector*)(*theSlowIntegralTable)[i];
                v->DumpValues();
           }
         }
           if (firstIntegralTable) {
           G4int PhysicsTableSize = firstIntegralTable->entries();
           G4PhysicsOrderedFreeVector *v;

           for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
           {
                v = (G4PhysicsOrderedFreeVector*)(*firstIntegralTable)[i];
                v->DumpValues();
           }
         }
           if (secondIntegralTable) {
           G4int PhysicsTableSize = secondIntegralTable->entries();
           G4PhysicsOrderedFreeVector *v;

           for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
           {
                v = (G4PhysicsOrderedFreeVector*)(*secondIntegralTable)[i];
                v->DumpValues();
           }
         }
           if (thirdIntegralTable) {
           G4int PhysicsTableSize = thirdIntegralTable->entries();
           G4PhysicsOrderedFreeVector *v;

           for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
           {
                v = (G4PhysicsOrderedFreeVector*)(*thirdIntegralTable)[i];
                v->DumpValues();
           }
        }*/
}

inline
G4double OMSimScintillation::single_exp(G4double t, G4double tau2)
{
         return std::exp(-1.0*t/tau2)/tau2;
}

inline
G4double OMSimScintillation::bi_exp(G4double t, G4double tau1, G4double tau2)
{
         return std::exp(-1.0*t/tau2)*(1-std::exp(-1.0*t/tau1))/tau2/tau2*(tau1+tau2);
}

#endif /* mdomScintillation_h */
