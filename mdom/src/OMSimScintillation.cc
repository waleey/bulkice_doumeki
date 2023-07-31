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
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        mdomScintillation.cc
// Description: RestDiscrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//              of energy deposited by particle type
//              Thanks to Zach Hartwig (Department of Nuclear
//              Science and Engineeering - MIT)
//              2010-09-22 by Peter Gumplinger
//              > scintillation rise time included, thanks to
//              > Martin Goettlich/DESY
//              2005-08-17 by Peter Gumplinger
//              > change variable name MeanNumPhotons -> MeanNumberOfPhotons
//              2005-07-28 by Peter Gumplinger
//              > add G4ProcessType to constructor
//              2004-08-05 by Peter Gumplinger
//              > changed StronglyForced back to Forced in GetMeanLifeTime
//              2002-11-21 by Peter Gumplinger
//              > change to use G4Poisson for small MeanNumberOfPhotons
//              2002-11-07 by Peter Gumplinger
//              > now allow for fast and slow scintillation component
//              2002-11-05 by Peter Gumplinger
//              > now use scintillation constants from G4Material
//              2002-05-09 by Peter Gumplinger
//              > use only the PostStepPoint location for the origin of
//                scintillation photons when energy is lost to the medium
//                by a neutral particle
//              2000-09-18 by Peter Gumplinger
//              > change: aSecondaryPosition=x0+rand*aStep.GetDeltaPosition();
//                        aSecondaryTrack->SetTouchable(0);
//              2001-09-17, migration of Materials to pure STL (mma)
//              2003-06-03, V.Ivanchenko fix compilation warnings
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTypes.hh"
#include "G4EmProcessSubType.hh"
//#include "mdomAnalysisManager.hh"
#include "OMSimScintillation.hh"
//extern MdomAnalysisManager gAnalysisManager;
/////////////////////////
// Class Implementation
/////////////////////////

//////////////
// Operators
//////////////

// mdomScintillation::operator=(const mdomScintillation &right)
// {
// }

/////////////////
// Constructors
/////////////////
extern G4double gElectronFactor;


OMSimScintillation::OMSimScintillation(const G4String& processName,
                                     G4ProcessType type)
: G4VRestDiscreteProcess(processName, type)
{       G4cout << gElectronFactor << G4endl;
    SetProcessSubType(fScintillation);

    fTrackSecondariesFirst = false;
    fFiniteRiseTime = false;

    YieldFactor = 1.0;
    ExcitationRatio = 1.0;

    scintillationByParticleType = false;

    theFastIntegralTable = NULL;
    theSlowIntegralTable = NULL;
    firstIntegralTable = NULL;
    secondIntegralTable = NULL;
    thirdIntegralTable = NULL;

    if (verboseLevel>0) {
        G4cout << GetProcessName() << " is created " << G4endl;
    }

    BuildThePhysicsTable();

    emSaturation = NULL;
}

////////////////
// Destructors
////////////////

OMSimScintillation::~OMSimScintillation()
{
    if (theFastIntegralTable != NULL) {
        theFastIntegralTable->clearAndDestroy();
        delete theFastIntegralTable;
    }
    if (theSlowIntegralTable != NULL) {
        theSlowIntegralTable->clearAndDestroy();
        delete theSlowIntegralTable;
    }
    if (firstIntegralTable != NULL) {
        firstIntegralTable->clearAndDestroy();
        delete firstIntegralTable;
    }
    if (secondIntegralTable != NULL) {
        secondIntegralTable->clearAndDestroy();
        delete secondIntegralTable;
    }
    if (thirdIntegralTable != NULL) {
        thirdIntegralTable->clearAndDestroy();
        delete thirdIntegralTable;
    }

}

////////////
// Methods
////////////

// AtRestDoIt
// ----------
//
G4VParticleChange*
OMSimScintillation::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine simply calls the equivalent PostStepDoIt since all the
// necessary information resides in aStep.GetTotalEnergyDeposit()

{
    return OMSimScintillation::PostStepDoIt(aTrack, aStep);
}

// PostStepDoIt
// -------------
//
G4VParticleChange*
OMSimScintillation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine is called for each tracking step of a charged particle
// in a scintillator. A Poisson/Gauss-distributed number of photons is
// generated according to the scintillation yield formula, distributed
// evenly along the track segment and uniformly into 4pi.

{
    aParticleChange.Initialize(aTrack);

    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    const G4Material* aMaterial = aTrack.GetMaterial();

    //G4cout << aTrack.GetDefinition()->GetParticleName() << G4endl;

    G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
    G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

    G4ThreeVector x0 = pPreStepPoint->GetPosition();
    G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
    G4double      t0 = pPreStepPoint->GetGlobalTime();

    G4double TotalEnergyDeposit = aStep.GetTotalEnergyDeposit();

    //G4cout << TotalEnergyDeposit << G4endl;

    G4MaterialPropertiesTable* aMaterialPropertiesTable =
    aMaterial->GetMaterialPropertiesTable();

    if (!aMaterialPropertiesTable)
        return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);



    G4MaterialPropertyVector* lScintSpectrum = aMaterialPropertiesTable->GetProperty("SCINTILLATIONSPECTRUM");
//     G4cout<< ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << G4endl;
//     G4cout << lScintSpectrum->Energy(0) << " " <<lScintSpectrum->Energy(2) << G4endl;
//     G4cout << lScintSpectrum->Value(lScintSpectrum->Energy(0)) << " " <<lScintSpectrum->Value(lScintSpectrum->Energy(2)) <<G4endl;
    //lScintSpectrum->DumpValues();
    G4int nscnt = 0;
    G4MaterialPropertyVector* lScintComponents = aMaterialPropertiesTable->GetProperty("FRACTIONLIFETIMES");

    if (!lScintSpectrum )
        return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

    if (lScintSpectrum) nscnt = lScintComponents->GetVectorLength();

    G4double ScintillationYield = 0.;

    ScintillationYield = aMaterialPropertiesTable->
    GetConstProperty("SCINTILLATIONYIELD");


    // Units: [# scintillation photons / MeV]
    G4ParticleDefinition *pDef = aParticle->GetDefinition();

    if (pDef==G4Electron::ElectronDefinition() || pDef==G4Gamma::GammaDefinition()){
        //G4cout << "I got the yield "<< ScintillationYield*gElectronFactor << G4endl ;
        ScintillationYield = ScintillationYield*gElectronFactor;
        ScintillationYield *= YieldFactor;

    }
    else{
        ScintillationYield *= YieldFactor;
    }


    G4double ResolutionScale    = aMaterialPropertiesTable->
    GetConstProperty("RESOLUTIONSCALE");

    G4double MeanNumberOfPhotons = ScintillationYield*TotalEnergyDeposit;


    G4int NumPhotons;

    if (MeanNumberOfPhotons > 10.)
    {
        G4double sigma = ResolutionScale * std::sqrt(MeanNumberOfPhotons);
        NumPhotons = G4int(G4RandGauss::shoot(MeanNumberOfPhotons,sigma)+0.5);

    }
    else
    {
        NumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));

    }

    if (NumPhotons <= 0)
    {
        // return unchanged particle and no secondaries

        aParticleChange.SetNumberOfSecondaries(0);

        return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

    ////////////////////////////////////////////////////////////////

    aParticleChange.SetNumberOfSecondaries(NumPhotons);

    if (fTrackSecondariesFirst) {
        if (aTrack.GetTrackStatus() == fAlive )
            aParticleChange.ProposeTrackStatus(fSuspend);
    }

    ////////////////////////////////////////////////////////////////

    G4int materialIndex = aMaterial->GetIndex();

    // Retrieve the Scintillation Integral for this material
    // new G4PhysicsOrderedFreeVector allocated to hold CII's



    std::vector<G4double> lFractions;
    std::vector<G4double> lTimes;

    G4double lFractionNormalization = 0;
    for (G4int scnt = 0; scnt < nscnt; scnt++) {
        G4double lCurrentFraction = lScintComponents->Energy(scnt); // don't get confused about the energy here, it is a fraction... I wanted to use the G4 vectors, and this is the nomenclature
        G4double lCurrentTime = lScintComponents->Value(lCurrentFraction);
        lFractionNormalization += lCurrentFraction;
        lTimes.push_back(lCurrentTime);
        lFractions.push_back(lCurrentFraction);
    }

    std::vector<G4int> Nums;
    G4int lTotalNr = 0;
        for (G4int scnt = 0; scnt < nscnt; scnt++) {
           lFractions[scnt] /= lFractionNormalization;
           Nums.push_back(G4int (round(lFractions[scnt]*NumPhotons)));
           lTotalNr+=Nums[scnt];
        }

    if (lTotalNr != NumPhotons){

        if (NumPhotons == 1){
            G4double random_nr = G4UniformRand();
            G4double lTempSum = 0;
            for (G4int scnt = 0; scnt < nscnt; scnt++) {
                lTempSum += lFractions[scnt];
                if (random_nr < lTempSum) {
                    Nums[scnt] = 1;
                    break;
                }
            }
        }

        else {
            G4double difference = lTotalNr - NumPhotons;
            std::vector<G4double> diff;
            for (G4int scnt = 0; scnt < nscnt; scnt++) {
                diff.push_back(Nums[scnt]-lFractions[scnt]*NumPhotons);
            }
            for (G4int i = 1; i <= abs(difference); i++){
                if (difference < 0){ // if i have rounded less photons than needed
                    G4int myindex = G4int(find(diff.begin(), diff.end(),  *std::min_element(diff.begin(),diff.end()))- diff.begin());
                    Nums[myindex] += 1;
                    diff[myindex] = 0;
                }
                else if (difference>0){// if i have rounded more photons than needed
                G4int myindex = G4int(find(diff.begin(), diff.end(),  *std::max_element(diff.begin(),diff.end()))- diff.begin());
                Nums[myindex] -= 1;
                diff[myindex] = 0;
                }
            }
        }
    }



    for (G4int scnt = 1; scnt <= nscnt; scnt++) {


        G4double ScintillationRiseTime = 0.*ns;

        G4int Num = Nums[scnt-1];

        G4double ScintillationTime   =  lTimes[scnt-1];
        G4PhysicsOrderedFreeVector* ScintillationIntegral = (G4PhysicsOrderedFreeVector*)((*firstIntegralTable)(materialIndex));

        if (!ScintillationIntegral) continue;

        // Max Scintillation Integral

        G4double CIImax = ScintillationIntegral->GetMaxValue();

        for (G4int i = 0; i < Num; i++) {

            // Determine photon energy

            G4double CIIvalue = G4UniformRand()*CIImax;
            G4double sampledEnergy =
            ScintillationIntegral->GetEnergy(CIIvalue);

            if (verboseLevel>1) {
                G4cout << "sampledEnergy = " << sampledEnergy << G4endl;
                G4cout << "CIIvalue =        " << CIIvalue << G4endl;
            }

            // Generate random photon direction

            G4double cost = 1. - 2.*G4UniformRand();
            G4double sint = std::sqrt((1.-cost)*(1.+cost));

            G4double phi = twopi*G4UniformRand();
            G4double sinp = std::sin(phi);
            G4double cosp = std::cos(phi);

            G4double px = sint*cosp;
            G4double py = sint*sinp;
            G4double pz = cost;

            // Create photon momentum direction vector

            G4ParticleMomentum photonMomentum(px, py, pz);

            // Determine polarization of new photon

            G4double sx = cost*cosp;
            G4double sy = cost*sinp;
            G4double sz = -sint;

            G4ThreeVector photonPolarization(sx, sy, sz);

            G4ThreeVector perp = photonMomentum.cross(photonPolarization);
            //G4cout << photonMomentum.x() << " " << photonMomentum.y() << " " << photonMomentum.z() << G4endl;
            phi = twopi*G4UniformRand();
            sinp = std::sin(phi);
            cosp = std::cos(phi);

            photonPolarization = cosp * photonPolarization + sinp * perp;

            photonPolarization = photonPolarization.unit();

            // Generate a new photon:

            G4DynamicParticle* aScintillationPhoton =
            new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),
                                  photonMomentum);
            aScintillationPhoton->SetPolarization
            (photonPolarization.x(),
             photonPolarization.y(),
             photonPolarization.z());

            aScintillationPhoton->SetKineticEnergy(sampledEnergy);

            // Generate new G4Track object:

            G4double rand;

            if (aParticle->GetDefinition()->GetPDGCharge() != 0) {
                rand = G4UniformRand();
            } else {
                rand = 1.0;
            }

            G4double delta = rand * aStep.GetStepLength();
            G4double deltaTime = delta /
            ((pPreStepPoint->GetVelocity()+
            pPostStepPoint->GetVelocity())/2.);

            // emission time distribution
            if (ScintillationRiseTime==0.0) {
                deltaTime = deltaTime -
                ScintillationTime * std::log( G4UniformRand() );
            } else {
                deltaTime = deltaTime +
                sample_time(ScintillationRiseTime, ScintillationTime);
            }

            G4double aSecondaryTime = t0 + deltaTime;

            //G4cout << deltaTime << " " << t0 << " " <<aSecondaryTime<<  G4endl;

            G4ThreeVector aSecondaryPosition =
            x0 + rand * aStep.GetDeltaPosition();

            G4Track* aSecondaryTrack =
            new G4Track(aScintillationPhoton,aSecondaryTime,aSecondaryPosition);

            aSecondaryTrack->SetTouchableHandle(
                aStep.GetPreStepPoint()->GetTouchableHandle());
            // aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);

            aSecondaryTrack->SetParentID(aTrack.GetTrackID());

            aParticleChange.AddSecondary(aSecondaryTrack);

        }
    }

    if (verboseLevel>0) {
        G4cout << "\n Exiting from mdomScintillation::DoIt -- NumberOfSecondaries = "
        << aParticleChange.GetNumberOfSecondaries() << G4endl;
    }

    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// BuildThePhysicsTable for the scintillation process
// --------------------------------------------------
//

void OMSimScintillation::BuildThePhysicsTable()
{       if (firstIntegralTable && secondIntegralTable && thirdIntegralTable ) return;
    if (theFastIntegralTable && theSlowIntegralTable) return;

    const G4MaterialTable* theMaterialTable =
    G4Material::GetMaterialTable();
    G4int numOfMaterials = G4Material::GetNumberOfMaterials();

    // create new physics table

    if(!theFastIntegralTable)theFastIntegralTable = new G4PhysicsTable(numOfMaterials);
    if(!theSlowIntegralTable)theSlowIntegralTable = new G4PhysicsTable(numOfMaterials);



    if(!firstIntegralTable)firstIntegralTable = new G4PhysicsTable(numOfMaterials);
    if(!secondIntegralTable)secondIntegralTable = new G4PhysicsTable(numOfMaterials);
    if(!thirdIntegralTable)thirdIntegralTable = new G4PhysicsTable(numOfMaterials);


    // loop for materials

    for (G4int i=0 ; i < numOfMaterials; i++)
    {
        G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
        new G4PhysicsOrderedFreeVector();
        G4PhysicsOrderedFreeVector* bPhysicsOrderedFreeVector =
        new G4PhysicsOrderedFreeVector();
        G4PhysicsOrderedFreeVector* cPhysicsOrderedFreeVector =
        new G4PhysicsOrderedFreeVector();
        G4PhysicsOrderedFreeVector* dPhysicsOrderedFreeVector =
        new G4PhysicsOrderedFreeVector();
        G4PhysicsOrderedFreeVector* ePhysicsOrderedFreeVector =
        new G4PhysicsOrderedFreeVector();

        // Retrieve vector of scintillation wavelength intensity for
        // the material from the material's optical properties table.

        G4Material* aMaterial = (*theMaterialTable)[i];

        G4MaterialPropertiesTable* aMaterialPropertiesTable =
        aMaterial->GetMaterialPropertiesTable();

        if (aMaterialPropertiesTable) {

            G4MaterialPropertyVector* theFastLightVector =
            aMaterialPropertiesTable->GetProperty("FASTCOMPONENT");
            if (theFastLightVector) {

                // Retrieve the first intensity point in vector
                // of (photon energy, intensity) pairs

                G4double currentIN = (*theFastLightVector)[0];

                if (currentIN >= 0.0) {

                    // Create first (photon energy, Scintillation
                    // Integral pair

                    G4double currentPM = theFastLightVector->Energy(0);

                    G4double currentCII = 0.0;

                    aPhysicsOrderedFreeVector->
                    InsertValues(currentPM , currentCII);

                    // Set previous values to current ones prior to loop

                    G4double prevPM  = currentPM;
                    G4double prevCII = currentCII;
                    G4double prevIN  = currentIN;

                    // loop over all (photon energy, intensity)
                    // pairs stored for this material

                    for (size_t ii = 1;
                         ii < theFastLightVector->GetVectorLength();
                    ++ii)
                         {
                             currentPM = theFastLightVector->Energy(ii);
                             currentIN = (*theFastLightVector)[ii];

                             currentCII = 0.5 * (prevIN + currentIN);

                             currentCII = prevCII +
                             (currentPM - prevPM) * currentCII;

                             aPhysicsOrderedFreeVector->
                             InsertValues(currentPM, currentCII);

                             prevPM  = currentPM;
                             prevCII = currentCII;
                             prevIN  = currentIN;
                         }

                }
            }

            G4MaterialPropertyVector* theSlowLightVector =
            aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");

            if (theSlowLightVector) {

                // Retrieve the first intensity point in vector
                // of (photon energy, intensity) pairs

                G4double currentIN = (*theSlowLightVector)[0];

                if (currentIN >= 0.0) {

                    // Create first (photon energy, Scintillation
                    // Integral pair

                    G4double currentPM = theSlowLightVector->Energy(0);

                    G4double currentCII = 0.0;

                    bPhysicsOrderedFreeVector->
                    InsertValues(currentPM , currentCII);

                    // Set previous values to current ones prior to loop

                    G4double prevPM  = currentPM;
                    G4double prevCII = currentCII;
                    G4double prevIN  = currentIN;

                    // loop over all (photon energy, intensity)
                    // pairs stored for this material

                    for (size_t ii = 1;
                         ii < theSlowLightVector->GetVectorLength();
                    ++ii)
                         {
                             currentPM = theSlowLightVector->Energy(ii);
                             currentIN = (*theSlowLightVector)[ii];

                             currentCII = 0.5 * (prevIN + currentIN);

                             currentCII = prevCII +
                             (currentPM - prevPM) * currentCII;

                             bPhysicsOrderedFreeVector->
                             InsertValues(currentPM, currentCII);

                             prevPM  = currentPM;
                             prevCII = currentCII;
                             prevIN  = currentIN;
                         }

                }
            }



            G4MaterialPropertyVector* firstLightVector =
            aMaterialPropertiesTable->GetProperty("FIRSTCOMPONENT");


            if (firstLightVector) {

                // Retrieve the first intensity point in vector
                // of (photon energy, intensity) pairs

                G4double currentIN = (*firstLightVector)[0];

                if (currentIN >= 0.0) {

                    // Create first (photon energy, Scintillation
                    // Integral pair

                    G4double currentPM = firstLightVector->Energy(0);

                    G4double currentCII = 0.0;

                    cPhysicsOrderedFreeVector->
                    InsertValues(currentPM , currentCII);

                    // Set previous values to current ones prior to loop

                    G4double prevPM  = currentPM;
                    G4double prevCII = currentCII;
                    G4double prevIN  = currentIN;

                    // loop over all (photon energy, intensity)
                    // pairs stored for this material

                    for (size_t ii = 1;
                         ii < firstLightVector->GetVectorLength();
                    ++ii)
                         {
                             currentPM = firstLightVector->Energy(ii);
                             currentIN = (*firstLightVector)[ii];

                             currentCII = 0.5 * (prevIN + currentIN);

                             currentCII = prevCII +
                             (currentPM - prevPM) * currentCII;

                             cPhysicsOrderedFreeVector->
                             InsertValues(currentPM, currentCII);

                             prevPM  = currentPM;
                             prevCII = currentCII;
                             prevIN  = currentIN;
                         }

                }
            }






            G4MaterialPropertyVector* secondLightVector =
            aMaterialPropertiesTable->GetProperty("SECONDCOMPONENT");


            if (secondLightVector) {

                // Retrieve the first intensity point in vector
                // of (photon energy, intensity) pairs

                G4double currentIN = (*secondLightVector)[0];

                if (currentIN >= 0.0) {

                    // Create first (photon energy, Scintillation
                    // Integral pair

                    G4double currentPM = secondLightVector->Energy(0);

                    G4double currentCII = 0.0;

                    dPhysicsOrderedFreeVector->
                    InsertValues(currentPM , currentCII);

                    // Set previous values to current ones prior to loop

                    G4double prevPM  = currentPM;
                    G4double prevCII = currentCII;
                    G4double prevIN  = currentIN;

                    // loop over all (photon energy, intensity)
                    // pairs stored for this material

                    for (size_t ii = 1;
                         ii < secondLightVector->GetVectorLength();
                    ++ii)
                         {
                             currentPM = secondLightVector->Energy(ii);
                             currentIN = (*secondLightVector)[ii];

                             currentCII = 0.5 * (prevIN + currentIN);

                             currentCII = prevCII +
                             (currentPM - prevPM) * currentCII;

                             dPhysicsOrderedFreeVector->
                             InsertValues(currentPM, currentCII);

                             prevPM  = currentPM;
                             prevCII = currentCII;
                             prevIN  = currentIN;
                         }

                }
            }

            G4MaterialPropertyVector* thirdLightVector =
            aMaterialPropertiesTable->GetProperty("THIRDCOMPONENT");


            if (thirdLightVector) {

                // Retrieve the first intensity point in vector
                // of (photon energy, intensity) pairs

                G4double currentIN = (*thirdLightVector)[0];

                if (currentIN >= 0.0) {

                    // Create first (photon energy, Scintillation
                    // Integral pair

                    G4double currentPM = thirdLightVector->Energy(0);

                    G4double currentCII = 0.0;

                    ePhysicsOrderedFreeVector->
                    InsertValues(currentPM , currentCII);

                    // Set previous values to current ones prior to loop

                    G4double prevPM  = currentPM;
                    G4double prevCII = currentCII;
                    G4double prevIN  = currentIN;

                    // loop over all (photon energy, intensity)
                    // pairs stored for this material

                    for (size_t ii = 1;
                         ii < thirdLightVector->GetVectorLength();
                    ++ii)
                         {
                             currentPM = thirdLightVector->Energy(ii);
                             currentIN = (*thirdLightVector)[ii];

                             currentCII = 0.5 * (prevIN + currentIN);

                             currentCII = prevCII +
                             (currentPM - prevPM) * currentCII;

                             ePhysicsOrderedFreeVector->
                             InsertValues(currentPM, currentCII);

                             prevPM  = currentPM;
                             prevCII = currentCII;
                             prevIN  = currentIN;
                         }

                }
            }


        }


        // The scintillation integral(s) for a given material
        // will be inserted in the table(s) according to the
        // position of the material in the material table.
        theFastIntegralTable->insertAt(i,aPhysicsOrderedFreeVector);
        theSlowIntegralTable->insertAt(i,bPhysicsOrderedFreeVector);
        firstIntegralTable->insertAt(i,cPhysicsOrderedFreeVector);
        secondIntegralTable->insertAt(i,dPhysicsOrderedFreeVector);
        thirdIntegralTable->insertAt(i,ePhysicsOrderedFreeVector);
    }





}

// Called by the user to set the scintillation yield as a function
// of energy deposited by particle type

void OMSimScintillation::SetScintillationByParticleType(const G4bool scintType)
{
    if (emSaturation) {
        G4Exception("mdomScintillation::SetScintillationByParticleType", "Scint02",
                    JustWarning, "Redefinition: Birks Saturation is replaced by ScintillationByParticleType!");
        RemoveSaturation();
    }
    scintillationByParticleType = scintType;
}

// GetMeanFreePath
// ---------------
//

G4double OMSimScintillation::GetMeanFreePath(const G4Track&,
                                            G4double ,
                                            G4ForceCondition* condition)
{
    *condition = StronglyForced;

    return DBL_MAX;

}

// GetMeanLifeTime
// ---------------
//

G4double OMSimScintillation::GetMeanLifeTime(const G4Track&,
                                            G4ForceCondition* condition)
{
    *condition = Forced;

    return DBL_MAX;

}

G4double OMSimScintillation::sample_time(G4double tau1, G4double tau2)
{
    // tau1: rise time and tau2: decay time

    while(1) {
        // two random numbers
        G4double ran1 = G4UniformRand();
        G4double ran2 = G4UniformRand();
        //
        // exponential distribution as envelope function: very efficient
        //
        G4double d = (tau1+tau2)/tau2;
        // make sure the envelope function is
        // always larger than the bi-exponential
        G4double t = -1.0*tau2*std::log(1-ran1);
        G4double gg = d*single_exp(t,tau2);
        if (ran2 <= bi_exp(t,tau1,tau2)/gg) return t;
    }
    return -1.0;
}
