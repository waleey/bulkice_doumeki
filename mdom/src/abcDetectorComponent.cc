

#include "abcDetectorComponent.hh"
#include "OMSimPMTConstruction.hh"
#include "OMSimLogger.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"

/**
 * Append one component to Components vector. 
 * @param pSolid Solid of component
 * @param pLogical Logical of component
 * @param pVector G4ThreeVector with position of component wrt 0
 * @param pRotation Rotation of component wrt 0
 * @param pName Name of component with which you can use the get method @see GetComponent()
 */
void abcDetectorComponent::AppendComponent(G4VSolid* pSolid, G4LogicalVolume* pLogical, G4ThreeVector pVector, G4RotationMatrix pRotation, G4String pName) {
    Components.push_back(new Component(
        pSolid,
        pLogical,
        pVector,
        pRotation,
        pName));
}

/**
 * Get a component from Components vector. 
 * @param pName Name of component which you wanna have
 */
Component abcDetectorComponent::GetComponent(G4String pName) {
    G4Transform3D lTrans;
    for (auto Component : Components) {
        if (Component->Name == pName) return *Component;
    }
    G4String mssg = pName+" not found in component list. This will probably throw a segmentation fault...";
    critical(mssg);
}


G4Transform3D abcDetectorComponent::GetNewPosition(G4ThreeVector pPosition, G4RotationMatrix pRotation, G4ThreeVector pObjectPosition, G4RotationMatrix pObjectRotation)
{
    return G4Transform3D(pObjectRotation.transform(pRotation), pPosition + pObjectPosition.transform(pRotation));
}


/**
 * Placement of the DetectorComponent. Each Component is placed in the same mother.
 * @param pPosition G4ThreeVector with position of the module (as in G4PVPlacement())
 * @param pRotation G4RotationMatrix with rotation of the module (as in G4PVPlacement())
 * @param pMother G4LogicalVolume where the module is going to be placed (as in G4PVPlacement())
 * @param pIncludeHarness bool Harness is placed if true
 * @param pNameExtension G4String name of the physical volume. You should not have two physicals with the same name
 */
void abcDetectorComponent::PlaceIt(G4ThreeVector pPosition, G4RotationMatrix pRotation, G4LogicalVolume*& pMother, G4String pNameExtension)
{
    mPlacedPositions.push_back(pPosition);
    mPlacedOrientations.push_back(pRotation);
    G4Transform3D lTrans;
    for (auto Component : Components) {
        lTrans = GetNewPosition(pPosition, pRotation, Component->Position, Component->Rotation);
        new G4PVPlacement(lTrans, Component->VLogical, Component->Name + pNameExtension, pMother, false, 0, mCheckOverlaps);
    }

}


/**
 * Integrate components of another abcDetectorComponent instance. You can translate/rotate before integrating for coordinate system matching.
 * @param pToIntegrate abcDetectorComponent instance whose components we want to integrate
 * @param pPosition G4ThreeVector with position of new detector component where it will be integrated
 * @param pRotation G4RotationMatrix with rotation of new detector component
 * @param pNameExtension G4String to extend original name of component
 */
void abcDetectorComponent::IntegrateDetectorComponent(abcDetectorComponent* pToIntegrate, G4ThreeVector pPosition, G4RotationMatrix pRotation, G4String pNameExtension)
{
    G4Transform3D lTrans;
    for (auto Component : pToIntegrate->Components) {
        lTrans = GetNewPosition(pPosition, pRotation, Component->Position, Component->Rotation);
        AppendComponent(Component->VSolid, Component->VLogical, lTrans.getTranslation(), lTrans.getRotation(), Component->Name+pNameExtension);
    }

}


G4SubtractionSolid* abcDetectorComponent::SubstractToVolume(G4VSolid* pInputVolume, G4ThreeVector pSubstractionPos, G4RotationMatrix pSubstractionRot, G4String pNewVolumeName)
{
    G4SubtractionSolid* lSubstractedVolume;
    G4Transform3D lTrans;
    G4int iCounter = 0;
    for (auto Component : Components) {
        lTrans = GetNewPosition(pSubstractionPos, pSubstractionRot, Component->Position, Component->Rotation);
        if (iCounter == 0) {
            lSubstractedVolume = new G4SubtractionSolid("SubstractedVolume", pInputVolume, Component->VSolid, lTrans);
        } else {
            lSubstractedVolume = new G4SubtractionSolid("SubstractedVolume", lSubstractedVolume, Component->VSolid, lTrans);
        }
        iCounter++;
    }
    return lSubstractedVolume;
}

