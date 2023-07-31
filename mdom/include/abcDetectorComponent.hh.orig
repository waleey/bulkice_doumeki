
#ifndef abcDetectorComponent_h
#define abcDetectorComponent_h 1

#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"
#include "OMSimInputData.hh"
#include "OMSimPMTConstruction.hh"
#include "G4MultiUnion.hh"
#include <sstream>

class abcDetectorComponent
{
public:
    abcDetectorComponent(){};
    virtual void Construction() = 0; // Abstract method you have to define in order to make a derived class from abcDetectorComponent

    OMSimInputData* mData; // Instance of OMSimInputdata, which should be started only once.
    bool mCheckOverlaps = true;
    
    struct Component{     // Struct of variables of component of a DetectorComponent
        G4VSolid* VSolid; //Solid Volume of component
        G4LogicalVolume* VLogical; //Logical Volume of component
        G4ThreeVector Position; //Position of component wrt to 0
        G4RotationMatrix Rotation; //Rotation of component wrt to 0
        G4String Name; // Name of Component
    };

    std::vector<G4ThreeVector> mPlacedPositions; //store the positions each time the components are placed
    std::vector<G4RotationMatrix> mPlacedOrientations; //store the orientations each time the components are placed
    std::vector<Component> Components;

    /**Methods in abcDetectorComponent.cc*/
    void AppendComponent(G4VSolid* pSolid, G4LogicalVolume* pLogical, G4ThreeVector pVector, G4RotationMatrix pRotation, G4String pName);
    Component GetComponent(G4String pName);
    G4Transform3D GetNewPosition(G4ThreeVector pPosition, G4RotationMatrix pRotation, G4ThreeVector pObjectPosition, G4RotationMatrix pObjectRotation);
    void IntegrateDetectorComponent(abcDetectorComponent* pToIntegrate, G4ThreeVector pPosition, G4RotationMatrix pRotation, G4String pNameExtension);
    void PlaceIt(G4ThreeVector pPosition, G4RotationMatrix pRotation, G4LogicalVolume*& pMother, G4String pNameExtension = "");
    G4SubtractionSolid* SubstractToVolume(G4VSolid* pInputVolume, G4ThreeVector pSubstractionPos, G4RotationMatrix pSubstractionRot, G4String pNewVolumeName);
    
protected:
    
    const G4VisAttributes* mGlassVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.8, 0.2));
    const G4VisAttributes* mGelVis = new G4VisAttributes(G4Colour(0.45, 0.5, 0.35, 0.2));
    const G4VisAttributes* mAluVis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.9, 1.0));
    const G4VisAttributes* mSteelVis = new G4VisAttributes(G4Colour(0.75, 0.75, 0.85, 1.0));
    const G4VisAttributes* mAbsorberVis = new G4VisAttributes(G4Colour(0.2, 0.2, 0.2, 1.0));
    const G4VisAttributes* mBoardVis = new G4VisAttributes(G4Colour(0, 1, 0, 1));
    const G4VisAttributes* mAirVis = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 1.0));
    const G4VisAttributes mInvisibleVis = G4VisAttributes::GetInvisible();
    const G4VisAttributes* mRedVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 1.0));
    const G4VisAttributes* mBlackVis = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0, 1.0));
    const G4VisAttributes* mLEDvis= new G4VisAttributes(G4Colour(0.2,0.6,0.8,0.5));	    

};


#endif

