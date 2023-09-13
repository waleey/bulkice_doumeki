/** @file OMSimDetectorConstruction.cc
 *  @brief User defined detector.
 *
 * You should define and construct the detector here...this template is an example for a single arg module.
 *
 *  @author Martin Unland - modified by Markus Dittmer
 *  @date October 2021 - modified at February 2022
 *  @Modified by Waly M Z Karim
 *  @version Geant4 10.7
 *
 *  @todo
 */

#include "OMSimDetectorConstruction.hh"

#include "OMSimRadioactivityData.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"




extern G4double gworldsize;
extern G4int gDOM;
extern G4bool gCADImport;
extern G4bool gPlaceHarness;

G4double OMSimRadioactivityData::fglassInRad = 0.0;
G4double OMSimRadioactivityData::fglassOutRad = 0.0;
G4int OMSimRadioactivityData::fomModel = 0;

OMSimDetectorConstruction::OMSimDetectorConstruction(G4int DOM, G4double worldSize)
    : mWorldSolid(0), mWorldLogical(0), mWorldPhysical(0), fDOM(DOM), fworldSize(worldSize)
{
    radData = new OMSimRadioactivityData();
}
OMSimDetectorConstruction::OMSimDetectorConstruction() : mWorldSolid(0), mWorldLogical(0), mWorldPhysical(0), fDOM(0), fworldSize(0)
{

}

OMSimDetectorConstruction::~OMSimDetectorConstruction()
{
    delete mWorldSolid;
    delete mWorldLogical;
    delete mWorldPhysical;
    delete mData;
    delete mPMTManager;
    delete fOuterSolid;
    delete fInnerSolid;
    delete fMDOM;
    delete fLOM16;
    delete fLOM18;
    delete fDEGG;
    delete fPDOM;
    delete fWOM;
    delete radData;
}

/**
 * Construct the world volume
 */

void OMSimDetectorConstruction::ConstructWorld()
{
    //ConstructWorldMat();

    mWorldSolid = new G4Box("World", fworldSize * m, fworldSize * m, fworldSize* m);
    //mWorldSolid = new G4Box("World", .5 * m, .5 * m, .5 * m);
    //mWorldLogical = new G4LogicalVolume(mWorldSolid, ice, "World_log", 0, 0, 0);
    //mWorldPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), mWorldLogical, "World_phys", 0, false, 0, true);
    mWorldLogical = new G4LogicalVolume(mWorldSolid, mData->GetMaterial("argWorld"), "World_log", 0, 0, 0);
    mWorldPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), mWorldLogical, "World_phys", 0, false, 0);
    G4VisAttributes* World_vis = new G4VisAttributes(G4Colour(0.45,0.5,0.35,0.2));
    //G4VisAttributes* World_vis = new G4VisAttributes(G4Colour(0,0,0,0));
    mWorldLogical -> SetVisAttributes(World_vis);

    G4cout << "::::::::::::::::::World Volume Constructed Sucessfully:::::::::::::::::::" << G4endl;

}

// DONE UP TO THIS POINT

/**
 * Construct all solids needed for your study. Call OMSimOMConstruction for simulations with optical modules
 * and OMSimPMTConstruction for simulations with only PMTs.
 * @return World physical for the main
 */
G4VPhysicalVolume *OMSimDetectorConstruction::Construct()
{

    OMSimRadioactivityData::fomModel = fDOM;
    mData = new OMSimInputData();
    //mData->SearchFolders("/home/waly/bulkice_doumeki/mdom/build/"); //Will change soon
    mData->SearchFolders("../build/"); // you need to change this path if data is saved somewhere else.

    ConstructWorld();

   if (fDOM == 0){ //Single PMT
        G4cout << "Constructing single PMT" << G4endl;
        mPMTManager = new OMSimPMTConstruction(mData);
        mPMTManager->SelectPMT("argPMT");

        G4RotationMatrix* lRot = new G4RotationMatrix();
        G4Transform3D lTransformers = G4Transform3D(*lRot, G4ThreeVector(0,0,0));
        mPMTManager->PlaceIt(lTransformers, mWorldLogical, "PMT");
    }
    else if (fDOM == 1){ //mDOM
        G4cout << "Constructing mDOM" << G4endl;
        fMDOM = new mDOM(mData,gPlaceHarness);
        glassOutRad = fMDOM -> getGlassOutRad();
        glassInRad = fMDOM -> getGlassInRad();
        OMSimRadioactivityData::SetGlassRad(glassOutRad, glassInRad);
        fMDOM->PlaceIt(G4ThreeVector(0, 0, 0), G4RotationMatrix(), mWorldLogical, "");
    }
    else if (fDOM == 2){ //PDOM
        G4cout << "Constructing PDOM" << G4endl;
        fPDOM = new pDOM(mData,gPlaceHarness);
        fPDOM->PlaceIt(G4ThreeVector(0, 0, 0), G4RotationMatrix(), mWorldLogical, "");
    }
    else if (fDOM == 3){ //LOM16
        G4cout << "Constructing LOM16" << G4endl;
        fLOM16 = new LOM16(mData,gPlaceHarness);
        fLOM16->PlaceIt(G4ThreeVector(0, 0, 0), G4RotationMatrix(), mWorldLogical, "");
    }
    else if (fDOM == 4){ //LOM18
        G4cout << "Constructing LOM18" << G4endl;
        fLOM18 = new LOM18(mData,gPlaceHarness);
        fLOM18->PlaceIt(G4ThreeVector(0, 0, 0), G4RotationMatrix(), mWorldLogical, "");
        G4cout << "::::::::::::::LOM18 successfully constructed::::::::::::" << G4endl;
    }
     else if (fDOM == 5){ //DEGG
        G4cout << "Constructing DEGG" << G4endl;
        fDEGG = new dEGG(mData,gPlaceHarness);
        fDEGG->PlaceIt(G4ThreeVector(0, 0, 0), G4RotationMatrix(), mWorldLogical, "");
        glassOutRad = 350;
        glassInRad = 50;
     }
     else if(fDOM == 6)
     {
        G4cout << "Constructing WOM" << G4endl;
        fWOM = new WOM(mWorldLogical, mData);
        fWOM -> PlaceIt();
     }
    else{ //Add your costume detector contruction here and call it with -m 6 (or greater)
        G4cout << "Constructing custome detector construction" << G4endl;
    }






    return mWorldPhysical;
}
G4ThreeVector OMSimDetectorConstruction::DrawFromVolume()
{
    switch(fDOM)
    {

        case 0:
            std::cerr << "Glass Solid Not Available for PMTs only! Aborting..." << std::endl;
            exit(0);
            break;
        case 1:
            fOuterSolid = fMDOM -> GetOuterSolid();
            fInnerSolid = fMDOM -> GetInnerSolid();
            break;
        case 2:
            fOuterSolid = fPDOM -> GetOuterSolid();
            fInnerSolid = fPDOM -> GetInnerSolid();
            break;
        case 3:
            fOuterSolid = fLOM16 -> GetOuterSolid();
            fInnerSolid = fLOM16 -> GetInnerSolid();
            break;
        case 4:
            fOuterSolid = fLOM18 -> GetOuterSolid();
            fInnerSolid = fLOM18 -> GetInnerSolid();
            break;
        case 5:
            fOuterSolid = fDEGG -> GetOuterSolid();
            fInnerSolid = fDEGG -> GetInnerSolid();
            break;
        default:
            std::cerr << "Invalid OM Model. Aborting..." << std::endl;
            exit(0);
    }

    G4double theta = radData -> RandomGen(0, M_PI);
    //G4double theta = M_PI / 2;
    G4double phi = radData -> RandomGen(0, M_PI * 2);

    G4int maxTries = 1000000;
    G4int iTry = 1;

    G4ThreeVector point;

    do
    {
        point.set(
            (radData -> RandomGen(glassInRad - 50, glassOutRad + 50) * sin(theta) * cos(phi)) * mm,
            (radData -> RandomGen(glassInRad - 50, glassOutRad + 50) * sin(theta) * sin(phi)) * mm,
            (radData -> RandomGen(glassInRad - 50, glassOutRad + 50) * cos(theta)) * mm
        );
    }while(!(fOuterSolid -> Inside(point) && !fInnerSolid -> Inside(point)) && ++iTry < maxTries);

    if(iTry == maxTries)
    {
        std::cout << "Unable to find a point inside the Volume! Check the volume params correctly! Aborting..." << std::endl;
        exit(0);
        point = G4ThreeVector(0,0,0);
        return point; //should not reach here.
    }
    else
    {
            std::cout << "x: " << point.x() / mm << std::endl
            << "y: " << point.y() / mm << std::endl
            << "z: " << point.z() / mm << std::endl;
        return point;
    }


}

