#include "WOM.hh"
/**
*Very rough draft of possible WOM Simulation
*Waly M Z Karim
*8/8/2023
**/
#include "WOMMaterial.hh"
#include "G4UserLimits.hh"

WOM::WOM(G4LogicalVolume* LogicMother, OMSimInputData* data) : fLogicMother(LogicMother), fData(data)
{
    std::cout << "******Constructing WOMs******" << std::endl;
    GetSharedData();

}
WOM::~WOM()
{
    //Deleting the solids
    delete fGlassSolid;
    delete fGelSolid;
    delete fWomTubeSolid;
    delete fWomInsideSolid;
    delete fPMTSolid1;

    //deleting the logicals
    delete fGlassLogical;
    delete fGelLogical;
    delete fWomTubeLogical;
    delete fWomTubeInsideLogical;
    delete fPMTBodyLogical1;
    delete fPMTBodyLogical2;
    delete fPMTCathodeLogical1;
    delete fPMTCathodeLogical2;
    delete fPMTCathode;
}
void WOM::PlaceIt()
{
    Construction();
}
void WOM::GetSharedData()
{
    fGlassTubeOuterRad = (173 / 2) *mm;
    fGlassTubeInnerRad = (145 / 2) *mm;
    fGlassCapOuterRad = (173 / 2) *mm;
    fGlassCapInnerRad = (145 / 2) *mm;
    fGlassTubeHalfLength = 550 * mm;

    fWOMPaintOuterRad = (117000 / 2) * micrometer;
    fWOMPaintThickness = 25;
    fWOMPaintInnerRad = (117000 / 2 - fWOMPaintThickness) *micrometer; //the paint is assumed to be 25 micrometer thick.

    fWomTubeOuterRad = (117000 / 2 - fWOMPaintThickness) *micrometer;
    fWomTubeInnerRad = (115 / 2) *mm;
    fWomTubeHalfLength = (760 / 2) *mm;

    fPMTHolderRad = (135 / 2) *mm;
    fPMTHolderZ = 19 *mm;
    fPMTHolderHalfLength = (38 / 2) *mm;

    fPMTConeHeight = 56.95509 *mm;
    fPMTConeCut = 25.5 *mm;
    fPMTConeRatio = 0.81863 * mm;
    fPMTConeZ = 63.5 *mm;

    fPMTTubeRad = (51.5 / 2) * mm;
    fPMTTubeZ = 122.0 * mm;
    fPMTTubeHalfLength = (66 / 2) *mm;

    fPMTBoardRad = (51.5 / 2) *mm;
    fPMTBoardZ = 171 * mm;
    fPMTBoardHalfLength = (32 / 2) *mm;

    fPMTCathodeRad = (135 / 2) *mm;
    fPMTCathodeHalfLength = 1 *mm;
    fPMTCathodeZ = fWomTubeHalfLength + 1 * mm;

    fPMTOffset = 2 *mm;

    fPMTGlobalZ = fWomTubeHalfLength;
;}
void WOM::GenerateLogicals()
{
    fGlassSolid = PressureVessel("glass", fGlassTubeOuterRad, fGlassCapOuterRad);
    fGelSolid = PressureVessel("Gel", fGlassTubeInnerRad, fGlassCapInnerRad); //Gel will contain the optical filling materials
    fWomTubeSolid = WOMTube("WOM", fWomTubeOuterRad);
    fWomInsideSolid = WOMTube("WOMInside", fWomTubeInnerRad);
    fPMTSolid1 = PMTConstruction();
    fPMTCathode = PMTCathodeConstruction();
    fWOMPaintSolid = WOMPaint("WOM");

    ConstructMaterial();

    fGlassLogical = new G4LogicalVolume(fGlassSolid, glassMaterial, "glassLogical");
    fGelLogical = new G4LogicalVolume(fGelSolid, fillerMaterial, "gelLogical");
    fWomTubeLogical = new G4LogicalVolume(fWomTubeSolid, tubeMaterial, "WOOMTubeLogical");
    fWomTubeInsideLogical = new G4LogicalVolume(fWomInsideSolid, tubeInsideMaterial, "WOMInsideLogical");
    fPMTBodyLogical1 = new G4LogicalVolume(fPMTSolid1, pmtAbsorber, "PMTBodyLogical1");
    fPMTBodyLogical2 = new G4LogicalVolume(fPMTSolid1, pmtAbsorber, "PMTBodyLogical2");
    fPMTCathodeLogical1 = new G4LogicalVolume(fPMTCathode, pmtPhotoCathode, "PMT_Cathode_Logical1");
    fPMTCathodeLogical2 = new G4LogicalVolume(fPMTCathode, pmtPhotoCathode, "PMT_Cathode_Logical2");
    fWOMPaintLogical = new G4LogicalVolume(fWOMPaintSolid, paintMaterial, "WOM_Paint_Logical");
    //fWOMPaintLogical = new G4LogicalVolume(fWOMPaintSolid, pmtPhotoCathode, "WOM_Paint_Logical"); //delete this later

    /*std::cout << "Reading out the volume materials!!" << std::endl
    << "Glass Logical: " << fGlassLogical -> GetMaterial() -> GetName() << std::endl
    << "Gel Logical: " << fGelLogical -> GetMaterial() -> GetName() << std::endl
    << "WOM Tube Logical: " << fWomTubeLogical -> GetMaterial() -> GetName() << std::endl
    << "WOM Tube Inside Logical: " << fWomTubeInsideLogical -> GetMaterial() -> GetName() << std::endl
    << "Paint Logical: " << fWOMPaintLogical -> GetMaterial() -> GetName() << std::endl;

    exit(0);*/
    /**
    *Adding a user limit for shorter absorption lengths
    **/

    G4double maxStepSize = 0.1 * micrometer;
    G4UserLimits* maxLimit = new G4UserLimits(maxStepSize);
    fWOMPaintLogical -> SetUserLimits(maxLimit);


   fGlassLogical -> SetVisAttributes(new G4VisAttributes(G4Colour::Blue()));
   fGelLogical -> SetVisAttributes(new G4VisAttributes(G4Colour::White()));
   fWomTubeLogical -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
   fWomTubeInsideLogical -> SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
   fPMTBodyLogical1 -> SetVisAttributes(new G4VisAttributes(G4Colour::Magenta()));
   fPMTBodyLogical2 -> SetVisAttributes(new G4VisAttributes(G4Colour::Magenta()));
   fPMTCathodeLogical1 -> SetVisAttributes(new G4VisAttributes(G4Colour::/*Red*/White()));
   fPMTCathodeLogical2 -> SetVisAttributes(new G4VisAttributes(G4Colour::/*Red*/White()));
   fWOMPaintLogical -> SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
}
void WOM::Construction()
{
    GenerateLogicals();
    //placing pressoure housing solids
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fGlassLogical, "glassPhysical", fLogicMother, false, 0, true);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fGelLogical, "gelPhysical", fGlassLogical, false, 0, true);
    //placing WOM Tube solids
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fWomTubeLogical, "WOMTubePhysical", fGelLogical, false, 0, true);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fWomTubeInsideLogical, "WOMTubeInsidePhysical", fWomTubeLogical, false, 0, true);
    //Placing PMT Solid
    new G4PVPlacement(0, G4ThreeVector(0, 0, fPMTGlobalZ), fPMTBodyLogical1, "PMT_Body_Physical1", fGelLogical, false, 0, true);
    //new G4PVPlacement(0, G4ThreeVector(0, 0, fPMTGlobalZ), fPMTBodyLogical1, "PMT_Body_Physical1", fLogicMother, false, 0, true); //delete this later!
    G4RotationMatrix* pmtGlobalRot = new G4RotationMatrix();
    pmtGlobalRot -> rotateY(CLHEP::pi);
    new G4PVPlacement(pmtGlobalRot, G4ThreeVector(0, 0, -fPMTGlobalZ), fPMTBodyLogical2, "PMT_Body_Physical2", fGelLogical, false, 0, true);
    //new G4PVPlacement(pmtGlobalRot, G4ThreeVector(0, 0, -fPMTGlobalZ), fPMTBodyLogical2, "PMT_Body_Physical2", fLogicMother, false, 0, true); //delete this later
    //Placing PMT Cathode Solid
    new G4PVPlacement(0, G4ThreeVector(0, 0, fPMTCathodeZ), fPMTCathodeLogical1, "PMT_Cathode_Physical1", fGelLogical, false, 0, true);
    //new G4PVPlacement(0, G4ThreeVector(0, 0, fPMTCathodeZ), fPMTCathodeLogical1, "PMT_Cathode_Physical1", fLogicMother, false, 0, true); //delete this later
    G4RotationMatrix* pmtCathodeRot = new G4RotationMatrix();
    pmtCathodeRot -> rotateY(CLHEP::pi);
    new G4PVPlacement(pmtCathodeRot, G4ThreeVector(0, 0, -fPMTCathodeZ), fPMTCathodeLogical2, "PMT_Cathode_Physical2", fGelLogical, false, 0, true);
    //new G4PVPlacement(pmtCathodeRot, G4ThreeVector(0, 0, -fPMTCathodeZ), fPMTCathodeLogical2, "PMT_Cathode_Physical2", fLogicMother, false, 0, true); //delete this later
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fWOMPaintLogical, "WOM_Paint_Logical", fGelLogical, false, 0, true);
    //new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fWOMPaintLogical, "WOM_Paint_Logical", fLogicMother, false, 0, true); //delete this later

}
G4MultiUnion* WOM::PressureVessel(const G4String& vesselName, G4double vesselTubeRad, G4double vesselCapRad) //for now it;s hardcoded but will change soon.
{
    G4MultiUnion* vesselSolid = new G4MultiUnion(vesselName);

//  Glass Tube Solid
    G4Tubs* vesselTubeSolid = new G4Tubs(vesselName + "TubeSolid", 0.0, vesselTubeRad, fGlassTubeHalfLength, 0., 2 * CLHEP::pi);
    G4ThreeVector vesselTubeCenter = G4ThreeVector(0., 0., 0.);
    G4RotationMatrix vesselTubeRot = G4RotationMatrix();
    G4Transform3D trVesselTube = G4Transform3D(vesselTubeRot, vesselTubeCenter);

//  Glass Cap Upper
    G4Sphere* vesselUpperCap = new G4Sphere(vesselName + "UpperCap", 0.0, vesselCapRad, 0.0, 2 * CLHEP::pi, 0.0, 0.5 * CLHEP::pi);
    G4ThreeVector vesselCapCenter1 = G4ThreeVector(0, 0, fGlassTubeHalfLength);
    G4RotationMatrix vesselCapUpperRot = G4RotationMatrix();
    G4Transform3D trVesselUpperCap = G4Transform3D(vesselCapUpperRot, vesselCapCenter1);

// Glass Cap Lower
    G4Sphere* vesselLowerCap = new G4Sphere(vesselName + "LowerCap", 0.0, vesselCapRad, 0.0, 2 * CLHEP::pi, 0.0, 0.5 * CLHEP::pi);
    G4ThreeVector vesselCapCenter2 = G4ThreeVector(0, 0, -fGlassTubeHalfLength);
    G4RotationMatrix vesselCapLowerRot = G4RotationMatrix();
    vesselCapLowerRot.rotateY(CLHEP::pi);
    G4Transform3D trVesselLowerCap = G4Transform3D(vesselCapLowerRot, vesselCapCenter2);

//  Adding all components to multiunion solid
    vesselSolid -> AddNode(*vesselTubeSolid, trVesselTube);
    vesselSolid -> AddNode(*vesselUpperCap, trVesselUpperCap);
    vesselSolid -> AddNode(*vesselLowerCap, trVesselLowerCap);

    vesselSolid -> Voxelize();
    return vesselSolid;
}
G4VSolid* WOM::WOMTube(const G4String& tubeName, G4double tubeRad)
{
    //WOm tube solid
    G4Tubs* vesselTubeSolid = new G4Tubs(tubeName + "TubeSolid", 0.0, tubeRad, fWomTubeHalfLength, 0., 2 * CLHEP::pi);

    return vesselTubeSolid;
}
G4VSolid* WOM::WOMPaint(const G4String& paintName)
{
    G4Tubs* paintTubeSolid = new G4Tubs(paintName + "PaintSolid", fWOMPaintInnerRad, fWOMPaintOuterRad, fWomTubeHalfLength, 0., 2 * CLHEP::pi);

    return paintTubeSolid;
}
G4VSolid* WOM::PMTConstruction()
{
    G4MultiUnion* pmtSolid = new G4MultiUnion("PMT_Body_Solid");

    //Holder
    G4Tubs* PMTHolderSolid = new G4Tubs("PMT_Holder_Solid", 0.0, fPMTHolderRad, fPMTHolderHalfLength, 0., 2 * CLHEP::pi);
    G4ThreeVector PMTHolderCenter = G4ThreeVector(0, 0, fPMTHolderZ + fPMTOffset);
    G4RotationMatrix PMTHolderRot = G4RotationMatrix();
    G4Transform3D trPMTHolder = G4Transform3D(PMTHolderRot, PMTHolderCenter);
    //Cone
    G4EllipticalCone* PMTConeSolid = new G4EllipticalCone("PMT_Cone_Solid", fPMTConeRatio, fPMTConeRatio, fPMTConeHeight, fPMTConeCut);
    G4ThreeVector PMTConeCenter = G4ThreeVector(0, 0, fPMTConeZ + fPMTOffset);
    G4RotationMatrix PMTConeRot = G4RotationMatrix();
    G4Transform3D trPMTCone = G4Transform3D(PMTConeRot, PMTConeCenter);
    //Tube
    G4Tubs* PMTTubeSolid = new G4Tubs("PMT_Tube_Solid", 0.0, fPMTTubeRad, fPMTTubeHalfLength, 0., 2 * CLHEP::pi);
    G4ThreeVector PMTTubeCenter = G4ThreeVector(0, 0, fPMTTubeZ + fPMTOffset);
    G4RotationMatrix PMTTubeRot = G4RotationMatrix();
    G4Transform3D trPMTTube = G4Transform3D(PMTTubeRot, PMTTubeCenter);
    //Board
    G4Tubs* PMTBoardSolid = new G4Tubs("PMT_Board_Solid", 0.0, fPMTBoardRad, fPMTBoardHalfLength, 0., 2 * CLHEP::pi);
    G4ThreeVector PMTBoardCenter = G4ThreeVector(0, 0, fPMTBoardZ + fPMTOffset);
    G4RotationMatrix PMTBoardRot = G4RotationMatrix();
    G4Transform3D trPMTBoard = G4Transform3D(PMTBoardRot, PMTBoardCenter);

    pmtSolid -> AddNode(*PMTHolderSolid, trPMTHolder);
    pmtSolid -> AddNode(*PMTConeSolid, trPMTCone);
    pmtSolid -> AddNode(*PMTTubeSolid, trPMTTube);
    pmtSolid -> AddNode(*PMTBoardSolid, trPMTBoard);

    pmtSolid -> Voxelize();
    return pmtSolid;
}
G4VSolid* WOM::PMTCathodeConstruction()
{
    G4Tubs* pmtCathode = new G4Tubs("PMT_Cathode_Solid", 0.0, fPMTHolderRad /*fWOMPaintOuterRad*/, fPMTCathodeHalfLength, 0., 2 * CLHEP::pi);
    return pmtCathode;
}
void WOM::SubtractHarnessPlug()
{

}
void WOM::SetPMTPosition()
{

}
void WOM::SelLEDPosition()
{

}
void WOM::ConstructMaterial()
{
    WOMMaterial* womMat = new WOMMaterial();

    glassMaterial = womMat -> GetGlassMaterial();
    fillerMaterial = womMat -> GetFillerMaterial();
    paintMaterial = womMat -> GetPaintMaterial();
    tubeMaterial = womMat -> GetTubeMaterial();
    tubeInsideMaterial = womMat -> GetTubeInsideMaterial();
    pmtAbsorber = fData -> GetMaterial("NoOptic_Absorber");
    pmtPhotoCathode = fData -> GetMaterial("RiAbs_Photocathode");

    G4NistManager* nist = G4NistManager::Instance();
    air = nist -> FindOrBuildMaterial("G4_AIR");

    delete womMat;
}
