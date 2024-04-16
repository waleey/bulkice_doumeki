/** @file OMSimDetectorConstruction.cc
 *  @brief User defined detector.
 *
 * You should define and construct the detector here...this template is an example for a single arg module.
 *
 *  @author Martin Unland - modified by Markus Dittmer
 *  @date October 2021 - modified at February 2022
 *  @Modified by Waly M Z Karim
 *  @version Geant4 11.1
 *
 *  @todo
 */

#include "OMSimDetectorConstruction.hh"

#include "OMSimRadioactivityData.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"



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
    /**
    *Defining the world params
    for GWSO Cylider
    **/
    G4double heightOffset = 0.05 * m;
    G4double outerRad = (5.2 / 2 ) * m;
    G4double innerRad = 0 * m;
    G4double diskInnerRad = 0.5 * m;
    G4double absLayerThickness = 1 * mm;
    G4double topHalfLength = (3.2 / 2) * m - heightOffset;
    G4double topZ = 0.4 * m;
    G4double bottomHalfLength = 0.4 * m;
    G4double bottomZ = -1.6 *m;
    G4double diskHalfLength = heightOffset / 2;
    G4double diskZ = - 1.175 * m;
    G4double reflHalfLength = 0.5 * mm;
    G4double reflZ = -2 * m - reflHalfLength;

    G4Material* worldMat = new G4Material("Air", density = 1.29 * mg/cm3, nelements = 2);
    G4Material* cylMat = new G4Material("Water", density = 1.0 * g/cm3, nelements = 2);
    ConstructWorldMat(cylMat, worldMat);
    G4Material* absorber = mData -> GetMaterial("NoOptic_Absorber");
    //G4Material* reflector = mData -> GetMaterial("NoOptic_Reflector");
    auto reflProperty = mData -> GetOpticalSurface("Refl_100polished");
    /**
    *Creating the world volume
    **/
    mWorldSolid = new G4Box("World", 20 * m, 20 * m, 20 * m);
    mWorldLogical = new G4LogicalVolume(mWorldSolid, worldMat, "World_log");
    mWorldPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), mWorldLogical, "World_phys", 0, false, 0);


    /**
    *Generating top cylinder
    **/
    G4VSolid* topSolid = new G4Tubs("Top_Cylinder", innerRad, outerRad, topHalfLength, 0, 2 * CLHEP::pi);
    G4LogicalVolume* topLogical = new G4LogicalVolume(topSolid, cylMat, "Top_Cylinder_log");
    G4VPhysicalVolume* topPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., topZ), topLogical, "Top_Cylinder_phys", mWorldLogical, false, 0, true);

    /**
    *Generating top absorption layer
    **/
    G4VSolid* topAbsSolid = new G4Tubs("Top_Abs_Cylinder", outerRad - absLayerThickness, outerRad, topHalfLength, 0, 2 * CLHEP::pi);
    G4LogicalVolume* topAbsLogical = new G4LogicalVolume(topAbsSolid, absorber, "Top_Abs_Cylinder_log");
    G4VPhysicalVolume* topAbsPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., 0), topAbsLogical, "Top_Abs_Cylinder_phys", topLogical, false, 0, true);


    /**
    *Generating bottom cylinder
    **/
    G4VSolid* bottomSolid = new G4Tubs("Bottom_Cylinder", innerRad, outerRad, bottomHalfLength, 0, 2 * CLHEP::pi);
    G4LogicalVolume* bottomLogical = new G4LogicalVolume(bottomSolid, cylMat, "Bottom_Cylinder_log");
    G4VPhysicalVolume* bottomPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., bottomZ), bottomLogical, "Bottom_Cylinder_phys", mWorldLogical, false, 0, true);

    /**
    * Generating bottom reflective layer


    G4VSolid* bottomReflSolid = new G4Tubs("Bottom_Refl_Cylinder", outerRad - absLayerThickness, outerRad, bottomHalfLength, 0, 2 * CLHEP::pi);
    G4LogicalVolume* bottomReflLogical = new G4LogicalVolume(bottomReflSolid, cylMat, "Bottom_Refl_Cylinder_log");
    G4VPhysicalVolume* bottomReflPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0.,0.), bottomReflLogical, "Bottom_Refl_Cylinder_phys", bottomLogical, false, 0, true);
    **/
    /**
    *Generating disk
    **/
    G4VSolid* diskSolid = new G4Tubs("Disk_Cylinder", diskInnerRad, outerRad, diskHalfLength, 0, 2 * CLHEP::pi);
    G4LogicalVolume* diskLogical = new G4LogicalVolume(diskSolid, absorber, "Disk_logical");
    G4VPhysicalVolume* diskPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, diskZ), diskLogical, "Disk_Physical", mWorldLogical, false, 0, true);

    /**
    *Generating reflecting disk

    G4VSolid* diskReflSolid = new G4Tubs("Disk_Refl_Cylinder", innerRad, outerRad, reflHalfLength, 0, 2 * CLHEP::pi);
    G4LogicalVolume* diskReflLogical = new G4LogicalVolume(diskReflSolid, reflector, "Disk_Refl_logical");
    G4VPhysicalVolume* diskReflPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, reflZ), diskReflLogical, "Disk_Refl_Physical", mWorldLogical, false, 0, true);
    **/

    /**
    *Generating a relfecting border surface
    **/
    auto reflSurface = new G4LogicalBorderSurface("Bottom_Reflective_sruface", bottomPhysical, mWorldPhysical, reflProperty);

    G4VisAttributes* Cyl_vis = new G4VisAttributes(G4Colour(0.45,0.5,0.85,0.2));
    G4VisAttributes* World_vis = new G4VisAttributes(G4Colour(0,0,0,0));


    mWorldLogical -> SetVisAttributes(World_vis);
    topLogical -> SetVisAttributes(Cyl_vis);
    bottomLogical -> SetVisAttributes(Cyl_vis);
    diskLogical -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));


    G4cout << "::::::::::::::::::World Volume Constructed Sucessfully:::::::::::::::::::" << G4endl;

}

// DONE UP TO THIS POINT
void OMSimDetectorConstruction::ConstructWorldMat(G4Material* cylMat, G4Material* worldMat)
{
/**
*This material properties
*are taken from Geant4 examples
*Geant4 OpNovice example
**/

 auto H = new G4Element("Hydrogen", "H", z = 1, a = 1.01 * g/mole);
 auto O = new G4Element("Oxygen", "O", z = 8, a = 16.00 * g/mole);
 auto N = new G4Element("Nitrogen", "N", z = 7, a = 14.01 * g / mole);
 worldMat -> AddElement(N, 70. * perCent);
 worldMat -> AddElement(O, 30. * perCent);

 cylMat -> AddElement(H, 2);
 cylMat -> AddElement(O, 1);


  // ------------ Generate & Add Material Properties Table ------------
  //
  std::vector<G4double> photonEnergy = {
    2.034 * eV, 2.068 * eV, 2.103 * eV, 2.139 * eV, 2.177 * eV, 2.216 * eV,
    2.256 * eV, 2.298 * eV, 2.341 * eV, 2.386 * eV, 2.433 * eV, 2.481 * eV,
    2.532 * eV, 2.585 * eV, 2.640 * eV, 2.697 * eV, 2.757 * eV, 2.820 * eV,
    2.885 * eV, 2.954 * eV, 3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV,
    3.353 * eV, 3.446 * eV, 3.545 * eV, 3.649 * eV, 3.760 * eV, 3.877 * eV,
    4.002 * eV, 4.136 * eV
  };

  // Water
  std::vector<G4double> refractiveIndex1 = {
    1.3435, 1.344,  1.3445, 1.345,  1.3455, 1.346,  1.3465, 1.347,
    1.3475, 1.348,  1.3485, 1.3492, 1.35,   1.3505, 1.351,  1.3518,
    1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,  1.3555, 1.356,
    1.3568, 1.3572, 1.358,  1.3585, 1.359,  1.3595, 1.36,   1.3608
  };
  std::vector<G4double> absorption = {
    3.448 * m,  4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m,
    15.152 * m, 17.241 * m, 18.868 * m, 20.000 * m, 26.316 * m, 35.714 * m,
    45.455 * m, 47.619 * m, 52.632 * m, 52.632 * m, 55.556 * m, 52.632 * m,
    52.632 * m, 47.619 * m, 45.455 * m, 41.667 * m, 37.037 * m, 33.333 * m,
    30.000 * m, 28.500 * m, 27.000 * m, 24.500 * m, 22.000 * m, 19.500 * m,
    17.500 * m, 14.500 * m
  };

  // Material properties can be added as arrays. However, in this case it is
  // up to the user to make sure both arrays have the same number of elements.
  G4double scintilFastArray[]{ 1.0, 1.0 };
  G4double energyArray[]{ 2.034 * eV, 4.136 * eV };
  G4int lenArray = 2;

  std::vector<G4double> scintilSlow = {
    0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
    7.00, 6.00, 4.00, 3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00, 4.00,
    5.00, 6.00, 7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 5.00, 4.00
  };

  auto myMPT1 = new G4MaterialPropertiesTable();

  // Values can be added to the material property table individually.
  // With this method, spline interpolation cannot be set. Arguments
  // createNewKey and spline both take their default values of false.
  // Need to specify the number of entries (1) in the arrays, as an argument
  // to AddProperty.
  G4int numEntries = 1;
  myMPT1->AddProperty("RINDEX", &photonEnergy[0], &refractiveIndex1[0],
                      numEntries);

  for(size_t i = 1; i < photonEnergy.size(); ++i)
  {
    myMPT1->AddEntry("RINDEX", photonEnergy[i], refractiveIndex1[i]);
  }

  // Check that group velocity is calculated from RINDEX
  if(myMPT1->GetProperty("RINDEX")->GetVectorLength() !=
     myMPT1->GetProperty("GROUPVEL")->GetVectorLength())
  {
    G4ExceptionDescription ed;
    ed << "Error calculating group velocities. Incorrect number of entries "
          "in group velocity material property vector.";
    G4Exception("OpNovice::OpNoviceDetectorConstruction", "OpNovice001",
                FatalException, ed);
  }

  // Adding a property from two std::vectors. Argument createNewKey is false
  // and spline is true.
  myMPT1->AddProperty("ABSLENGTH", photonEnergy, absorption, false, true);

  // Adding a property using a C-style array.
  // Spline interpolation isn't used for scintillation.
  // Arguments spline and createNewKey both take default value false.
  myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", energyArray, scintilFastArray,
                      lenArray);

  myMPT1->AddProperty("SCINTILLATIONCOMPONENT2", photonEnergy, scintilSlow,
                      false, true);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD", 50. / MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE", 1.0);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1. * ns);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 10. * ns);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD1", 0.8);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD2", 0.2);
  std::vector<G4double> energy_water = {
    1.56962 * eV, 1.58974 * eV, 1.61039 * eV, 1.63157 * eV, 1.65333 * eV,
    1.67567 * eV, 1.69863 * eV, 1.72222 * eV, 1.74647 * eV, 1.77142 * eV,
    1.7971 * eV,  1.82352 * eV, 1.85074 * eV, 1.87878 * eV, 1.90769 * eV,
    1.93749 * eV, 1.96825 * eV, 1.99999 * eV, 2.03278 * eV, 2.06666 * eV,
    2.10169 * eV, 2.13793 * eV, 2.17543 * eV, 2.21428 * eV, 2.25454 * eV,
    2.29629 * eV, 2.33962 * eV, 2.38461 * eV, 2.43137 * eV, 2.47999 * eV,
    2.53061 * eV, 2.58333 * eV, 2.63829 * eV, 2.69565 * eV, 2.75555 * eV,
    2.81817 * eV, 2.88371 * eV, 2.95237 * eV, 3.02438 * eV, 3.09999 * eV,
    3.17948 * eV, 3.26315 * eV, 3.35134 * eV, 3.44444 * eV, 3.54285 * eV,
    3.64705 * eV, 3.75757 * eV, 3.87499 * eV, 3.99999 * eV, 4.13332 * eV,
    4.27585 * eV, 4.42856 * eV, 4.59258 * eV, 4.76922 * eV, 4.95999 * eV,
    5.16665 * eV, 5.39129 * eV, 5.63635 * eV, 5.90475 * eV, 6.19998 * eV
  };

  // Rayleigh scattering length is calculated by G4OpRayleigh

  // Mie: assume 100 times larger than the rayleigh scattering
  std::vector<G4double> mie_water = {
    167024.4 * m, 158726.7 * m, 150742 * m,   143062.5 * m, 135680.2 * m,
    128587.4 * m, 121776.3 * m, 115239.5 * m, 108969.5 * m, 102958.8 * m,
    97200.35 * m, 91686.86 * m, 86411.33 * m, 81366.79 * m, 76546.42 * m,
    71943.46 * m, 67551.29 * m, 63363.36 * m, 59373.25 * m, 55574.61 * m,
    51961.24 * m, 48527.00 * m, 45265.87 * m, 42171.94 * m, 39239.39 * m,
    36462.50 * m, 33835.68 * m, 31353.41 * m, 29010.30 * m, 26801.03 * m,
    24720.42 * m, 22763.36 * m, 20924.88 * m, 19200.07 * m, 17584.16 * m,
    16072.45 * m, 14660.38 * m, 13343.46 * m, 12117.33 * m, 10977.70 * m,
    9920.416 * m, 8941.407 * m, 8036.711 * m, 7202.470 * m, 6434.927 * m,
    5730.429 * m, 5085.425 * m, 4496.467 * m, 3960.210 * m, 3473.413 * m,
    3032.937 * m, 2635.746 * m, 2278.907 * m, 1959.588 * m, 1675.064 * m,
    1422.710 * m, 1200.004 * m, 1004.528 * m, 833.9666 * m, 686.1063 * m
  };

  // Mie: gforward, gbackward, forward backward ratio
  G4double mie_water_const[3] = { 0.99, 0.99, 0.8 };

  myMPT1->AddProperty("MIEHG", energy_water, mie_water, false, true);
  myMPT1->AddConstProperty("MIEHG_FORWARD", mie_water_const[0]);
  myMPT1->AddConstProperty("MIEHG_BACKWARD", mie_water_const[1]);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO", mie_water_const[2]);

  G4cout << "Water G4MaterialPropertiesTable:" << G4endl;
  myMPT1->DumpTable();

  cylMat->SetMaterialPropertiesTable(myMPT1);

  // Set the Birks Constant for the Water scintillator
  cylMat->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  // Air
  std::vector<G4double> refractiveIndex2 = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                             1.0, 1.0, 1.0, 1.0 };

  auto myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", photonEnergy, refractiveIndex2);

  G4cout << "Air G4MaterialPropertiesTable:" << G4endl;
  myMPT2->DumpTable();

  worldMat->SetMaterialPropertiesTable(myMPT2);
}

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
    G4double diskZ = - 1.175 * m;

    ConstructWorld();
    G4ThreeVector modulePos = G4ThreeVector(0, 0, diskZ);

   if (fDOM == 0){ //Single PMT
        G4cout << "Constructing single PMT" << G4endl;
        mPMTManager = new OMSimPMTConstruction(mData);
        mPMTManager->SelectPMT("argPMT");

        G4RotationMatrix* lRot = new G4RotationMatrix();
        G4Transform3D lTransformers = G4Transform3D(*lRot, modulePos);
        mPMTManager->PlaceIt(lTransformers, mWorldLogical, "PMT");
    }
    else if (fDOM == 1){ //mDOM
        G4cout << "Constructing mDOM" << G4endl;
        fMDOM = new mDOM(mData,gPlaceHarness);
        glassOutRad = fMDOM -> getGlassOutRad();
        glassInRad = fMDOM -> getGlassInRad();
        OMSimRadioactivityData::SetGlassRad(glassOutRad, glassInRad);

        //"Rotating the module for GWSO"
        G4RotationMatrix* mRot = new G4RotationMatrix();
        //mRot -> rotateY(CLHEP::pi / 2);
        fMDOM->PlaceIt(modulePos, *mRot, mWorldLogical, "");
        delete mRot;
    }
    else if (fDOM == 2){ //PDOM
        G4cout << "Constructing PDOM" << G4endl;
        fPDOM = new pDOM(mData,gPlaceHarness);
        //"Rotating the module for GWSO"
        G4RotationMatrix* mRot = new G4RotationMatrix();
        //mRot -> rotateY(CLHEP::pi / 2);
        fPDOM->PlaceIt(modulePos, *mRot, mWorldLogical, "");
        delete mRot;
    }
    else if (fDOM == 3){ //LOM16
        G4cout << "Constructing LOM16" << G4endl;
        fLOM16 = new LOM16(mData,gPlaceHarness);
        G4RotationMatrix* mRot = new G4RotationMatrix();
        // -> rotateY(90 * deg);
        fLOM16->PlaceIt(modulePos, *mRot, mWorldLogical, "");
        delete mRot;
    }
    else if (fDOM == 4){ //LOM18
        G4cout << "Constructing LOM18" << G4endl;
        fLOM18 = new LOM18(mData,gPlaceHarness);
        //"Rotating the module for GWSO"
        G4RotationMatrix* mRot = new G4RotationMatrix();
        //mRot -> rotateY(90 * deg);
        fLOM18->PlaceIt(modulePos, *mRot, mWorldLogical, ""); //temporarily z = -1
        G4cout << "::::::::::::::LOM18 successfully constructed::::::::::::" << G4endl;
        delete mRot;
    }
     else if (fDOM == 5){ //DEGG
        G4cout << "Constructing DEGG" << G4endl;
        fDEGG = new dEGG(mData,gPlaceHarness);
        //"Rotating the module for GWSO"
        G4RotationMatrix* mRot = new G4RotationMatrix();
        //mRot -> rotateY(CLHEP::pi / 2);
        fDEGG->PlaceIt(modulePos, *mRot, mWorldLogical, "");
        glassOutRad = 350;
        glassInRad = 50;
        delete mRot;
     }
     else if(fDOM == 6)
     {
        G4cout << "Constructing WOM" << G4endl;
        fWOM = new WOM(mWorldLogical, mData);
        G4RotationMatrix* mRot = new G4RotationMatrix();
        //mRot -> rotateY(CLHEP::pi / 2);
        fWOM -> PlaceIt(*mRot, modulePos);
        delete mRot;
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

