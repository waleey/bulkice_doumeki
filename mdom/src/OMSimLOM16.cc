/** @file OMSimLOM16.cc
 *  @brief Construction of LOM16.
 *
 *  @author Javi Vara & Markus Dittmer
 *  @date February 2022
 *
 *  @version Geant4 10.7
 */

#include "OMSimLOM16.hh"
#include "abcDetectorComponent.hh"
#include <dirent.h>
#include <stdexcept>
#include <cstdlib>

#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalTube.hh"
#include "G4Sphere.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "G4Polycone.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include <G4UnitsTable.hh>
#include "G4VisAttributes.hh"
#include "G4Torus.hh"
#include "G4EllipticalCone.hh"

#include "CADMesh.hh" //CAD import from .obj files (tesselated -> not for optical stuff)



extern G4int gDOM;
extern G4bool gharness_ropes;
extern G4bool gVisual;
extern G4double gmdomseparation;
extern G4int gn_mDOMs;
extern G4bool gCADImport;

//main
//ToDo:
//Implement Harness
//Redo these imports: mTotalLenght = mData->GetValue("pmt_Hamamatsu_4inch", "jOuterShape.jTotalLenght");
LOM16::LOM16(OMSimInputData* pData, G4bool pPlaceHarness) {
    mData = pData;
    mPMTManager = new OMSimPMTConstruction(mData);
    mPMTManager->SelectPMT("argPMT");
    //mPMTManager->SimulateInternalReflections();
    GetSharedData();

    mPlaceHarness = pPlaceHarness;
    if (mPlaceHarness){
        //mHarness = new mDOMHarness(this, mData);
        //IntegrateDetectorComponent(mHarness, G4ThreeVector(0,0,0), G4RotationMatrix(), "");
        G4cout << "LOM16 harness not implemented yet" << G4endl;
    }

    Construction();
}

//Parameters from json files
void LOM16::GetSharedData() {
    //Module Parameters
    mGlassOutRad = mData->GetValue(mDataKey, "jGlassOutRad"); // outer radius of galss cylinder (pressure vessel)
    mCylHigh = mData->GetValue(mDataKey, "jCylHigh");         // height of cylindrical part of glass half-vessel
    mGlassThick = mData->GetValue(mDataKey, "jGlassThick");                     // maximum Glass thickness
    mGelThicknessFrontPMT = mData->GetValue(mDataKey, "jGelThicknessFrontPMT"); // distance between inner glass surface and tip of PMTs
    mThetaPolar = mData->GetValue(mDataKey, "jThetaPolar"); //theta angle polar pmts
    mThetaEquatorial = mData->GetValue(mDataKey, "jThetaEquatorial"); //theta angle equatorial pmts
    mCylinderAngle = mData->GetValue(mDataKey, "jCylinderAngle");  // Deviation angle of cylindrical part of the pressure vessel
    mNrPolarPMTs = mData->GetValue(mDataKey,"jNrPolarPMTs");
    mNrEqPMTs = mData->GetValue(mDataKey,"jNrEqPMTs");
    mEqTiltAngle = mData->GetValue(mDataKey,"jEqTiltAngle"); //tilt angle of gel pad axis in respect to PMT axis
    mPolEqPMTPhiPhase = mData->GetValue(mDataKey,"jPolEqPMTPhiPhase"); //rotation of equatorial PMTs in respect to polar PMTs
    mPolPadOpeningAngle = mData->GetValue(mDataKey,"jPolPadOpeningAngle"); //
    mEqPadOpeningAngle = mData->GetValue(mDataKey,"jEqPadOpeningAngle"); //
    
    mGlassInRad = mGlassOutRad - mGlassThick;
    mTotalNrPMTs = (mNrPolarPMTs + mNrEqPMTs) * 2;

    //PMT parameters
    mPMToffset = mPMTManager->GetDistancePMTCenterToPMTtip();
    mMaxPMTRadius = mPMTManager->GetMaxPMTMaxRadius() + 2 * mm;

    mTotalLenght = mData->GetValue("pmt_Hamamatsu_4inch", "jOuterShape.jTotalLenght");
    mOutRad = mData->GetValue("pmt_Hamamatsu_4inch", "jOuterShape.jOutRad");
    mSpherePos_y = mData->GetValue("pmt_Hamamatsu_4inch", "jOuterShape.jSpherePos_y");
    mEllipsePos_y = mData->GetValue("pmt_Hamamatsu_4inch", "jOuterShape.jEllipsePos_y");
}


//Placement function
void LOM16::Construction()
{
    //Create pressure vessel and inner volume
    G4UnionSolid* lGlassSolid = PressureVessel(mGlassOutRad, "Glass");
    G4UnionSolid* lGelSolid = PressureVessel(mGlassInRad, "Gel"); // Fill entire vessel with gel as logical volume (not placed) for intersectionsolids with gelpads
    
    //Set positions and rotations of PMTs and gelpads
    SetPMTAndGelpadPositions();
   
    //Logicals
    G4LogicalVolume* lGlassLogical = new G4LogicalVolume(lGlassSolid, mData->GetMaterial("argVesselGlass")," Glass_log"); //Vessel
    G4LogicalVolume* lInnerVolumeLogical = new G4LogicalVolume(lGelSolid, mData->GetMaterial("Ri_Air"), "Inner volume logical"); //Inner volume of vessel (mothervolume of all internal components)  
    CreateGelpadLogicalVolumes(lGelSolid); //logicalvolumes of all gelpads saved globally to be placed below
    G4LogicalVolume* lEquatorbandLogical = CreateEquatorBand(mGlassOutRad);

    //Placements 
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lInnerVolumeLogical, "Gel_physical", lGlassLogical, false, 0); //Innervolume (mother volume for all components)
    PlacePMTs(lInnerVolumeLogical);
    PlaceGelpads(lInnerVolumeLogical);
    if (gCADImport) PlaceCADSupportStructure(lInnerVolumeLogical); 
    new G4PVPlacement( 0 , G4ThreeVector(0, 0, 0) , lEquatorbandLogical, "Equatorband" , lInnerVolumeLogical, false, 0); 
    
    

    // ------------------ Add outer shape solid to MultiUnion in case you need substraction -------------------------------------------
    //Each Component needs to be appended to be places in abcDetectorComponent. Everything is placed in the gel which is placed in the glass which is the mother volume. This is the reason why not everything is appended on its own
    AppendComponent(lGlassSolid, lGlassLogical, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "LOM16");

    // ---------------- visualisation attributes --------------------------------------------------------------------------------
    lGlassLogical->SetVisAttributes(mGlassVis);
    lInnerVolumeLogical->SetVisAttributes(mInvisibleVis); //Material defined as Ri_Air
    for (int i = 0; i <= mTotalNrPMTs-1; i++) {
        mGelPad_logical[i]->SetVisAttributes(mGelVis); //mGelVis
    }
    lEquatorbandLogical->SetVisAttributes(mAbsorberVis);
    //mSupportStructureLogical->SetVisAttributes(mAluVis);

}


// ---------------- Component functions --------------------------------------------------------------------------------

G4UnionSolid* LOM16::PressureVessel(const G4double pOutRad, G4String pSuffix)
{
    G4Ellipsoid* lTopSolid = new G4Ellipsoid("SphereTop solid" + pSuffix, pOutRad, pOutRad, pOutRad, -5 * mm, pOutRad + 5 * mm);
    G4Ellipsoid* lBottomSolid = new G4Ellipsoid("SphereBottom solid" + pSuffix, pOutRad, pOutRad, pOutRad, -(pOutRad + 5 * mm), 5 * mm);

    G4double zCorners[] = { mCylHigh * 1.001, mCylHigh, 0, -mCylHigh, -mCylHigh * 1.001 };
    G4double rCorners[] = { 0, pOutRad, pOutRad + mCylHigh * sin(mCylinderAngle), pOutRad, 0 };
    G4Polycone* lCylinderSolid = new G4Polycone("Cylinder solid" + pSuffix, 0, 2 * CLHEP::pi, 5, rCorners, zCorners);

    G4UnionSolid* lTempUnion = new G4UnionSolid("temp" + pSuffix, lCylinderSolid, lTopSolid, 0, G4ThreeVector(0, 0, mCylHigh));
    G4UnionSolid* lUnionSolid = new G4UnionSolid("OM body" + pSuffix, lTempUnion, lBottomSolid, 0, G4ThreeVector(0, 0, -mCylHigh));
    return lUnionSolid;
}


G4LogicalVolume* LOM16::CreateEquatorBand(const G4double pOutRad)
{
    G4double lWidth = 45*mm; //Total width (both halves)
    G4double lThickness = 1*mm; //Thickness since its a 3D object

    G4Box* lCuttingBox = new G4Box("Cutter",pOutRad+10*mm,pOutRad+10*mm,lWidth/2.);
    G4UnionSolid* lOuter = PressureVessel(pOutRad + lThickness, "BandOuter");
    G4UnionSolid* lInner = PressureVessel(pOutRad , "BandInner");

    G4IntersectionSolid* lIntersectionSolid = new G4IntersectionSolid("BandThicknessBody" , lCuttingBox, lOuter, 0, G4ThreeVector(0, 0, 0));
    G4SubtractionSolid* lEquatorbandSolid = new G4SubtractionSolid("Equatorband_solid" , lIntersectionSolid, lInner, 0, G4ThreeVector(0, 0, 0));
    
    G4LogicalVolume* lEquatorbandLogical = new G4LogicalVolume(lEquatorbandSolid, mData->GetMaterial("NoOptic_Absorber"),"Equatorband_log"); 
    return lEquatorbandLogical;
}


// ---------------- Module specific funtions below --------------------------------------------------------------------------------

//Places a shrunken vessel in the centre acting as an absorber
void LOM16::PlaceDummySupportStructure(G4LogicalVolume* lInnerVolumeLogical)
{
    G4double shrinkfactor = 0.6;
    G4UnionSolid* SupportStructureSolid = PressureVessel(mGlassOutRad*shrinkfactor, "SimpleSupportStructure");
    mSupportStructureLogical = new G4LogicalVolume(SupportStructureSolid, mData->GetMaterial("NoOptic_Absorber")," Simple support structure"); //Vessel
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), mSupportStructureLogical, "Simple support structure", lInnerVolumeLogical, false, 0); //Innervolume (mother volume for all components)
    mSupportStructureLogical->SetVisAttributes(mAluVis);
}

//Places internal components based on a CAD file. OBJ filetype, modified by python script, tesselated via mesh.cc / Documentation here: https://github.com/christopherpoole/CADMesh
//To do: 
//Recreate CAD obj files (no SetScale & more variety to choose from)
//Each component has its own label in visualizer -> just one ... another function for penetrator, dummy main boards, ...
void LOM16::PlaceCADSupportStructure(G4LogicalVolume* lInnerVolumeLogical)
{
    //select file
    std::stringstream CADfile;
    CADfile.str(""); 
    //CADfile << "LOM_Internal.obj";
    CADfile << "LOM_Internal_WithPen_WithCali.obj";
    G4cout <<  "using the following CAD file for support structure: "  << CADfile.str()  << G4endl;

    //load mesh
    auto mesh = CADMesh::TessellatedMesh::FromOBJ("../data/CADmeshes/LOM16/" + CADfile.str() );
    G4ThreeVector CADoffset = G4ThreeVector(68.248*mm, 0, -124.218*mm); //measured from CAD file since origin =!= Module origin
    mesh->SetOffset(CADoffset);
   // mesh->SetScale(10); //did a mistake...this LOM_Internal file needs cm -> mm -> x10

    // Place all of the meshes it can find in the file as solids individually.
    for (auto solid : mesh->GetSolids())
    { 
        mSupportStructureLogical  = new G4LogicalVolume( solid , mData->GetMaterial("NoOptic_Absorber") , "logical" , 0, 0, 0);
        mSupportStructureLogical->SetVisAttributes(mAluVis);
        new G4PVPlacement( 0 , G4ThreeVector(0, 0, 0) , mSupportStructureLogical, "Support structure" , lInnerVolumeLogical, false, 0);
    }
}


//ToDo:
//sin(90 +- ...) -> cos(...)
void LOM16::SetPMTAndGelpadPositions()
{
    G4double Z_center_module_tobottomPMT_polar = 70.9*mm; //measured z-offset from vessel origin from CAD file
    G4double Z_center_module_tobottomPMT_equatorial = 25.4*mm; //measured z-offset from vessel origin from CAD file
    jGelPadDZ   = 30* mm;  // semiaxis (along pmt axis) of gelpads ... simply needs to be larger then 5mm (+ some more for tilted pads)...could be 100

    //calculate distances
    const G4double lFrontToEllipse_y = mOutRad + mSpherePos_y - mEllipsePos_y; //Center of PMT (PMT solid in Geant4) to tip
    const G4double lZ_centertobottomPMT= mTotalLenght-lFrontToEllipse_y; //Distance from bottom base of PMT to center of the Geant4 solid
    
    //calculate PMT and pads positions in the usual 4 for loop way 
    for (int i = 0; i <= mTotalNrPMTs-1; i++) {

        //upper polar
        if (i>=0 && i<=mNrPolarPMTs-1){ 
            lPMT_theta=mThetaPolar; 
            lPMT_phi=mPolEqPMTPhiPhase+i*90.0*deg;

            PMT_rho = (147.7406-4.9)*mm;  // For Position of PMT to Vessel wall (147.7406). 4.9mm is distance from PMT photocathode to vessel inner surface
            GelPad_rho = PMT_rho + 2*jGelPadDZ;
            
            PMT_z = Z_center_module_tobottomPMT_polar+lZ_centertobottomPMT*sin(90*deg-lPMT_theta);
            GelPad_z = PMT_z+2*jGelPadDZ*sin(90*deg-lPMT_theta);
        }

        //upper equatorial
        if (i>=mNrPolarPMTs && i<=mNrPolarPMTs + mNrEqPMTs -1){ 
            lPMT_theta= mThetaEquatorial; 
            lPMT_phi=i*90.0*deg; 

            PMT_rho = (120.1640-5)*mm; // For Position of PMT to Vessel wall (147.7406). 4.9mm is distance from PMT photocathode to vessel inner surface
            GelPad_rho = PMT_rho + 2*jGelPadDZ;

            PMT_z=Z_center_module_tobottomPMT_equatorial+lZ_centertobottomPMT*sin(90*deg-lPMT_theta);
        }

        //lower equatorial
        if (i>=mNrPolarPMTs + mNrEqPMTs && i<= mNrPolarPMTs + mNrEqPMTs + mNrEqPMTs -1){ 
            lPMT_theta= 180*deg - mThetaEquatorial; //118*deg; // 118.5*deg;
            lPMT_phi=i*90.0*deg;

            PMT_rho = (120.1640-5)*mm; // For Position of PMT to Vessel wall (147.7406). 4.9mm is distance from PMT photocathode to vessel inner surface
            GelPad_rho = PMT_rho + 2*jGelPadDZ ;

            PMT_z=(-Z_center_module_tobottomPMT_equatorial)-lZ_centertobottomPMT*sin(90*deg-(180*deg-lPMT_theta)); //rodo to cos ...
        }

        //lower polar
        if (i>=mTotalNrPMTs-mNrPolarPMTs && i<=mTotalNrPMTs-1){ 
            lPMT_theta=180*deg - mThetaPolar; //144*deg; //152*deg;
            lPMT_phi=mPolEqPMTPhiPhase +((i)*90.0)*deg;

            PMT_rho = (147.7406-4.9)*mm;  // For Position of PMT to Vessel wall (147.7406). 4.9mm is distance from PMT photocathode to vessel inner surface
            GelPad_rho = PMT_rho + 2*jGelPadDZ;

            PMT_z=(-Z_center_module_tobottomPMT_polar)-lZ_centertobottomPMT*sin(90*deg-(180*deg-lPMT_theta));
            GelPad_z = PMT_z-2*jGelPadDZ*sin(90*deg-(180*deg-lPMT_theta));
        }

    //PMTs
    lPMT_x = (PMT_rho*mm)* sin(lPMT_theta) * cos(lPMT_phi);
    lPMT_y = (PMT_rho*mm)* sin(lPMT_theta) * sin(lPMT_phi);

    //Gelpads
    GelPad_x = GelPad_rho * sin(lPMT_theta) * cos(lPMT_phi);
    GelPad_y = GelPad_rho * sin(lPMT_theta) * sin(lPMT_phi);    

    //save positions in arrays
    mPMTPositions.push_back( G4ThreeVector(lPMT_x,lPMT_y,PMT_z) );
    mGelpadPositions.push_back( G4ThreeVector(GelPad_x,GelPad_y,GelPad_z) );

    //save angles in arrays
    mPMT_theta.push_back(lPMT_theta);
    mPMT_phi.push_back(lPMT_phi);
    }
}    

//Todo
//4*jGelPadDZ -> 2 2*jGelPadDZ -> 1 ...  not needed if jGelPadDZ is long enough
//tra and transformers declaration uniformely.
//rename some stuff for clarity
void LOM16::CreateGelpadLogicalVolumes(G4UnionSolid* lGelSolid) 
{
    //getting the PMT solid
    G4UnionSolid* lPMTsolid = mPMTManager->GetPMTSolid();

    //Definition of helper volumes
    G4Cons* lGelPadBasicSolid;
    G4EllipticalCone* lGelPadBasicSolid_eq;
    G4IntersectionSolid* Cut_cone;
    G4LogicalVolume* GelPad_logical;
    G4SubtractionSolid* Cut_cone_final;

    // Definition of semiaxes in (elliptical section) cone for titled gel pads
    G4double dx=std::cos(mEqTiltAngle)/(2*(1+std::sin(mEqTiltAngle)*std::tan(mEqPadOpeningAngle)))*mMaxPMTRadius*2;  //semiaxis y at -ztop
    G4double ztop=2*jGelPadDZ;
    G4double dy= mMaxPMTRadius;  //semiaxis x at -ztop
    G4double Dy=dy+2*ztop*std::tan(mEqPadOpeningAngle);  //semiaxis x at +ztop
    G4double Dx=dx*Dy/dy;  //     Dx/dx=Dy/dy always ;  //semiaxis y at -ztop
    G4double xsemiaxis=(Dx-dx)/(2*ztop);  // Best way it can be defined
    G4double ysemiaxis=(Dy-dy)/(2*ztop); // Best way it can be defined
    G4double zmax= (Dx+dx)/(2*xsemiaxis); // Best way it can be defined

    // For the placement of tilted  equatorial pads
    G4double dz=-dx*std::sin(mEqTiltAngle);
    G4double dY3=mMaxPMTRadius-(ztop*std::sin(mEqTiltAngle)+dx*std::cos(mEqTiltAngle));

    //create logical volume for each gelpad    
    for(int k=0 ; k<=mTotalNrPMTs-1 ; k++){
            G4Transform3D* tra;
            G4Transform3D* tra2;
            
            converter.str("");
            converter2.str("");
            converter << "GelPad_" << k << "_solid";
            converter2 << "Gelpad_final" << k << "_logical";
            
            // polar gel pads
            if(k<=mNrPolarPMTs-1 or k>=mTotalNrPMTs-mNrPolarPMTs){ 
                lGelPadBasicSolid = new G4Cons("GelPadBasic", 0,mMaxPMTRadius , 0, mMaxPMTRadius  + 4*jGelPadDZ*tan(mPolPadOpeningAngle), 2*jGelPadDZ, 0, 2*CLHEP::pi);
                
                //rotation and position of gelpad
                G4RotationMatrix* rotation = new G4RotationMatrix();
                rotation->rotateY(mPMT_theta[k]);
                rotation->rotateZ(mPMT_phi[k]);

                tra = new G4Transform3D(*rotation, G4ThreeVector(mGelpadPositions[k]));
                G4Transform3D transformers = G4Transform3D(*rotation, G4ThreeVector(mPMTPositions[k]));

                //creating volumes ... basic cone, subtract PMT, logical volume of gelpad
                Cut_cone = new G4IntersectionSolid(converter.str(), lGelSolid, lGelPadBasicSolid, *tra);
                Cut_cone_final = new G4SubtractionSolid(converter.str(), Cut_cone, lPMTsolid, transformers);
                GelPad_logical = new G4LogicalVolume(Cut_cone_final, mData->GetMaterial("argGel"), converter2.str());
            };

            //upper equatorial
            if(k>=mNrPolarPMTs && k<=mNrPolarPMTs+mNrEqPMTs-1){
                //the tilted pad is a cone with elliptical section. Moreover, there is an union with a disc called "Tube". This union does not affect the geometry, but it is used to shift the center of the geant solid so it can be placed in the same way as the PMT. In addition to that, Cut_tube is used to do a better substraction of the residual part that lies below the photocatode. 
                lGelPadBasicSolid_eq =  new G4EllipticalCone("cone",xsemiaxis,ysemiaxis,zmax,ztop);

                //rotation and position of tilted gelpad
                G4RotationMatrix* rotation = new G4RotationMatrix();
                rotation->rotateY((180*deg+mEqTiltAngle));
                tra2 = new G4Transform3D(*rotation, G4ThreeVector(-dY3,0,ztop*std::cos(mEqTiltAngle)+dz));

                //place tilted cone
                G4Tubs *Tube = new G4Tubs("tube",0,mMaxPMTRadius,0.1*mm,0.,2*CLHEP::pi);
                G4UnionSolid *Tilted_cone = new G4UnionSolid("union",Tube,  lGelPadBasicSolid_eq,*tra2);

                //For subtraction of gelpad reaching below the photocathode
                rotation = new G4RotationMatrix();
                tra = new G4Transform3D(*rotation, G4ThreeVector(0,0,-30*mm)); //30 -> thickness/width of subtraction box ... just has to be large eough to cut overlapping gelpads under photocathode

                //rotation and position of PMT
                G4RotationMatrix* rot = new G4RotationMatrix();
                rot->rotateY(mPMT_theta[k]);
                rot->rotateZ(mPMT_phi[k]);
                G4Transform3D transformers = G4Transform3D(*rot, G4ThreeVector(mPMTPositions[k]));
                
                //creating volumes ... basic cone, tilted cone, subtract PMT, logical volume of gelpad
                G4Tubs* Cut_tube = new G4Tubs("tub",0,2*mMaxPMTRadius,30*mm,0,2*CLHEP::pi);
                G4SubtractionSolid* Tilted_final = new G4SubtractionSolid("cut",Tilted_cone,Cut_tube,*tra);
                Cut_cone = new G4IntersectionSolid(converter.str(), lGelSolid, Tilted_final,  transformers);
                Cut_cone_final = new G4SubtractionSolid(converter.str(), Cut_cone, lPMTsolid, transformers);
                GelPad_logical = new G4LogicalVolume(Cut_cone_final, mData->GetMaterial("argGel"), converter2.str());
            };  

            //lower equatorial
            if(k >=mNrPolarPMTs+mNrEqPMTs && k<=mTotalNrPMTs-mNrPolarPMTs-1){
                lGelPadBasicSolid_eq =  new G4EllipticalCone("cone",xsemiaxis,ysemiaxis,zmax,ztop);

                //rotation and position of tilted gelpad
                G4RotationMatrix* rotation = new G4RotationMatrix();
                rotation->rotateY((180*deg-mEqTiltAngle));
                tra2 = new G4Transform3D(*rotation, G4ThreeVector(dY3,0,ztop*std::cos(mEqTiltAngle)+dz));

                //place tilted cone
                G4Tubs *Tube= new G4Tubs("tube",0,mMaxPMTRadius,0.1*mm,0.,2*CLHEP::pi);
                G4UnionSolid *Tilted_cone=new G4UnionSolid("union",Tube,  lGelPadBasicSolid_eq,*tra2);
                
                //For subtraction of gelpad reaching below the photocathode
                rotation = new G4RotationMatrix();
                tra = new G4Transform3D(*rotation, G4ThreeVector(0,0,-30*mm));

                //rotation and position of PMT
                G4RotationMatrix* rot = new G4RotationMatrix();
                rot->rotateY(mPMT_theta[k]);
                rot->rotateZ(mPMT_phi[k]);
                G4Transform3D transformers = G4Transform3D(*rot, G4ThreeVector(mPMTPositions[k]));
                
                //creating volumes ... basic cone, tilted cone, subtract PMT, logical volume of gelpad
                G4Tubs* Cut_tube=new G4Tubs("tub",0,2*mMaxPMTRadius,30*mm,0,2*CLHEP::pi);
                G4SubtractionSolid* Tilted_final = new G4SubtractionSolid("cut",Tilted_cone,Cut_tube,*tra);
                Cut_cone = new G4IntersectionSolid(converter.str(), lGelSolid, Tilted_final,  transformers);
                Cut_cone_final = new G4SubtractionSolid(converter.str(), Cut_cone, lPMTsolid, transformers);
                GelPad_logical = new G4LogicalVolume(Cut_cone_final, mData->GetMaterial("argGel"), converter2.str());
            };   

    //save logicalvolume of gelpads in array
    mGelPad_logical.push_back( GelPad_logical ); //
    }
}


void LOM16::PlacePMTs(G4LogicalVolume* lInnerVolumeLogical)
{
    for(int k=0 ; k<=mTotalNrPMTs-1 ; k++){
        converter.str("");
        converter << k << "_physical";

        lRot = new G4RotationMatrix();
        lRot->rotateY(mPMT_theta[k]);
        lRot->rotateZ(mPMT_phi[k]);
        lTransformers = G4Transform3D(*lRot, G4ThreeVector(mPMTPositions[k]));

        mPMTManager->PlaceIt(lTransformers, lInnerVolumeLogical, converter.str());
    }

}

void LOM16::PlaceGelpads(G4LogicalVolume* lInnerVolumeLogical)
{
    for(int k=0 ; k<=mTotalNrPMTs-1 ; k++){
        converter.str("");
        converter << "GelPad_" << k;

        lRot = new G4RotationMatrix();
        lRot->rotateY(mPMT_theta[k]);
        lRot->rotateZ(mPMT_phi[k]);
        lTransformers = G4Transform3D(*lRot, G4ThreeVector(mPMTPositions[k]));

        new G4PVPlacement(0, G4ThreeVector(0,0,0), mGelPad_logical[k], converter.str(), lInnerVolumeLogical, false, 0);
    }
}
