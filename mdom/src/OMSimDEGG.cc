/** @file OMSimDEGG.cc
 *  @brief Construction of DEGG.
 *
 *  @author , modified by Berit Schlüter
 *  @date February 2022
 *
 *  @version Geant4 10.7
 *
 */


#include "OMSimDEGG.hh"
#include "OMSimDEGGHarness.hh" 
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
#include "CADMesh.hh" 


extern G4bool gCADImport;

dEGG::dEGG(OMSimInputData* pData, G4bool pPlaceHarness){
   mData = pData;
   mPMTManager = new OMSimPMTConstruction(mData);
   mPMTManager->SelectPMT("pmt_dEGG");

   if (pPlaceHarness){
      mHarness = new dEGGHarness(this,mData);
      IntegrateDetectorComponent(mHarness, G4ThreeVector(0,0,0), G4RotationMatrix(), "");
   }

   GetSharedData();
   Construction();
}


/**
 * @brief Get the Values from the json file om_DEGG.dat. If you want to change any values, change them there.
 * 
 */
void dEGG::GetSharedData(){
   // Outer shape of the vessel
   mOutSegments1 = mData->GetValue(mDataKey, "jOutSegments1");
   mOutSphereRadiusMax = mData->GetValue(mDataKey, "jOutSphereRadiusMax");
   mOutSphereDtheta = mData->GetValue(mDataKey, "jOutSphereDtheta");
   mOutTransformZ = mData->GetValue(mDataKey, "jOutTransformZ");
   mOutTorusRadius1 = mData->GetValue(mDataKey, "jOutTorusRadius1"); // radius of small spindle torus sphere
   mOutCenterOfTorusRadius1 = mData->GetValue(mDataKey, "jOutCenterOfTorusRadius1"); // distance from center of torus to z-axis
   mOutSegments2 = mData->GetValue(mDataKey, "jOutSegments2");
   mOutTorusRadius2 = mData->GetValue(mDataKey, "jOutTorusRadius2"); // radius of large spindle torus sphere
   mOutCenterOfTorusRadius2 = mData->GetValue(mDataKey, "jOutCenterOfTorusRadius2"); // distance from center of torus to z-axis
   mOutCenterOfTorusZ2 = mData->GetValue(mDataKey, "jOutCenterOfTorusZ2");
   mOutTorusZmin2 = mData->GetValue(mDataKey, "jOutTorusZmin2"); //minimum z shift from z=0 in positive z direction
   mOutTorusZmax2 = mData->GetValue(mDataKey, "jOutTorusZmax2"); //maximum z shift from z=0 in positive z direction
   mOutTorusZ0 = mData->GetValue(mDataKey, "jOutTorusZ0");
   mOutTorusTransformZ = mData->GetValue(mDataKey, "jOutTorusTransformZ");

   // Inner shape of the vessel
   mInnSegments1 = mData->GetValue(mDataKey, "jInnSegments1");
   mInnSphereRmax = mData->GetValue(mDataKey, "jInnSphereRmax");
   mInnSphereDtheta = mData->GetValue(mDataKey, "jInnSphereDtheta");
   mInnTransformZ = mData->GetValue(mDataKey, "jInnTransformZ");
   mInnTorusRadius1 = mData->GetValue(mDataKey, "jInnTorusRadius1"); // radius of small spindle torus sphere
   mInnCenterOfTorusRadius1 = mData->GetValue(mDataKey, "jInnCenterOfTorusRadius1"); // distance from center of torus to z-axis
   mInnSegments2 = mData->GetValue(mDataKey, "jInnSegments2");
   mInnTorusRadius2 = mData->GetValue(mDataKey, "jInnTorusRadius2"); // radius of large spindle torus sphere
   mInnCenterOfTorusRadius2 = mData->GetValue(mDataKey, "jInnCenterOfTorusRadius2"); // distance from center of torus to z-axis
   mInnCenterOfTorusZ2 = mData->GetValue(mDataKey, "jInnCenterOfTorusZ2");
   mInnTorusZmin2 = mData->GetValue(mDataKey, "jInnTorusZmin2"); //minimum z shift from z=0 in positive z direction
   mInnTorusZmax2 = mData->GetValue(mDataKey, "jInnTorusZmax2"); //maximum z shift from z=0 in positive z direction
   mInnTorusZ0 = mData->GetValue(mDataKey, "jInnTorusZ0");
   mInnTorusTransformZ = mData->GetValue(mDataKey, "jInnTorusTransformZ");

   // PMT distance to vessel & height of Gel
   mPmtDistance = mData->GetValue(mDataKey,"jPmtDistance");
   mGelHeight = mData->GetValue(mDataKey,"jGelHeight");
   mGelOffset = mData->GetValue(mDataKey,"jGelOffset");

   // Mainboard and HV Board Sizes and Position 
   mMainBoardRmin = mData->GetValue(mDataKey,"jMainBoardRmin");
   mMainBoardRmax =mData->GetValue(mDataKey,"jMainBoardRmax");
   mMainBoardDz = mData->GetValue(mDataKey,"jMainBoardDz");
   mMainBoardPosition = mData->GetValue(mDataKey,"jMainBoardPosition");
    
   mHVBoardRmin = mData->GetValue(mDataKey,"jHVBoardRmin");
   mHVBoardRmax =mData->GetValue(mDataKey,"jHVBoardRmax");
   mHVBoardDz = mData->GetValue(mDataKey,"jHVBoardDz");
   mHVBoardPosition = mData->GetValue(mDataKey,"jHVBoardPosition");
}
/**
 * @brief Construction of the whole DEGG. If you want to change any component, you have to change it at the specific function.
 * 
 */
void dEGG::Construction(){
   //Create pressure vessel and inner volume
   G4VSolid *lOuterGlass = CreateEggSolid(mOutSegments1,mOutSphereRadiusMax,mOutSphereDtheta,mOutTransformZ,mOutTorusRadius1,mOutCenterOfTorusRadius1,mOutSegments2,mOutTorusRadius2,mOutCenterOfTorusRadius2,mOutCenterOfTorusZ2,mOutTorusZmin2,mOutTorusZmax2,mOutTorusZ0,mOutTorusTransformZ);
   G4VSolid *lInnerGlass = CreateEggSolid(mInnSegments1,mInnSphereRmax,mInnSphereDtheta,mInnTransformZ,mInnTorusRadius1,mInnCenterOfTorusRadius1,mInnSegments2,mInnTorusRadius2,mInnCenterOfTorusRadius2,mInnCenterOfTorusZ2,mInnTorusZmin2,mInnTorusZmax2,mInnTorusZ0,mInnTorusTransformZ);
   
   //Logicals
   G4LogicalVolume* PDOM_Glass_logical = new G4LogicalVolume (lOuterGlass, mData->GetMaterial("argVesselGlass"), "Glass_phys");
   lInnerVolumeLogical = new G4LogicalVolume (lInnerGlass, mData->GetMaterial("Ri_Air"), "Glass inside");
   
   //Placements
   G4RotationMatrix *rot = new G4RotationMatrix();
   new G4PVPlacement (rot, G4ThreeVector(0,0,0), lInnerVolumeLogical, "VacuumGlass", PDOM_Glass_logical, false, 0, true);
   PlaceGel(); // placed in lInnerVolumeLogical
   PlacePMT(); // placed in the logical of the gel (see PlaceGel)
   if (gCADImport) PlaceCADSupportStructure(lInnerVolumeLogical);
   //if (gCADImport) PlaceCADPenetrator(lInnerVolumeLogical);

   // ------------------ Add outer shape solid to MultiUnion in case you need substraction -------------------------------------------
   //Each Component needs to be appended to be places in abcDetectorComponent. Everything is placed in the InnerVolume which is placed in the glass which is the mother volume. This is the reason why not everything is appended on its own
   AppendComponent(lOuterGlass, PDOM_Glass_logical, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "dEGG");
   
   // ---------------- visualisation attributes --------------------------------------------------------------------------------
   PDOM_Glass_logical->SetVisAttributes(mGlassVis);
   lInnerVolumeLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
   
}

/**
 * @brief Construction the Gel of the DEGG.
 * 
 */
void dEGG::PlaceGel(){
   // get Gel
   G4RotationMatrix *rot = new G4RotationMatrix();
   G4RotationMatrix *rot1 = new G4RotationMatrix();
   rot1->rotateY(180*deg);

   G4VSolid *degg_inner = Egg_Inner(mInnSegments1,mInnSphereRmax,mInnSphereDtheta,mInnTransformZ,mInnTorusRadius1,mInnCenterOfTorusRadius1,mInnSegments2,mInnTorusRadius2,mInnCenterOfTorusRadius2,mInnCenterOfTorusZ2,mInnTorusZmin2,mInnTorusZmax2,mInnTorusZ0,mInnTorusTransformZ);
   G4Box* box = new G4Box("Box",mOutSphereRadiusMax,mOutSphereRadiusMax, mGelHeight );
   G4ThreeVector box_trans(0,0,mGelHeight*2+mGelOffset);
   G4IntersectionSolid *intersection = new G4IntersectionSolid("EggInner * Box",degg_inner,box ,rot,box_trans);
   G4IntersectionSolid *intersection1 = new G4IntersectionSolid("EggInner * Box",degg_inner,box ,rot1,box_trans);
   
   //logicals
   lgelsolid = new G4LogicalVolume (intersection, mData->GetMaterial("argGel"), "intersection logical");
   lgelsolid1 = new G4LogicalVolume (intersection1, mData->GetMaterial("argGel"), "intersection logical");
   
   //placements
   new G4PVPlacement (rot, G4ThreeVector(0,0,0), lgelsolid, "Gel", lInnerVolumeLogical, false, 0, true);
   new G4PVPlacement (rot1, G4ThreeVector(0,0,0), lgelsolid1, "Gel1", lInnerVolumeLogical, false, 0, true);
   lgelsolid->SetVisAttributes(mGelVis);
   lgelsolid1->SetVisAttributes(mGelVis);
}
/**
 * @brief Construction the PMT of the DEGG. PMTs are placed in the logical of the Gel (see PlaceGel()).
 * 
 */
void dEGG::PlacePMT(){
   G4RotationMatrix *rot = new G4RotationMatrix();
   G4RotationMatrix *rot1 = new G4RotationMatrix();
   rot1->rotateY(180*deg);
   mPMTManager->PlaceIt(G4ThreeVector(0, 0, mPmtDistance), rot, lgelsolid, "PMT1");
   mPMTManager->PlaceIt(G4ThreeVector(0, 0, mPmtDistance), rot, lgelsolid1, "PMT2");
}
/**
 * @brief Placement of the SupportStructure (from CAD)
 * 
 */
void dEGG::PlaceCADSupportStructure(G4LogicalVolume* lInnerVolumeLogical)
{
   //select file
   std::stringstream CADfile;
   CADfile.str(""); 
   //CADfile << "Internal_Everything.obj";
   CADfile << "Internal_Everything_NoMainboard.obj";
   G4cout <<  "using the following CAD file for support structure: "  << CADfile.str()  << G4endl;

   //load mesh
   auto mesh = CADMesh::TessellatedMesh::FromOBJ("../data/CADmeshes/DEGG/" + CADfile.str() );
   G4ThreeVector CADoffset = G4ThreeVector(-427.6845*mm, 318.6396*mm, 152.89*mm); //measured from CAD file since origin =!= Module origin
   mesh->SetOffset(CADoffset);

   // Place all of the meshes it can find in the file as solids individually.
   for (auto solid : mesh->GetSolids())
   { 
      G4LogicalVolume* mSupportStructureLogical  = new G4LogicalVolume( solid , mData->GetMaterial("NoOptic_Absorber") , "logical" , 0, 0, 0);
      mSupportStructureLogical->SetVisAttributes(mAluVis);
      new G4PVPlacement( 0 , G4ThreeVector(0, 0, 0) , mSupportStructureLogical, "Support structure" , lInnerVolumeLogical, false, 0);
   }
}

/**
 * Placement of the CreateEggSolid. 
 * @param segments_1 G4int 
 * @param sphere_rmax G4double outer radius sphere
 * @param sphere_dtheta G4double delta theta angle of the segment 
 * @param sphere_transform_z G4double shift of sphere in z direction
 * @param torus1_r G4double radius of small spindle torus sphere 
 * @param centeroftorus1_r G4double distance from center of torus 1 to z-axis 
 * @param segments_2 G4int 
 * @param torus2_r G4double radius of large spindle torus sphere
 * @param centeroftorus2_r G4double distance from center of torus2_r to z-axis (signed)
 * @param centeroftorus2_z G4double distance from center of torus2_r to z-axis (signed)
 * @param torus2_zmin G4double minimum z shift from z=0 in positive z direction
 * @param torus2_zmax G4double maximum z shift from z=0 in positive z direction
 * @param torus2_z0 G4double
 * @param torus1_transform_z G4double
 * @return return the outer or inner shape of the glass vessel
 * 
 */
G4VSolid* dEGG::CreateEggSolid(G4int segments_1,
                           G4double sphere_rmax,
                           G4double sphere_dtheta,
                           G4double sphere_transform_z,
                           G4double torus1_r,
                           G4double centeroftorus1_r,
                           G4int segments_2,
                           G4double torus2_r,
                           G4double centeroftorus2_r,
                           G4double centeroftorus2_z,
                           G4double torus2_zmin,
                           G4double torus2_zmax,
                           G4double torus2_z0,
                           G4double torus1_transform_z)
{
   G4double rmin = 0;
   G4double rmax = sphere_rmax;
   G4double sphi = 0. *degree;
   G4double dphi = 2*M_PI;
   G4double dtheta = sphere_dtheta;
   G4double stheta = 0. *degree;
    
   // Create Egg sphere 
   G4Sphere* sphere1 = new G4Sphere("sphere", rmin, rmax, sphi, dphi, stheta, dtheta);
   G4VSolid* sphere = sphere1;
   G4ThreeVector centerOfSphereUp(0, 0, sphere_transform_z);

   //Torus Part 1
   //buidling small polycones, define full size of torus
   G4double step = torus1_r / segments_1;
   G4double torus1_zmax = torus1_r;
   std::vector<G4double> tempZ_1, tempOuter_1;

   G4double r;
   G4double torus_relative_zmax = torus1_zmax;
   for (G4int j=0; j<=segments_1; ++j) {
      tempZ_1.push_back(torus1_zmax - j*step);
      r = sqrt((torus1_r + torus_relative_zmax - j*step)*(torus1_r - torus_relative_zmax + j*step));
      tempOuter_1.push_back(centeroftorus1_r + r);
   }
   G4double rInner[segments_1+1], rOuter[segments_1+1], zPlane[segments_1+1];
   for (G4int i=0; i<=segments_1; ++i) {
      zPlane[i] = tempZ_1[i];
      rInner[i] = 0.;
      rOuter[i] = tempOuter_1[i];
   }

   G4Polycone * torus11 = new G4Polycone("torus1", 0, 2*M_PI, segments_1+1, zPlane, rInner, rOuter);
   G4VSolid * torus1 = torus11;

   //Torus Part 2
   // building large sphere revolution

   G4double zminrelative = torus2_zmin - centeroftorus2_z; //minimum z shift from center of torus in positive z direction
   G4double zmaxrelative = torus2_zmax - centeroftorus2_z; //maximum z shift from center of torus in positive z direction
   std::vector<G4double> tempZ2, tempOuter2;
   G4double step2 = (zmaxrelative-zminrelative)/(segments_2-1);

   G4double r2;
   G4double torus_relative_zmax2 = torus2_zmax - centeroftorus2_z;
   for (G4int j=0; j<=segments_2-1; ++j) {
      tempZ2.push_back(torus2_zmax - j*step2);
      r2 = sqrt((torus2_r + torus_relative_zmax2 - j*step2)*(torus2_r - torus_relative_zmax2 + j*step2));
      tempOuter2.push_back(centeroftorus2_r + r2);
   }

   tempZ2.push_back(0.);
   tempOuter2.push_back(torus2_z0);

   G4double rInner2[segments_2+1], rOuter2[segments_2+1], zPlane2[segments_2+1];
   for (G4int i=0; i<=segments_2; i++) {
      rInner2[i] = 0;
      rOuter2[i] = tempOuter2[i];
      zPlane2[i] = tempZ2[i];
      }

   G4Polycone * torus21 = new G4Polycone("polycone2", 0, 2*M_PI, segments_2+1, zPlane2, rInner2, rOuter2);
   G4VSolid * torus2 = torus21;

   //Create Vessel

   G4ThreeVector centerOfPolycone1(0, 0, torus1_transform_z); 

   G4UnionSolid *solid1= new G4UnionSolid("solid1", torus2, torus1, 0, centerOfPolycone1);
    
   G4UnionSolid *solid= new G4UnionSolid("solid", solid1, sphere, 0, centerOfSphereUp);
    
   G4VSolid * deggup = solid;
   G4VSolid * deggdown = solid;

   G4RotationMatrix *rot = new G4RotationMatrix();
   rot->rotateY(180.0*deg);
   G4UnionSolid * degg = new G4UnionSolid("degg", deggup, deggdown, rot, G4ThreeVector());

   return degg;
}


G4VSolid* dEGG::Egg_Inner(G4int segments_1,
                            G4double sphere_rmax,
                            G4double sphere_dtheta,
                            G4double sphere_transform_z,
                            G4double torus1_r,
                            G4double centeroftorus1_r,
                            G4int segments_2,
                            G4double torus2_r,
                            G4double centeroftorus2_r,
                            G4double centeroftorus2_z,
                            G4double torus2_zmin,
                            G4double torus2_zmax,
                            G4double torus2_z0,
                            G4double torus1_transform_z)
{   G4double rmin = 0;
    G4double rmax = sphere_rmax;
    G4double sphi = 0. *degree;
    G4double dphi = 2*M_PI;
    G4double dtheta = sphere_dtheta;
    G4double stheta = 0. *degree;
    
    // Create Egg sphere for top part
    G4Sphere* sphere1 = new G4Sphere("sphere", rmin, rmax, sphi, dphi, stheta, dtheta);
    G4VSolid* sphere = sphere1;
    
    G4ThreeVector centerOfSphereUp(0, 0, sphere_transform_z);

    //Torus Part 1
    //buidling small polycones, define full size of torus
    G4double step = torus1_r / segments_1;
    G4double torus1_zmax = torus1_r;

    std::vector<G4double> tempZ_1, tempOuter_1;

    G4double r;
    G4double torus_relative_zmax = torus1_zmax;
    for (G4int j=0; j<=segments_1; ++j) {
       tempZ_1.push_back(torus1_zmax - j*step);
       r = sqrt((torus1_r + torus_relative_zmax - j*step)*(torus1_r - torus_relative_zmax + j*step));
       tempOuter_1.push_back(centeroftorus1_r + r);
    }
    G4double rInner[segments_1+1], rOuter[segments_1+1], zPlane[segments_1+1];
    for (G4int i=0; i<=segments_1; ++i) {
       zPlane[i] = tempZ_1[i];
       rInner[i] = 0.;
       rOuter[i] = tempOuter_1[i];
    }

    G4Polycone * torus11 = new G4Polycone("torus1", 0, 2*M_PI, segments_1+1, zPlane, rInner, rOuter);
    G4VSolid * torus1 = torus11;

    //Torus Part 2
    // building large sphere revolution

    G4double zminrelative = torus2_zmin - centeroftorus2_z; //minimum z shift from center of torus in positive z direction
    G4double zmaxrelative = torus2_zmax - centeroftorus2_z; //maximum z shift from center of torus in positive z direction

    std::vector<G4double> tempZ2, tempOuter2;
    G4double step2 = (zmaxrelative-zminrelative)/(segments_2-1);

    G4double r2;
    G4double torus_relative_zmax2 = torus2_zmax - centeroftorus2_z;
    for (G4int j=0; j<=segments_2-1; ++j) {
       tempZ2.push_back(torus2_zmax - j*step);
       // r = sqrt(torus_r*torus_r - (torus_relative_zmax-j*step)*(torus_relative_zmax-j*step))
       r2 = sqrt((torus2_r + torus_relative_zmax2 - j*step2)*(torus2_r - torus_relative_zmax2 + j*step2));
       tempOuter2.push_back(centeroftorus2_r + r2);
    }


    tempZ2.push_back(0.);
    tempOuter2.push_back(torus2_z0);

    G4double rInner2[segments_2+1], rOuter2[segments_2+1], zPlane2[segments_2+1];
    for (G4int i=0; i<=segments_2; i++) {
       rInner2[i] = 0;
       rOuter2[i] = tempOuter2[i];
       zPlane2[i] = tempZ2[i];
       //std::cout<<"EggTorus2 "<<i<<" "<<zPlane2[i]<<" "<<rInner2[i]<<" "<<rOuter2[i]<<std::endl;
    }
    G4Polycone * torus21 = new G4Polycone("polycone2", 0, 2*M_PI, segments_2+1, zPlane2, rInner2, rOuter2);
    G4VSolid * torus2 = torus21;

    G4ThreeVector centerOfPolycone1(0, 0, torus1_transform_z); 

    G4UnionSolid *solid1= new G4UnionSolid("solid1", torus2, torus1, 0, centerOfPolycone1);
    
    G4UnionSolid *solid= new G4UnionSolid("solid", solid1, sphere, 0, centerOfSphereUp);
    
    return solid;

 }