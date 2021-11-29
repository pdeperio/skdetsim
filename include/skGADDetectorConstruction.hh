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
#ifndef skGADDetectorConstruction_H
#define skGADDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Box;
class G4Trd;
class G4Tubs;
class skGADPMTSD;
class skGADHPDSD;
class skGADHQESD;
class skGADScintSD;
class G4Sphere;

#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"

class skGADDetectorConstruction : public G4VUserDetectorConstruction
{

  public:

    skGADDetectorConstruction();
    ~skGADDetectorConstruction();

    G4VPhysicalVolume* Construct();



  //  G4ThreeVector GetPosOfPMT( G4int ipmt );

  G4Material* GetWater() {return Water;};
  G4Material* GetGdWater() {return GdWater;};

private:

  void DefineMaterials();
  G4VPhysicalVolume* ConstructDetector();

  G4bool updated;
   
  G4Tubs* experimentalHall_box;
  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;

  G4Sphere* waterSphere_box;
  G4LogicalVolume* waterSphere_log;
  G4VPhysicalVolume* waterSphere_phys;

  G4Tubs*     gdTube_box;
  G4LogicalVolume*     gdTube_log;
  G4VPhysicalVolume*     gdTube_phys;

  G4Box*     BGO_box;
  G4LogicalVolume*     BGO_log;
  G4VPhysicalVolume*     BGO_phys;

  G4Trd*     INSPC_box;
  G4LogicalVolume*     INSPC_log;
  G4VPhysicalVolume*     INSPC_phys;


  G4Tubs*     AmBe_box;
  G4LogicalVolume*     AmBe_log;
  G4VPhysicalVolume*     AmBe_phys;

  G4VPhysicalVolume* world_sphere_phys;
  //Materials & Elements
  G4Element* N;
  G4Element* O;
  G4Element* C;
  G4Element* H;
  G4Element* B;
  G4Element* Si;
  G4Element* S;
  G4Element* Na;
  G4Element* Al;
  G4Element* Cl;
  G4Element* Fe;
  G4Element* Co;
  G4Element* Gd;

  G4Element* Ge;
  G4Element* Bi;

  G4Material* BGO;
  G4Material* GdCl3;
  G4Material* Gd2S3O12;
  G4Material* Water;
  G4Material* GdWater;
  G4Material* GdWater2005;
  G4Material* GdWater2000;
  G4Material* GdWater250;
  G4Material* GdWater235;
  G4Material* GdWater230;
  G4Material* GdWater120;
  G4Material* GdWater115;
  G4Material* GdWater90;
  G4Material* Air;
  G4Material* Vacuum;
  G4Material* Iron;
  G4Material* Aluminum;
  G4Material* Glass;
  G4Material* SiO2;
  G4Material* B2O3;
  G4Material* Na2O;
  G4Material* Al2O3;
  G4Material* Pstyrene;
  G4Material* PMMA;
  G4Material* Pethylene;
  G4Material* fPethylene; 
  G4Material* PET; 
  G4Material* sSteel; 

  //Geometry
  G4double scint_diameter;
  G4double scint_height;
  G4double bsheet_thickness;
  G4double d_mtl;
  G4double diameter_pmt;
  G4double spacing_pmt;
  G4double curveture_radius_pmt;
  G4double glass_thickness_pmt;
  G4double diameter_hpd;
  G4double curveture_radius_hpd;
  G4double glass_thickness_hpd;
  G4double diameter_hqe;
  G4double curveture_radius_hqe;
  G4double glass_thickness_hqe;

  G4double inner_curve_r_acryl;
  G4double thickness_acryl;
  G4double position_offset_acryl;

  //  G4int nfibers;
  static G4bool sphereOn;
  G4double refl;
  //  G4bool WLSslab;
  G4bool mainVolume;
  //  G4double slab_z;

  G4MaterialPropertiesTable* WATER_mt;
  G4MaterialPropertiesTable* MPTPStyrene;
  G4MaterialPropertiesTable* BGO_mt;

  void MaxNoOfPMTs();

  G4int AmBe_type;
  G4int Water_status;

};

#endif




