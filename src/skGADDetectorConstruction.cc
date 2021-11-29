//
// ********************************************************************
// * License and Disclaimer						      *
// *									      *
// * The	Geant4 software	is	copyright of the Copyright Holders of *
// * the Geant4 Collaboration.	It is provided	under	the terms and *
// * conditions of the Geant4 Software License,	included in the file  *
// * LICENSE and available at	http://cern.ch/geant4/license .	These *
// * include a list of copyright holders.*
// * *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work	make any representation or	warranty, express or implied, *
// * regarding this	software system or assume any liability for its   *
// * use.	Please see the license in the file	LICENSE	and URL above *
// * for the full disclaimer and the limitation of liability.		  *
// * *
// * This code implementation is the result of	the	scientific and    *
// * technical work of the GEANT4 collaboration.			     	  *
// * By using,	copying,	modifying or	distributing the software (or *
// * any work based	on the software)	you	agree	to acknowledge its *
// * use	in	resulting	scientific	publications,	and indicate your *
// * acceptance of all terms of the Geant4 Software license.   *
// ********************************************************************
//

#include "skGADDetectorConstruction.hh"

#include "G4String.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4UImanager.hh"

using namespace CLHEP;

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
	skGADDetectorConstruction::skGADDetectorConstruction()
: WATER_mt(NULL), MPTPStyrene(NULL), BGO_mt(NULL)
{


}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
skGADDetectorConstruction::~skGADDetectorConstruction(){
}



//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void skGADDetectorConstruction::DefineMaterials(){
  G4double a;	// atomic mass
  G4double z;	// atomic number
  G4double density;
  
  G4int iz, n;
  G4int ncomponents, natoms;
  G4double abundance, fractionmass;
  
  G4int polyPMMA = 1;
  G4int nC_PMMA  = 3+2*polyPMMA;
  G4int nH_PMMA  = 6+2*polyPMMA;
  
  G4int polyeth = 1;
  G4int nC_eth  = 2*polyeth;
  G4int nH_eth  = 4*polyeth;
  
  //***Elements
  H  = new G4Element("Hydrogen",  "H",  z=1.,  a=1.01*g/mole);
  O  = new G4Element("Oxygen" ,   "O",  z=8.,  a=16.00*g/mole);

  Water = new G4Material("Water_", density= 1.0*g/cm3, ncomponents=2);
  Water->AddElement(H, natoms=2);
  Water->AddElement(O, natoms=1);
  
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4VPhysicalVolume* skGADDetectorConstruction::Construct(){
  DefineMaterials();
  return ConstructDetector();
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4VPhysicalVolume* skGADDetectorConstruction::ConstructDetector()
{
  G4double expHall_sphi = 0 *rad;
  G4double expHall_dphi = 2.0*M_PI *rad;
  
  // Create experimental hall
  experimentalHall_box
    = new G4Tubs("expHall_box", 0, 20.*m, 20.*m, expHall_sphi, expHall_dphi);
  experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box, Water, "expHall_log", 0,0,0);
  experimentalHall_phys
    = new G4PVPlacement(0, G4ThreeVector(), experimentalHall_log, "expHall",0,false,0);
  experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  return experimentalHall_phys;
}

