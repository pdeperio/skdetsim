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
#include "skGADPrimaryGeneratorAction.hh"
#include "G4ReactionProductVector.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GeneralParticleSource.hh" // 0413
#include "G4Gamma.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

#include "globals.hh"
#include "Randomize.hh"
#include <iostream>
#include <CLHEP/Units/SystemOfUnits.h>

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
skGADPrimaryGeneratorAction::skGADPrimaryGeneratorAction(){
	
	   G4int n_particle = 1;
	   particleGun = new G4ParticleGun(n_particle);

	   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	   

	//nigenFlag = false;
	//create a messenger for this class
	//gunMessenger = new skGADPrimaryGeneratorMessenger(this);


	
	G4String particleName;
	particleGun->SetParticleDefinition(particleTable->
	FindParticle(particleName="neutron"));

	//-thremal, default
	particleGun->SetParticleMomentum(2.*CLHEP::MeV);
	particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0 *CLHEP::cm));
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
	

	//particleGun = new G4GeneralParticleSource (); //0413

}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
skGADPrimaryGeneratorAction::~skGADPrimaryGeneratorAction(){
	delete particleGun;
}

G4ParticleGun* skGADPrimaryGeneratorAction::GetParticleGun()
{
  return particleGun;
} 

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void skGADPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){

  //G4double position = -0.5*(myDetector->GetWorldFullLength());
  //particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,1.*cm));
  
  particleGun->GeneratePrimaryVertex(anEvent);

	//-org  particleGun->GeneratePrimaryVertex(anEvent);
	/*
	   G4ThreeVector position(0.,0.,0.);
	   G4ThreeVector direction(0.,0.,0.);

	//-ak-0430----------------
	evID = anEvent->GetEventID();
	*/
	/*
	   ReadGeneratorFile();
	   skGADGeneParams* skGadGeneParams = skGADGeneParams::getInstance();

	   if( evID > skGadGeneParams->getNumOfParamBins() ) return;

	   G4double neutronEne = skGadGeneParams->getGeneParamsNeutronEne( evID );
	   G4double neutronDirX = skGadGeneParams->getGeneParamsNeutronDir( evID ).x();
	   G4double neutronDirY = skGadGeneParams->getGeneParamsNeutronDir( evID ).y();
	   G4double neutronDirZ = skGadGeneParams->getGeneParamsNeutronDir( evID ).z();

	   G4double neutronPosX = skGadGeneParams->getGeneParamsVertexPos( evID ).x();
	   G4double neutronPosY = skGadGeneParams->getGeneParamsVertexPos( evID ).y();
	   G4double neutronPosZ = skGadGeneParams->getGeneParamsVertexPos( evID ).z();

	   if (rndmFlag == "file"){

	   particleGun->SetParticleEnergy( neutronEne );

	   position.setX( neutronPosX );
	   position.setY( neutronPosY );
	   position.setZ( neutronPosZ );
	   particleGun -> SetParticlePosition(position);

	   direction.setX( neutronDirX );
	   direction.setY( neutronDirY );
	   direction.setZ( neutronDirZ );
	   particleGun -> SetParticleMomentumDirection(direction);

	   }
	   */
	//-ak-0430----------------

	/*
	   if (rndmFlag == "on"){
	   skGADDetectorConstruction* constructor = new skGADDetectorConstruction();
	   G4double scint_diameter = constructor->GetScintDiameter();
	   G4double scint_height = constructor->GetScintHeight();

	//  Random generation of gamma source point inside the Iodium core ...
	G4double x,y,z;

	do{
	x = (G4UniformRand()-0.5)*scint_diameter;
	y = (G4UniformRand()-0.5)*scint_diameter;
	}while( (x*x+y*y) > (scint_diameter*scint_diameter/4) );

	z = (G4UniformRand()-0.5)*scint_height ;

	position.setX(x);
	position.setY(y);
	position.setZ(z);
	particleGun -> SetParticlePosition(position);

	// Random generation of the impulse direction of primary particles ...
	G4double a,b,c;
	G4double n;
	do{
	a = (G4UniformRand()-0.5)/0.5;
	b = (G4UniformRand()-0.5)/0.5;
	c = (G4UniformRand()-0.5)/0.5;
	n = a*a+b*b+c*c;
	}while(n > 1 || n == 0.0);
	n = std::sqrt(n);
	a /= n;
	b /= n;
	c /= n;

	direction.setX(a);
	direction.setY(b);
	direction.setZ(c);
	particleGun -> SetParticleMomentumDirection(direction);

	}
	*/



	/*if( nigenFlag == true )
	{
		SetUpNickelParticleGun(anEvent);
	}
	else
	{
		particleGun->GeneratePrimaryVertex(anEvent);
	}*/

	/*
	   neutronEne = particleGun->GetParticleMomentum(); //-ak0412

	   G4cout <<  " gun-pos1 " << neutronEne /MeV << " " 
	   << direction.x() /cm <<" "<< direction.y() /cm << " " << direction.z() /cm << " " 
	   << position.x() /cm <<" "<< position.y() /cm << " " << position.z() /cm << " " 
	   << G4endl;      

	//-ak-0413 generate Rn
	G4ParticleDefinition* particleDef = particleGunGPS->GetParticleDefinition();
	G4cout <<  " gun-pos2 " 
	<< particleDef->GetParticleName() << " " 
	//	 << neutronEne /MeV << " " 
	<< particleGunGPS->GetParticlePosition().x() /cm <<" "
	<< particleGunGPS->GetParticlePosition().y() /cm <<" "
	<< particleGunGPS->GetParticlePosition().z() /cm <<" "
	//  	 << position.x() /cm <<" "<< position.y() /cm << " " << position.z() /cm << " " 
	<< G4endl;      

*/

}


