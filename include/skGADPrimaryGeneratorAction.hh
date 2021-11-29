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

#ifndef skGADPrimaryGeneratorAction_h
#define skGADPrimaryGeneratorAction_h 1

#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class skGADPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
	public:
		skGADPrimaryGeneratorAction();
		~skGADPrimaryGeneratorAction();

	public:
		G4String      rndmFlag;
		G4double      rndmRadius;
		G4bool        nigenFlag;
		G4ThreeVector nigenPos;

		void GeneratePrimaries(G4Event* anEvent);
		void SetRndmFlag(G4String val)   { rndmFlag   = val; };
		void SetRndmRadius(G4double val) { rndmRadius = val; };
		void SetNiGenFlag(G4bool val)    { nigenFlag  = val; };
		void SetNiPos(G4ThreeVector val) { nigenPos   = val; };
                G4ParticleGun* GetParticleGun();

	private:
		G4ParticleGun* particleGun;
		//G4GeneralParticleSource*      particleGun;

		G4int iBin;
		G4double neutrinoEne, positronEne, neutronEne;
		G4double neutrinoDirX, positronDirX, neutronDirX;
		G4double neutrinoDirY, positronDirY, neutronDirY;
		G4double neutrinoDirZ, positronDirZ, neutronDirZ;
		G4double vertexPosX, vertexPosY, vertexPosZ; // position of input neutrino


		G4int evID;
};

#endif
