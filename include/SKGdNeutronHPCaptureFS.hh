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
//
//
#ifndef SKGdNeutronHPCaptureFS_h
#define SKGdNeutronHPCaptureFS_h 1

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4ReactionProductVector.hh"
#include "G4NeutronHPNames.hh"
#include "G4NeutronHPPhotonDist.hh"
#include "G4NeutronHPEnAngCorrelation.hh"
#include "SKGdCaptureGammas.hh"
#include "GdCaptureGammas.hh"
#include "G4Allocator.hh"
class SKGdNeutronHPCaptureFS : public G4NeutronHPFinalState
{
  public:
  
  SKGdNeutronHPCaptureFS()
  {
    hasXsec = false; 
    hasExactMF6 = false;
    targetMass = 0;
    G4cout<<"Initial FS"<<G4endl;
    //theFinalGammas = new SKGdCaptureGammas(); 
  }
  
  ~SKGdNeutronHPCaptureFS()
  {
  }
  
  void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aFSType);
  void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aFSType, G4ParticleDefinition* )
  {
    return Init(A, Z, M, dirName, aFSType);
  }

  G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack);
  G4NeutronHPFinalState * New() 
  {
   SKGdNeutronHPCaptureFS * theNew = new SKGdNeutronHPCaptureFS;
   return theNew;
  }
  
  private:
  
  G4double targetMass;
  
  G4NeutronHPPhotonDist theFinalStatePhotons;

  SKGdCaptureGammas theFinalGammas;

   G4NeutronHPEnAngCorrelation theMF6FinalState;
   G4bool hasExactMF6;
  
  G4NeutronHPNames theNames;
  
//  G4double theCurrentA;
//  G4double theCurrentZ;
};
#endif
