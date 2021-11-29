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
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4String.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"
//#include "Random.hh"
//#include "Randomize.hh"

//---------------------------
// Parameterisation manager:
//---------------------------
#include "G4GlobalFastSimulationManager.hh"

//------------
// Geometries:
//------------

//#ifdef G4VIS_USE
//#include "G4VisExecutive.hh"
//#endif

#include "skGADDetectorConstruction.hh"
#include "skGADPrimaryGeneratorAction.hh"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4HadronCaptureProcess.hh"
//#include "myG4NeutronHPorLCapture.hh"
#include "G4NeutronHPorLCapture.hh"
#include "GdNeutronHPCapture.hh"
//#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4CrossSectionDataStore.hh"

#include "G4NeutronInelasticProcess.hh"
//#include "myG4NeutronHPorLEInelastic.hh"
#include "G4NeutronHPorLEInelastic.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"

#include "G4HadronElasticProcess.hh"
//#include "myG4NeutronHPorLElastic.hh"
#include "G4NeutronHPorLElastic.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"

#include "G4HadronicAbsorptionBertini.hh"
#include "G4MuonMinusCapture.hh"

#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Step.hh"

#ifndef userInterface2_h
#define userInterface2_h 1

using namespace CLHEP;

struct Secondary_Vector
{
  G4int evId;
  G4double productionTime; //ns
  G4double productionPosX;
  G4double productionPosY;
  G4double productionPosZ;
  G4double momentumX;
  G4double momentumY;
  G4double momentumZ;
  G4double totalEnergy;
  G4long pdgID;
};

class userInterface2
{
public:
  //userInterface2();
  userInterface2( G4int , G4int , G4int , G4float );
  ~userInterface2();
  
  void Initialize();
  void Finalize();

  // [Comments by Y. Ashida]
  //G4double GetCapXsec(G4double nEnergy);
  //G4double GetElaXsec(G4double nEnergy);
  //G4double GetInElaXsec(G4double nEnergy);
  // for Gadolinium
  G4double GetGdCapXsec(G4double nEnergy);
  G4double GetGdElaXsec(G4double nEnergy);
  G4double GetGdInElaXsec(G4double nEnergy);
  // for Hydrogen
  G4double GetHCapXsec(G4double nEnergy);
  G4double GetHElaXsec(G4double nEnergy);
  G4double GetHInElaXsec(G4double nEnergy);
  // for Oxygen
  G4double GetOCapXsec(G4double nEnergy);
  G4double GetOElaXsec(G4double nEnergy);
  G4double GetOInElaXsec(G4double nEnergy);

  //std::vector< Secondary_Vector > GetCap2nd(G4double px,G4double py,G4double pz);
  //std::vector< Secondary_Vector > GetEla2nd(G4double px,G4double py,G4double pz);
  //std::vector< Secondary_Vector > GetInEla2nd(G4double px,G4double py,G4double pz);
  std::vector< Secondary_Vector > Get2nd(G4double px,G4double py,G4double pz,G4int process);
  G4HadFinalState * CheckResult(const G4HadProjectile & aPro,
		const G4Nucleus & aNucleus,
		G4HadFinalState * result,G4HadronicInteraction * hadModel);

  std::vector< Secondary_Vector > DoPionMinusCapture(G4double, G4double, G4double);
  std::vector< Secondary_Vector > DoMuonMinusCapture(G4double, G4double, G4double);
private:
  skGADDetectorConstruction* det;
  G4RunManager* runManager;
  G4UImanager* UI;

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
  G4Material* GdSulfateOctahydrate;
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

  G4HadronCaptureProcess* theCaptureProcess;
  //myG4NeutronHPorLCapture * theLENeutronCaptureModel;  // to load Gd only 
  GdNeutronHPCapture * theLENeutronCaptureModel;
  //G4NeutronHPorLCapture * theLENeutronCaptureModel;
  G4NeutronHPCaptureData * theNeutronData3;
  G4CrossSectionDataStore* capProcessXsecDataStore;

  G4HadronElasticProcess* theNeutronElasticProcess;
  //myG4NeutronHPorLElastic * theElasticNeutron;  // to load Gd only
  G4NeutronHPorLElastic * theElasticNeutron;
  //G4NeutronHPElastic * theElasticNeutron;
  G4NeutronHPElasticData * theNeutronData;
  G4CrossSectionDataStore* elaProcessXsecDataStore;

  G4NeutronInelasticProcess* theInelasticProcess;
  //myG4NeutronHPorLEInelastic* theInelasticModel;  // to load Gd only
  G4NeutronHPorLEInelastic* theInelasticModel;
  //G4NeutronHPInelastic* theInelasticModel;
  G4NeutronHPInelasticData * theNeutronData1;
  G4CrossSectionDataStore* inelaProcessXsecDataStore;

  G4ParticleDefinition* particle;
  G4DynamicParticle* neutronParticle;
  G4Step* aStep;

  G4double nTotalAtomsPerVolume;
  G4double nHAtomsPerVolume;
  G4double nOAtomsPerVolume;
  G4double nGdAtomsPerVolume;
  //
  G4double rateHtoTotal;  
  G4double rateOtoTotal;  
  G4double rateGdtoTotal;  


  //G4VisManager* visManager;
  G4Step							*aStepMuMi=NULL;
  G4Step							*aStepPiMi=NULL;

  G4Track* aTrack;
  
  G4Step *PiMiStep;
	G4StepPoint* PiMiStepPoint;
  G4Track* aTrackPiMi;
  G4DynamicParticle			*piminusParticle2;

  G4DynamicParticle			*muminusParticle;
	G4Step *MuMiStep;
	G4StepPoint* MuMiStepPoint;
  G4Track* aTrackMuMi;

  G4HadFinalState* result;
};

#endif

