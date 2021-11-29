

///////////////////////////////////////////////////////////////////////////////
//                    Spectrum of radiative neutron capture by Gadolinium            
//                                    version 1.0.0                                
//                                    (Sep.09.2005)                               
                            
//                Author : karim.zbiri@subatech.in2p3.fr                  

//Modified class from original G4NeutronHPCaptureFS class to deexcite and
//add correctly the secondary to the hadronic final state

// Karim Zbiri, Aug, 2005
///////////////////////////////////////////////////////////////////////////////



#ifndef GdNeutronHPCaptureFS_h
#define GdNeutronHPCaptureFS_h 1

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4ReactionProductVector.hh"
#include "G4NeutronHPNames.hh"
#include "G4NeutronHPPhotonDist.hh"
#include "G4Nucleus.hh"
#include "G4Fragment.hh"
#include "G4Allocator.hh"
#include "GdCaptureGammas.hh"
#include "SKGdCaptureGammas.hh"

#include "ANNRIGd_GdNCaptureGammaGenerator.hh"
#include "ANNRIGd_GeneratorConfigurator.hh"
#include "ANNRIGd_OutputConverter.hh"

namespace ANNRIGdGammaSpecModel { class ANNRIGd_GdNCaptureGammaGenerator; }

class GdNeutronHPCaptureFS : public G4NeutronHPFinalState
{
  public:
  
  GdNeutronHPCaptureFS();
  
  ~GdNeutronHPCaptureFS();
  
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aFSType)
  { return Init(A, Z, dirName, aFSType); }
  void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aFSType, G4ParticleDefinition* )
  { return Init(A, Z, dirName, aFSType); }


  G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack);
  G4NeutronHPFinalState * New() 
  {
   GdNeutronHPCaptureFS * theNew = new GdNeutronHPCaptureFS;
   return theNew;
  }

  void   UpdateNucleus( const G4Fragment* , G4double );
  void InitANNRIGdGenerator();
  G4ReactionProductVector* GenerateWithANNRIGdGenerator();
  
  private:

  G4int Gd_CAPTURE;
  G4int Gd_CASCADE;
  G4int fMOEDEL;
  G4Fragment * nucleus;

  G4DynamicParticle * theTwo ;
  //G4ReactionProduct theTarget; 
  //G4Nucleus aNucleus;
  //G4ReactionProduct theNeutron;

  G4double targetMass;
  
  GdCaptureGammas theGdFinalGammas;
  SKGdCaptureGammas theSKGdFinalGammas;

  G4NeutronHPPhotonDist theFinalStatePhotons;
  G4NeutronHPNames theNames;
  
  G4double theCurrentA;
  G4double theCurrentZ;

  private:
  static ANNRIGdGammaSpecModel::ANNRIGd_GdNCaptureGammaGenerator* ANNRIGammaGen;
  G4String Gd155_ROOTFile;
  G4String Gd157_ROOTFile;

};
#endif
