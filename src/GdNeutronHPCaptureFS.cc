

///////////////////////////////////////////////////////////////////////////////
//                   Spectrum of radiative neutron capture by Gadolinium           
//                                    version 1.0.0                                
//                                    (Sep.09.2005)                               

//                Author : karim.zbiri@subatech.in2p3.fr                  

//Modified class from original G4NeutronHPCaptureFS class to deexcite and
//add correctly the secondary to the hadronic final state

// Karim Zbiri, Aug, 2005
///////////////////////////////////////////////////////////////////////////////



#include "GdNeutronHPCaptureFS.hh"
#include "G4Gamma.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4PhotonEvaporation.hh"
#include "G4Fragment.hh"
#include "G4ParticleTable.hh" 
#include "G4IonTable.hh" 
#include "G4NeutronHPDataUsed.hh"
#include "userInterface2.hh"


ANNRIGdGammaSpecModel::ANNRIGd_GdNCaptureGammaGenerator* GdNeutronHPCaptureFS::ANNRIGammaGen = 0;

extern G4int fMODEL;

GdNeutronHPCaptureFS::GdNeutronHPCaptureFS()
:nucleus(0), theTwo(0), targetMass(0) 
{
  hasXsec = false; 
  targetMass = 0;

  if(fMODEL == 3) InitANNRIGdGenerator();
}

GdNeutronHPCaptureFS::~GdNeutronHPCaptureFS()
{}


G4HadFinalState * GdNeutronHPCaptureFS::ApplyYourself(const G4HadProjectile & theTrack)
{
  G4int i;
  if ( theResult.Get() == NULL ) theResult.Put( new G4HadFinalState );
  theResult.Get()->Clear();
  //theResult.Clear();
  // prepare neutron
  G4double eKinetic = theTrack.GetKineticEnergy();
  const G4HadProjectile *incidentParticle = &theTrack;
  G4ReactionProduct theNeutron = const_cast<G4ParticleDefinition *>(incidentParticle->GetDefinition()) ;
  theNeutron.SetMomentum( incidentParticle->Get4Momentum().vect() );
  theNeutron.SetKineticEnergy( eKinetic );

  // prepare target
  G4ReactionProduct theTarget; 
  G4Nucleus aNucleus;
  G4double eps = 0.0001;
  if(targetMass<500*CLHEP::MeV)
    targetMass = ( G4NucleiProperties::GetNuclearMass(static_cast<G4int>(theBaseA+eps),static_cast<G4int>(theBaseZ+eps))) /
      G4Neutron::Neutron()->GetPDGMass();
  G4ThreeVector neutronVelocity = 1./G4Neutron::Neutron()->GetPDGMass()*theNeutron.GetMomentum();
  G4double temperature = theTrack.GetMaterial()->GetTemperature();
  theTarget = aNucleus.GetBiasedThermalNucleus(targetMass, neutronVelocity, temperature);
  theTarget.SetDefinition( G4IonTable::GetIonTable()->GetIon( G4int(theBaseZ), G4int(theBaseA) , 0.0 ) );  //TESTPHP

  // go to nucleus rest system
  theNeutron.Lorentz(theNeutron, -1*theTarget);
  eKinetic = theNeutron.GetKineticEnergy();

  // dice the photons
  G4ReactionProductVector * thePhotons = NULL;    
  if(fMODEL==1){
    thePhotons = theGdFinalGammas.GetGammas();
  }else if(fMODEL==2){
    thePhotons = theSKGdFinalGammas.GetGammas();
  }else if(fMODEL==3){
    thePhotons = GenerateWithANNRIGdGenerator();
  }else{
    G4cout<<" This model number (="<<fMODEL<<") is not supported now"<<G4endl;
    G4cout<<" Please select in 0-3 "<<G4endl;
  }

  // update the nucleus
  /*G4ThreeVector aCMSMomentum = theNeutron.GetMomentum()+theTarget.GetMomentum();
  G4LorentzVector p4(aCMSMomentum, theTarget.GetTotalEnergy() + theNeutron.GetTotalEnergy());
  nucleus = new G4Fragment(static_cast<G4int>(theBaseA+1), static_cast<G4int>(theBaseZ) ,p4);
  */
  G4int nPhotons = 0;
  if(thePhotons!=NULL) nPhotons=thePhotons->size();

  if ( DoNotAdjustFinalState() ) {
//Make at least one photon  
//101203 TK
    if ( nPhotons == 0 )
    {
       G4ReactionProduct * theOne = new G4ReactionProduct;
       theOne->SetDefinition( G4Gamma::Gamma() ); 
       G4double theta = pi*G4UniformRand();
       G4double phi = twopi*G4UniformRand();
       G4double sinth = std::sin(theta);
       G4ThreeVector direction( sinth*std::cos(phi), sinth*std::sin(phi), std::cos(theta) );
       theOne->SetMomentum( direction ) ;
       thePhotons->push_back(theOne);
       nPhotons++; // 0 -> 1
    }
//One photon case: energy set to Q-value 
//101203 TK
    //if ( nPhotons == 1 )
    if ( nPhotons == 1 && thePhotons->operator[](0)->GetDefinition()->GetBaryonNumber() == 0 )
    {
       G4ThreeVector direction = thePhotons->operator[](0)->GetMomentum().unit();

       G4double Q = G4IonTable::GetIonTable()->GetIonMass(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA), 0) + G4Neutron::Neutron()->GetPDGMass()
         - G4IonTable::GetIonTable()->GetIonMass(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA+1), 0);

       thePhotons->operator[](0)->SetMomentum( Q*direction );
    } 
//
  }

  /*for(i=0;i<nPhotons;i++)
  {
    G4Fragment * theOne;
    G4ThreeVector pGamma(thePhotons->operator[](i)->GetMomentum());
    G4LorentzVector gamma(pGamma,thePhotons->operator[](i)->GetTotalEnergy());
    theOne= new G4Fragment(gamma,G4Gamma::GammaDefinition());
    UpdateNucleus(theOne,thePhotons->operator[](i)->GetTotalEnergy());
  }
  */

  /*theTwo = new G4DynamicParticle;
  theTwo->SetDefinition(G4IonTable::GetIonTable()
      ->GetIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA+1), 0));
  theTwo->SetMomentum(nucleus->GetMomentum());
  */
  // add them to the final state

  //G4int nParticles = nPhotons;
  //if(1==nPhotons) nParticles = 2;

  // back to lab system
  for(i=0; i<nPhotons; i++)
  {
    thePhotons->operator[](i)->Lorentz(*(thePhotons->operator[](i)), theTarget);
  }
  // Recoil, if only one gamma
  if (1==nPhotons && thePhotons->operator[](0)->GetDefinition()->GetBaryonNumber() == 0 )
  {
    G4DynamicParticle * theOne = new G4DynamicParticle;
    //G4ParticleDefinition * aRecoil = G4IonTable::GetIonTable()
    //  ->FindIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA+1), 0, static_cast<G4int>(theBaseZ));
    G4ParticleDefinition * aRecoil = G4IonTable::GetIonTable()
                                        ->GetIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA+1), 0);
    theOne->SetDefinition(aRecoil);
    // Now energy; 
    // Can be done slightly better @
    G4ThreeVector aMomentum =  theTrack.Get4Momentum().vect()
                              +theTarget.GetMomentum()
                              -thePhotons->operator[](0)->GetMomentum();

    /*G4ThreeVector theMomUnit = aMomentum.unit();
    G4double aKinEnergy =  theTrack.GetKineticEnergy()
      +theTarget.GetKineticEnergy(); // gammas come from Q-value
    G4double theResMass = aRecoil->GetPDGMass();
    G4double theResE = aRecoil->GetPDGMass()+aKinEnergy;
    G4double theAbsMom = std::sqrt(theResE*theResE - theResMass*theResMass);
    G4ThreeVector theMomentum = theAbsMom*theMomUnit;
    theOne->SetMomentum(theMomentum);*/
    theOne->SetMomentum(aMomentum);
    theResult.Get()->AddSecondary(theOne);
  }
  // Now fill in the gammas.

  for(i=0; i<nPhotons; i++)
  {
    // back to lab system
    G4DynamicParticle * theOne = new G4DynamicParticle;
    theOne->SetDefinition(thePhotons->operator[](i)->GetDefinition());
    theOne->SetMomentum(thePhotons->operator[](i)->GetMomentum());
    theResult.Get()->AddSecondary(theOne);
    delete thePhotons->operator[](i);
  }
  delete thePhotons; 

  //ADD deexcited nucleus

  //theResult.Get()->AddSecondary(theTwo);

//101203TK
  G4bool residual = false;
  G4ParticleDefinition * aRecoil = G4IonTable::GetIonTable()
                                   ->GetIon(static_cast<G4int>(theBaseZ), static_cast<G4int>(theBaseA+1), 0);
  for ( G4int j = 0 ; j != theResult.Get()->GetNumberOfSecondaries() ; j++ )
    {
       if ( theResult.Get()->GetSecondary(j)->GetParticle()->GetDefinition() == aRecoil ) residual = true;
    }
  if ( residual == false )
  {
       G4int nNonZero = 0;
       G4LorentzVector p_photons(0,0,0,0);
       for ( G4int j = 0 ; j != theResult.Get()->GetNumberOfSecondaries() ; j++ )
       {
          p_photons += theResult.Get()->GetSecondary(j)->GetParticle()->Get4Momentum();
          // To many 0 momentum photons -> Check PhotonDist 
          if ( theResult.Get()->GetSecondary(j)->GetParticle()->Get4Momentum().e() > 0 ) nNonZero++;
       }

       // Can we include kinetic energy here?
       G4double deltaE = ( theTrack.Get4Momentum().e() + theTarget.GetTotalEnergy() )
                       - ( p_photons.e() + aRecoil->GetPDGMass() );

//Add photons
       if ( nPhotons - nNonZero > 0 ) 
       {
              //G4cout << "TKDB G4ParticleHPCaptureFS::ApplyYourself we will create additional " << nPhotons - nNonZero << " photons" << G4endl;
          std::vector<G4double> vRand;
          vRand.push_back( 0.0 );
          for ( G4int j = 0 ; j != nPhotons - nNonZero - 1 ; j++ )
          { 
             vRand.push_back( G4UniformRand() );
          }
          vRand.push_back( 1.0 );
          std::sort( vRand.begin(), vRand.end() );

          std::vector<G4double> vEPhoton;
          for ( G4int j = 0 ; j < (G4int)vRand.size() - 1 ; j++ )
          {
             vEPhoton.push_back( deltaE * ( vRand[j+1] - vRand[j] ) );
          }
          std::sort( vEPhoton.begin(), vEPhoton.end() );

          for ( G4int j = 0 ; j < nPhotons - nNonZero - 1 ; j++ )
          {
             //Isotopic in LAB OK?
             G4double theta = pi*G4UniformRand();
             G4double phi = twopi*G4UniformRand();
             G4double sinth = std::sin(theta);
             G4double en = vEPhoton[j];
             G4ThreeVector tempVector(en*sinth*std::cos(phi), en*sinth*std::sin(phi), en*std::cos(theta) );
              
             p_photons += G4LorentzVector ( tempVector, tempVector.mag() );
             G4DynamicParticle * theOne = new G4DynamicParticle;
             theOne->SetDefinition( G4Gamma::Gamma() );
             theOne->SetMomentum( tempVector );
             theResult.Get()->AddSecondary(theOne);
          }

//        Add last photon 
          G4DynamicParticle * theOne = new G4DynamicParticle;
          theOne->SetDefinition( G4Gamma::Gamma() );
//        For better momentum conservation 
          G4ThreeVector lastPhoton = -p_photons.vect().unit()*vEPhoton.back();
          p_photons += G4LorentzVector( lastPhoton , lastPhoton.mag() );
          theOne->SetMomentum( lastPhoton );
          theResult.Get()->AddSecondary(theOne);
       }

//Add residual 
       G4DynamicParticle * theOne = new G4DynamicParticle;
       G4ThreeVector aMomentum = theTrack.Get4Momentum().vect() + theTarget.GetMomentum()
			       - p_photons.vect();
       theOne->SetDefinition(aRecoil);
       theOne->SetMomentum( aMomentum );
       theResult.Get()->AddSecondary(theOne);

  }
//101203TK END
  // clean up the primary neutron
  theResult.Get()->SetStatusChange(stopAndKill);
  return theResult.Get();
}



void GdNeutronHPCaptureFS::UpdateNucleus( const G4Fragment* gamma , G4double eGamma )
{

  G4LorentzVector p4Gamma = gamma->GetMomentum();
  G4ThreeVector pGamma(p4Gamma.vect());

  G4LorentzVector p4Nucleus(nucleus->GetMomentum() );

  static G4ParticleTable * partTab = G4ParticleTable::GetParticleTable();
  G4IonTable * 	ionTable = partTab->GetIonTable();
  G4double m1 = ionTable->GetIonMass(static_cast<G4int>(nucleus->GetZ()),static_cast<G4int>(nucleus->GetA()));
  //G4double m1 = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(static_cast<G4int>(nucleus->GetZ()),
  //									       static_cast<G4int>(nucleus->GetA()));
  G4double m2 = nucleus->GetZ() *  G4Proton::Proton()->GetPDGMass() + 
    (nucleus->GetA()- nucleus->GetZ())*G4Neutron::Neutron()->GetPDGMass();

  //G4cout<<"m1="<<m1<<" m2="<<m2<<G4endl;

  G4double Mass = std::min(m1,m2);

  G4double newExcitation = p4Nucleus.mag() - Mass - eGamma;

  // G4cout<<" Egamma = "<<eGamma<<G4endl;
  //   G4cout<<" NEW EXCITATION ENERGY = "<< newExcitation <<G4endl;

  if(newExcitation < 0)
    newExcitation = 0;

  G4ThreeVector p3Residual(p4Nucleus.vect() - pGamma);
  G4double newEnergy = std::sqrt(p3Residual * p3Residual +
      (Mass + newExcitation) * (Mass + newExcitation));
  G4LorentzVector p4Residual(p3Residual, newEnergy);

  // Update excited nucleus parameters

  nucleus->SetMomentum(p4Residual);

  //  G4cout<<"nucleus excitation energy = "<<nucleus->GetExcitationEnergy() <<G4endl;  

  return;
}

void GdNeutronHPCaptureFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
  G4String tString = "/FS/";
  G4bool dbool;
  G4NeutronHPDataUsed aFile = theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), dirName, tString, dbool);
  G4String filename = aFile.GetName();
  theBaseA = A;
  theBaseZ = G4int(Z+.5);
  if(!dbool || ( Z<2.5 && ( std::abs(theBaseZ - Z)>0.0001 || std::abs(theBaseA - A)>0.0001)))
  {
    hasAnyData = false;
    hasFSData = false; 
    hasXsec = false;
    return;
  }
  std::ifstream theData(filename, std::ios::in);

  hasFSData = theFinalStatePhotons.InitMean(theData); 
  if(hasFSData)
  {
    targetMass = theFinalStatePhotons.GetTargetMass();
    theFinalStatePhotons.InitAngular(theData); 
    theFinalStatePhotons.InitEnergies(theData); 
  }
  theData.close();
}


void GdNeutronHPCaptureFS::InitANNRIGdGenerator()
{
  Gd_CAPTURE = 1;
  Gd_CASCADE = 1;
  if(!getenv("SKDETSIMDIR"))
          throw G4HadronicException(__FILE__, __LINE__, "Please setenv SKDETSIMDIR to point to the skdetsim directory.");
  G4String dirName = getenv("SKDETSIMDIR");

  Gd155_ROOTFile = dirName + "/156GdContTbl__E1SLO4__HFB.root";
  Gd157_ROOTFile = dirName + "/158GdContTbl__E1SLO4__HFB.root";

	if(ANNRIGammaGen) delete ANNRIGammaGen;
	ANNRIGammaGen = new ANNRIGdGammaSpecModel::ANNRIGd_GdNCaptureGammaGenerator();

	ANNRIGdGammaSpecModel::ANNRIGd_GeneratorConfigurator::Configure(*ANNRIGammaGen,
			        Gd_CAPTURE, Gd_CASCADE, Gd155_ROOTFile, Gd157_ROOTFile);
}


G4ReactionProductVector* GdNeutronHPCaptureFS::GenerateWithANNRIGdGenerator()
{
	G4ReactionProductVector* theProducts = 0;

	if(ANNRIGammaGen) {
		ANNRIGdGammaSpecModel::ReactionProductVector products;
		if(Gd_CAPTURE == 1)                          // nat. Gd
			products = ANNRIGammaGen->Generate_NatGd();
		else if(Gd_CAPTURE == 2 and Gd_CASCADE == 1)  // 157Gd, discrete + continuum
			products = ANNRIGammaGen->Generate_158Gd();
		else if(Gd_CAPTURE == 2 and Gd_CASCADE == 2)  // 157Gd, discrete
			products = ANNRIGammaGen->Generate_158Gd_Discrete();
		else if(Gd_CAPTURE == 2 and Gd_CASCADE == 3)  // 157Gd, continuum
			products = ANNRIGammaGen->Generate_158Gd_Continuum();
		else if(Gd_CAPTURE == 3 and Gd_CASCADE == 1)  // 155Gd, discrete + continuum
			products = ANNRIGammaGen->Generate_156Gd();
		else if(Gd_CAPTURE == 3 and Gd_CASCADE == 2)  // 155Gd, discrete
			products = ANNRIGammaGen->Generate_156Gd_Discrete();
		else if(Gd_CAPTURE == 3 and Gd_CASCADE == 3)  // 155Gd, continuum
			products = ANNRIGammaGen->Generate_156Gd_Continuum();

		theProducts = ANNRIGdGammaSpecModel::ANNRIGd_OutputConverter::ConvertToG4(products);
	}
	else {
		//Printing.Error(File_Name,"AnnriGammaGen");
		theProducts = new G4ReactionProductVector();
	}

	return theProducts;
}
