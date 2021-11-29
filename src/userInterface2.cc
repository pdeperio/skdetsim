#include "userInterface2.hh"
#include "G4NucleiProperties.hh"
#include "ExN01PhysicsList.hh"
#include "G4ParticleHPManager.hh"

G4int NumGamma;
G4double GammaEnergies[15];
G4int fMODEL;
/*****************************************************************************/
 userInterface2::userInterface2( G4int iseed , G4int gd_water , G4int gd_model , G4float gd_concentration )
/*****************************************************************************/
{
//----------------------------------------------------------------------------
//  Random seed setting in GEANT4 side
//
  CLHEP::HepRandom::setTheEngine( new CLHEP::HepJamesRandom );
  CLHEP::HepRandom::getTheEngine()->setSeed( iseed, 0 );
  G4cout << " userInterface2::userInterface2" << G4endl;
  G4cout << "   - Initialized CLHEP::HepJamesRandom with " << iseed << G4endl;


//----------------------------------------------------------------------------
//  Most basic setup is needed. 
//
  runManager = new G4RunManager;
  //
  G4cout<<" <Detector Construction> "<<G4endl;
  runManager->SetUserInitialization( det = new skGADDetectorConstruction );
  G4cout<<" <Physics List> "<<G4endl;
  runManager->SetUserInitialization( new ExN01PhysicsList() );
  G4cout<<" <Primary Generator Action> "<<G4endl;
  runManager->SetUserAction( new skGADPrimaryGeneratorAction );
  runManager->Initialize();


//----------------------------------------------------------------------------
//  Define detector materials (actually, only GdWater is necessary).
//
  G4cout<<"\t Material Definition \t#########"<<G4endl;
  // basic parameters 
  //
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4int iz, n;
  G4int ncomponents, natoms;
  G4double abundance, fractionmass;
  // 
  G4int polyPMMA = 1;
  G4int nC_PMMA  = 3 + 2 * polyPMMA;
  G4int nH_PMMA  = 6 + 2 * polyPMMA;
  G4int polyeth = 1;
  G4int nC_eth  = 2 * polyeth;
  G4int nH_eth  = 4 * polyeth;

  // elements 
  //
  H  = new G4Element( "Hydrogen",  "H",  z = 1.,  a = 1.01 * g/mole );
  B  = new G4Element( "Boron",     "B",  z = 5.,  a = 10.81 * g/mole );
  C  = new G4Element( "Carbon",    "C",  z = 6.,  a = 12.01 * g/mole );
  N  = new G4Element( "Nitrogen",  "N",  z = 7.,  a = 14.01 * g/mole );
  O  = new G4Element( "Oxygen" ,   "O",  z = 8.,  a = 16.00 * g/mole );
  Na = new G4Element( "Sodium" ,   "Na", z = 11., a = 22.99 * g/mole );
  Al = new G4Element( "Aluminum" , "Al", z = 13., a = 26.98 * g/mole );
  Si = new G4Element( "Silicon",   "Si", z = 14., a = 28.09 * g/mole ); 
  S  = new G4Element( "Sulfur",    "S",  z = 16., a = 32.07 * g/mole );
  Cl = new G4Element( "Chlorine",  "Cl", z = 17., a = 35.453 * g/mole );
  Co = new G4Element( "Cobalt",    "Co", z = 27., a = 58.9332 * g/mole ); 
  Bi = new G4Element( "Bismuth",   "Bi", z = 83., a = 208.980 * g/mole ); 
  // (Fe)
  G4Isotope* Fe54 = new G4Isotope( "Fe54", iz = 26, n = 54, a = 53.9396105 * g/mole );
  G4Isotope* Fe56 = new G4Isotope( "Fe56", iz = 26, n = 56, a = 55.9349375 * g/mole );
  G4Isotope* Fe57 = new G4Isotope( "Fe57", iz = 26, n = 57, a = 56.9353940 * g/mole );
  G4Isotope* Fe58 = new G4Isotope( "Fe58", iz = 26, n = 58, a = 57.9332756 * g/mole );
  Fe = new G4Element( "Iron", "Fe", ncomponents = 4 );
  Fe->AddIsotope( Fe54, abundance = 5.8 * perCent );
  Fe->AddIsotope( Fe56, abundance = 91.72 * perCent );
  Fe->AddIsotope( Fe57, abundance = 2.2 * perCent );
  Fe->AddIsotope( Fe58, abundance = 0.28 * perCent );
  // (Ge)
  G4Isotope* Ge70 = new G4Isotope( "Ge70", iz = 32, n = 70, a = 69.9242474 * g/mole );
  G4Isotope* Ge72 = new G4Isotope( "Ge72", iz = 32, n = 72, a = 71.9220758 * g/mole );
  G4Isotope* Ge73 = new G4Isotope( "Ge73", iz = 32, n = 73, a = 72.9234589 * g/mole );
  G4Isotope* Ge74 = new G4Isotope( "Ge74", iz = 32, n = 74, a = 73.9211778 * g/mole );
  G4Isotope* Ge76 = new G4Isotope( "Ge76", iz = 32, n = 76, a = 75.9214026 * g/mole );
  Ge = new G4Element( "Germanium", "Ge", ncomponents = 5 ); 
  Ge->AddIsotope( Ge70, abundance = 20.38 * perCent );
  Ge->AddIsotope( Ge72, abundance = 27.31 * perCent );
  Ge->AddIsotope( Ge73, abundance = 7.76 * perCent );
  Ge->AddIsotope( Ge74, abundance = 36.72 * perCent );
  Ge->AddIsotope( Ge76, abundance = 7.83 * perCent );
  // (Gd)
  G4Isotope* Gd152 = new G4Isotope( "Gd152", iz = 64, n = 152, a = 151.919786 * g/mole );
  G4Isotope* Gd154 = new G4Isotope( "Gd154", iz = 64, n = 154, a = 153.920861 * g/mole );
  G4Isotope* Gd155 = new G4Isotope( "Gd155", iz = 64, n = 155, a = 154.922618 * g/mole );
  G4Isotope* Gd156 = new G4Isotope( "Gd156", iz = 64, n = 156, a = 155.922118 * g/mole );
  G4Isotope* Gd157 = new G4Isotope( "Gd157", iz = 64, n = 157, a = 156.923956 * g/mole );
  G4Isotope* Gd158 = new G4Isotope( "Gd158", iz = 64, n = 158, a = 157.924099 * g/mole );
  G4Isotope* Gd160 = new G4Isotope( "Gd160", iz = 64, n = 160, a = 159.927049 * g/mole );
  Gd = new G4Element( "Gadolinium", "Gd", ncomponents = 7 );
  Gd->AddIsotope( Gd152, abundance = 0.20 * perCent );
  Gd->AddIsotope( Gd154, abundance = 2.18 * perCent );
  Gd->AddIsotope( Gd155, abundance = 14.80 * perCent );
  Gd->AddIsotope( Gd156, abundance = 20.47 * perCent );
  Gd->AddIsotope( Gd157, abundance = 15.65 * perCent );
  Gd->AddIsotope( Gd158, abundance = 24.84 * perCent );
  Gd->AddIsotope( Gd160, abundance = 21.86 * perCent );
  
  // materials 
  //
  GdCl3 = new G4Material( "GdCl3", density = 2.424 * g/cm3, ncomponents = 4 );
  GdCl3->AddElement( Gd, natoms = 1 );
  GdCl3->AddElement( Cl, natoms = 3 );
  GdCl3->AddElement( H,  natoms = 12 );
  GdCl3->AddElement( O,  natoms = 6 );
  // 
  G4double density_gdsulfate = 3.01; 
  GdSulfateOctahydrate = new G4Material("GdSulfateOctahydrate", density= density_gdsulfate * g/cm3, ncomponents=4);
  GdSulfateOctahydrate->AddElement( Gd, natoms = 2 );
  GdSulfateOctahydrate->AddElement( S,  natoms = 3 );
  GdSulfateOctahydrate->AddElement( H,  natoms = 16 );
  GdSulfateOctahydrate->AddElement( O,  natoms = 20 );
  //
  Water = new G4Material( "Water", density = 1.0 * g/cm3, ncomponents = 2 );
  Water->AddElement( H, natoms = 2 );
  Water->AddElement( O, natoms = 1 );
  // 
  BGO = new G4Material( "BGO", density = 7.13 * g/cm3, ncomponents = 3 );
  BGO->AddElement( Bi, natoms = 4 );
  BGO->AddElement( Ge, natoms = 3 );
  BGO->AddElement( O,  natoms = 12 );

  // (note when testing Water<->GdWater need to change water_mt also ...) 
//     M.Harada modified for selecting Gd concentration(2020 Jul. 1st)
	G4double AAS = gd_concentration;
  //G4double AAS = 0.2178;
  //G4double AAS = 0.554963;
  G4double conc = AAS; // * 0.6538;
  GdWater = new G4Material( "GdWater", density = (1.00 / (1 - conc * perCent)) * g/cm3, ncomponents = 2 );
  //GdWater = new G4Material( "GdWater", density = 1.0138 * g/cm3, ncomponents = 2 );
  //GdWater = new G4Material( "GdWater", density = 1.01742 * g/cm3, ncomponents = 2 );
  GdWater->AddMaterial( Water,    fractionmass = ( 100. - conc ) * perCent );
  GdWater->AddMaterial( GdSulfateOctahydrate, fractionmass = conc * perCent );


//----------------------------------------------------------------------------
//  Calculate each material rate to the total amount. 
//
  G4int nElements = GdWater->GetNumberOfElements();
  const G4double* nAtomsPerVolume = GdWater->GetVecNbOfAtomsPerVolume();
  nTotalAtomsPerVolume = 0;
  nHAtomsPerVolume = 0;
  nOAtomsPerVolume = 0;
  nGdAtomsPerVolume = 0;
  rateHtoTotal = 0; 
  rateOtoTotal = 0; 
  rateGdtoTotal = 0; 

  for ( G4int i = 0; i < nElements; ++i ) 
  {
    G4cout << ( *GdWater->GetElementVector() )[ i ]->GetName() 
           << " : " << nAtomsPerVolume[ i ] * ( cm3/mm3 ) << " [/cm^3]" 
	   << G4endl;
    //
    nTotalAtomsPerVolume = nTotalAtomsPerVolume + nAtomsPerVolume[ i ]; 
    if ( ( *GdWater->GetElementVector() )[ i ]->GetName() == "Hydrogen" ) 
      nHAtomsPerVolume  = nAtomsPerVolume[ i ];
    if ( ( *GdWater->GetElementVector() )[ i ]->GetName() == "Oxygen" ) 
      nOAtomsPerVolume  = nAtomsPerVolume[ i ];
    if ( ( *GdWater->GetElementVector() )[ i ]->GetName() == "Gadolinium" ) 
      nGdAtomsPerVolume = nAtomsPerVolume[ i ];
  }
  //std::cout << "nHAtomsPerVolume = "  << nHAtomsPerVolume  << std::endl; 
  //std::cout << "nOAtomsPerVolume = "  << nOAtomsPerVolume  << std::endl; 
  //std::cout << "nGdAtomsPerVolume = " << nGdAtomsPerVolume << std::endl; 
  //
  rateHtoTotal  = nHAtomsPerVolume/nTotalAtomsPerVolume; 
  rateOtoTotal  = nOAtomsPerVolume/nTotalAtomsPerVolume; 
  rateGdtoTotal = nGdAtomsPerVolume/nTotalAtomsPerVolume; 
  //std::cout << "rateHtoTotal = "  << rateHtoTotal  << std::endl; 
  //std::cout << "rateOtoTotal = "  << rateOtoTotal  << std::endl; 
  //std::cout << "rateGdtoTotal = " << rateGdtoTotal << std::endl; 
  
  // area expansion [/mm^3] --> [/cm^3] (note that this does not mean unit conversion.)
  nHAtomsPerVolume  = nHAtomsPerVolume * ( cm3/mm3 ); 
  nOAtomsPerVolume  = nOAtomsPerVolume * ( cm3/mm3 ); 
  nGdAtomsPerVolume = nGdAtomsPerVolume * ( cm3/mm3 ); 

  runManager->GeometryHasBeenModified();
//----------------------------------------------------------------------------
//  Load GEANT4 neutron interaction model class. 
//
  particle = G4Neutron::NeutronDefinition();

  // capture 
  //
  fMODEL = gd_model;
  theCaptureProcess = new G4HadronCaptureProcess;
  // (for loading only Gd)
  //theLENeutronCaptureModel = new myG4NeutronHPorLCapture;
  //theLENeutronCaptureModel = new G4NeutronHPorLCapture;  
  theLENeutronCaptureModel = new GdNeutronHPCapture;
  //theLENeutronCaptureModel = new G4NeutronHPCapture;  
  theCaptureProcess->RegisterMe( theLENeutronCaptureModel );
  theNeutronData3 = new G4NeutronHPCaptureData;
  theCaptureProcess->AddDataSet( theNeutronData3 );
  capProcessXsecDataStore = theCaptureProcess->GetCrossSectionDataStore();
  capProcessXsecDataStore->AddDataSet(new G4ParticleHPCaptureData);
  capProcessXsecDataStore->BuildPhysicsTable(*particle);
  capProcessXsecDataStore->DumpPhysicsTable(*particle);
  

  // elastic 
  //
  theNeutronElasticProcess = new G4HadronElasticProcess;
  // (for loading only Gd)
  //theElasticNeutron = new myG4NeutronHPorLElastic(); 
  theElasticNeutron = new G4NeutronHPorLElastic();
  //theElasticNeutron = new G4NeutronHPElastic();//theElasticNeutron->BuildPhysicsTable(*particle);
  theNeutronElasticProcess->RegisterMe( theElasticNeutron );
  theNeutronData = new G4NeutronHPElasticData();
  theNeutronElasticProcess->AddDataSet( theNeutronData );
  elaProcessXsecDataStore = theNeutronElasticProcess->GetCrossSectionDataStore();
  elaProcessXsecDataStore->AddDataSet(new G4ParticleHPElasticData);
  elaProcessXsecDataStore->BuildPhysicsTable(*particle);
  elaProcessXsecDataStore->DumpPhysicsTable(*particle);
  
  
  // inelastic
  //
  theInelasticProcess = new G4NeutronInelasticProcess( "inelastic" );
  // (for loading only Gd)
  //theInelasticModel = new myG4NeutronHPorLEInelastic(); 
  theInelasticModel = new G4NeutronHPorLEInelastic();
  //theInelasticModel = new G4NeutronHPInelastic();//theInelasticModel->BuildPhysicsTable(*particle);
  theInelasticProcess->RegisterMe( theInelasticModel );
  theNeutronData1 = new G4NeutronHPInelasticData;
  theInelasticProcess->AddDataSet( theNeutronData1 );
  inelaProcessXsecDataStore = theInelasticProcess->GetCrossSectionDataStore();
  inelaProcessXsecDataStore->AddDataSet(new G4ParticleHPInelasticData);
  inelaProcessXsecDataStore->BuildPhysicsTable(*particle);
  inelaProcessXsecDataStore->DumpPhysicsTable(*particle);


//----------------------------------------------------------------------------
//  Define some parameters for later use. 
//
  G4ThreeVector aMomentumDirection( 0, 0, 1 );
  G4ThreeVector aPosition( 0, 0, 0 );
  //particle = G4Neutron::NeutronDefinition();
  neutronParticle = new G4DynamicParticle( particle, aMomentumDirection, 0.02 * eV );
  aStep = new G4Step();
  G4StepPoint* aStepPoint = new G4StepPoint();
	
  // for pure water simulation 
  // for Gd water simulation
  if(gd_water) aStepPoint->SetMaterial( GdWater );
  else aStepPoint->SetMaterial( Water );
  aStep->SetPreStepPoint( aStepPoint );

  aTrack = new G4Track( neutronParticle, 0., aPosition );

	G4ParticleDefinition		*piminus=G4PionMinus::PionMinusDefinition();
	piminusParticle2=new G4DynamicParticle(piminus, aMomentumDirection, 0.*eV);
	PiMiStep			=new G4Step();
	PiMiStepPoint 	=new G4StepPoint();
  if(gd_water) PiMiStepPoint->SetMaterial( GdWater );
  else PiMiStepPoint->SetMaterial( Water );
	PiMiStep					->SetPreStepPoint(PiMiStepPoint);
	aTrackPiMi 		=new G4Track(piminusParticle2, 0., aPosition);
	aTrackPiMi					->SetStep(PiMiStep);

	G4ParticleDefinition		*muminus=G4MuonMinus::MuonMinusDefinition();
	muminusParticle=new G4DynamicParticle(muminus, aMomentumDirection, 0.*eV);
	MuMiStep			=new G4Step();
	MuMiStepPoint 	=new G4StepPoint();
  if(gd_water) MuMiStepPoint->SetMaterial( GdWater );
  else MuMiStepPoint->SetMaterial( Water );
	MuMiStep					->SetPreStepPoint(MuMiStepPoint);
	aTrackMuMi 		=new G4Track(muminusParticle, 0., aPosition);
	aTrackMuMi					->SetStep(MuMiStep);
}


/*****************************************************************************/
 userInterface2::~userInterface2()
/*****************************************************************************/
{
}


/*****************************************************************************/
 void userInterface2::Initialize()
/*****************************************************************************/
{
}


/*****************************************************************************/
 void userInterface2::Finalize()
/*****************************************************************************/
{
  delete runManager;
}


//----------------------------------------------------------------------------
//  Get H cross sections 
//   - GetCrossSection(): [cm^2] 
//   - nHAtomsPerVolume: [/mm^3] (note that the value is area-expanded to /cm^3.)
//
//  You are always getting H and O xsec from GdWater
//  Any difference bewteen GdWater and Water?
/*****************************************************************************/
 G4double userInterface2::GetHCapXsec( G4double nEnergy )
/*****************************************************************************/
{
  neutronParticle->SetKineticEnergy( nEnergy * GeV );
  return capProcessXsecDataStore->GetCrossSection( neutronParticle, H, GdWater )/cm2 
                                                 * nHAtomsPerVolume/( 1./mm3 );
}

/*****************************************************************************/
 G4double userInterface2::GetHElaXsec( G4double nEnergy )
/*****************************************************************************/
{
  neutronParticle->SetKineticEnergy( nEnergy * GeV );
  return elaProcessXsecDataStore->GetCrossSection( neutronParticle, H, GdWater )/cm2 
                                                 * nHAtomsPerVolume/( 1./mm3 ); 
}

/*****************************************************************************/
 G4double userInterface2::GetHInElaXsec( G4double nEnergy )
/*****************************************************************************/
{
  neutronParticle->SetKineticEnergy( nEnergy * GeV );
  return inelaProcessXsecDataStore->GetCrossSection( neutronParticle, H, GdWater )/cm2 
                                                   * nHAtomsPerVolume/( 1./mm3 ); 
}


//----------------------------------------------------------------------------
//  Get O cross sections 
//   - GetCrossSection(): [cm^2] 
//   - nOAtomsPerVolume: [/mm^3] (note that the value is area-expanded to /cm^3.)
//
/*****************************************************************************/
 G4double userInterface2::GetOCapXsec( G4double nEnergy )
/*****************************************************************************/
{
  neutronParticle->SetKineticEnergy( nEnergy * GeV );
  return capProcessXsecDataStore->GetCrossSection( neutronParticle, O, GdWater )/cm2 
                                                 * nOAtomsPerVolume/( 1./mm3 ); 
}

/*****************************************************************************/
 G4double userInterface2::GetOElaXsec( G4double nEnergy )
/*****************************************************************************/
{
  neutronParticle->SetKineticEnergy( nEnergy * GeV );
  return elaProcessXsecDataStore->GetCrossSection( neutronParticle, O, GdWater )/cm2 
                                                 * nOAtomsPerVolume/( 1./mm3 );
}

/*****************************************************************************/
 G4double userInterface2::GetOInElaXsec( G4double nEnergy )
/*****************************************************************************/
{
  neutronParticle->SetKineticEnergy( nEnergy * GeV );
  return inelaProcessXsecDataStore->GetCrossSection( neutronParticle, O, GdWater )/cm2 
                                                   * nOAtomsPerVolume/( 1./mm3 );
}


//----------------------------------------------------------------------------
//  Get Gd cross sections 
//   - GetCrossSection(): [cm^2] 
//   - nGdAtomsPerVolume: [/mm^3] (note that the value is area-expanded to /cm^3)
//
/*****************************************************************************/
 G4double userInterface2::GetGdCapXsec( G4double nEnergy )
/*****************************************************************************/
{
  neutronParticle->SetKineticEnergy( nEnergy * GeV );
  return capProcessXsecDataStore->GetCrossSection( neutronParticle, Gd, GdWater )/cm2 
                                                 * nGdAtomsPerVolume/( 1./mm3 ); 
}

/*****************************************************************************/
 G4double userInterface2::GetGdElaXsec( G4double nEnergy )
/*****************************************************************************/
{
  neutronParticle->SetKineticEnergy( nEnergy * GeV );
  return elaProcessXsecDataStore->GetCrossSection( neutronParticle, Gd, GdWater )/cm2 
                                                 * nGdAtomsPerVolume/( 1./mm3 ); 
}

/*****************************************************************************/
 G4double userInterface2::GetGdInElaXsec( G4double nEnergy )
/*****************************************************************************/
{
  neutronParticle->SetKineticEnergy( nEnergy * GeV );
  return inelaProcessXsecDataStore->GetCrossSection( neutronParticle, Gd, GdWater )/cm2 
                                                   * nGdAtomsPerVolume/( 1./mm3 ); 
}

//G4Track* aTrack;
//----------------------------------------------------------------------------
//  Get secondary particles
//   - processNumber = { 10, 11, 12, 20, 21, 22 } set in gtneut.F
//
/*****************************************************************************/
 std::vector< Secondary_Vector > userInterface2::Get2nd( G4double px, 
                                                         G4double py, 
							 G4double pz, 
							 G4int process )
/*****************************************************************************/
{
  // Select target nucleus among the material based on cross sections. 
  //
  G4ThreeVector aMomentum( px * GeV, py * GeV, pz * GeV );
  neutronParticle->SetMomentum( aMomentum );
  //G4Element* anElement = nullptr;
  G4Nucleus targetNucleus;
  //
  //std::cout << " processNumber (G4) = " << process << std::endl; 
  //std::cout << "TEST (anElement) : " << anElement << std::endl; 
  
  // Get the detailed interaction information.  
  //
  G4HadProjectile thePro;
  //G4ThreeVector aPosition( 0, 0, 0 );
  //aTrack = new G4Track( neutronParticle, 0., aPosition );
  G4double mass=neutronParticle->GetDefinition()->GetPDGMass();
  G4double KE=neutronParticle->GetKineticEnergy();
  aTrack->SetMomentumDirection(aMomentum.unit());
  aTrack->SetKineticEnergy(KE);
  aTrack->SetStep( aStep );
  thePro.Initialise( *aTrack );
  G4int reentryCount = 0;
  result = nullptr;

  do
  {
    try
    { 
      // Get the interaction results.
      //
      if ( process == 10 )       // neutron capture by H or O
        result = theLENeutronCaptureModel->ApplyYourself( thePro, targetNucleus );
      //else if ( process == 11 )  // neutron elastic scattering on H or O 
        //result = theElasticNeutron->ApplyYourself( thePro, targetNucleus );
      else if ( process == 12 )  // neutron inelastic scattering on H or O 
        result = theInelasticModel->ApplyYourself( thePro, targetNucleus );
      else if ( process == 20 ){  // neutron capture by H or O or Gd 
        result = theLENeutronCaptureModel->ApplyYourself( thePro, targetNucleus );
			}
      else if ( process == 11 || process == 21 ){  // neutron elastic scattering on H or O or Gd
        // NOTE:  Very low energy scatters were causing numerical (FPE) errors
        //        in earlier releases; these limits have not been changed since.
        if (KE <= 1.e-6*eV) {
          G4cout<<"GEANT4: Neutron KE="<<KE/eV<<"eV too low to be handled by Elastic process."<<G4endl;
          std::vector< Secondary_Vector > secondaryVector;
          Secondary_Vector secVec;
          secVec.evId = 0;
          secVec.productionTime = 0;  // [ns]
          secVec.productionPosX = 0;
          secVec.productionPosY = 0;
          secVec.productionPosZ = 0;
          secVec.momentumX = px;
          secVec.momentumY = py;
          secVec.momentumZ = pz;
          secVec.totalEnergy = neutronParticle->GetTotalEnergy()/GeV;
          secVec.pdgID = aTrack->GetParticleDefinition()->GetPDGEncoding();
          secondaryVector.push_back( secVec );
          free(result);
          return secondaryVector;
        }
        result = theElasticNeutron->ApplyYourself( thePro, targetNucleus );
      }
      else                       // neutron inelastic scattering on H or O or Gd
        result = theInelasticModel->ApplyYourself( thePro, targetNucleus );
      //
      ++reentryCount;
    }
    catch ( G4HadronicException aR )
    {
      G4ExceptionDescription ed;
      aR.Report( ed );
      ed << "Call for " << theLENeutronCaptureModel->GetModelName() << G4endl;
      ed //<< "Target element " << anElement->GetName() 
         << "  Z = "
         << targetNucleus.GetZ_asInt()
         << "  A = " << targetNucleus.GetA_asInt() << G4endl;
      ed << " ApplyYourself failed" << G4endl;
      G4Exception( "G4HadronicProcess::PostStepDoIt", "had006", FatalException, ed );
    }
 

    // Check the result for catastrophic energy non-conservation.
    //
    result = CheckResult( thePro, targetNucleus, result, theLENeutronCaptureModel );
    // 
    if ( reentryCount > 100 ) 
    {
      G4ExceptionDescription ed;
      ed << "Call for " << theLENeutronCaptureModel->GetModelName() << G4endl;
      ed //<< "Target element " << anElement->GetName() 
         << "  Z = "
         << targetNucleus.GetZ_asInt()
         << "  A = " << targetNucleus.GetA_asInt() << G4endl;
      //theCaptureProcess->DumpState( aTrack, "ApplyYourself", ed );
      ed << " ApplyYourself does not completed after 100 attempts" << G4endl;
      G4Exception( "G4HadronicProcess::PostStepDoIt", "had006", FatalException, ed );
    }
  }
  while ( !result );

  if (process == 11 || process==21){
    //elastic process needs a different treatment
    G4ThreeVector indir = aTrack->GetMomentumDirection();
    G4double phi = CLHEP::twopi*G4UniformRand();
    G4ThreeVector it(0., 0., 1.);
    G4ThreeVector outdir = result->GetMomentumChange();
    std::vector< Secondary_Vector > secondaryVector;
    G4double efinal = result->GetEnergyChange();
    if(efinal < 0.0) { efinal = 0.0; }
    if(result->GetNumberOfSecondaries() > 0){
        G4DynamicParticle* p = result->GetSecondary(0)->GetParticle();
        G4ThreeVector pdir = p->GetMomentumDirection();
        pdir.rotate(phi, it);
        pdir.rotateUz(indir);
        p->SetMomentumDirection(pdir);
        Secondary_Vector secVec;
        secVec.evId = 0;
        secVec.productionTime = 0;  // [ns]
        secVec.productionPosX = 0;
        secVec.productionPosY = 0;
        secVec.productionPosZ = 0;
        secVec.momentumX = p->GetMomentum().x()/GeV;
        secVec.momentumY = p->GetMomentum().y()/GeV;
        secVec.momentumZ = p->GetMomentum().z()/GeV;
        secVec.totalEnergy = p->GetTotalEnergy()/GeV;
        secVec.pdgID = p->GetParticleDefinition()->GetPDGEncoding();
        secondaryVector.push_back( secVec );
        delete p;
    }
    if(efinal > 0.0) {
        outdir.rotate(phi, it);
        outdir.rotateUz(indir);
        G4double newE = efinal + mass;
        G4double newP = std::sqrt( efinal * ( efinal + 2 * mass ) );
        Secondary_Vector secVec;
        secVec.evId = 0;
        secVec.productionTime = 0;  // [ns]
        secVec.productionPosX = 0;
        secVec.productionPosY = 0;
        secVec.productionPosZ = 0;
        secVec.momentumX = newP*outdir.x()/GeV;
        secVec.momentumY = newP*outdir.y()/GeV;
        secVec.momentumZ = newP*outdir.z()/GeV;
        secVec.totalEnergy = newE/GeV;
        secVec.pdgID = aTrack->GetParticleDefinition()->GetPDGEncoding();
        secondaryVector.push_back( secVec );
    }
    result->Clear();
    result = nullptr;
    free(result);
    return secondaryVector;
  }

  // Output secondaries for capture/inelastic processes
  result->SetTrafoToLab( thePro.GetTrafoToLab() );
  G4double rotation = CLHEP::twopi * G4UniformRand();
  G4ThreeVector it( 0., 0., 1. );


  // Put the secondaries into a vector output.
  //
  G4int nSec = result->GetNumberOfSecondaries();
  std::vector< Secondary_Vector > secondaryVector;
  if ( nSec > 0 ) 
  {
    //double totene = 0.;//test 
    for ( G4int i = 0; i < nSec; ++i ) 
    {
      G4LorentzVector theM = result->GetSecondary( i )->GetParticle()->Get4Momentum();
      theM.rotate( rotation, it );
      theM *= result->GetTrafoToLab();
      result->GetSecondary( i )->GetParticle()->Set4Momentum( theM );
 
      G4double time = result->GetSecondary( i )->GetTime();
      G4DynamicParticle* secParticle = result->GetSecondary( i )->GetParticle();
      G4ThreeVector secMomentum = secParticle->GetMomentum();
      Secondary_Vector secVec;
      secVec.evId = 0;
      secVec.productionTime = 0;  // [ns]
      secVec.productionPosX = 0;
      secVec.productionPosY = 0;
      secVec.productionPosZ = 0;
      secVec.momentumX = secMomentum.x()/GeV;
      secVec.momentumY = secMomentum.y()/GeV;
      secVec.momentumZ = secMomentum.z()/GeV;
      secVec.totalEnergy = secParticle->GetTotalEnergy()/GeV;
      secVec.pdgID = secParticle->GetPDGcode();
			//if(process == 20)std::cout<<" pdgID : "<<secVec.pdgID<<" secVec.totalEnergy : "<<secVec.totalEnergy<<std::endl;
      //if (process == 20&&secParticle->GetPDGcode()==22) totene +=secVec.totalEnergy;
      secondaryVector.push_back( secVec );
      delete secParticle;
    }
    //if(process==20) std::cout<<" Total gamma energy : "<<totene <<std::endl;
  }
  G4double efinal = result->GetEnergyChange();
  if ( efinal < 0.0 ) { efinal = 0.0; }

  // Output the neutron too if it is not killed.
  //
  if ( result->GetStatusChange() != stopAndKill && efinal != 0 ) 
  {
    G4double newE = efinal + mass;
    G4double newP = std::sqrt( efinal * ( efinal + 2 * mass ) );
    G4ThreeVector newPV = newP * result->GetMomentumChange();
    G4LorentzVector newP4( newE, newPV );
    newP4.rotate( rotation, it );
    newP4 *= result->GetTrafoToLab();
    G4ThreeVector secMomentum = newP4.vect();
    newE = newP4.e() - mass;
    
    if ( newE <= 0.0 ) 
    {
      G4ExceptionDescription ed;
      //DumpState( aT, "Primary has zero energy after interaction", ed );
      G4Exception( "G4HadronicProcess::FillResults", "had011", JustWarning, ed );
      newE = 0.0;
    }
    
    Secondary_Vector secVec;
    secVec.evId = 0;
    secVec.productionTime = 0;  // [ns]
    secVec.productionPosX = 0;
    secVec.productionPosY = 0;
    secVec.productionPosZ = 0;
    secVec.momentumX = secMomentum.x()/GeV;
    secVec.momentumY = secMomentum.y()/GeV;
    secVec.momentumZ = secMomentum.z()/GeV;
    secVec.totalEnergy = newP4.e()/GeV;
    secVec.pdgID = aTrack->GetParticleDefinition()->GetPDGEncoding();
    secondaryVector.push_back( secVec );
  }
  result->Clear();
  result = nullptr;
  free(result);
  return secondaryVector;
}


/*****************************************************************************/
 G4HadFinalState* userInterface2::CheckResult( const G4HadProjectile &aPro,
		                               const G4Nucleus &aNucleus,
		                               G4HadFinalState* aresult,
					       G4HadronicInteraction* hadModel )
/*****************************************************************************/
{
  // check for catastrophic energy non-conservation, to re-sample the interaction
  // 
  G4HadronicInteraction* theModel = hadModel;
  G4double nuclearMass( 0 );
  if ( theModel ) 
  {
     // Compute final-state total energy
     //
     G4double finalE( 0. );

     G4int nSec = aresult->GetNumberOfSecondaries();
     nuclearMass = G4NucleiProperties::GetNuclearMass( aNucleus.GetA_asInt(),
                                                       aNucleus.GetZ_asInt() );

     if ( aresult->GetStatusChange() != stopAndKill ) 
     {
       // Interaction didn't complete, returned "do nothing" state          => reset nucleus
       //  or  the primary survived the interaction (e.g. electro-nuclear ) => keep  nucleus
       finalE = aresult->GetLocalEnergyDeposit() 
              + aPro.GetDefinition()->GetPDGMass() 
	      + aresult->GetEnergyChange();
       
       if ( nSec == 0 ) {
         // Since there are no secondaries, there is no recoil nucleus.
         // To check energy balance we must neglect the initial nucleus too.
         nuclearMass = 0.0;
       }
     }
     
     for ( G4int i = 0; i < nSec; i++ ) {
       finalE += aresult->GetSecondary( i )->GetParticle()->GetTotalEnergy();
     }
     G4double deltaE = nuclearMass + aPro.GetTotalEnergy() - finalE;

     std::pair< G4double, G4double > checkLevels = theModel->GetFatalEnergyCheckLevels();  // [relative, absolute]
     if ( std::abs( deltaE ) > checkLevels.second 
       && std::abs( deltaE ) > checkLevels.first * aPro.GetKineticEnergy() )
     {
        // do not delete result, this is a pointer to a data member;
        //
	      aresult->Clear();
        aresult = nullptr;
        G4ExceptionDescription desc;
        desc << "Warning: Bad energy non-conservation detected, will "
             //<< (epReportLevel<0 ? "abort the event" : "re-sample the interaction") << G4endl
             << " Process / Model: " <<  "G4HadronicCaptureProcess" << " / " << theModel->GetModelName() << G4endl
             << " Primary: " << aPro.GetDefinition()->GetParticleName()
             << " (" << aPro.GetDefinition()->GetPDGEncoding() << "),"
             << " E= " <<  aPro.Get4Momentum().e()
             << ", target nucleus (" << aNucleus.GetZ_asInt() << ","<< aNucleus.GetA_asInt() << ")" << G4endl
             << " E(initial - final) = " << deltaE << " MeV." << G4endl;
        //G4Exception("G4HadronicProcess:CheckResult()", "had012", epReportLevel<0 ? EventMustBeAborted : JustWarning,desc);
     }
  }

  return aresult;
} 	

////////////////////
std::vector< Secondary_Vector > userInterface2::DoPionMinusCapture(G4double px,
																   G4double py,
																   G4double pz)

{
	G4ThreeVector 				aMomentum(px*GeV, py*GeV, pz*GeV);
  piminusParticle2->SetMomentum( aMomentum );
  G4double mass=piminusParticle2->GetDefinition()->GetPDGMass();
  G4double KE=piminusParticle2->GetKineticEnergy();
  aTrackPiMi->SetMomentumDirection(aMomentum.unit());
  aTrackPiMi->SetKineticEnergy(KE);
	
	aTrackPiMi					->SetStep(PiMiStep);

	G4HadronicAbsorptionBertini *thePiMiCap=new G4HadronicAbsorptionBertini(G4PionMinus::PionMinusDefinition());
	G4VParticleChange *vpc=thePiMiCap->AtRestDoIt(*aTrackPiMi, *aStepPiMi);
  	std::vector< Secondary_Vector > secondaryVector;
	for(G4int iscnd=0; iscnd<vpc->GetNumberOfSecondaries(); iscnd++)
	{
		G4String 			parName="";
		G4double			totEne=-1.;
		G4int	 			pdgCode	=-1;
		G4ThreeVector		momScnd(0., 0., 0.);

		parName		=vpc->GetSecondary( iscnd )->GetDefinition()->GetParticleName();
		totEne		=vpc->GetSecondary( iscnd )->GetTotalEnergy();
		pdgCode		=vpc->GetSecondary( iscnd )->GetDynamicParticle()->GetPDGcode();
		momScnd		=vpc->GetSecondary( iscnd )->GetMomentum();

		Secondary_Vector secVec;
		secVec.evId				=0; // Actually, they are obsolete.
		secVec.productionTime	=0; // 
		secVec.productionPosX	=0; //
		secVec.productionPosY	=0; //
		secVec.productionPosZ	=0; //
		secVec.momentumX		=momScnd.x()/GeV;
		secVec.momentumY		=momScnd.y()/GeV;
		secVec.momentumZ		=momScnd.z()/GeV;
		secVec.totalEnergy		=totEne/GeV;
		secVec.pdgID			=pdgCode;
		secondaryVector.push_back(secVec);
	}
	vpc->Clear();
  vpc = nullptr;
  free(vpc);
	delete thePiMiCap;
  return secondaryVector;
}

std::vector< Secondary_Vector > userInterface2::DoMuonMinusCapture(G4double px,
																   G4double py,
																   G4double pz)
{
	G4ThreeVector 				aMomentum(px*GeV, py*GeV, pz*GeV);
  muminusParticle->SetMomentum( aMomentum );
  G4double mass=muminusParticle->GetDefinition()->GetPDGMass();
  G4double KE=muminusParticle->GetKineticEnergy();
  aTrackMuMi->SetMomentumDirection(aMomentum.unit());
  aTrackMuMi->SetKineticEnergy(KE);

	G4MuonMinusCapture *theMuMiCap =new G4MuonMinusCapture();
	G4VParticleChange *vpc=theMuMiCap->AtRestDoIt(*aTrackMuMi, *aStepMuMi);
	
  std::vector< Secondary_Vector > secondaryVector;
	for(G4int iscnd=0; iscnd<vpc->GetNumberOfSecondaries(); iscnd++)
	{
		G4String 			parName="";
		G4double			totEne=-1.;
		G4int	 			pdgCode	=-1;
		G4ThreeVector		momScnd(0., 0., 0.);

		parName		=vpc->GetSecondary( iscnd )->GetDefinition()->GetParticleName();
		totEne		=vpc->GetSecondary( iscnd )->GetTotalEnergy();
		pdgCode		=vpc->GetSecondary( iscnd )->GetDynamicParticle()->GetPDGcode();
		momScnd		=vpc->GetSecondary( iscnd )->GetMomentum();

		Secondary_Vector secVec;
		secVec.evId				=0; // Actually, they are obsolete.
		secVec.productionTime	=0; // 
		secVec.productionPosX	=0; //
		secVec.productionPosY	=0; //
		secVec.productionPosZ	=0; //
		secVec.momentumX		=momScnd.x()/GeV;
		secVec.momentumY		=momScnd.y()/GeV;
		secVec.momentumZ		=momScnd.z()/GeV;
		secVec.totalEnergy		=totEne/GeV;
		secVec.pdgID			=pdgCode;
		secondaryVector.push_back(secVec);
	}
	vpc->Clear();
  vpc = nullptr;
  free(vpc);
	delete theMuMiCap;
 	return secondaryVector;
}
