
///////////////////////////////////////////////////////////////////////////////
//                   Spectrum of radiative neutron capture by Gadolinium            
//                                    version 1.0.0                                
//                                    (Sep.09.2005)                               
                            
//Modified class from original G4NeutronHPCapture class to include 
//the gammas spectrum of radiative neutron capture by Gadolinium.

// Karim Zbiri, April, 2005
///////////////////////////////////////////////////////////////////////////////


#include "GdNeutronHPCapture.hh"
/////////#include "OtherHPCaptureFS.hh"
#include "G4NeutronHPCaptureFS.hh"
#include "G4NeutronHPDeExGammas.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "userInterface2.hh"
#include "G4NeutronHPManager.hh"
#include "GdNeutronHPCaptureFS.hh"

extern G4int fMODEL;

  GdNeutronHPCapture::GdNeutronHPCapture()
   :G4HadronicInteraction("GdNeutronHPCapture")
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 20.*CLHEP::MeV );
    G4cout << "Capture : start of construction!!!!!!!!"<<G4endl;
    if(!getenv("G4NEUTRONHPDATA")) 
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv NeutronHPCrossSections to point to the neutron cross-section files.");
    dirName = getenv("G4NEUTRONHPDATA");
    G4String tString = "/Capture/";
    dirName = dirName + tString;
    //dirName = "/home/takaaki/geant4.9.6.p02/G4NDL4.2/Capture/";
    numEle = G4Element::GetNumberOfElements();
   G4cout << "+++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
   //    G4cout <<"Disname="<<dirName<<" numEle="<<numEle<<G4endl;
  
    //theCapture = new G4NeutronHPChannel[numEle];
    //    G4cout <<"G4NeutronHPChannel constructed"<<G4endl;
    G4NeutronHPCaptureFS * theFS = new G4NeutronHPCaptureFS; 
    //    OtherHPCaptureFS * theFS = new OtherHPCaptureFS; 
    GdNeutronHPCaptureFS * theGdFS = new GdNeutronHPCaptureFS;

    G4cout<<" Register : "<<numEle<<G4endl;
    for (G4int i=0; i<numEle; i++)
    {
      G4cout<<" \t"<<i<<"  "<<(*(G4Element::GetElementTable()))[i]->GetZ() <<G4endl;
      if((*(G4Element::GetElementTable()))[i]->GetZ() != 64)
      {
        theCapture.push_back( new G4NeutronHPChannel );
        (*theCapture[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
        (*theCapture[i]).Register(theFS);
      }
      else {
        theCapture.push_back( new G4NeutronHPChannel );
        (*theCapture[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
        if( fMODEL!=0&&fMODEL!=4 )(*theCapture[i]).Register(theGdFS);
        else  (*theCapture[i]).Register(theFS);
        // theCapture[i].Register(theFS);
      }
    }
    delete theFS;
    delete theGdFS;
//    G4cout << "-------------------------------------------------"<<G4endl;
//    G4cout << "Leaving GdNeutronHPCapture::GdNeutronHPCapture"<<G4endl;
  }
  
  GdNeutronHPCapture::~GdNeutronHPCapture()
  {
//    delete [] theCapture;
//    G4cout << "Leaving GdNeutronHPCapture::~GdNeutronHPCapture"<<G4endl;
      for ( std::vector<G4NeutronHPChannel*>::iterator 
            ite = theCapture.begin() ; ite != theCapture.end() ; ite++ )
      {
         delete *ite;
      }
      theCapture.clear();
  }
  
  #include "G4NeutronHPThermalBoost.hh"
  G4HadFinalState * GdNeutronHPCapture::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus)
  {

     //if ( numEle < (G4int)G4Element::GetNumberOfElements() ) addChannelForNewElement();
    G4NeutronHPManager::GetInstance()->OpenReactionWhiteBoard();
    if(getenv("NeutronHPCapture")) G4cout <<" ####### GdNeutronHPCapture called"<<G4endl;
    const G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = theMaterial->GetNumberOfElements();
    G4int index = theMaterial->GetElement(0)->GetIndex();
    if(n!=1)
    {
      G4double* xSec = new G4double[n];
      G4double sum=0;
      G4int i;
      const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
      G4double rWeight;    
      G4NeutronHPThermalBoost aThermalE;
      for (i=0; i<n; i++)
      {
        index = theMaterial->GetElement(i)->GetIndex();
        rWeight = NumAtomsPerVolume[i];
        xSec[i] = (*theCapture[index]).GetXsec(aThermalE.GetThermalEnergy(aTrack,
  		                                                     theMaterial->GetElement(i),
  								     theMaterial->GetTemperature()));
        xSec[i] *= rWeight;
        sum+=xSec[i];
      }
      G4double random = G4UniformRand();
      G4double running = 0;
      for (i=0; i<n; i++)
      {
        running += xSec[i];
        index = theMaterial->GetElement(i)->GetIndex();
        //if(random<=running/sum) break;
        if( sum == 0 || random <= running/sum ) break;
      }
      if(i==n) i=std::max(0, n-1);
      delete [] xSec;
    }
    //G4cout <<" ####### Loop for ele xsec"<<G4endl;
    //return theCapture[index].ApplyYourself(aTrack);
    G4HadFinalState* result = (*theCapture[index]).ApplyYourself(aTrack);
    //G4cout <<" ####### (*theCapture[index]).ApplyYourself(aTrack);"<<G4endl;
    //Overwrite target parameters
    aNucleus.SetParameters(G4NeutronHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA(),G4NeutronHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargZ());
    //G4cout <<" ####### aNucleus.SetParameters"<<G4endl;
    const G4Element* target_element = (*G4Element::GetElementTable())[index];
    const G4Isotope* target_isotope=NULL;
    G4int iele = target_element->GetNumberOfIsotopes();
    for ( G4int j = 0 ; j != iele ; j++ ) { 
       target_isotope=target_element->GetIsotope( j );
       if ( target_isotope->GetN() == G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA() ) break; 
    }
     aNucleus.SetIsotope( target_isotope );
     G4NeutronHPManager::GetInstance()->CloseReactionWhiteBoard();
     return result;
  }

  const std::pair<G4double, G4double> GdNeutronHPCapture::GetFatalEnergyCheckLevels() const
  {
     //return std::pair<G4double, G4double>(10*perCent,10*GeV);
     return std::pair<G4double, G4double>(10*CLHEP::perCent,DBL_MAX);
  }
 
  void GdNeutronHPCapture::addChannelForNewElement()
  {/*
     G4NeutronHPCaptureFS* theFS = new G4NeutronHPCaptureFS;
     for ( G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; i++ ) 
     {
        G4cout << "G4NeutronHPCapture Prepairing Data for the new element of " << (*(G4Element::GetElementTable()))[i]->GetName() << G4endl;
        theCapture.push_back( new G4NeutronHPChannel );
        (*theCapture[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
        (*theCapture[i]).Register(theFS);
     }
     delete theFS;
     numEle = (G4int)G4Element::GetNumberOfElements();*/
  }
