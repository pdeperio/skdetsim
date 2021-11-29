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
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
// 080319 Compilation warnings - gcc-4.3.0 fix by T. Koi
//

// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPorLElastic.hh"
#include "G4SystemOfUnits.hh"
#include "G4NeutronHPElasticFS.hh"
#include "G4ParticleHPManager.hh"
#include "G4Threading.hh"

G4NeutronHPorLElastic::G4NeutronHPorLElastic()
  :G4HadronicInteraction("NeutronHPElastic")
  ,theElastic(NULL)
  ,numEle(0)
{
   overrideSuspension = false;
   G4NeutronHPElasticFS * theFS = new G4NeutronHPElasticFS;
   if(!getenv("G4NEUTRONHPDATA")) 
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
   dirName = getenv("G4NEUTRONHPDATA");
   G4String tString = "/Elastic/";
   dirName = dirName + tString;
//    G4cout <<"G4NeutronHPorLElastic::G4NeutronHPorLElastic testit "<<dirName<<G4endl;
   numEle = G4Element::GetNumberOfElements();
   theElastic = new std::vector<G4ParticleHPChannel*>;
   unavailable_elements.clear();
   for (G4int i=0; i<numEle; i++)
   {
      theElastic->push_back( new G4ParticleHPChannel );
      ((*theElastic)[i])->Init((*(G4Element::GetElementTable()))[i], dirName);
      //theElastic[i].Init((*(G4Element::GetElementTable()))[i], dirName);
      //while(!theElastic[i].Register(theFS));
      try { while(!((*theElastic)[i])->Register(theFS)) ; }
      catch ( G4HadronicException )
      {
          unavailable_elements.insert ( (*(G4Element::GetElementTable()))[i]->GetName() ); 
      }
   }
   delete theFS;
   SetMinEnergy(0.*eV);
   SetMaxEnergy(20.*MeV);
   if ( unavailable_elements.size() > 0 ) 
   {
      std::set< G4String>::iterator it;
      G4cout << "HP Elastic data are not available for thess  elements "<< G4endl;
      for ( it = unavailable_elements.begin() ; it != unavailable_elements.end() ; it++ )
         G4cout << *it << G4endl;
      G4cout << "Low Energy Parameterization Models will be used."<< G4endl;
   }

   createXSectionDataSet();
   //G4ParticleHPManager::GetInstance()->OpenReactionWhiteBoard();
}
  
G4NeutronHPorLElastic::~G4NeutronHPorLElastic()
{
   //delete [] theElastic;
   if ( theElastic != NULL ) {
            for ( std::vector<G4ParticleHPChannel*>::iterator
                it = theElastic->begin() ; it != theElastic->end() ; it++ ) {
                delete *it;
            }
            theElastic->clear();
        }
   delete theDataSet; 
   //G4ParticleHPManager::GetInstance()->CloseReactionWhiteBoard();
}
  
#include "G4NeutronHPThermalBoost.hh"
  
G4HadFinalState * G4NeutronHPorLElastic::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus)
{
   G4ParticleHPManager::GetInstance()->OpenReactionWhiteBoard();
   const G4Material * theMaterial = aTrack.GetMaterial();
   G4int n = theMaterial->GetNumberOfElements();
   G4int index = theMaterial->GetElement(0)->GetIndex();
   if(n!=1)
   {
      G4int i;
      G4double* xSec = new G4double[n];
      G4double sum=0;
      const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
      G4double rWeight;    
      G4NeutronHPThermalBoost aThermalE;
      for (i=0; i<n; i++)
      {
        index = theMaterial->GetElement(i)->GetIndex();
        rWeight = NumAtomsPerVolume[i];
        //G4double x = aThermalE.GetThermalEnergy(aTrack, theMaterial->GetElement(i), theMaterial->GetTemperature());

        //xSec[i] = theElastic[index].GetXsec(aThermalE.GetThermalEnergy(aTrack,
        //	                               theMaterial->GetElement(i),
        //				       theMaterial->GetTemperature()));
        //xSec[i] = theElastic[index].GetXsec(x);
        xSec[i] = ((*theElastic)[index])->GetXsec(aThermalE.GetThermalEnergy(aTrack,
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
        if(sum == 0 ||random<=running/sum) break;
      }
      delete [] xSec;
      // it is element-wise initialised.
    }
    //G4HadFinalState* finalState = theElastic[index].ApplyYourself(aTrack);
    G4HadFinalState* finalState = ((*theElastic)[index])->ApplyYourself(aTrack);
    if (overrideSuspension) finalState->SetStatusChange(isAlive);

   //Overwrite target parameters
    aNucleus.SetParameters(G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA(),G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargZ());
    const G4Element* target_element = (*G4Element::GetElementTable())[index];
    const G4Isotope* target_isotope=NULL;
    G4int iele = target_element->GetNumberOfIsotopes();
    for ( G4int j = 0 ; j != iele ; j++ ) { 
       target_isotope=target_element->GetIsotope( j );
       if ( target_isotope->GetN() == G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA() ) break; 
    }
    //G4cout << "Target Material of this reaction is " << theMaterial->GetName() << G4endl;
    //G4cout << "Target Element of this reaction is " << target_element->GetName() << G4endl;
    //G4cout << "Target Isotope of this reaction is " << target_isotope->GetName() << G4endl;
    aNucleus.SetIsotope( target_isotope );
    
    G4ParticleHPManager::GetInstance()->CloseReactionWhiteBoard();

    return finalState;
}



G4bool G4NeutronHPorLElastic::IsThisElementOK( G4String name )
{
   if ( unavailable_elements.find( name ) == unavailable_elements.end() ) 
      return TRUE;
   else 
      return FALSE; 
}



void G4NeutronHPorLElastic::createXSectionDataSet()
{
   //theDataSet = new G4NeutronHPorLElasticData ( theElastic , &unavailable_elements );
   //theDataSet = new G4NeutronHPElasticData ( /*theElastic , &unavailable_elements*/ );
}

const std::pair<G4double, G4double> G4NeutronHPorLElastic::GetFatalEnergyCheckLevels() const
{
   //return std::pair<G4double, G4double>(10*perCent,10*GeV);
   return std::pair<G4double, G4double>(10*perCent,DBL_MAX);
}

G4int G4NeutronHPorLElastic::GetVerboseLevel() const 
{
   return G4ParticleHPManager::GetInstance()->GetVerboseLevel();
}
void G4NeutronHPorLElastic::SetVerboseLevel( G4int newValue ) 
{
   G4ParticleHPManager::GetInstance()->SetVerboseLevel(newValue);
}

void G4NeutronHPorLElastic::BuildPhysicsTable(const G4ParticleDefinition&)
{

  /* G4ParticleHPManager* hpmanager = G4ParticleHPManager::GetInstance();

   theElastic = hpmanager->GetElasticFinalStates();

   if ( G4Threading::IsMasterThread() ) {

      if ( theElastic == NULL ) theElastic = new std::vector<G4ParticleHPChannel*>;

      if ( numEle == (G4int)G4Element::GetNumberOfElements() ) return;

      if ( theElastic->size() == G4Element::GetNumberOfElements() ) {
         numEle = G4Element::GetNumberOfElements();
         return;
      }

      G4NeutronHPElasticFS * theFS = new G4NeutronHPElasticFS;
      if(!getenv("G4NEUTRONHPDATA")) 
         throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
      dirName = getenv("G4NEUTRONHPDATA");
      G4String tString = "/Elastic";
      dirName = dirName + tString;
      for ( G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; i++ ) {
         if((*(G4Element::GetElementTable()))[i]->GetZ() != 64) continue;//only register Gd
         theElastic->push_back( new G4ParticleHPChannel );
         ((*theElastic)[i])->Init((*(G4Element::GetElementTable()))[i], dirName);
         //while(!((*theElastic)[i])->Register(theFS)) ;
         ((*theElastic)[i])->Register(theFS) ;
      }
      delete theFS;
      hpmanager->RegisterElasticFinalStates( theElastic );

   }
   numEle = 1;//G4Element::GetNumberOfElements();*/
}
void G4NeutronHPorLElastic::ModelDescription(std::ostream& outFile) const
{
   outFile << "High Precision model based on Evaluated Nuclear Data Files (ENDF) for inelastic reaction of neutrons below 20MeV\n";
}
