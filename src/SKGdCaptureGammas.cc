//-------------------------------------------------------------------
// Emission gamma energy spectrum.  
// Based on experimental data, L.V. Groshev et al., Nucl. Data Tab. A5(1968) 1
// The sampling method is derived from DYB Geant3 simulation program.
//------------------------------------------------------------------
// Author: Liang Zhan, 2006/01/27
// Modified: bv@bnl.gov, 2008/4/16 for DetSim
//------------------------------------------------------------------

#include "SKGdCaptureGammas.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream>
#include <sstream> 
#include "globals.hh"
//#include <TRandom.h>
#include <vector>
#include <iostream>
#include "G4Gamma.hh"
#include "G4Poisson.hh"
//#include "TROOT.h"
//#include "TFile.h"
//#include "TH1F.h"
//#include "TMath.h"

/////////////////////////////////////////////////////////////////////

using namespace std;
//using namespace CLHEP;

SKGdCaptureGammas::SKGdCaptureGammas () 
{

  //input gamma spectrum
  //this is just the Gd-157 spectrum from G4NDL for testing purpose, can be replaced by other spectrum
  //data format is in gamma energy(MeV):rate
  //G4String specfilename = "gdcap_gamma.dat";
  if(!getenv("SKDETSIMDIR")) 
    throw G4HadronicException(__FILE__, __LINE__, "Please setenv SKDETSIMDIR to point to the skdetsim directory.");

  G4String dirName = getenv("SKDETSIMDIR");
  G4String specfilename = "/gdcap_gamma.dat";
  specfilename=dirName+specfilename;

  std::ifstream theData(specfilename, std::ios::in);
  if (theData.good() == true) {
    G4cout<<"Read in gdcap_gamma.dat"<<G4endl;
    theData >> numOfPairs;
    G4double x,y;
    total = 0;
    for (G4int i=0; i<numOfPairs; i++) {
      theData >> x >> y;
      gammaSpectrum.push_back(x/1e6);
      gammaSpectrum.push_back(y);
      total += y;
    }
  } else G4cout << "Can not open gdcap_gamma.dat. "<< G4endl;


  // Gamma model prepared by Takatomi Yano, Dec, 2014
  //
  //When you present using this model, please reference these two papers:
  //1. http://www.sciencedirect.com/science/article/pii/S0168900216306477
  //2. hagiwara et al. will be published soon as KMI2017 proceedings
  //

  //   Prepare a table and kopecky_int(200,200)
  // i : E_res-E_gam
  // j : E_gam

  std::string Gd155_File = "/Gd155.dat";
  std::string Gd157_File = "/Gd157.dat";
  Gd155_File=dirName.data()+Gd155_File;
  Gd157_File=dirName.data()+Gd157_File;

  //ifstream ifs155("Gd155.dat");
  //ifstream ifs157("Gd157.dat");

  std::ifstream ifs155(Gd155_File.c_str());
  std::ifstream ifs157(Gd157_File.c_str());

  int nbin=200;
  double prob_int;
  for (int i=0;i<nbin;i++){
    for (int j=0;j<nbin;j++){
      if(ifs155>>prob_int){
        kopecky155_int[i][j]=prob_int;
      }else {
        printf("error while reading Gd155.dat");
      }

      if(ifs157>>prob_int){
        kopecky157_int[i][j]=prob_int;
      }else {
        printf("error while reading Gd157.dat");
      }
    }
  }

}

SKGdCaptureGammas::~SKGdCaptureGammas () 
{
}

G4ReactionProductVector* SKGdCaptureGammas::GetGammas ()
{
  G4ReactionProductVector* theGammas = new G4ReactionProductVector;

  // These 2 variates are modified from Daya Bay code
  //std::vector<double> energy = GetEnergy();
  //std::vector<double> energy = GetEnergy2();

  // Gamma model prepared by Takatomi Yano
  std::vector<double> energy = GetEnergy_ggarnet();

  for(unsigned int i=0; i<energy.size(); i++) {

    G4ReactionProduct* theOne = new G4ReactionProduct;
    theOne->SetDefinition( G4Gamma::Gamma() );

    // Get the gammas direction. 
    // Isotropic emission.
    G4double costheta = 2.*G4UniformRand()-1;
    G4double theta = acos(costheta);
    G4double phi = CLHEP::twopi*G4UniformRand();
    G4double sinth = sin(theta);
    theOne->SetTotalEnergy( energy[i] );
    G4ThreeVector temp(energy[i]*sinth*cos(phi), 
        energy[i]*sinth*sin(phi),
        energy[i]*cos(theta) );
    theOne->SetMomentum( temp ) ;
    theGammas->push_back(theOne);
  }
  return theGammas;
}

G4double SKGdCaptureGammas::getRandomCLHEP()
{
  // return a random number distributed according the histogram bin contents.

  G4int i;
  G4double sum = 0;
  G4double rand = G4UniformRand()*total;

  for (i=0; i<numOfPairs; i++) {
    sum += gammaSpectrum[2*i+1];
    if (sum >= rand) break;
  }
  G4double energy;
  energy = gammaSpectrum[2*i];
  //if (sum > rand) 
  energy += (gammaSpectrum[2*i+2]-gammaSpectrum[2*i])*(sum-rand)/(gammaSpectrum[2*i+1]);
  return energy;

}


// The total energy of emission gammas is determinate.
// The gammas number and energy can be modified to make the sampling spectrum 
// to be close to the experimental spectrum as best as one can.
// The number of th gammas is about 3.

std::vector<double>  SKGdCaptureGammas::GetEnergy ()
{
  std::vector<double> Energy;
  G4double TotalEnergy;
  if (G4UniformRand()<0.815)
    TotalEnergy = 7.9373; // captured by Gd157
  else
    TotalEnergy = 8.53; // captured by Gd155


  //resample:

  G4double energy_sum = 0.0;
  G4int    gamma_num = 0;
  G4double gamma_energy[100]={0};
  G4double energy_tmp;

sample:
  G4double energy = getRandomCLHEP();

  energy_sum = energy_sum+energy;
  gamma_num = gamma_num+1;
  gamma_energy[gamma_num] = energy;

  if(energy_sum>TotalEnergy) {
    if (gamma_num==2) {
      gamma_energy[gamma_num]=TotalEnergy-gamma_energy[gamma_num-1];
      goto goon;
    }
    energy_tmp = gamma_energy[gamma_num-1]+energy-(energy_sum-TotalEnergy);
    if(energy_tmp>6.58) {
      /*energy_sum = energy_sum-energy;
        gamma_num = gamma_num-1;
        goto sample;*/
      gamma_energy[gamma_num]-=energy_sum-TotalEnergy;
      goto goon;
    }
    gamma_num = gamma_num-1;
    gamma_energy[gamma_num]=energy_tmp;
    goto goon;
  }
  if ((TotalEnergy-energy_sum)<0.4) {
    if(gamma_energy[gamma_num]+TotalEnergy-energy_sum>6.58){
      gamma_num = gamma_num+1;
      gamma_energy[gamma_num] = TotalEnergy-energy_sum;
    } else gamma_energy[gamma_num] += TotalEnergy-energy_sum;
    goto goon;
  }
  goto sample;

goon:

  // set the sequence
  for(int i=gamma_num;i>1;i--) {
    for(int j=1;j<i;j++) {
      if(gamma_energy[j]<gamma_energy[j+1]) {
        energy_tmp = gamma_energy[j+1];
        gamma_energy[j+1] = gamma_energy[j];
        gamma_energy[j] = energy_tmp;
      }
    }
  }

ok:

  for(int k=1;k<=gamma_num;k++)
    Energy.push_back(gamma_energy[k]); 

  return Energy;
}

std::vector<double>  SKGdCaptureGammas::GetEnergy2 ()
{
  std::vector<double> Energy;
  G4double TotalEnergy;
  if (G4UniformRand()<0.815)
    TotalEnergy = 7.9373; // captured by Gd157
  else
    TotalEnergy = 8.53; // captured by Gd155
  G4double energy_sum = 0.0;
  G4int    gamma_num = 0;
  G4double gamma_energy[100]={0};
  G4double energy_tmp;

sample:


  //Use the CLHEP random engine to sample the distribution instead
  //Jianglai 10/01/2006
  //G4double energy=capspec->GetRandom();// sample from the experimental spectrum.
  G4double energy = getRandomCLHEP();
  //G4cout<<energy_sum<<G4endl;
  energy_sum = energy_sum+energy;
  gamma_num = gamma_num+1;
  gamma_energy[gamma_num] = energy;

  if(energy_sum>TotalEnergy) {
    if((energy_sum-TotalEnergy)<0.8) {
      energy_tmp = gamma_energy[gamma_num-1]+energy-(energy_sum-TotalEnergy);
      if(energy_tmp>7.4) {
        if (gamma_energy[gamma_num-1]>7.4) {
          gamma_energy[gamma_num] = energy-(energy_sum-TotalEnergy);
          goto goon;
        }
        energy_sum = energy_sum-energy;
        gamma_num = gamma_num-1;
        goto sample;
      }
      gamma_num = gamma_num-1;
      gamma_energy[gamma_num]=energy_tmp;
      goto goon;
    }
    energy_sum = energy_sum-energy;
    gamma_num = gamma_num-1;
    goto sample;
  }
  goto sample;

goon:


  // set the sequence
  for(int i=gamma_num;i>1;i--) {
    for(int j=1;j<i;j++) {
      if(gamma_energy[j]<gamma_energy[j+1]) {
        energy_tmp = gamma_energy[j+1];
        gamma_energy[j+1] = gamma_energy[j];
        gamma_energy[j] = energy_tmp;
      }
    }
  }
  if(gamma_num<=3)
    goto  ok;

  if(gamma_energy[gamma_num]<0.8) {
    energy_tmp = gamma_energy[gamma_num]+gamma_energy[gamma_num-1];
    if(energy_tmp>7.4)
      goto ok;
    gamma_energy[gamma_num-1]=energy_tmp;
    gamma_num=gamma_num-1;
  }

ok:

  for(int k=1;k<=gamma_num;k++)
    Energy.push_back(gamma_energy[k]); 

  return Energy;
}

//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::GetEnergy_ggarnet ()
  //////////////////////////////////////////////
{
  std::vector<double> Eg;
  double DefElem;
  DefElem = G4UniformRand()*48.764;	//14.8% (abundance of Gd155) * 60900 (cross-section of Gd155) + 15.65% (abundance of Gd157) * 254000 (cross-section of Gd157) = 4.87642e+6
  if(DefElem < 39.751){
    Eg = CapGad157(); 			//15.65% (abundance of Gd157) * 254000 (cross-section of Gd157) = 3.9751e+6
  }else{
    Eg = CapGad155();
  }
  return Eg;
}
////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::CapGad155()
  ////////////////////////////////////////////
{
  //  gammas from GAD155
  //  either 2 gammas
  //  either a continuum
  //  total energy = 8.46 MeV
  double Nlevel = G4UniformRand();
  double a=0;
  std::vector<double> Egammas;

  a=0.14;	//Default Continuum is 86% and Discrete is 14%

  double ProbD1 = 0.2134*a;	//6.429MeV+...
  double ProbD2 = 0.3775*a;	//6.766MeV+...
  double ProbD3 = 0.2519*a;	//7.288MeV+...
  double ProbD4 = 0.1437*a;	//7.382MeV+...
  double ProbD5 = 0.0135*a;	//8.447MeV+...
  double ProbC  = 1.-(ProbD1+ProbD2+ProbD3+ProbD4+ProbD5);
  Nlevel = Nlevel * ( ProbC + ProbD1 + ProbD2 + ProbD3 + ProbD4 + ProbD5 );

  if(Nlevel < ProbC ) Egammas=continuum155();
  else if(Nlevel < ProbC + ProbD1) Egammas=casca1Gd155(); 
  else if(Nlevel < ProbC + ProbD1 + ProbD2)  Egammas=casca2Gd155(); 
  else if(Nlevel < ProbC + ProbD1 + ProbD2 + ProbD3)  Egammas=casca3Gd155(); 
  else if(Nlevel < ProbC + ProbD1 + ProbD2 + ProbD3 + ProbD4)  Egammas=casca4Gd155(); 
  else if(Nlevel < ProbC + ProbD1 + ProbD2 + ProbD3 + ProbD4 + ProbD5)  Egammas=casca5Gd155(); 

  return Egammas;
}

////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::CapGad157()
  ////////////////////////////////////////////
{
  //  gammas from GAD157
  //                   either 2 gammas
  //                   either a continuum
  //                   total energy = 7.87 MeV
  //                   either 1 gamma 

  double Nlevel = G4UniformRand();
  double b=0;
  std::vector<double> Egammas;

  b=0.029015;	//Default Continuum is 97.0985% and Discrete is 2.9015% (calculated by T.Tanaka at 2017.03.03)

  double ProbD1 = 0.0008*b;	//7.937MeV+...
  double ProbD2 = 0.0033*b;	//7.857MeV+...
  double ProbD3 = 0.0029*b;	//6.960MeV+...
  double ProbD4 = 0.0183*b;	//6.914MeV+...
  double ProbD5 = 0.3516*b;	//6.750MeV+...
  double ProbD6 = 0.0276*b;	//6.672MeV+...
  double ProbD7 = 0.0472*b;	//6.420MeV+...
  double ProbD8 = 0.0227*b;	//6.006MeV+...
  double ProbD9 = 0.1703*b;	//5.903MeV+...
  double ProbD10 = 0.0286*b;	//5.784MeV+...
  double ProbD11 = 0.0896*b;	//5.677MeV+...
  double ProbD12 = 0.0931*b;	//5.610MeV+...
  double ProbD13 = 0.0341*b;	//5.543MeV+...
  double ProbD14 = 0.0234*b;	//5.436MeV+...
  double ProbD15 = 0.0865*b;	//5.179MeV+...
  double ProbC  = 1.-(ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7+ProbD8+ProbD9+ProbD10+ProbD11+ProbD12+ProbD13+ProbD14+ProbD15);
  Nlevel = Nlevel*(ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7+ProbD8+ProbD9+ProbD10+ProbD11+ProbD12+ProbD13+ProbD14+ProbD15);

  if(Nlevel<ProbC) Egammas = continuum157();
  else if(Nlevel<ProbC+ProbD1) Egammas = casca1Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2) Egammas = casca2Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3) Egammas = casca3Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4) Egammas = casca4Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5) Egammas = casca5Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6) Egammas = casca6Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7) Egammas = casca7Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7+ProbD8) Egammas = casca8Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7+ProbD8+ProbD9) Egammas = casca9Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7+ProbD8+ProbD9+ProbD10) Egammas = casca10Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7+ProbD8+ProbD9+ProbD10+ProbD11) Egammas = casca11Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7+ProbD8+ProbD9+ProbD10+ProbD11+ProbD12) Egammas = casca12Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7+ProbD8+ProbD9+ProbD10+ProbD11+ProbD12+ProbD13) Egammas = casca13Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7+ProbD8+ProbD9+ProbD10+ProbD11+ProbD12+ProbD13+ProbD14) Egammas = casca14Gd157(); 
  else if(Nlevel<ProbC+ProbD1+ProbD2+ProbD3+ProbD4+ProbD5+ProbD6+ProbD7+ProbD8+ProbD9+ProbD10+ProbD11+ProbD12+ProbD13+ProbD14+ProbD15) Egammas = casca15Gd157(); 

  return Egammas;
}

//====================================//
///////////155Gd cascade////////////////
//====================================//
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca1Gd155()
  //////////////////////////////////////////////
{
  std::vector<double> energy;
  double cas155_1 = G4UniformRand();
  double cas155_11 = 0.618;
  double cas155_12 = 0.382;

  energy.push_back(6.429);
  if(cas155_1 < cas155_11){
    energy.push_back(2.018);
    energy.push_back(0.089);
  }else if(cas155_1 < cas155_11+cas155_12){
    energy.push_back(1.819);
    energy.push_back(0.200);
    energy.push_back(0.089);
  }
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca2Gd155()
  //////////////////////////////////////////////
{
  std::vector<double> energy;
  double cas155_2 = G4UniformRand();
  double cas155_21 = 0.535;
  double cas155_22 = 0.465;

  energy.push_back(6.766);
  if(cas155_2 < cas155_21){
    energy.push_back(1.682);
    energy.push_back(0.089);
  }else if(cas155_2 < cas155_21+cas155_22){
    energy.push_back(1.158);
    energy.push_back(0.200);
    energy.push_back(0.089);
  }
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca3Gd155()
  //////////////////////////////////////////////
{
  vector<double> energy;
  double cas155_3 = G4UniformRand();
  double cas155_31 = 0.833;
  double cas155_32 = 0.167;

  energy.push_back(7.288);
  if(cas155_3 < cas155_31){
    energy.push_back(1.158);
    energy.push_back(0.089);
  }else if(cas155_3 < cas155_31+cas155_32){
    energy.push_back(0.960);
    energy.push_back(0.200);
    energy.push_back(0.089);
  }
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca4Gd155()
  //////////////////////////////////////////////
{
  vector<double> energy;
  double cas155_4 = G4UniformRand();
  double cas155_41 = 0.518;
  double cas155_42 = 0.482;

  energy.push_back(7.382);
  if(cas155_4 < cas155_41){
    energy.push_back(1.155);
  }else if(cas155_4 < cas155_41+cas155_42){
    energy.push_back(1.066);
    energy.push_back(0.089);
  }
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca5Gd155()
  //////////////////////////////////////////////
{
  vector<double> energy;
  energy.push_back(8.447);
  energy.push_back(0.089);
  return energy;
}

//====================================//
////////////////Gd157 cascade///////////
//====================================//
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca1Gd157()
  //////////////////////////////////////////////
{
  ////7937keV (Q-value)/////
  vector<double> energy;
  energy.push_back(7.937);
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca2Gd157()
  //////////////////////////////////////////////
{
  ////7857keV + 80keV /////
  vector<double> energy;
  energy.push_back(7.857);
  energy.push_back(0.080);
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca3Gd157()
  //////////////////////////////////////////////
{
  ////6960keV  /////
  vector<double> energy;
  energy.push_back(6.960);
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca4Gd157()
  //////////////////////////////////////////////
{
  ////6914keV + 944keV + 80keV /////
  vector<double> energy;
  energy.push_back(6.914);
  energy.push_back(0.944);
  energy.push_back(0.080);
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca5Gd157()
  //////////////////////////////////////////////
{
  ////6750keV + 1187 or 1107keV /////
  vector<double> energy;
  double cas157_3 = G4UniformRand();
  double cas157_31 = 0.501;
  double cas157_32 = 0.499;

  energy.push_back(6.750);
  if(cas157_3 < cas157_31){
    energy.push_back(1.187);
  }else if(cas157_3 < cas157_31+cas157_32){
    energy.push_back(1.107);
    energy.push_back(0.080);
  }
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca6Gd157()
  //////////////////////////////////////////////
{
  ////6672keV + .../////
  vector<double> energy;
  double cas157_4 = G4UniformRand();
  double cas157_41 = 0.848;
  double cas157_42 = 0.152;

  energy.push_back(6.672);
  if(cas157_4 < cas157_41){
    energy.push_back(1.186);
  }else if(cas157_4 < cas157_41+cas157_42){
    energy.push_back(1.004);
    energy.push_back(0.182);
    energy.push_back(0.080);
  }
  return energy;
}
///////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca7Gd157()
  //////////////////////////////////////////////
{
  ////6420keV + .... /////
  vector<double> energy;
  double cas157_5 = G4UniformRand();
  double cas157_51 = 0.381;
  double cas157_52 = 0.415;
  double cas157_53 = 0.204;

  energy.push_back(6.420);
  if(cas157_5 < cas157_51){
    energy.push_back(1.517);
  }else if(cas157_5 < cas157_51+cas157_52){
    energy.push_back(1.438);
    energy.push_back(0.080);
  }else if(cas157_5 < cas157_51+cas157_52+cas157_53){
    energy.push_back(1.256);
    energy.push_back(0.182);
    energy.push_back(0.080);
  }
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca8Gd157()
  //////////////////////////////////////////////
{
  ////6000.5keV (6006,5995 double peaks)/////
  vector<double> energy;
  energy.push_back(6.001);
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca9Gd157()
  //////////////////////////////////////////////
{
  ////5903keV + .... /////
  vector<double> energy;
  double cas157_7 = G4UniformRand();
  double cas157_7_1 = G4UniformRand();
  double cas157_7_2 = G4UniformRand();
  double cas157_71 = 0.392;
  double cas157_72 = 0.254;
  double cas157_73 = 0.222;
  double cas157_74 = 0.132;

  double cas157_73_1 = 0.848;
  double cas157_73_2 = 0.152;
  double cas157_74_1 = 0.755;
  double cas157_74_2 = 0.245;

  energy.push_back(5.903);
  if(cas157_7 < cas157_71){
    energy.push_back(1.010);
    energy.push_back(0.944);
    energy.push_back(0.080);
  }else if(cas157_7 < cas157_71+cas157_72){
    energy.push_back(0.875);
    energy.push_back(0.898);
    energy.push_back(0.182);
    energy.push_back(0.080);
  }else if(cas157_7 < cas157_71+cas157_72+cas157_73){
    energy.push_back(0.769);
    if(cas157_7_1 < cas157_73_1){
      energy.push_back(1.186);
      energy.push_back(0.080);
    }else if(cas157_7 < cas157_73_1+cas157_73_2){
      energy.push_back(1.004);
      energy.push_back(0.182);
      energy.push_back(0.080);
    }
  }else if(cas157_7 < cas157_71+cas157_72+cas157_73+cas157_74){
    energy.push_back(0.676);
    if(cas157_7_2 < cas157_74_1){
      energy.push_back(1.097);
      energy.push_back(0.182);
      energy.push_back(0.080);
    }else if(cas157_7 < cas157_74_1+cas157_74_2){
      energy.push_back(1.279);
      energy.push_back(0.080);
    }
  }
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca10Gd157()
  //////////////////////////////////////////////
{
  ////5784keV /////
  vector<double> energy;
  energy.push_back(5.784);
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca11Gd157()
  //////////////////////////////////////////////
{
  ////5669keV (5677,5661 double peaks)/////
  vector<double> energy;
  energy.push_back(5.669);
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca12Gd157()
  //////////////////////////////////////////////
{
  ////5595keV (5610,5593,5582 triple peaks)/////
  vector<double> energy;
  energy.push_back(5.595);
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca13Gd157()
  //////////////////////////////////////////////
{
  ////5543keV /////
  vector<double> energy;
  energy.push_back(5.543);
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca14Gd157()
  //////////////////////////////////////////////
{
  ////5436keV /////
  vector<double> energy;
  energy.push_back(5.436);
  return energy;
}
//////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::casca15Gd157()
  //////////////////////////////////////////////
{
  ////5167keV (5179,5155 double peaks)/////
  vector<double> energy;
  energy.push_back(5.167);
  return energy;
}

///////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::continuum155()
  ///////////////////////////////////////////////
{
  // continuum part of gadolinium
  //cross sections
  vector<double> energy;

  double BindingEnergy = 8.46; // MeV, for 157Gd
  int nbin = 200; // constant for reading table
  double E_max = BindingEnergy; //constant for reading table, maximum value of E

  double E_res = BindingEnergy-0.0001;
  double E_gam = 0;
  int ithbin = 0;
  //	int maxbin=0;
  int ybin = 0;
  int ngam = 0;
  double random = 0;

  while(E_res>0.20){
    random = G4UniformRand();
    E_gam = random*E_res;
    ybin = (int)((BindingEnergy-E_res)/E_max*nbin);
    ithbin = (int)(E_gam/E_max*nbin);

    random=G4UniformRand();
    if (kopecky155_int[ithbin][ybin]>random){
      energy.push_back(E_gam);
      E_res=E_res-E_gam;
      ngam=ngam+1;
    }
  }
  energy.push_back(E_res);
  return energy;
}

///////////////////////////////////////////////
std::vector<double>  SKGdCaptureGammas::continuum157()
  ///////////////////////////////////////////////
{
  // continuum part of gadolinium
  //cross sections
  vector<double> energy;

  double BindingEnergy = 7.87; // MeV, for 157Gd
  int nbin = 200; // constant for reading table
  double E_max = BindingEnergy; //constant for reading table, maximum value of E

  double E_res = BindingEnergy-0.0001;
  double E_gam = 0;
  int ithbin = 0;
  //	int maxbin=0;
  int ybin = 0;
  int ngam = 0;
  double random = 0;

  double E182 = G4UniformRand();
  double add182 = 0.153;

  while(E_res>0.20){
    random = G4UniformRand();
    E_gam = random*E_res;
    ybin = (int)((BindingEnergy-E_res)/E_max*nbin);
    ithbin = (int)(E_gam/E_max*nbin);

    random = G4UniformRand();
    if (kopecky157_int[ithbin][ybin]>random){
      energy.push_back(E_gam);
      E_res = E_res-E_gam;
      ngam = ngam+1;
    }
  }

  if(E182<add182){
    energy.push_back(0.182);
    E_res = E_res-0.182;
    energy.push_back(E_res);
  }else{
    energy.push_back(E_res);
  }
  return energy;
}
