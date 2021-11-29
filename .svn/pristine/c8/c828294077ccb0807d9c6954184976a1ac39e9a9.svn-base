//-------------------------------------------------------
// Dice the emission gamma energy and direction based on experimental data.
//-------------------------------------------------------
// Author: Liang Zhan, 2006/01/27
//-------------------------------------------------------

#ifndef SKGdCaptureGammas_hh
#define SKGdCaptureGammas_hh

#include "G4ReactionProductVector.hh"
#include <vector>
#include "G4Allocator.hh"
//class TFile;
//class TH1F;

/////////////////////////////////////////////////////////

struct pair{
  G4double x; 
  G4double y;
};

class SKGdCaptureGammas
{
public:
   
    SKGdCaptureGammas();
    ~SKGdCaptureGammas();
    G4ReactionProductVector * GetGammas ();
    std::vector<double>  GetEnergy ();
    std::vector<double>  GetEnergy2 ();
    std::vector<double>  GetEnergy_ggarnet ();

		std::vector<double> CapGad155();
		std::vector<double> CapGad157();
		std::vector<double> continuum155();
		std::vector<double> continuum157();

		std::vector<double> casca1Gd155();
		std::vector<double> casca2Gd155();
		std::vector<double> casca3Gd155();
		std::vector<double> casca4Gd155();
		std::vector<double> casca5Gd155();

		std::vector<double> casca1Gd157();
		std::vector<double> casca2Gd157();
		std::vector<double> casca3Gd157();
		std::vector<double> casca4Gd157();
		std::vector<double> casca5Gd157();
		std::vector<double> casca6Gd157();
		std::vector<double> casca7Gd157();
		std::vector<double> casca8Gd157();
		std::vector<double> casca9Gd157();
		std::vector<double> casca10Gd157();
		std::vector<double> casca11Gd157();
		std::vector<double> casca12Gd157();
		std::vector<double> casca13Gd157();
		std::vector<double> casca14Gd157();
		std::vector<double> casca15Gd157();

    G4double getRandomCLHEP();
    double  kopecky155_int[200][200];
    double  kopecky157_int[200][200];

private:
   
    //TFile * specfile;
    //TH1F  * capspec; // the experimental gamma spectrum.
    G4int numOfPairs;
    G4double total;
    std::vector<G4double>  gammaSpectrum;
    std::vector<G4double>  unusedEnergy;
};

#endif
