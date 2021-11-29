/* Modified from sk1-4 version (Roger?) by Baran (July 2019)
*
* Reasons for modification: 
* 	-Increased type of PMTs: Installiation of >200 PMTs that are different from previous
* 	classification of "old" and "new". The information about which PMT is which is read from
* 	connection table for sk5
* 	-New saturation measurement with laser data for 538 PMTs, other PMTs saturation parameter is
* 	assigned from the distribution of the same type of PMTs everytime open function is called.
* 
* 
* - 2020-01-24 Baran: Remove ConnectionTableReader usage and read table directly, 
*   read (for memory usage concerns) 
*/


// root 
#include "TH1D.h"
#include "TFile.h" 
#include "TTree.h"
#include "TString.h"

//c++
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>

// skofl
#include "skheadC.h"
//#include "ConnectionTable.h"

#define odmax 1885

int local_cable;
int photons;
float q;

TFile * _file;
TFile * _out = 0;
TFile * _file_dark;
TFile * _file_dark_cur;
TH1D  * _old ;
TH1D  * _new ;

TH1D  * _pmt[ odmax ];
TH1D  * _od_dark[5]; // sk 1 to 5
float odDarkTree[odmax];
TTree * _od_dark_tree; // For reading darkr dark rates
TH1D  * _od_sat;

const int nType = 6;
TH1D * _sat_byType[nType]; // for sk5 only
TH1D * _pe_byType[nType]; // for sk5 only
TH1D * _pc2pe; // for sk5 only

TTree * tree;

//ConnectionTable *conTable; // to obtain PMT types when necessary
int pmtTypes[odmax];

TString typeString[nType] = {"all","IMB","KEK","ipmu","UCI","rec"}; // pmt types in od
const int OD_PMT_OFFSET = 20000;

// obtaines required constants
// delete the 1st input after sk5 geometry is added
extern "C" void open_root_(int *inSkGeom, int *darkRun)
{
	int skGeom = *inSkGeom;
	bool flag = false;
  bool typeFlag = false;
  _old = 0;
  _new = 0;
 
	//conTable = new ConnectionTable(skGeom); // do not use this one
		
	//only used in sk-5, read connection table manually, not using ConnectionTableReader.cc
	std::stringstream conFilename;
	conFilename << getenv("SKOFL_ROOT") << "/const/connection.super.sk-5.dat";
	std::ifstream conFile( conFilename.str().c_str() );
	while( conFile.good() ) {
		std::string line;
		int pmtNo, pmtType;
		getline(conFile,line);
		if(line.find("#",0) != std::string::npos || line.empty()) continue; // comment
		// Read only relevant inf	
		pmtNo = atoi(line.substr(3,5).c_str()); 
		pmtType = atoi(line.substr(83,2).c_str()); 

		// This still uses 1.5 times more memory (better than initial 8 x) but is cleaner...
		/*
		TString lineTs = (TString) line;
		TObjString *temp1, *temp2;
		TObjArray * tokens = lineTs.Tokenize(" ");
		temp1 = (TObjString*) tokens->At(0); 
		pmtNo = temp1->GetString().Atoi();
		temp2 = (TObjString*) tokens->At(14); 
		pmtType = temp2->GetString().Atoi();
		tokens->Clear();
		tokens = NULL;
		*/
		
		// Initialize Pmt type array
		if(pmtNo>OD_PMT_OFFSET && pmtNo<=OD_PMT_OFFSET+odmax) {
			//std::cout << "Read Cable: " << pmtNo << " Type: " << pmtType << std::endl;
			pmtTypes[pmtNo-1-OD_PMT_OFFSET] = pmtType;
		}
	}
	conFile.close(); // close the input file
	
	if(skGeom>=5)	{ // for sk5 and maybe later...
		_file = new TFile("$SKOFL_ROOT/const/od.1pe.sk5.root");
  	_od_sat = (TH1D*) _file->Get("satValHist"); // sat param of each pmt (if <0 then not known)
  	_pc2pe = (TH1D*) _file->Get("pc2peHist"); //pc2pe value of each pmt
  	for (int i = 0; i < nType; i++) {
  		_sat_byType[i] = (TH1D*) _file->Get("satDist_"+typeString[i]);	
  		_pe_byType[i] = (TH1D*) _file->Get("total1pe_"+typeString[i]);	
  		if(_sat_byType[i] == 0 || _pe_byType[i] == 0) typeFlag = true; // warning if no histogram 
  	}
		// Assign appropriate random value to saturation parameters of non-saturated PMTs using
		// distributions, I find it appropriate to do this everytime we call this function rather
		// than doing it only once for all simulations. We do not know the values, and that
		// uncertainty should carry into the simulation I think, if I am wrong you are free to
		// change this...
		
		for (int i = 1; i <= _od_sat->GetNbinsX(); i++) {
			if(_od_sat->GetBinContent(i) < 0) {
				//int tableType = atoi(conTable->GetFlag(OD_PMT_OFFSET+i)); // should be commented out
				int tableType = pmtTypes[i-1];
				// Order: imb, rec (various versions), kek, ipmu, uci, no pmt or some other...
				// If the tableType==<number> looks weird see connection table for types
				if(tableType == 13) {
					if(_sat_byType[1]->Integral() > 0) _od_sat->SetBinContent(i,_sat_byType[1]->GetRandom());
					else _od_sat->SetBinContent(i,_sat_byType[0]->GetRandom());
				}
				else if(tableType>10 && tableType<17 && tableType != 13) {
					if(_sat_byType[5]->Integral() > 0)_od_sat->SetBinContent(i,_sat_byType[5]->GetRandom());
					else _od_sat->SetBinContent(i,_sat_byType[0]->GetRandom());
				}
				else if(tableType == 20) {
					if(_sat_byType[2]->Integral() > 0)_od_sat->SetBinContent(i,_sat_byType[2]->GetRandom());
					else _od_sat->SetBinContent(i,_sat_byType[0]->GetRandom());
				}
				else if(tableType == 18) {
					if(_sat_byType[3]->Integral() > 0)_od_sat->SetBinContent(i,_sat_byType[3]->GetRandom());
					else _od_sat->SetBinContent(i,_sat_byType[0]->GetRandom());
				}
				else if(tableType == 19) {
					if(_sat_byType[4]->Integral() > 0)_od_sat->SetBinContent(i,_sat_byType[4]->GetRandom());
					else _od_sat->SetBinContent(i,_sat_byType[0]->GetRandom());
				}
				else _od_sat->SetBinContent(i,_sat_byType[0]->GetRandom());
				if(_od_sat->GetBinContent(i)<0) _od_sat->SetBinContent(i,0); //make sure nothing remained <0
			} // end of if saturation parameter not measured assign random
		} // end of loop over saturation bins
	} // end of sk5 or after if

	else { // sk4 and before
		_file = new TFile("$SKOFL_ROOT/const/od.1pe.root");
  	_old  = (TH1D*) _file->Get("old1pe");
  	_new  = (TH1D*)	_file->Get("new1pe");
  	_od_sat = (TH1D*) _file->Get("pmt_saturation");
	}

	// Get dark rates
  _file_dark = new TFile("$SKOFL_ROOT/const/od.dark.root");
  _file_dark->cd();
  _od_dark[0] = (TH1D*) _file_dark->Get("od_dark_sk1");
  _od_dark[1] = (TH1D*) _file_dark->Get("od_dark_sk2");
  _od_dark[2] = (TH1D*) _file_dark->Get("od_dark_sk3");
  _od_dark[3] = (TH1D*) _file_dark->Get("od_dark_sk4");
  _od_dark[4] = (TH1D*) _file_dark->Get("od_dark_sk5");


	if(*darkRun > 0) {
		TString darkFileName = Form("$SKOFL_ROOT/const/darkr/darkr.0%d.root",*darkRun);
		_file_dark_cur = new TFile(darkFileName,"READ"); 
		_od_dark_tree = (TTree*) _file_dark_cur->Get("skdark"); 
		_od_dark_tree->SetBranchAddress("dark_rate_od",odDarkTree);
		_od_dark_tree->GetEntry(0);
	}


  for( int i = 0 ; i < 5 ; i++ )
    if( _od_dark[i] == 0 ) flag = true; // for error message
 

  //  for tube by tube spe , same in sk 1-5
  std::stringstream ss;
  _file->cd();
	for( int i = 1 ; i <= odmax ; i++ )
  {
    ss.str(""); ss << "pmt1pe_" << i ;
    _pmt[i-1] = (TH1D*) _file->Get( ss.str().c_str() ); 
		if( _pmt[i-1] == 0 ) flag =true;
  }

	// Error message in various conditions resulting in non accessing data
  if( (skGeom<5 && (! _file || ! _old || !_new || flag || ! _od_sat)) ||
  		(skGeom>=5 && (! _file || typeFlag || flag || ! _od_sat || !_pc2pe)) )
  {
     printf(" Error od single pe histograms do not exist. \n" );
     printf(" Please check that SKOFL_ROOT is set properly and \n" );
     printf(" if problems persist please contact expert: rvw or Baran \n" );
     printf("Sk Geom:%d, peFile:%d, sat:%d, pmt1pe:%d, old:%d, new:%d ,typeDist:%d, pc2pe:%d \n"
     		 ,skGeom,!_file,!_od_sat,flag,!_old,!_new,typeFlag,!_pc2pe);
     abort();
  } 
	printf("od info loaded for skGeom: %d \n",skGeom);
}


///
///  These should all be read in at the begining
///  so we assume existence and simple delete
/// 
extern "C" void close_root_(int * inSkGeom, int * darkRun)
{	
	int skGeom = *inSkGeom; // temporary until sk geometry 5 is added to skdetsim
	if(skGeom>=5) {
		for (int i = 0; i < nType; i++) {
			delete _sat_byType[i];
			delete _pe_byType[i];
		}
		delete _pc2pe;
	}
	else { // for sk4 and before
 		delete _old;
  	delete _new;
  }

	// Common for all, delete pointers
  
  for( unsigned i = 1 ; i <= odmax ; i++ )
    delete _pmt[i-1];

  for( unsigned sk = 0 ; sk < 5 ; sk++ ) // THIS WAS 4, terrible terrible memory leak...
    delete  _od_dark[sk];

	delete _od_dark_tree;

  delete _od_sat;

  _file->Close();
  _file_dark->Close();

  delete( _file );
  delete( _file_dark );

  if(*darkRun>0) _file_dark_cur->Close();
	delete( _file_dark_cur);

	//delete conTable; // should be commented out
}

// call for sk5 and after only

// Gives a pe, from pe distribution of all same type of pmts
// Example, the PMT in cable 1749 is a KEK, so this will return a random pe from KEK spe dist.
// instead of just this PMT's spe distribution 
// cable should be [1, 1885]
extern "C" void od_type1pe_( float * pe, int * cable )	{
	//int tableType = atoi(conTable->GetFlag(OD_PMT_OFFSET+*cable)); // should be commented out
	int tableType = pmtTypes[*cable-1];
	
	// Order: imb, rec (various versions), kek, ipmu, uci, no pmt or some other...
	// If the tableType==<number> looks weird see connection table for types
	if(tableType == 13) *pe = _pe_byType[1]->GetRandom();
	else if(tableType>10 && tableType<17 && tableType != 13) *pe = _pe_byType[5]->GetRandom();
	else if(tableType == 20) *pe = _pe_byType[2]->GetRandom();
	else if(tableType == 18) *pe = _pe_byType[3]->GetRandom();
	else if(tableType == 19) *pe = _pe_byType[4]->GetRandom();
	else if(tableType == 99) *pe = 0; // No PMT
	else *pe = _pe_byType[0]->GetRandom();
  *pe = (float) *pe;
	//cout << "Cable: " << *cable << " Type: " << tableType << " pe: " << *pe << endl; // To test
}

// call for sk5 and after only
// get pc2pe value (used as calibration constant) for a given cable
// can be used instead of 1 hard coded parameter in ds_anti_sk4.F for all pmts
extern "C" void od_pc2pe_( float * pc2pe, int * cable )	{
  // root histogram binning is [1, NbinsX] so do adjustment
  *pc2pe = (float ) _pc2pe->GetBinContent( *cable );
}

// call for sk4 and before only
extern "C" void od_old1pe_( float * pe )
{
    *pe = (float ) _old->GetRandom();
}

// call for sk4 and before only
extern "C" void od_new1pe_( float * pe )
{
    *pe = (float ) _new->GetRandom();
}


//cable should be [1, 1885]
extern "C" void od_pmt1pe_( float * pe , int * cable )
{
    int integral = _pmt[ *cable -1 ]->Integral() ; 
    // suppress root errors form empty histograms
    if( integral <1000 )
    { 
			//Get From Type
      od_type1pe_(pe,cable);
      return;
      // will be sqaushed by discriminator
      //*pe = -1.0;
      // return;
    }

    *pe = (float ) _pmt[ *cable -1 ]->GetRandom();
}


// cable should be [1, 1885]
// skx should be [1,4] -> Update: should work for 5 too (Baran)
extern "C" void od_pmt_dark_( float * rate, int *cable, int * skx, int *darkRun )
{
	if(*darkRun<=0) {
  	// root histogram binning is [1, NbinsX] so do not adjustment
  	*rate = (float ) _od_dark[*skx-1]->GetBinContent( *cable );
	}
	else { // from darkr
		*rate = (float) odDarkTree[*cable-1]*1000;	// check x1000 multiplier
	}
	//Test
	//if(*cable==647) cout << "Rate of OD channel 647: " << *rate << endl;
}

//cable should be [1, 1885]
extern "C" void od_pmt_sat_( float * sat, int *cable )
{

  // root histogram binning is [1, NbinsX] so do adjustment
  *sat = (float ) _od_sat->GetBinContent( *cable );

}

// These 3 functions for outputing
extern "C" void open_root_output_( int * instance )
{

  char name[100]; 

  if( ! _out ) 
  {
     sprintf(name, "$odtune/out/root/sim/treeout.%d.root", (*instance) );	
     _out = new TFile( name, "recreate" );
     tree = new TTree("otree","tuning stuff" );
     tree->Branch("cable"   , &local_cable , "cable/I"   );
     tree->Branch("photons" , &photons     , "photons/I" );
     tree->Branch("q"       , &q           , "q/F" );
  }

}



extern "C" void output_fill_( int * cable, int * n , float * lq )
{

   q = *lq; 
   local_cable = *cable;
   photons = *n ;
   _out->cd();
   tree->Fill();
}

extern "C" void close_root_output_(  )
{
   _out->cd();
   tree->Write();

   _out->Close();
   delete _out;

}




