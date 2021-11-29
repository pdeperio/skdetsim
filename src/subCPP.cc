//----------------------------------------------------------------------------
//  Fortran-C++ interface to utilize GEANT4 functions. 
//

#include "userInterface2.hh"
#include "Randomize.hh"
#include <iostream>

typedef double secvec[ 100 ][ 9 ];

extern "C" 
{
  void g4_initialize_( int&, int&, int &, float & );
  void g4_finalize_();
  void g4_getrndm_( double& ); 
  void g4_getgdxsec_( float&, double&, double&, double& );
  void g4_gethxsec_( float&, double&, double&, double& );
  void g4_getoxsec_( float&, double&, double&, double& );
  void g4_run_( float&, float&, float&, float&, float&, float&, secvec&, int& );
  void g4_pimicap_(float*, float*, float*, int&, secvec&); 
  void g4_mumicap_(float*, float*, float*, int&, secvec&);
}

userInterface2* testUI ;


//----------------------------------------------------------------------------
//  Initialize & Finalize GEANT4. 
//
/*****************************************************************************/
 void g4_initialize_( int& iseed, int& gd_water, int& gd_model, float& gd_concentration )
/*****************************************************************************/
{
  std::cout << "(TEST) GEANT4 initialize, iseed : "<<iseed<<", GD_WTR : "<<gd_water<<", gd_conc : "<<gd_concentration << std::endl;  
  testUI = new userInterface2( iseed, gd_water, gd_model, gd_concentration );
}

/*****************************************************************************/
 void g4_finalize_()
/*****************************************************************************/
{
  testUI->Finalize();
  std::cout << "GEANT4 finalize" << std::endl; 
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
//  Get a uniform random number later used to determine interaction type. 
//
/*****************************************************************************/
 void g4_getrndm_( double& aRndm ) 
/*****************************************************************************/
{
  aRndm = G4UniformRand();
  //std::cout << "TEST (random): " << aRndm << std::endl;
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
//  Get cross sections from GEANT4. 
//   - g4_gethxsec()  for hydrogen 
//   - g4_getoxsec()  for oxygen 
//   - g4_getgdxsec() for gadolinium
//
//  Unit: 
//   - nEnergy [GeV] : neutron kinetic energy 
//   - cap/ela/inelaXsec [/cm] : cross section [cm^2] * density [/cm^3]
//
/*****************************************************************************/
 void g4_getgdxsec_( float& nEnergy, 
                     double& capXsec, double& elaXsec, double& inelaXsec ) 
/*****************************************************************************/
{
  capXsec = testUI->GetGdCapXsec( nEnergy );
  elaXsec = testUI->GetGdElaXsec( nEnergy );
  inelaXsec = testUI->GetGdInElaXsec( nEnergy );
//  std::cout << "TEST (Gd): " << nEnergy << " " 
//            << capXsec << " " << elaXsec << " " << inelaXsec << std::endl;  
}

/*****************************************************************************/
 void g4_gethxsec_( float& nEnergy, 
                    double& capXsec, double& elaXsec, double& inelaXsec ) 
/*****************************************************************************/
{
  capXsec = testUI->GetHCapXsec( nEnergy );
  elaXsec = testUI->GetHElaXsec( nEnergy );
  inelaXsec = testUI->GetHInElaXsec( nEnergy );
//  std::cout << "TEST (H): " << nEnergy << " " 
//            << capXsec << " " << elaXsec << " " << inelaXsec << std::endl;  
}

/*****************************************************************************/
 void g4_getoxsec_( float& nEnergy, 
                    double& capXsec, double& elaXsec, double& inelaXsec ) 
/*****************************************************************************/
{
  capXsec = testUI->GetOCapXsec( nEnergy );
  elaXsec = testUI->GetOElaXsec( nEnergy );
  inelaXsec = testUI->GetOInElaXsec( nEnergy );
//  std::cout << "TEST (O): " << nEnergy << " " 
//            << capXsec << " " << elaXsec << " " << inelaXsec << std::endl;  
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
//  Get secondary products from the selected interaction. 
//
/*****************************************************************************/
 void g4_run_( float& posx, float& posy, float& posz, 
               float& px,   float& py,   float& pz, 
	       secvec& secVec, int& process )
/*****************************************************************************/
{
  std::vector< Secondary_Vector > secondaryVector;
  secondaryVector = testUI->Get2nd( px, py, pz, process );

  int nSecondary = secondaryVector.size();
  int iterMax = nSecondary;
  if ( iterMax > 100 ) iterMax = 100; 

	//if(process==20)std::cout<<"Process : "<<process<<std::endl;
  for ( int iter = 0; iter < iterMax; iter++ ) 
  {
    Secondary_Vector vec = secondaryVector[ nSecondary - 1 - iter ];
    
    // obtain shorter pdgID for ions
    if ( vec.pdgID > 1000000000 ) vec.pdgID -= 1000000000;
    
    // get kinematics 
		//if(process == 20)std::cout<<iter<<"  "<<vec.pdgID<<"  "<<vec.totalEnergy<<std::endl;

    secVec[ iter ][ 0 ] = vec.pdgID;
    secVec[ iter ][ 1 ] = vec.productionTime;
    secVec[ iter ][ 2 ] = vec.productionPosX;
    secVec[ iter ][ 3 ] = vec.productionPosY;
    secVec[ iter ][ 4 ] = vec.productionPosZ;
    secVec[ iter ][ 5 ] = vec.momentumX;
    secVec[ iter ][ 6 ] = vec.momentumY;
    secVec[ iter ][ 7 ] = vec.momentumZ;
    secVec[ iter ][ 8 ] = vec.totalEnergy;
  }
  std::vector< Secondary_Vector >().swap(secondaryVector);
}
//----------------------------------------------------------------------------

int getPDG2G3ID(int pdg)
{
//
// Convert PDG particle ID to Geant3 one
//
    const int PDG_NEUTRON	=2112;		
    const int PDG_PROTON	=2212;
    const int PDG_DEUTERON	=10020; 
    const int PDG_TRITON	=10030;	
    const int PDG_HE3		=20030;
    const int PDG_ALPHA		=20040;
    const int PDG_GAMMA		=22;
    
    const int G3_NEUTRON	=13;		
    const int G3_PROTON		=14;
    const int G3_DEUTERON	=45;
    const int G3_TRITON		=46;
    const int G3_HE3		=49;
    const int G3_ALPHA		=47;
    const int G3_GAMMA		=1;
    
    const int NUCLEI		=-1;
	int g3=-99;

    if( pdg>1000000000 )
	{
		pdg-=1000000000;
	}

	switch ( pdg )
	{
		case PDG_NEUTRON : 
			g3	=G3_NEUTRON;	break;
		case PDG_PROTON :
			g3	=G3_PROTON;		break;
		case PDG_DEUTERON :
			g3	=G3_DEUTERON; 	break;
		case PDG_TRITON :
			g3	=G3_TRITON;		break;
		case PDG_HE3 :
			g3	=G3_HE3;		break;
		case PDG_ALPHA :
			g3	=G3_ALPHA;		break;
		case PDG_GAMMA :
			g3	=G3_GAMMA;		break;
		default :
			g3	=NUCLEI;		break;
	}
	return g3;
}

//
// Run Geant4 pi- capture on 16O 
//
void g4_pimicap_(float *px,
			     float *py,
				 float *pz,
				 int &nscnd,
				 secvec &secVec)
{
	// Do pi- capture and then get resultand secondaries
	std::vector< Secondary_Vector > secondaryVector;
	secondaryVector	=testUI->DoPionMinusCapture((G4double)*px, (G4double)*py, (G4double)*pz);

	//
	// Only the below particles are brought back into SKDETSIM
	// Basically, residual nuclei are ignored.
	//
	// --- Supported particles ---
	// 	1. neutron
	// 	2. proton
	// 	3. deutron
	// 	4. triton
	// 	5. He3
	// 	6. alpha
	// 	7. gamma

	int ntot	=-1,
		nselec	=-1,
		g3id	=-99;
	ntot	=secondaryVector.size();
	nselec	=0;

	for(int i=0; i<ntot; i++)
	{
		Secondary_Vector vec=secondaryVector[i];

		g3id	=-99;
		g3id	=getPDG2G3ID( vec.pdgID );
		if( g3id>=0 )
		{
			secVec[nselec][0]	=g3id;
			secVec[nselec][1]	=vec.productionTime;
			secVec[nselec][2]	=vec.productionPosX;
			secVec[nselec][3]	=vec.productionPosY;
			secVec[nselec][4]	=vec.productionPosZ;
			secVec[nselec][5]	=vec.momentumX;
			secVec[nselec][6]	=vec.momentumY;
			secVec[nselec][7]	=vec.momentumZ;
			secVec[nselec][8]	=vec.totalEnergy;
			nselec				+=1;
		}
	}
	nscnd=nselec;
}

void g4_mumicap_(float *px,
			     float *py,
				 float *pz,
				 int &nscnd,
				 secvec &secVec)
{
static const int	G3_NEUTRON			=13;		

	// Do mu- capture and then get resultant particles
	// 
	// Since both of mu- decay and mu- capture are managed with a single function, 
    // ** DoMuonMinusCapture ** is repeated until its result become a mu- capture.
	int  nLoop=0,
		 Nnuebar=0,
		 Nnumu=0;
	bool isDec=true;
	std::vector<Secondary_Vector> secondaryVector;
	while ( isDec )
	{
		Nnuebar	=0;
		Nnumu	=0;
		secondaryVector.clear();
		secondaryVector=testUI->DoMuonMinusCapture((G4double)*px, (G4double)*py, (G4double)*pz);

		for(int i=0; i<(int)secondaryVector.size(); i++)
		{
			if( secondaryVector[i].pdgID==-12 ){ Nnuebar+=1; }
			if( secondaryVector[i].pdgID==14  ){ Nnumu  +=1; }
		}

		if( nLoop>1000 )
		{
			std::cout<<" More than 1000 trials of userInterface2::DoMuonMinusCapture " <<std::endl;
			std::cout<<" -> EXIT " <<std::endl;
			exit( -1 );
		}

		if( Nnuebar>=1 && Nnumu>=1 ){ isDec=true;  		   } // This call resulted in calling decay
		else						{ isDec=false;	break; } // capture!
		nLoop+=1;
	}	

	// Only neutrons are brought back into SKDETSIM, 
    // in order not to affect existing gamma ray emission
    int ntot	=-1,
	nselec	=-1,
	g3id	=-99;
	ntot	=secondaryVector.size();
	nselec	=0;
	for(int i=0; i<ntot; i++)
	{
		Secondary_Vector vec=secondaryVector[i];

		g3id	=-99;
		g3id	=getPDG2G3ID( vec.pdgID );
		if( g3id==G3_NEUTRON )
		{
			secVec[nselec][0]	=g3id;
			secVec[nselec][1]	=vec.productionTime;
			secVec[nselec][2]	=vec.productionPosX;
			secVec[nselec][3]	=vec.productionPosY;
			secVec[nselec][4]	=vec.productionPosZ;
			secVec[nselec][5]	=vec.momentumX;
			secVec[nselec][6]	=vec.momentumY;
			secVec[nselec][7]	=vec.momentumZ;
			secVec[nselec][8]	=vec.totalEnergy;
			nselec				+=1;
		}
	}
	nscnd=nselec;
}
