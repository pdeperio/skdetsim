#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "skqtfitscat.h"
#include "vcvrtxC.h"
#include "vcworkC.h"
#include "skheadC.h"
#include "skparmC.h"
#include "geopmtC.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"
#include "TVector3.h"

//skqtfitscat* scattable;

TFile* outfile;
TTree* sttree;
struct hitphtinfo {
  int iPMT;
  int iSTRY;
};
struct ScatVars {
  float zs;
  float rs;
  float ast;
  float ct;
  float phi;
  float zt;
  float rt;
  float ctcor;
  float phicor;
} scatvars;
float srcpos[3],srcdir[3];
int ihPMT,isct;

std::vector<TVector3> photoninipos;
//std::vector<int> hitpmtlist;
std::vector<hitphtinfo> hphtlist;

extern"C" {
  void initscattable_();
  void recordphotoninipos_(float* x, float* y, float* z);
  void addtoscattable_(int* nopmt, int* ISTRY);
  void fillscattable_();
  void writescattable_();
}

void initscattable_(){
  //  scattable = new skqtfitscat();
  outfile = new TFile("sttree.root","RECREATE");
  sttree = new TTree("sttree","Source-Tube Tree");
  sttree->Branch("scatvars",&scatvars.zs,"zs/F:rs:ast:ct:phi:zt:rt:ctcor:phicor");
  sttree->Branch("srcpos",srcpos,"srcpos[3]/F");
  sttree->Branch("srcdir",srcdir,"srcdir[3]/F");
  sttree->Branch("ihPMT",&ihPMT,"ihPMT/I");
  sttree->Branch("isct",&isct,"isct/I");
//  sttree->Branch("zs",&zs,"zs/F");
//  sttree->Branch("rs",&rs,"rs/F");
//  sttree->Branch("ast",&ast,"ast/F");
//  sttree->Branch("ct",&ct,"ct/F");
//  sttree->Branch("phi",&phi,"phi/F");
//  sttree->Branch("zt",&zt,"zt/F");
//  sttree->Branch("rt",&rt,"rt/F");
}

void recordphotoninipos_(float* x, float* y, float* z){
  photoninipos.push_back(TVector3(*x,*y,*z));
}

void addtoscattable_(int* nopmt, int* ISTRY) {
  hitphtinfo tmpphinfo;
  tmpphinfo.iPMT=*nopmt;
  tmpphinfo.iSTRY=*ISTRY;
  
  hphtlist.push_back(tmpphinfo);
}

void fillscattable_() {
  //  scattable->AddToScatTable(*nopmt);

  //  std::cout << "# photons in, # hit tubes = " << photoninipos.size() << ", " << hitpmtlist.size() << std::endl;

  TVector3 avgpos(0,0,0);
  for (std::vector<TVector3>::iterator iphotpos = photoninipos.begin(); iphotpos != photoninipos.end(); iphotpos++) {
    avgpos += *iphotpos;
  }
  avgpos *= 1./photoninipos.size();
  photoninipos.clear();
  
  double tmpmom=sqrt(vcwork_.pvc[0][0]*vcwork_.pvc[0][0]+vcwork_.pvc[0][1]*vcwork_.pvc[0][1]+vcwork_.pvc[0][2]*vcwork_.pvc[0][2]);
  
  for (int j=0; j<3; j++) {
    srcpos[j]=vcwork_.posivc[0][j];//vertex
    srcdir[j]=vcwork_.pvc[0][j]/tmpmom;//direction
  }

  TVector3 sourcepos(srcpos[0],srcpos[1],srcpos[2]);
  TVector3 sourcedir(srcdir[0],srcdir[1],srcdir[2]);
  TVector3 sourcecordir = (avgpos - sourcepos).Unit();

  for (std::vector<hitphtinfo>::iterator ipht = hphtlist.begin(); ipht != hphtlist.end(); ipht++) {
    ihPMT=(*ipht).iPMT-1;
    TVector3 tubepos(geopmt_.xyzpm[ihPMT][0],geopmt_.xyzpm[ihPMT][1],geopmt_.xyzpm[ihPMT][2]);

    scatvars.zs = sourcepos.Z();
    scatvars.rs = sourcepos.Perp();
    scatvars.ast = tubepos.DeltaPhi(sourcepos);
    scatvars.ct = srcdir[2];
    scatvars.phi = sourcedir.DeltaPhi(tubepos-sourcepos);
    scatvars.zt = tubepos.Z();
    scatvars.rt = tubepos.Perp();
    
    isct = (*ipht).iSTRY;

    sttree->Fill();

  }

  hphtlist.clear();

}

void writescattable_(){

//  scattable->WriteScatTable();
//  delete scattable;
  sttree->Write();
  outfile->Close();

}

skqtfitscat::skqtfitscat() {

    pi = 3.1415926535897932384626433832795;

    //    nzsbins = 71;
    nzsbins = 1;
    minzs = -1750.;
    maxzs = 1750.;

    nrsbins = 34;
    minrs = 0.;
    maxrs = 1650.;

    nastbins = 32;
    minast = -pi;
    maxast = pi;

    nctbins = 32;
    minct = -1.;
    maxct = 1.;

    nphibins = 32;
    minphi = -pi;
    maxphi = pi;

    nztbins = 71;
    minzt = -1750.;
    maxzt = 1750.;

    nrtbins = 34;
    minrt = 0.;
    maxrt = 1650.;

    side_scattable_size = nztbins*nzsbins*nrsbins*nastbins*nctbins*nphibins;
    end_scattable_size = nrtbins*nzsbins*nrsbins*nastbins*nctbins*nphibins;

    side_scattable = new int[side_scattable_size];
    end_scattable = new int[end_scattable_size];

}

skqtfitscat::~skqtfitscat() {

    delete[] side_scattable;
    delete[] end_scattable;

}

void skqtfitscat::AddToScatTable(int nopmt) {

  double tubex = geopmt_.xyzpm[nopmt-1][0];
  double tubey = geopmt_.xyzpm[nopmt-1][1];
  double tubez = geopmt_.xyzpm[nopmt-1][2];
  double tubercyl = sqrt(tubex*tubex+tubey*tubey);

  double sourcex = vcvrtx_.pvtxvc[0][0];
  double sourcey = vcvrtx_.pvtxvc[0][1];
  double sourcez = vcvrtx_.pvtxvc[0][2];

  std::cout << "tubepos (" << tubex << "," << tubey << "," << tubez << "), sourcepos (" << sourcex << "," << sourcey << "," << sourcez << ")" << std::endl;
  std::cout << "sourcepos2 = (" << vcwork_.posvc[0] << "," << vcwork_.posvc[1] << "," << vcwork_.posvc[2] << ")" << std::endl;

}

void skqtfitscat::WriteScatTable() {

  std::cout << "side_scattable_size, end_scattable_size = " << side_scattable_size << ", " << end_scattable_size << std::endl;

  std::cout << "opening output file" << std::endl;
  FILE* outfile = fopen("scattables.bin", "rb");
  std::cout << "writing end table" << std::endl;
  fwrite(end_scattable, sizeof(end_scattable[0]), end_scattable_size, outfile);
  std::cout << "writing side table" << std::endl;
  fwrite(side_scattable, sizeof(side_scattable[0]), side_scattable_size, outfile);
  std::cout << "closing file" << std::endl;
  fclose(outfile);

}
