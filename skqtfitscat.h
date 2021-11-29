#ifndef skqtfitscat_hxx_seen
#define skqtfitscat_hxx_seen

class skqtfitscat {

 public:

  float pi;
  
  int nzsbins;
  float minzs;
  float maxzs;

  int nrsbins;
  float minrs;
  float maxrs;

  int nastbins;
  float minast;
  float maxast;

  int nctbins;
  float minct;
  float maxct;

  int nphibins;
  float minphi;
  float maxphi;

  int nztbins;
  float minzt;
  float maxzt;

  int nrtbins;
  float minrt;
  float maxrt;

  int side_scattable_size;
  int* side_scattable;
  int end_scattable_size;
  int* end_scattable;

  skqtfitscat();
  ~skqtfitscat();

  int GetIndexSide(int izt, int izs, int irs, int iast, int ict, int iphi){
    return izt + nztbins*(izs + nzsbins*(irs + nrsbins*(iast + nastbins*(ict + nctbins*(iphi)))));
  }

  int GetIndexEnd(int irt, int izs, int irs, int iast, int ict, int iphi){
    return irt + nrtbins*(izs + nzsbins*(irs + nrsbins*(iast + nastbins*(ict + nctbins*(iphi)))));
  }

  void AddToScatTable(int nopmt);
  void WriteScatTable();

};

#endif
