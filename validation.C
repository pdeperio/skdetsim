//
//  ROOT script to make TQ distribution from validation output
//
//  How to run:
//   %  root -l -q -b validation.C\(\"validation.dat\"\)
//
//

void validation( char* inputfile = "validation.dat" ) {

  cout << "Input file is " << inputfile << endl ; 

  std::ifstream fin( inputfile ) ; 

  std::string line ;
  TString sline ;
  Int_t iseq, n=0 ;
  Double_t t, q, qsum[2];
  TString s;  

  TH1F *h[10] ; 
  h[0] = new TH1F("h_0", "qisk", 100, 0., 10.);
  h[1] = new TH1F("h_1", "qsum0", 100, 0., 5000.);
  h[2] = new TH1F("h_2", "qsum1", 100, 0., 300.);
  h[3] = new TH1F("h_3", "tisk", 100, -20000., 50000.);
  h[4] = new TH1F("h_4", "tisk", 100, 500., 1500.);

  while ( !fin.eof() ) {

    std::getline( fin, line );

    if ( line.size() < 5 ) {
      continue ; 
    }
    if ( line.substr(1,5).compare("Event") == 0 ) {
      //cout << line << endl ; 
      qsum[0] = 0 ; 
      qsum[1] = 0 ; 
      continue ; 
    }
    if ( line.substr(1,2).compare("ID") == 0 ) {
      //cout << line << endl ; 
      continue ; 
    }
    if ( line.substr(1,2).compare("OD") == 0 ) {
      //cout << line << endl ; 
      h[1]->Fill( qsum[0] ) ; 
      h[2]->Fill( qsum[1] ) ; 
      continue ; 
    }

    n ++ ; 
    sline = (const char*)line.c_str() ; 

    s = sline(0,7) ; 
    iseq = s.Atoi() ;
    s = sline(7,8) ; 
    q    = s.Atof() ;
    s = sline(15,9) ; 
    t    = s.Atof() ;
    h[0]->Fill( q ) ; 
    h[3]->Fill( t ) ; 
    h[4]->Fill( t ) ; 
    qsum[0] += q ; 
    if ( t>900. && t<1200 )  qsum[1] += q ; 
    //cout << iseq << " " << q << " " << t << endl ;
    //printf("%d %f %f\n", iseq,q,t);


  }


  TFile *fout = new TFile("validation.root","recreate");
  for ( Int_t i=0 ; i<5 ; i++ ) {
    h[i]->Write();
  }

  cout << "h_2 Entries " << h[2]->GetEntries() << endl ;
  cout << "h_2 mean    " << h[2]->GetMean() << endl ;
  cout << "h_2 RMS     " << h[2]->GetRMS() << endl ;

}
