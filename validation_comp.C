
void validation_comp() {

  TString inputroot1 = "validation_r29076_20k_g77.root" ;
  TString inputroot2 = "validation_r29076_20k.root" ;

  cout << "Input file 1 " << inputroot1 << endl ; 
  cout << "Input file 2 " << inputroot2 << endl ; 

  TFile *f1 = new TFile(inputroot1);
  TFile *f2 = new TFile(inputroot2);

  TH1F *h1[10] ; 
  TH1F *h2[10] ; 

  for (Int_t i=0 ; i<5 ; i++ ) {
    cout << TString().Format("h_%d",i).Data() << endl ; 
    h1[i] = (TH1F*)f1->Get( TString().Format("h_%d",i).Data() ) ; 
    h2[i] = (TH1F*)f2->Get( TString().Format("h_%d",i).Data() ) ; 
  }

  TCanvas *c1 = new TCanvas();
  c1->SetGridx( 0 );
  c1->SetGridy( 0 );
  c1->SetFillColor( 10 );
  c1->SetLeftMargin( 0.15 );
  c1->SetRightMargin( 0.15 );
  c1->SetTopMargin( 0.15 );
  c1->SetBottomMargin( 0.15 );
  c1->SetTickx( 0 );
  c1->SetTicky( 0 );

  for ( Int_t i=0 ; i<5 ; i++ ) {
    h1[i]->Draw("");
    h2[i]->SetLineColor(kBlue);
    h2[i]->Draw("same");
    c1->Print( TString().Format("validation_h%d.png",i) );
  }



}
