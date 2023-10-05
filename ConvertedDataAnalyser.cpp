#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "Inclusion.h"

#include <vector>

void offlineQA() //(int runN, int partN)
{

  gStyle->SetTitleOffset(0.95, "X");
  gStyle->SetTitleOffset(0.95, "Y");
  gStyle->SetTitleSize(0.055, "XY");
  gStyle->SetLabelSize(0.05, "XY");
  gStyle->SetHistLineWidth(2);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetNdivisions(510, "XY");
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);   
  gROOT->ForceStyle();

  auto *histo = new TH2F("histo", "; # scifi hits; # mu hits", 100, 0, 100, 100, 0, 100);
  
  auto *fEventTree = new TChain("rawConv");
  char *outputName;
  
  fEventTree->Add("root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/run_100678/sndsw_raw-0000.root");

  auto mu_hits = new TClonesArray("MuFilterHit");
  fEventTree->SetBranchAddress("Digi_MuFilterHits", &mu_hits);
  auto sf_hits = new TClonesArray("sndScifiHit");
  fEventTree->SetBranchAddress("Digi_ScifiHits", &sf_hits);
  auto header = new SNDLHCEventHeader();
  fEventTree->SetBranchAddress("EventHeader.", &header);
  
  //cout<<fEventTree->GetEntries()<<endl;
  // Loop over events
  int iMax = fEventTree->GetEntries();
  for ( int m =0; m < iMax; m++ ){ // iMax=fEventTree->GetEntries()
    if (m % 100 == 0) std::cout << "Processing event: " << m  << std::endl; //<< '\r'
    fEventTree->GetEntry(m);
    int sf_max=sf_hits->GetEntries();
    int mu_max=mu_hits->GetEntries();
    if (sf_max > 2 && mu_max> 2) histo->Fill(sf_max, mu_max);

    if (m%100 == 0 ) cout << header->GetEventNumber() << endl;
    
    /*cout << "*****************" << endl;
    for (int i=0 ; i<sf_max; i++) {
        auto sf_hit = (sndScifiHit*) sf_hits->At(i);
        sf_hit->Print();

    }
    for (int i=0 ; i<mu_max; i++) {
        auto mu_hit = (MuFilterHit*) mu_hits->At(i);
        mu_hit->Print();

    }*/
  }

  auto *c = new TCanvas("c", "c",  500, 500, 500, 500);
  histo->Draw("box");
}
