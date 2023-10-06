#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "Inclusion.h"

const double FREQ{160.316E6};
const double TDC2ns = 1E9/FREQ;
const int NSIPM{8};
const int NSIDE{2};
const int NsidesNch{16};

void definePlots( std::map<std::string, TH1*> &m_plots, std::vector<std::string> &tags) {
  m_plots["histo"] = new TH2F("histo", "; # scifi hits; # mu hits", 100, 0, 100, 100, 0, 100);

  //events characteristics plots
  m_plots["EventTimeDistribution"] = new TH1D("EventTimeDistribution", "EventTimeDistribution; time (ns); events", 6E9, 0, 6E9); 

  //basic quantities plots for scifi and mu 
  for (auto tag : tags) {
    const auto t{tag.c_str()};
    m_plots[Form("%s_times", t)] = new TH1F(Form("%s_times", t), Form("%s_times; time (ns) ; entries", t), 150, -5, 145);
    m_plots[Form("%s_signals", t)] = new TH1F(Form("%s_signals", t), Form("%s_signals; qdc? ; entries", t), 100, -30, 80);
    m_plots[Form("%s_channels", t)] = new TH1F(Form("%s_channels", t), Form("%s_channels; n channel; entries", t), 100, 0, 100);
    m_plots[Form("%s_nSIPM", t)] = new TH1F(Form("%s_nSIPM", t), Form("%s_nSIPM; n fired SIPM ?; entries", t), 1000, 0, 1000);
    m_plots[Form("%s_charge", t)] = new TH1F(Form("%s_charge", t), Form("%s_charge; qdc; entries", t), 200, 0, 200);
    m_plots[Form("%s_station", t)] = new TH1F(Form("%s_station", t), Form("%s_station; station ; entries", t), 6, -0.5, 5.5);
  }
}

void runAnalysis() //(int runN, int partN)
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

  
  std::map<std::string, TH1*> plots;
  std::vector<std::string> tags;
  tags.push_back("Scifi");
  tags.push_back("MuFilter");

  definePlots(plots, tags);
  
  auto *fEventTree = new TChain("rawConv");
  TFile outputFile("output.root", "RECREATE"); 
  outputFile.cd();

  fEventTree->Add("root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/run_100678/sndsw_raw-0000.root");

  auto mu_hits = new TClonesArray("MuFilterHit");
  fEventTree->SetBranchAddress("Digi_MuFilterHits", &mu_hits);
  auto sf_hits = new TClonesArray("sndScifiHit");
  fEventTree->SetBranchAddress("Digi_ScifiHits", &sf_hits);
  auto header = new SNDLHCEventHeader();
  fEventTree->SetBranchAddress("EventHeader.", &header);
  
  // Loop over events
  int iMax = fEventTree->GetEntries();
  for ( int m =0; m < iMax; m++ ){ // iMax=fEventTree->GetEntries()
    if (m % 100 == 0) std::cout << "Processing event: " << m << '\r' << std::flush;
    fEventTree->GetEntry(m);
    int sf_max=sf_hits->GetEntries();
    int mu_max=mu_hits->GetEntries();

    //remove almost empty events
    if (sf_max < 2 || mu_max < 2) continue;

    plots["EventTimeDistribution"]->Fill(header->GetEventTime());
    plots["histo"]->Fill(sf_max, mu_max);

    for (int i=0 ; i<sf_max; i++) {
        auto t = tags[0].c_str();
        auto sf_hit = (sndScifiHit*) sf_hits->At(i);
        
        int station = sf_hit->GetStation();
        plots[Form("%s_station", t)]->Fill(station);


        // for each event->loop over all SIPMs in one side
        for (int j = 0;  j<NsidesNch; ++j){
          int channel = sf_hit->Getchannel(j);
          if (channel > -1) plots[Form("%s_channels", t)]->Fill(channel);
          float charge = sf_hit->GetSignal(j);
          if (charge > -1) plots[Form("%s_charge", t)]->Fill(charge);
          float time = sf_hit->GetTime(j)*TDC2ns;
          plots[Form("%s_times", t)] ->Fill(time);
        }

    }

    for (int i=0 ; i<mu_max; i++) {
        auto t = tags[1].c_str();
        auto mu_hit = (MuFilterHit*) mu_hits->At(i);
        for (int j = 0; j<NsidesNch; ++j){
          int channel = mu_hit->Getchannel(j);
          if (channel > -1) plots[Form("%s_channels", t)]->Fill(channel);
          float charge = mu_hit->GetSignal(j);
          if (charge > -1) plots[Form("%s_charge", t)]->Fill(charge);
        }
        int station = mu_hit->GetPlane();
        plots[Form("%s_station", t)]->Fill(station);
        float time = mu_hit->GetTime()*TDC2ns;
        plots[Form("%s_times", t)] ->Fill(time);
    }
  }


  auto *c = new TCanvas("c", "c",  500, 500, 500, 500);
  c->Divide(2, 2);
  c->cd(1);
  plots[Form("Scifi_channels")]->Draw();
  c->cd(2);
  plots[Form("MuFilter_channels")]->Draw();
  c->cd(3);
  plots[Form("Scifi_charge")]->Draw();
  c->cd(4);
  plots[Form("MuFilter_charge")]->Draw();

  auto *c1 = new TCanvas("c1", "c1", 500, 500, 500, 500);
  c1->Divide(2, 2);
  c1->cd(1);
  plots["Scifi_station"]->Draw();
  c1->cd(2);
  plots[Form("MuFilter_station")]->Draw();
  c1->cd(3);
  plots[Form("Scifi_times")]->Draw();
  c1->cd(4);
  plots[Form("MuFilter_times")]->Draw();

  auto *c2 = new TCanvas("c2", "c2", 500, 500, 500, 500);
  plots["EventTimeDistribution"]->Draw();
  
  std::cout << "\nDone" << std::endl;
  outputFile.Write();
  outputFile.Close();
}
