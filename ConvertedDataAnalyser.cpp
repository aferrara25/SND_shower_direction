#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "Inclusion.h"

const double FREQ{160.316E6};
const double TDC2ns = 1E9/FREQ;
const int NSIPM{8};
const int NSIDE{2};
const int NsidesNch{16};
const int TOFPETperBOARD{8};
const int SCIFIMINHITS{2};
const int MUMINHITS{2};

namespace TB{
  const int SCIFISTATION{4};  
  const int MUSTATION{5};
  const int NWALLS{3};
  const int SCIFITHRESHOLD{56};
}

namespace TI18{
  const int SCIFISTATION{5};
  const int MUSTATION{8};
  const int NWALLS{5};
  const int SCIFITHRESHOLD{40};
}

void definePlots( std::map<std::string, TH1*> &m_plots, std::map<std::string, double> &m_counters, std::vector<std::string> &tags) {
  m_plots["histo"] = new TH2F("histo", "; # scifi hits; # mu hits", 100, 0, 100, 100, 0, 100);

  //events characteristics plots
  //m_plots["EventTimeDistribution"] = new TH1D("EventTimeDistribution", "EventTimeDistribution; time (ns); events", 6E9, 0, 6E9); 
  m_plots["ShowerStart"] = new TH1D("ShowerStart", "ShowerStart; station; entries", 8, -2, 5);
  m_plots["ShowerStartProbability"] = new TH1D("ShowerStartProbability", "ShowerStartProbability; station; probability", 6, -.5, 5.5);

  // together x and y
  for (int st = 0; st < TB::SCIFISTATION; ++st){
    m_plots[Form("TimeResidualStation_%d", st)] = new TH1D (Form("TimeResidualStation_%d", st), Form("TimeResidualStation_%d;t_x - t_y;entries", st), 200, -100, 100);
  }
  // separately x and y 
  for (int st = 0; st < 2*TB::SCIFISTATION; ++st){
    m_plots[Form("HitsperStation_%d", st)] = new TH1D (Form("HitsperStation_%d", st), Form("HitsperStation_%d;station;entries", st), 500, 0, 500);
  }
 
  //basic quantities plots for scifi and mu 
  for (auto tag : tags) {
    const auto t{tag.c_str()};
    m_plots[Form("%s_times", t)] = new TH1F(Form("%s_times", t), Form("%s_times; time (ns) ; entries", t), 150, -5, 145);
    m_plots[Form("%s_signals", t)] = new TH1F(Form("%s_signals", t), Form("%s_signals; qdc? ; entries", t), 100, -30, 80);
    m_plots[Form("%s_channels", t)] = new TH1F(Form("%s_channels", t), Form("%s_channels; n channel; entries", t), 100, 0, 100);
    m_plots[Form("%s_nSIPM", t)] = new TH1F(Form("%s_nSIPM", t), Form("%s_nSIPM; n fired SIPM ?; entries", t), 1000, 0, 1000);
    m_plots[Form("%s_charge", t)] = new TH1F(Form("%s_charge", t), Form("%s_charge; qdc; entries", t), 200, 0, 200);
    m_plots[Form("%s_station", t)] = new TH1F(Form("%s_station", t), Form("%s_station; station ; entries", t), 6, -0.5, 5.5);
    m_plots[Form("%s_tofpet", t)] = new TH1F(Form("%s_tofpet", t), Form("%s_tofpet; tofpet number; entries", t), 50, 0, 50);
    

    //2D plots
    m_plots[Form("%s_channelPerTofpet", t)] = new TH2F(Form("%s_channelPerTofpet", t), Form("%s_channelPerTofpet; #tofpet, #channel", t), 50, 0, 50, 100, 0, 100);
  }
}

int showerStartWall( std::string dataType, std::array<int, 2*TB::SCIFISTATION> &hits, int hitThr) {
  //GetTofpetID
  if (dataType == "TB"){
    for (int i = 0; i < TB::SCIFISTATION; ++i){
      //std::cout << hits[i] << "   ||    " << hits[i+4] << std::endl; 
      if (hits[i] > hitThr && hits[i+4] > hitThr){ 
        return i;} // the Fe wall number is the same as the ScifiStation in front of it
    }
  }
  return -1;
}

void runAnalysis() //(int runN, int partN)
{
  auto start = std::chrono::system_clock::now();
  auto now = std::chrono::system_clock::to_time_t(start);
  std::cout << "Start: " << std::ctime(&now)  << "\n" <<std::flush;

  // ##################### Read file #####################
  int runNumber{100633};  
  // 100633: pion 140 GeV 3 walls file
  // 100635: pion 180 GeV 3 walls file
  // 100637: pion 240 GeV 3 walls file
  // 100639: pion 300 GeV 3 walls file

  auto *fEventTree = new TChain("rawConv");
  for (int i = 0; i<3; ++i){
    fEventTree->Add(Form("root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/run_%d/sndsw_raw-%04d.root", runNumber, i)); 
  }

  TFile outputFile(Form("outputRun_%d.root", runNumber), "RECREATE"); 
  outputFile.cd();

  std::map<std::string, double> counters;
  std::map<std::string, TH1*> plots;
  std::vector<std::string> tags;
  tags.push_back("Scifi");
  tags.push_back("MuFilter");


  definePlots(plots, counters, tags);
  
  // ##################### Read hits from Scifi and Mufilter  #####################

  auto mu_hits = new TClonesArray("MuFilterHit");
  fEventTree->SetBranchAddress("Digi_MuFilterHits", &mu_hits);
  auto sf_hits = new TClonesArray("sndScifiHit");
  fEventTree->SetBranchAddress("Digi_ScifiHits", &sf_hits);
  auto header = new SNDLHCEventHeader();
  fEventTree->SetBranchAddress("EventHeader.", &header);
  
  // Loop over events
  int iMax = fEventTree->GetEntries();
  for ( int m =0; m < iMax; m++ ){ 
    if (m % 100 == 0) std::cout << "Processing event: " << m << '\r' << std::flush;
    fEventTree->GetEntry(m);
    int sf_max=sf_hits->GetEntries();
    int mu_max=mu_hits->GetEntries();

    //remove almost empty events
    if (sf_max < SCIFIMINHITS || mu_max < MUMINHITS) continue;

    std::array<int, 2*TB::SCIFISTATION> hitsPerStation = {0};     //1x, 2x, 3x, 4x, 1y, 2y, 3y, 4y

    plots["histo"]->Fill(sf_max, mu_max);

    for (int i=0 ; i<sf_max; i++) {

        auto t = tags[0].c_str();
        auto sf_hit = (sndScifiHit*) sf_hits->At(i);

        //plot some basic info
        int station = sf_hit->GetStation();
        plots[Form("%s_station", t)]->Fill(station);
        

        if (sf_hit->isVertical()) {
          int bin = 4+(station-1);
          hitsPerStation[bin] +=1 ;
          plots[Form("HitsperStation_%d", bin)]->Fill(hitsPerStation[bin]);
        }
        else {
          int bin = station-1;
          hitsPerStation[bin] +=1;
          plots[Form("HitsperStation_%d", bin)]->Fill(hitsPerStation[bin]);
        }

    }
    //check where the shower starts
    int showerStart= showerStartWall("TB", hitsPerStation, TB::SCIFITHRESHOLD);
    plots["ShowerStart"]->Fill(showerStart);
    plots["ShowerStartProbability"]->Fill(showerStart);

    for (int i=0 ; i<mu_max; i++) {
        auto t = tags[1].c_str();
        auto mu_hit = (MuFilterHit*) mu_hits->At(i);

/*        REDO AFTER SIMONA COMMENT
          for (int j = 0; j<NsidesNch; ++j){
          int channel = mu_hit->Getchannel(j);
          if (channel > -1) plots[Form("%s_channels", t)]->Fill(channel);
          float charge = mu_hit->GetSignal(j);
          if (charge > -1) plots[Form("%s_charge", t)]->Fill(charge);
        }*/

        //plot some basic info
        int station = mu_hit->GetPlane();
        plots[Form("%s_station", t)]->Fill(station);
        float time = mu_hit->GetTime()*TDC2ns;
        plots[Form("%s_times", t)] ->Fill(time);
    }
  }

  auto stop = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = stop-start;
  auto end = std::chrono::system_clock::to_time_t(stop);
  std::cout << "\nDone: " << std::ctime(&end)  << std::endl;

  // Write ratio of showering event in different stations
  std::cout << "Shower start ratio for station 2 (" << plots["ShowerStart"]->GetBinContent(5) << ") and station 1 (" << plots["ShowerStart"]->GetBinContent(4) 
            << ") is " << plots["ShowerStart"]->GetBinContent(5)/plots["ShowerStart"]->GetBinContent(4) << std::endl;

  std::cout << "Shower start ratio for station 3 (" << plots["ShowerStart"]->GetBinContent(6) << ") and station 2 (" << plots["ShowerStart"]->GetBinContent(5) 
            << ") is " << plots["ShowerStart"]->GetBinContent(6)/plots["ShowerStart"]->GetBinContent(5) << std::endl;
  
  // ##################### Write results to file #####################
  outputFile.Write();
  outputFile.Close();
}