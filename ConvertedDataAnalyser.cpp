#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "Inclusion.h"
#include "Scifihits.h"

const double FREQ{160.316E6};
const double TDC2ns = 1E9/FREQ;
const int NSIPM{8};
const int NSIDE{2};
const int NsidesNch{16};
const int TOFPETperBOARD{8};
const int TOFPETCHANNELS{64};


struct cfg
{
  int SCIFISTATION{-1};
  int MUSTATION{-1};
  int NWALLS{-1};
  int SCIFITHRESHOLD{-1};
  int SCIFIMAXGAP{-1};
  int SCIFISIDECUT{-1};
  int SCIFIMINHITS{999};
  int MUMINHITS{999};
  int BOARDPERSTATION{-1};

  //geometry parameters
  double SCIFIDIM{-1};
  
  const char *INFILENAME;
  const char *OUTFILENAME;
};

cfg setCfg( bool istb ) {
  cfg config;
  if (istb) {
    config.SCIFISTATION = 4;
    config.MUSTATION = 5;
    config.NWALLS = 3;
    config.SCIFITHRESHOLD = 56;
    config.SCIFIMAXGAP = 999;
    config.SCIFISIDECUT = 0;
    config.SCIFIMINHITS = 10;
    config.MUMINHITS = 2;
    config.BOARDPERSTATION = 1;
    config.SCIFIDIM = 13;
    config.INFILENAME = "root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/";
    config.OUTFILENAME = "TB_output";
  }
  else {
    config.SCIFISTATION = 5;
    config.MUSTATION = 8;
    config.NWALLS = 5;
    config.SCIFITHRESHOLD = 56;
    config.SCIFIMINHITS = 2;
    config.MUMINHITS = 2;
    config.BOARDPERSTATION = 3;
    config.INFILENAME = "root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/physics/2023/";
    config.OUTFILENAME = "TI18_output";
  }
  return config;
}


void definePlots( cfg configuration, std::map<std::string, TH1*> &m_plots, std::map<std::string, double> &m_counters, std::vector<std::string> &tags) {
  m_plots["histo"] = new TH2F("histo", "; # scifi hits; # mu hits", 100, 0, 100, 100, 0, 100);

  // events characteristics plots
  m_plots["ShowerStart"] = new TH1D("ShowerStart", "ShowerStart; station; entries", 8, -2, 5);
  // together x and y
  for (int st = 0; st < configuration.SCIFISTATION; ++st){
    m_plots[Form("HitDistribution_st%d", st)] = new TH2D (Form("HitDistribution_st%d", st), Form("HitDistribution_st%d; n hit %dX; n hit %dY", st, st, st), 275, 0, 550, 275, 0, 550);
  }
  // separately x and y 
  for (int st = 0; st < 2*configuration.SCIFISTATION; ++st){
    m_plots[Form("HitsperStation_%d", st)] = new TH1D (Form("HitsperStation_%d", st), Form("HitsperStation_%dX;station;entries", st+1), 500, 0, 500);
    if (st>3) {
      m_plots[Form("HitsperStation_%d", st)]->SetTitle(Form("HitsperStation_%dY", st-3));
    }

  }
 
  //basic quantities plots for scifi and mu 
  for (auto tag : tags) {
    const auto t{tag.c_str()};
    m_plots[Form("%s_times", t)] = new TH1D(Form("%s_times", t), Form("%s_times; time (ns) ; entries", t), 150, -5, 145);
    m_plots[Form("%s_signals", t)] = new TH1D(Form("%s_signals", t), Form("%s_signals; qdc? ; entries", t), 100, -30, 80);
    m_plots[Form("%s_channels", t)] = new TH1D(Form("%s_channels", t), Form("%s_channels; n channel; entries", t), 100, 0, 100);
  m_plots[Form("%s_charge", t)] = new TH1D(Form("%s_charge", t), Form("%s_charge; qdc?; entries", t), 200, 0, 200);
    m_plots[Form("%s_station", t)] = new TH1D(Form("%s_station", t), Form("%s_station; station ; entries", t), 6, -0.5, 5.5);
    m_plots[Form("%s_tofpet", t)] = new TH1I(Form("%s_tofpet", t), Form("%s_tofpet; tofpet number; entries", t), 10, 0, 10);
    m_plots[Form("%s_boardID", t)] = new TH1I(Form("%s_boardID", t), Form("%s_boardID; board number; entries", t), 100, 0, 100);

    //for stations
    for (int st = 0; st < configuration.SCIFISTATION; ++st){
      m_plots[Form("%s_Xposition_st%d", t, st)] = new TH1D(Form("%s_Xposition_st%d", t, st), Form("%s_Xposition_st%d; x (cm); entries", t, st), configuration.SCIFIDIM+2, -1, configuration.SCIFIDIM+1);
      m_plots[Form("%s_Yposition_st%d", t, st)] = new TH1D(Form("%s_Yposition_st%d", t, st), Form("%s_Yposition_st%d; x (cm); entries", t, st), configuration.SCIFIDIM+2, -1, configuration.SCIFIDIM+1);
    }
  }
}

void runAnalysis(int runNumber, int nFiles, bool isTB) //(int runN, int partN)
{

  auto start = std::chrono::system_clock::now();
  auto now = std::chrono::system_clock::to_time_t(start);
  std::cout << "Start: " << std::ctime(&now)  << "\n" <<std::flush;

  // ##################### Set right parameters for data type (TB/TI18) #####################
  cfg configuration = setCfg(isTB);
  // ##################### Read file #####################
  //int runNumber{100633};  
  // 100633: pion 140 GeV 3 walls file
  // 100635: pion 180 GeV 3 walls file
  // 100637: pion 240 GeV 3 walls file
  // 100639: pion 300 GeV 3 walls file

  auto *fEventTree = new TChain("rawConv");
  //for (int i = 0; i<3; ++i){
  for (int i = 0; i<nFiles; ++i){
    fEventTree->Add(Form("%srun_%06d/sndsw_raw-%04d.root", configuration.INFILENAME, runNumber, i)); 
  }

  TFile outputFile(Form("%sRun_%d.root", configuration.OUTFILENAME, runNumber), "RECREATE"); 
  outputFile.cd();

  std::map<std::string, double> counters;
  std::map<std::string, TH1*> plots;
  std::vector<std::string> tags;
  tags.push_back("Scifi");
  //tags.push_back("MuFilter");


  definePlots(configuration, plots, counters, tags);
  
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
    if (sf_max < configuration.SCIFIMINHITS || mu_max < configuration.MUMINHITS) continue;

    std::vector<Scifihits> Hits;
    for (int i = 0; i < configuration.SCIFISTATION; ++i) Hits.emplace_back(Scifihits(configuration.BOARDPERSTATION*TOFPETperBOARD*TOFPETCHANNELS));

    plots["histo"]->Fill(sf_max, mu_max);

    for (int i=0 ; i<sf_max; i++) {

        auto t = tags[0].c_str();
        auto sf_hit = (sndScifiHit*) sf_hits->At(i);

        //plot some basic info
        int station = sf_hit->GetStation();
        plots[Form("%s_station", t)]->Fill(station);
        
        double time = sf_hit->GetTime(0)*TDC2ns;
        plots[Form("%s_times", t)]->Fill(time);

        int tofpet = sf_hit->GetTofpetID(0);
        plots[Form("%s_tofpet", t)]->Fill(tofpet);

        int boardID = sf_hit->GetBoardID(0);
        plots[Form("%s_boardID", t)]->Fill(boardID);

        double signal = sf_hit->GetSignal(0);
        plots[Form("%s_signals", t)]->Fill(signal);

        int channel = sf_hit->Getchannel(0);
        plots[Form("%s_channels", t)]->Fill(channel);

        bool vertical = sf_hit->isVertical();

        double pos = (64*tofpet+63-channel)*0.025;

        if (vertical) plots[Form("%s_Yposition_st%d", t, station-1)]->Fill(pos); 
        else plots[Form("%s_Xposition_st%d", t, station-1)]->Fill(pos);

        Hits[station-1].addHit(tofpet, channel, vertical);

    }

    //check where the shower starts
    std::vector<bool> isShowering;

    for (int st = 0;  st < configuration.SCIFISTATION; ++st) {
      plots[Form("HitDistribution_st%d", st)]->Fill(Hits[st].getXhits(), Hits[st].getYhits());
      plots[Form("HitsperStation_%d", st)]->Fill(Hits[st].getXhits());
      plots[Form("HitsperStation_%d", st+4)]->Fill(Hits[st].getYhits());

      isShowering.emplace_back(Hits[st].checkShower(configuration.SCIFITHRESHOLD, configuration.SCIFIMAXGAP, configuration.SCIFISIDECUT));
    
      if (Hits[st].checkMultipleHits()) std::cout << "Same channel firing more than once in st " << st << std::endl;
    }
    
    for (int i = 0; i < configuration.SCIFISTATION; ++i) {
      if (isShowering[i]) {
        plots["ShowerStart"]->Fill(i);
        break;
      }
      if (i == (configuration.SCIFISTATION-1)) plots["ShowerStart"]->Fill(-1);
    }


/*    for (int i=0 ; i<mu_max; i++) {
        auto t = tags[1].c_str();
        auto mu_hit = (MuFilterHit*) mu_hits->At(i);
        //plot some basic info
        int station = mu_hit->GetPlane();
        plots[Form("%s_station", t)]->Fill(station);
        float time = mu_hit->GetTime()*TDC2ns;
        plots[Form("%s_times", t)] ->Fill(time);
    }*/
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