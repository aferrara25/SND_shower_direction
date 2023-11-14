#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "Inclusion.h"
#include "SciFiPlaneView.h"

const double FREQ{160.316E6};
const double TDC2ns = 1E9/FREQ;
const int NSIPM{8};
const int NSIDE{2};
const int NsidesNch{16};
const int TOFPETperBOARD{8};
const int TOFPETCHANNELS{64};

cfg setCfg( bool istb ) {
  cfg config;
  if (istb) {
    config.SCIFISTATION = 4;
    config.MUSTATION = 5;
    config.NWALLS = 3;
    config.SCIFITHRESHOLD = 35;
    config.SCIFIMAXGAP = 5;
    config.SCIFISIDECUT = 90;
    config.SCIFIMINHITS = 10;
    config.MUMINHITS = 5;
    config.BOARDPERSTATION = 1;
    config.TIMECUT = 1;
    config.SCIFIDIM = 13;
    config.INFILENAME = "root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/";
    config.OUTFILENAME = "output/TB_output";
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
    config.OUTFILENAME = "output/TI18_output";
  }
  return config;
}



void definePlots( cfg configuration, std::map<std::string, TH1*> &m_plots, std::map<std::string, double> &m_counters, std::vector<std::string> &tags) {
  m_plots["histo"] = new TH2F("histo", "; # scifi hits; # mu hits", 100, 0, 100, 100, 0, 100);

  // events characteristics plots
  m_plots["ShowerStart"] = new TH1D("ShowerStart", "ShowerStart; station; entries", 7, -2, 5);
  m_plots["TimeCut_ShowerStart"] = new TH1D("TimeCut_ShowerStart", "TimeCut_ShowerStart; station; entries", 7, -2, 5);
  // together x and y
  // SPOSTA IN %s
  for (int st = 0; st < configuration.SCIFISTATION; ++st){
    m_plots[Form("HitDistribution_st%d", st)] = new TH2D (Form("HitDistribution_st%d", st), Form("HitDistribution_st%d; n hit %dX; n hit %dY", st+1, st+1, st+1), 275, 0, 550, 275, 0, 550);
    m_plots[Form("TimeCut_HitDistribution_st%d", st)] = new TH2D (Form("TimeCut_HitDistribution_st%d", st), Form("TimeCut_HitDistribution_st%d; n hit %dX; n hit %dY", st+1, st+1, st+1), 275, 0, 550, 275, 0, 550);
  }
  // separately x and y 
  for (int st = 0; st < 2*configuration.SCIFISTATION; ++st){
    m_plots[Form("HitsperStation_%d", st)] = new TH1D (Form("HitsperStation_%d", st), Form("HitsperStation_%dX;station;entries", st+1), 500, 0, 500);
    m_plots[Form("TimeCut_HitsperStation_%d", st)] = new TH1D (Form("TimeCut_HitsperStation_%d", st), Form("TimeCut_HitsperStation_%dX;station;entries", st+1), 500, 0, 500);
    if (st>3) {
      m_plots[Form("HitsperStation_%d", st)]->SetTitle(Form("HitsperStation_%dY", st-3));
      m_plots[Form("TimeCut_HitsperStation_%d", st)]->SetTitle(Form("TimeCut_HitsperStation_%dY", st-3));
    }

  }
 
  //basic quantities plots for scifi and mu 
  for (auto tag : tags) {
    
    //SCIFI plots
    if (tag == "Scifi"){
      const auto t{tag.c_str()};
      m_plots[Form("%s_times", t)] = new TH1D(Form("%s_times", t), Form("%s_times; time (ns) ; entries", t), 150, -5, 145);
      m_plots[Form("%s_signals", t)] = new TH1D(Form("%s_signals", t), Form("%s_signals; qdc? ; entries", t), 100, -30, 80);
      m_plots[Form("%s_channels", t)] = new TH1D(Form("%s_channels", t), Form("%s_channels; n channel; entries", t), 100, 0, 100);
      m_plots[Form("%s_charge", t)] = new TH1D(Form("%s_charge", t), Form("%s_charge; qdc?; entries", t), 200, 0, 200);
      m_plots[Form("%s_station", t)] = new TH1D(Form("%s_station", t), Form("%s_station; station ; entries", t), 6, -0.5, 5.5);
      m_plots[Form("%s_tofpet", t)] = new TH1D(Form("%s_tofpet", t), Form("%s_tofpet; tofpet number; entries", t), 10, 0, 10);
      m_plots[Form("%s_boardID", t)] = new TH1D(Form("%s_boardID", t), Form("%s_boardID; board number; entries", t), 100, 0, 100);
      
      for (int st = 0; st < configuration.SCIFISTATION; ++st){
        m_plots[Form("%s_Xposition_st%d", t, st)] = new TH1D(Form("%s_Xposition_st%d", t, st), Form("%s_Xposition_st%d; x (cm); entries", t, st), configuration.SCIFIDIM+2, -1, configuration.SCIFIDIM+1);
        m_plots[Form("%s_Yposition_st%d", t, st)] = new TH1D(Form("%s_Yposition_st%d", t, st), Form("%s_Yposition_st%d; x (cm); entries", t, st), configuration.SCIFIDIM+2, -1, configuration.SCIFIDIM+1);
        m_plots[Form("%s_TimeCut_signals_st%d", t, st)] = new TH1D(Form("%s_TimeCut_signals_st%d", t, st), Form("%s_TimeCut_signals_st%d; qdc? ; entries", t, st+1), 100, -30, 80);
        m_plots[Form("%s_signals_st%d", t, st)] = new TH1D(Form("%s_signals_st%d", t, st), Form("%s_signals_st%d; qdc? ; entries", t, st+1), 100, -30, 80);
        m_plots[Form("%s_tofpet_st%d", t, st)] = new TH1D(Form("%s_tofpet_st%d", t, st), Form("%s_tofpet_st%d; tofpet number; entries", t, st+1), 10, 0, 10);
        m_plots[Form("%s_TimeCut_Xposition_st%d", t, st)] = new TH1D(Form("%s_TimeCut_Xposition_st%d", t, st), Form("%s_TimeCut_Xposition_st%d; n channel; entries", t, st+1), configuration.SCIFIDIM+2, -1, configuration.SCIFIDIM+1);
        m_plots[Form("%s_TimeCut_Yposition_st%d", t, st)] = new TH1D(Form("%s_TimeCut_Yposition_st%d", t, st), Form("%s_TimeCut_Yposition_st%d; n channel; entries", t, st+1), configuration.SCIFIDIM+2, -1, configuration.SCIFIDIM+1);
        m_plots[Form("%s_Centroid_Position_st%d", t, st)] = new TH2D(Form("%s_Centroid_Position_st%d", t, st), Form("%s_Centroid_Position_st%d; channel x; channel y", t, st+1), 513, -0.5, 512.5, 513, -0.5, 512.5);
      }

      for (int st = 0; st < 2*configuration.SCIFISTATION; ++st){
        m_plots[Form("%s_channels_st%d", t, st)] = new TH1D(Form("%s_channels_st%d", t, st), Form("%s_channels_st%dX; n channel; entries", t ,st+1), 514, -0.5, 512.5);
        m_plots[Form("%s_TimeCut_channels_st%d", t, st)] = new TH1D(Form("%s_TimeCut_channels_st%d", t, st), Form("%s_TimeCut_channels_st%dX; n channel; entries", t ,st+1), 513, -0.5, 512.5);   
        m_plots[Form("%s_Shower_TimeCut_channels_st%d", t, st)] = new TH1D(Form("%s_Shower_TimeCut_channels_st%d", t, st), Form("%s_Shower_TimeCut_channels_st%dX; n channel; entries", t ,st+1), 513
        , -0.5, 512.5);          
        if (st>3) {
          m_plots[Form("%s_channels_st%d", t, st)]->SetTitle(Form("%s_channels_st%dY", t, st-3));
          m_plots[Form("%s_TimeCut_channels_st%d", t, st)]->SetTitle(Form("%s_TimeCut_channels_st%dY", t, st-3));    
          m_plots[Form("%s_Shower_TimeCut_channels_st%d", t, st)]->SetTitle(Form("%s_Shower_TimeCut_channels_st%dY", t, st-3));
        }  
      }
    }

    //Muon plots
    else if ( tag == "MuFilter") {

      const auto t{tag.c_str()};
      m_plots[Form("%s_times", t)] = new TH1D(Form("%s_times", t), Form("%s_times; time (ns) ; entries", t), 150, -5, 145);
      m_plots[Form("%s_signals", t)] = new TH1D(Form("%s_signals", t), Form("%s_signals; qdc? ; entries", t), 100, -30, 80);
      m_plots[Form("%s_channels", t)] = new TH1D(Form("%s_channels", t), Form("%s_channels; n channel; entries", t), 100, 0, 100);
      m_plots[Form("%s_charge", t)] = new TH1D(Form("%s_charge", t), Form("%s_charge; qdc?; entries", t), 200, 0, 200);
      m_plots[Form("%s_station", t)] = new TH1D(Form("%s_station", t), Form("%s_station; station ; entries", t), 6, -0.5, 5.5);
      m_plots[Form("%s_tofpet", t)] = new TH1D(Form("%s_tofpet", t), Form("%s_tofpet; tofpet number; entries", t), 10, 0, 10);
      m_plots[Form("%s_boardID", t)] = new TH1D(Form("%s_boardID", t), Form("%s_boardID; board number; entries", t), 100, 0, 100);
      
      for (int st = 0; st < configuration.MUSTATION; ++st){

      m_plots[Form("%s_signals_st%d", t, st)] = new TH1D(Form("%s_signals_st%d", t, st), Form("%s_signals_st%d; qdc? ; entries", t, st+1), 100, -30, 1000);

      }
    }
  }
}

std::vector<SciFiPlaneView> fillSciFi(cfg configuration, TClonesArray *sf_hits){

  std::vector<SciFiPlaneView> scifi_planes;

  int begin{0};
  int count{0};

  int n_sf_hits{sf_hits->GetEntries()};
  //std::cout<<n_sf_hits<<std::endl;

  for (int st{1}; st <= configuration.SCIFISTATION; ++st) {
    begin = count;
    while (count < n_sf_hits &&
           st == static_cast<sndScifiHit *>(sf_hits->At(count))->GetStation()) {
      ++count;
    }

    auto plane = SciFiPlaneView(configuration, sf_hits, begin, count, st);
    plane.fillQDC();
    /*if (sf_hits->GetEntries()>0) {
        for (int j{0}; j<512; j++) {
            //std::cout<<plane.qdc.x[j]<<"\t";
        }
        std::cout<<std::endl;
    }*/
    scifi_planes.emplace_back(plane);
  }

  return scifi_planes;

}

int checkShower(std::vector<SciFiPlaneView> scifi_planes ) {
  //find start of shower
  for (auto &plane : scifi_planes) {
    if (plane.sizes().x > plane.getConfig().SCIFITHRESHOLD && plane.sizes().y > plane.getConfig().SCIFITHRESHOLD) return plane.getStation(); 
  }
  return -1;
}

std::array<int, 2> findCentroid(SciFiPlaneView plane, int windowSize) {
  std::array<int,2> centroid = {-1};
  double maxSignal{-1};
  auto config = plane.getConfig();
  for (int index{0}; index <2; ++index) {
    for (int i{0}; i < config.BOARDPERSTATION*TOFPETperBOARD*TOFPETCHANNELS -windowSize; ++i) {
  
      double signalSum{0};

      for (int j{i}; j < (i+windowSize); ++j) {
        int den = windowSize;
        double signal;
        if (index == 0) signal = plane.qdc.x[j];
        else signal = plane.qdc.y[j];
        
        if ( signal != -999 ){
          signalSum += signal;
        } else {
          den -=1;
        }

        if (j == i+windowSize-1) {
          double ratio =  signalSum/den;
          if ( maxSignal < ratio ) {
            maxSignal = ratio;
            centroid[index] = i;
          }
        }
      }
    }
  }

  return centroid;
}

// timecut -> vector scifiplaneview time cut 
// nel file skimmato ho comunque sempre solo un evento in stazione 1-> leggo tempo di quello (parametro è vector scifiplaneview con tutti hit) e butto via i fuori tempo

// 

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
  tags.push_back("MuFilter");


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
  for ( int m = 0; m < iMax; m++ ){ 
    if (m % 100 == 0) std::cout << "Processing event: " << m << '\r' << std::flush;
    //if (m >1000) break;
    fEventTree->GetEntry(m);


    int sf_max=sf_hits->GetEntries();
    int mu_max=mu_hits->GetEntries();

    if (sf_max < 15 && mu_max < 3 ) continue;
    
    auto scifi_planes = fillSciFi(configuration, sf_hits);
    int showerStart = checkShower(scifi_planes);
    plots["ShowerStart"]->Fill(showerStart);
    for (auto plane : scifi_planes){
      auto centroid = findCentroid(plane, 6);
      plots[Form("%s_Centroid_Position_st%d", tags[0].c_str(), plane.getStation()-1)]->Fill(centroid[0], centroid[1]);
      //std::cout << "In station: " << plane.getStation() << "  " << centroid[0] << "   ||    "<< centroid[1] << std::endl;
    }

/*
    //################ SCIFI HITS ##################

    for (int i=0 ; i<sf_max; i++) {

        auto t = tags[0].c_str();
        auto sf_hit = (sndScifiHit*) sf_hits->At(i);

        //plot some basic info
        bool vertical = sf_hit->isVertical();

        int station = sf_hit->GetStation();

        double time = sf_hit->GetTime(0);
       // double time = clocktime;//*TDC2ns;

        int tofpet = sf_hit->GetTofpetID(0);

        int boardID = sf_hit->GetBoardID(0);

        double signal = sf_hit->GetSignal(0);
        
        int channel = sf_hit->Getchannel(0);
        
        double pos = (64*tofpet+63-channel)*0.025;
   
    }

    //################ MU FILTER HITS ##################
    for (int i=0 ; i<mu_max; i++) {
        auto t = tags[1].c_str();
        auto mu_hit = (MuFilterHit*) mu_hits->At(i);

        int boardID = mu_hit->GetBoardID(0);

        int station = mu_hit->GetPlane();

        bool vertical = mu_hit->isVertical();

        for (int sipm = 0; sipm < NsidesNch; ++sipm ){
          double time = mu_hit->GetTime(sipm); //*TDC2ns;
        
          if (vertical && sipm > 0) continue;             //DS vertical has only 1 sipm to read
          int tofpet = mu_hit->GetTofpetID(sipm);
        
          double signal = mu_hit->GetSignal(sipm);
        
          int channel = mu_hit->Getchannel(sipm);
              
        }
      }
  */  
    }

  auto stop = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = stop-start;
  auto end = std::chrono::system_clock::to_time_t(stop);
  std::cout << "\nDone: " << std::ctime(&end)  << std::endl;

  // ##################### Write results to file #####################
  outputFile.Write();
  outputFile.Close();
}