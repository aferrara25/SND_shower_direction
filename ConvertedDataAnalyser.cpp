#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "Inclusion.h"
#include "SciFiPlaneView.h"

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
  int nChannels = configuration.BOARDPERSTATION*TOFPETperBOARD*TOFPETCHANNELS; //512 for Test Beam
  
  for (auto tag : tags) {

    const auto t{tag.c_str()};
    //plot per event
    m_plots[Form("%s_ShowerStart", t)] = new TH1D(Form("%s_ShowerStart", t), Form("%s_ShowerStart; station; entries", t), 5, 0.5, 5.5);
    m_plots[Form("%s_Times", t)] = new TH1D(Form("%s_Times", t), Form("%s_Times; time (ns) ; entries", t), 150, -5, 145);
    m_plots[Form("%s_Station", t)] = new TH1D(Form("%s_Station", t), Form("%s_Station; station ; entries", t), 6, -0.5, 5.5);
    
    //plot per station
    for (int st = 1; st < configuration.SCIFISTATION+1; ++st){
      m_plots[Form("%s_HitsperStation_st%dX", t, st)] = new TH1D (Form("%s_HitsperStation_%dX", t, st), Form("%s_HitsperStation_%dX;n hit in event;entries", t, st), nChannels, 0, nChannels);
      m_plots[Form("%s_HitsperStation_st%dY", t, st)] = new TH1D (Form("%s_HitsperStation_%dY", t, st), Form("%s_HitsperStation_%dY;n hit in event;entries", t, st), nChannels, 0, nChannels);
      m_plots[Form("%s_Position_st%dX", t, st)] = new TH1D(Form("%s_Position_st%dX", t, st), Form("%s_Position_st%dX; x (cm); entries", t, st), 135, -.5, 13);
      m_plots[Form("%s_Position_st%dY", t, st)] = new TH1D(Form("%s_Position_st%dY", t, st), Form("%s_Position_st%dY; y (cm); entries", t, st), 135, -.5, 13);
      m_plots[Form("%s_Signals_st%dX", t, st)] = new TH1D(Form("%s_Signals_st%dX", t, st), Form("%s_Signals_st%dX; qdc (a.u.) ; entries", t, st), 100, -30, 80);
      m_plots[Form("%s_Signals_st%dY", t, st)] = new TH1D(Form("%s_Signals_st%dY", t, st), Form("%s_Signals_st%dY; qdc (a.u.) ; entries", t, st), 100, -30, 80);
      m_plots[Form("%s_Tofpet_st%dX", t, st)] = new TH1D(Form("%s_Tofpet_st%dX", t, st), Form("%s_Tofpet_st%dX; tofpet number; entries", t, st), 10, 0, 10);
      m_plots[Form("%s_Tofpet_st%dY", t, st)] = new TH1D(Form("%s_Tofpet_st%dY", t, st), Form("%s_Tofpet_st%dY; tofpet number; entries", t, st), 10, 0, 10);
      m_plots[Form("%s_Centroid_Position_st%d", t, st)] = new TH2D(Form("%s_Centroid_Position_st%d", t, st), Form("%s_Centroid_Position_st%d; x (cm); y (cm)", t, st), nChannels+1, -0.5*.025, (nChannels+0.5)*.025, nChannels+1, -0.5*.025, (nChannels+0.5)*.025);
      m_plots[Form("%s_HitDistribution_st%d", t, st)] = new TH2D (Form("%s_HitDistribution_st%d", t, st), Form("%s_HitDistribution_st%d; n hit %dX; n hit %dY", t,  st, st, st), nChannels, 0, nChannels, nChannels, 0, nChannels);
    }
  }
}


void fillPlots (std::vector<SciFiPlaneView> detector, std::map<std::string, TH1*> &plots, std::string &t) {
  for (auto plane : detector){
    auto centroid = plane.getCentroid();
    plots[Form("%s_Centroid_Position_st%d", t.c_str(), plane.getStation())]->Fill(centroid[0], centroid[1]);
    
    int nhitsX = std::count_if(plane.qdc.x.begin(), plane.qdc.x.end(), [] (double t) {return t > DEFAULT;});
    int nhitsY = std::count_if(plane.qdc.y.begin(), plane.qdc.y.end(), [] (double t) {return t > DEFAULT;});

    plots[Form("%s_Station", t.c_str())]->Fill(plane.getStation(), nhitsX + nhitsY);
    plots[Form("%s_HitsperStation_st%dX", t.c_str(), plane.getStation())]->Fill(nhitsX);
    plots[Form("%s_HitsperStation_st%dY", t.c_str(), plane.getStation())]->Fill(nhitsY);
    plots[Form("%s_HitDistribution_st%d", t.c_str(), plane.getStation())]->Fill(nhitsX, nhitsY);
    for (int i{0}; i<plane.getConfig().BOARDPERSTATION*TOFPETperBOARD*TOFPETCHANNELS; ++i) {
      if (plane.qdc.x[i] != DEFAULT) {
        plots[Form("%s_Signals_st%dX", t.c_str(), plane.getStation())]->Fill(plane.qdc.x[i]);
        plots[Form("%s_Times", t.c_str())]->Fill(plane.hitTimestamps.x[i]);
        plots[Form("%s_Position_st%dX", t.c_str(), plane.getStation())]->Fill(i*0.025);
        plots[Form("%s_Tofpet_st%dX", t.c_str(), plane.getStation())]->Fill(static_cast<int>(i/64));
      }
      if (plane.qdc.y[i] != DEFAULT) {
        plots[Form("%s_Signals_st%dY", t.c_str(), plane.getStation())]->Fill(plane.qdc.y[i]);
        plots[Form("%s_Times", t.c_str())]->Fill(plane.hitTimestamps.y[i]);
        plots[Form("%s_Position_st%dY", t.c_str(), plane.getStation())]->Fill(i*0.025);
        plots[Form("%s_Tofpet_st%dY", t.c_str(), plane.getStation())]->Fill(static_cast<int>(i/64));
      }
    }
  }
}

std::vector<SciFiPlaneView> fillSciFi(cfg configuration, TClonesArray *sf_hits){

  std::vector<SciFiPlaneView> scifi_planes;

  int begin{0};
  int count{0};

  int n_sf_hits{sf_hits->GetEntries()};

  for (int st{1}; st <= configuration.SCIFISTATION; ++st) {
    begin = count;
    while (count < n_sf_hits &&
           st == static_cast<sndScifiHit *>(sf_hits->At(count))->GetStation()) {
      ++count;
    }

    auto plane = SciFiPlaneView(configuration, sf_hits, begin, count, st);
    plane.fillQDC();
    plane.fillTimestamps();
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

// timecut -> vector scifiplaneview time cut 
// nel file skimmato ho comunque sempre solo un evento in stazione 1-> leggo tempo di quello (parametro è vector scifiplaneview con tutti hit) e butto via i fuori tempo

bool hitCut (std::vector<SciFiPlaneView> detector){
  for (auto plane : detector){
    if (plane.getStation() == 1 && (plane.sizes().x != 1 || plane.sizes().y !=1 ) ) return false;
    else if (plane.getStation() > 1){
      int thr = plane.getConfig().SCIFITHRESHOLD;
      if (plane.sizes().x > thr && plane.sizes().y > thr) return true;
    }
  }
  // if I finish the for loop without finding a plane above threshold return false
  return false;
}

void timeCut (std::vector<SciFiPlaneView> detector) {
  double referenceTime{-1};
  for (auto plane : detector) {
    cfg config = plane.getConfig();
    int station = plane.getStation();
    std::array<double, 512> timeArrayX = plane.hitTimestamps.x;
    std::array<double, 512> timeArrayY = plane.hitTimestamps.y;
    if (station == 1) {
      if ( std::count_if(timeArrayX.begin(), timeArrayX.end(), [] (double t) {return t > DEFAULT;} ) == 1 &&
           std::count_if(timeArrayY.begin(), timeArrayY.end(), [] (double t) {return t > DEFAULT;} ) == 1) {
            double timeX = *std::max_element(timeArrayX.begin(), timeArrayX.end());
            double timeY = *std::max_element(timeArrayY.begin(), timeArrayY.end());
            if (std::abs(timeX - timeY) < 1) referenceTime = timeX;
           }
    }
    else {
      for (int i{0}; i < 512; ++i) {
        if ( (timeArrayX[i] - referenceTime) > config.TIMECUT) plane.resetHit(false, i); //lo cancello 
        if ( (timeArrayY[i] - referenceTime) > config.TIMECUT) plane.resetHit(true, i); //lo cancello 

      }
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
  // 100631: pion 100 GeV 3 walls file
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
  tags.push_back("NoCut");
  tags.push_back("Cut");
  tags.push_back("Cluster");


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
    //if (m >10000) break;
    fEventTree->GetEntry(m);


    int sf_max=sf_hits->GetEntries();
    int mu_max=mu_hits->GetEntries();

    if (sf_max < 15 && mu_max < 3 ) continue;
    
    auto scifi_planes = fillSciFi(configuration, sf_hits);

    //Before cut
    int showerStart = checkShower(scifi_planes);
    plots[Form("%s_ShowerStart", tags[0].c_str())]->Fill(showerStart);
    for (auto plane : scifi_planes)  plane.findCentroid(6);
    fillPlots(scifi_planes, plots, tags[0]);

    //After cut
    if ( !hitCut(scifi_planes) ) continue;
    timeCut(scifi_planes);
    showerStart = checkShower(scifi_planes);

    plots[Form("%s_ShowerStart", tags[1].c_str())]->Fill(showerStart);
    for (auto plane : scifi_planes) plane.findCentroid(6);
    fillPlots(scifi_planes, plots, tags[1]);
    
    //After cluster
    for (auto plane : scifi_planes){
      plane.findCluster();
      plane.findCentroid(6);
    }
    showerStart = checkShower(scifi_planes);
    plots[Form("%s_ShowerStart", tags[2].c_str())]->Fill(showerStart);
    fillPlots(scifi_planes, plots, tags[2]);
    
  }
  
  auto stop = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = stop-start;
  auto end = std::chrono::system_clock::to_time_t(stop);
  std::cout << "\nDone: " << std::ctime(&end)  << std::endl;

  // ##################### Write results to file #####################
  outputFile.Write();
  outputFile.Close();
}