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
    config.SCIFIMAXGAP = 2;
    config.SCIFISIDECUT = 90;
    config.SCIFIMINHITS = 10;
    config.MUMINHITS = 5;
    config.BOARDPERSTATION = 1;
    config.TIMECUT = 1;
    config.MUTIMECUT = 3;
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



void definePlots( cfg configuration, std::map<std::string, TH1*> &plots, std::map<std::string, double> &m_counters, std::vector<std::string> &tags) {
  int nChannels = configuration.BOARDPERSTATION*TOFPETperBOARD*TOFPETCHANNELS; //512 for Test Beam
  
  for (auto tag : tags) {

    const auto t{tag.c_str()};
    //plot per event
    plots[Form("%s_ShowerStart", t)] = new TH1D(Form("%s_ShowerStart", t), Form("%s_ShowerStart; station; entries", t), 5, 0.5, 5.5);
    plots[Form("%s_Times", t)] = new TH1D(Form("%s_Times", t), Form("%s_Times; time (clk cycles) ; entries", t), 60, -5, 25);
    plots[Form("%s_Station", t)] = new TH1D(Form("%s_Station", t), Form("%s_Station; station ; entries", t), 6, -0.5, 5.5);
    plots[Form("%s_QDCUS_vs_QDCScifi", t)] = new TH2D(Form("%s_QDCUS_vs_QDCScifi", t), Form("%s_QDCUS_vs_QDCScifi; US qdc; SciFi qdc;", t), 250, 0, 25000, 250, 0, 9000);

    //plot per station
    for (int st = 1; st < configuration.SCIFISTATION+1; ++st){
      plots[Form("%s_ClusterSize_st%dX", t, st)] = new TH1D (Form("%s_ClusterSize_st%dX", t, st), Form("%s_ClusterSize_st%dX; cluster size; entries", t, st), nChannels, 0, nChannels);
      plots[Form("%s_ClusterSize_st%dY", t, st)] = new TH1D (Form("%s_ClusterSize_st%dY", t, st), Form("%s_ClusterSize_st%dY; cluster size; entries", t, st), nChannels, 0, nChannels);
      plots[Form("%s_HitsperStation_st%dX", t, st)] = new TH1D (Form("%s_HitsperStation_%dX", t, st), Form("%s_HitsperStation_%dX;n hit in event;entries", t, st), nChannels, 0, nChannels);
      plots[Form("%s_HitsperStation_st%dY", t, st)] = new TH1D (Form("%s_HitsperStation_%dY", t, st), Form("%s_HitsperStation_%dY;n hit in event;entries", t, st), nChannels, 0, nChannels);
      plots[Form("%s_Position_st%dX", t, st)] = new TH1D(Form("%s_Position_st%dX", t, st), Form("%s_Position_st%dX; x (cm); entries", t, st), 135, -.5, 13);
      plots[Form("%s_Position_st%dY", t, st)] = new TH1D(Form("%s_Position_st%dY", t, st), Form("%s_Position_st%dY; y (cm); entries", t, st), 135, -.5, 13);
      plots[Form("%s_Signals_st%dX", t, st)] = new TH1D(Form("%s_Signals_st%dX", t, st), Form("%s_Signals_st%dX; qdc (a.u.) ; entries", t, st), 100, -30, 80);
      plots[Form("%s_Signals_st%dY", t, st)] = new TH1D(Form("%s_Signals_st%dY", t, st), Form("%s_Signals_st%dY; qdc (a.u.) ; entries", t, st), 100, -30, 80);
      plots[Form("%s_Tofpet_st%dX", t, st)] = new TH1D(Form("%s_Tofpet_st%dX", t, st), Form("%s_Tofpet_st%dX; tofpet number; entries", t, st), 10, 0, 10);
      plots[Form("%s_Tofpet_st%dY", t, st)] = new TH1D(Form("%s_Tofpet_st%dY", t, st), Form("%s_Tofpet_st%dY; tofpet number; entries", t, st), 10, 0, 10);
      plots[Form("%s_Centroid_Position_st%d", t, st)] = new TH2D(Form("%s_Centroid_Position_st%d", t, st), Form("%s_Centroid_Position_st%d; x (cm); y (cm)", t, st), nChannels+1, -0.5*.025, (nChannels+0.5)*.025, nChannels+1, -0.5*.025, (nChannels+0.5)*.025);
      plots[Form("%s_HitDistribution_st%d", t, st)] = new TH2D (Form("%s_HitDistribution_st%d", t, st), Form("%s_HitDistribution_st%d; n hit %dX; n hit %dY", t,  st, st, st), nChannels, 0, nChannels, nChannels, 0, nChannels);
      plots[Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d", t, st)] = new TH2D(Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d", t, st), Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d; US qdc; SciFi qdc;", t, st), 250, 0, 25000, 250, 0, 9000);
    }
    for (int st = 2; st < configuration.SCIFISTATION+1; ++st){
      plots[Form("%s_Centroid_Residuals_st%dX", t, st)] = new TH1D (Form("%s_Centroid_Residuals_st%dX", t, st), Form("%s_Centroid_Residuals_st%dX; x-x_ref (cm);entries", t, st), 26*100, -13, 13);
      plots[Form("%s_Centroid_Residuals_st%dY", t, st)] = new TH1D (Form("%s_Centroid_Residuals_st%dY", t, st), Form("%s_Centroid_Residuals_st%dY; y-y_ref (cm);entries", t, st), 26*100, -13, 13);
    }
    for (int start_st = 1; start_st < configuration.SCIFISTATION+1; ++start_st){
      plots[Form("%s_Shower_SciFi_QDC_shStart%d", t, start_st)] = new TH1D (Form("%s_Shower_SciFi_QDC_shStart%d", t, start_st), Form("%s_Shower_SciFi_QDC_shStart%d; qdc sum (a.u.);entries", t, start_st), 1000, -100, 8000);
      for (int st = start_st; st < configuration.SCIFISTATION+1; ++st){
        plots[Form("%s_Hits_st%d_start%d", t, st, start_st)] = new TH1D (Form("%s_Hits_st%d_start%d", t, st, start_st), Form("%s_Hits_st%d_start%d;n hit in event;entries", t, st, start_st), nChannels, 0, 2*nChannels);
      }
      for (int st = start_st+1; st < configuration.SCIFISTATION+1; ++st){
        plots[Form("%s_Hits_st%d_vs_start%d", t, st, start_st)] = new TH2D (Form("%s_Hits_st%d_vs_start%d", t, st, start_st), Form("%s_Hits_st%d_vs_start%d;n hit in st%d;n hit in st%d", t, st, start_st, st, start_st), nChannels, 0, 2*nChannels, nChannels, 0, 2*nChannels);
      }
    }
    // US planes
    for (int pl{1}; pl < 6; ++pl) {
      plots[Form("%s_US_Timestamps_pl%ds", t, pl)] = new TH1D (Form("%s_US_Timestamps_pl%ds", t, pl), Form("%s_US_Timestamps_pl%ds;timestamps (clk cycles);entries", t, pl), 100, -2, 15);
      plots[Form("%s_US_Timestamps_pl%dl", t, pl)] = new TH1D (Form("%s_US_Timestamps_pl%dl", t, pl), Form("%s_US_Timestamps_pl%dl;timestamps (clk cycles);entries", t, pl), 100, -2, 15);
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
    scifi_planes.emplace_back(plane);
  }

  return scifi_planes;

}

std::vector<USPlaneView> fillUS(cfg configuration, TClonesArray *mufi_hits){

  std::vector<USPlaneView> us_planes;

  int begin{0};
  int count{0};

  int n_mufi_hits{mufi_hits->GetEntries()};
  //skip veto/beam monitor
  while (count < n_mufi_hits &&
        static_cast<MuFilterHit *>(mufi_hits->At(count))->GetSystem() != 2) {
    ++count;
  }
  //plane count starts from 0
  for (int pl{0}; pl < 5; ++pl) {
    begin = count;
    while (count < n_mufi_hits &&
           pl == static_cast<MuFilterHit *>(mufi_hits->At(count))->GetPlane() &&
           static_cast<MuFilterHit *>(mufi_hits->At(count))->GetSystem() == 2) { //stop before DS
      ++count;
    }
    auto plane = USPlaneView(configuration, mufi_hits, begin, count, pl+1);
    us_planes.emplace_back(plane);
  }

  return us_planes;

}

int checkShower(std::vector<SciFiPlaneView> scifi_planes ) {
  //find start of shower
  for (auto &plane : scifi_planes) {
    if (plane.infoCluster()) return plane.getStation(); 
  }
  return -1;
}

// timecut -> vector scifiplaneview time cut 
// nel file skimmato ho comunque sempre solo un evento in stazione 1-> leggo tempo di quello (parametro è vector scifiplaneview con tutti hit) e butto via i fuori tempo

bool hitCut (std::vector<SciFiPlaneView> &detector){
  for (auto &plane : detector){
    if (plane.getStation() == 1 && (plane.sizes().x != 1 || plane.sizes().y !=1 ) ) return false;
    else if (plane.getStation() > 1){
      int thr = plane.getConfig().SCIFITHRESHOLD;
      if (plane.sizes().x > thr && plane.sizes().y > thr) return true;
    }
  }
  // if I finish the for loop without finding a plane above threshold return false
  return false;
}

void timeCut (std::vector<SciFiPlaneView> &Scifi, std::vector<USPlaneView> US ) {
  double referenceTime{-1};

  for (auto &plane : Scifi) {
    int station = plane.getStation();
    auto time = plane.getTime();
    //first look for the reference time
    if (station == 1) {
      if ( plane.sizes().x == 1 && plane.sizes().y == 1 ) {
            double timeX = *std::max_element(time.x.begin(), time.x.end());
            double timeY = *std::max_element(time.y.begin(), time.y.end());
            if (std::abs(timeX - timeY) < 1) referenceTime = timeX;
           }
    }
    else if (station > 1){
      if (referenceTime == -1) continue;
      plane.timeCut(referenceTime);
    }
  }
  for (auto &plane : US) {
    if (referenceTime == -1) continue;
    plane.timeCut(referenceTime);
  }
}

void fillPlots (std::vector<SciFiPlaneView> &Scifi_detector, std::vector<USPlaneView> US, std::map<std::string, TH1*> &plots, std::string &t, int shStart) {
  int showerHits{0};
  double ScifiQDCSum{0}, partialScifiQDCSum{0};
  double USQDCSum{0};
  double Small_USQDCSum{0}, Large_USQDCSum{0};
  auto refCentroid{Scifi_detector[0].getCentroid()};
  int showerStart = checkShower(Scifi_detector);
  
  for (auto plane : Scifi_detector){
    //const int nchannel{plane.getConfig().BOARDPERSTATION*TOFPETperBOARD*TOFPETCHANNELS};
    auto centroid{plane.getCentroid()};
    auto station{plane.getStation()};
    plots[Form("%s_Centroid_Position_st%d", t.c_str(), station)]->Fill(centroid.x, centroid.y);
    if (station > 1) {
      if (refCentroid.x != -1 && centroid.x != -1) {
        plots[Form("%s_Centroid_Residuals_st%dX", t.c_str(), station)]->Fill(centroid.x - refCentroid.x);
      }
      if (refCentroid.y != -1 && centroid.y != -1) {
        plots[Form("%s_Centroid_Residuals_st%dY", t.c_str(), station)]->Fill(centroid.y - refCentroid.y);
      }
    } 
    
    int nhitsX{plane.sizes().x};
    int nhitsY{plane.sizes().y};
    std::array<double, 512> timeX{plane.getTime().x};
    std::array<double, 512> timeY{plane.getTime().y};
    std::array<double, 512> qdcX{plane.getQDC().x};
    std::array<double, 512> qdcY{plane.getQDC().y};
    auto sumScifi = plane.getTotQDC();
    auto clusterSize{plane.getClusterSize()};

    ScifiQDCSum += sumScifi.x;
    ScifiQDCSum += sumScifi.y;
    if (station >= showerStart) {
      partialScifiQDCSum += sumScifi.x;
      partialScifiQDCSum += sumScifi.y;
    } 

    plots[Form("%s_Station", t.c_str())]->Fill(station, nhitsX + nhitsY);
    plots[Form("%s_HitsperStation_st%dX", t.c_str(), station)]->Fill(nhitsX);
    plots[Form("%s_HitsperStation_st%dY", t.c_str(), station)]->Fill(nhitsY);
    plots[Form("%s_HitDistribution_st%d", t.c_str(), station)]->Fill(nhitsX, nhitsY);
    plots[Form("%s_ClusterSize_st%dX", t.c_str(), station)]->Fill(clusterSize.x);
    plots[Form("%s_ClusterSize_st%dY", t.c_str(), station)]->Fill(clusterSize.y);
    
    if (station == shStart) {
      showerHits = nhitsX + nhitsY;
    }
    if (shStart != -1 && station >= shStart) {
      plots[Form("%s_Hits_st%d_start%d", t.c_str(), station, shStart)]->Fill(nhitsX + nhitsY);
    }
    if (shStart != -1 && station > shStart) {
      plots[Form("%s_Hits_st%d_vs_start%d", t.c_str(), station, shStart)]->Fill(nhitsX + nhitsY, showerHits);
    }
    for (int i{0}; i<512; ++i) {
      if (qdcX[i] != DEFAULT) {
        plots[Form("%s_Signals_st%dX", t.c_str(), station)]->Fill(qdcX[i]);
        plots[Form("%s_Times", t.c_str())]->Fill(timeX[i]);
        plots[Form("%s_Position_st%dX", t.c_str(), station)]->Fill(i*0.025);
        plots[Form("%s_Tofpet_st%dX", t.c_str(), station)]->Fill(static_cast<int>(i/64));
        if (station == 1 && nhitsX == 1 && nhitsY == 1){
          refCentroid.x = i*0.025;
        }
      }
      if (qdcY[i] != DEFAULT) {
        plots[Form("%s_Signals_st%dY", t.c_str(), station)]->Fill(qdcY[i]);
        plots[Form("%s_Times", t.c_str())]->Fill(timeY[i]);
        plots[Form("%s_Position_st%dY", t.c_str(), station)]->Fill(i*0.025);
        plots[Form("%s_Tofpet_st%dY", t.c_str(), station)]->Fill(static_cast<int>(i/64));
        if (station == 1 && nhitsX == 1 && nhitsY == 1){
          refCentroid.y = i*0.025;
        }
      }
    }
  }

  for (auto &plane : US){
    auto sumQDC = plane.getTotQDC();
    auto timesUS = plane.getTime();
    int pl = plane.getStation();
    Small_USQDCSum += sumQDC.s;
    Large_USQDCSum += sumQDC.l;
    USQDCSum += (sumQDC.s + sumQDC.l);
    for (int i{0}; i<NCHANNELS; ++i) {
      if (timesUS[i] != DEFAULT) {
        if ((i%16)%8==2 || (i%16)%8==5) {
          plots[Form("%s_US_Timestamps_pl%ds", t.c_str(), pl)]->Fill(timesUS[i]);
        }
        else {
          plots[Form("%s_US_Timestamps_pl%dl", t.c_str(), pl)]->Fill(timesUS[i]);
        }
      }
    }    
  }
  plots[Form("%s_QDCUS_vs_QDCScifi", t.c_str())]->Fill(Large_USQDCSum, ScifiQDCSum);  // only large?
  if (showerStart > 0) {
    plots[Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d", t.c_str(), showerStart)]->Fill(Large_USQDCSum, partialScifiQDCSum); // only large?
    plots[Form("%s_Shower_SciFi_QDC_shStart%d", t.c_str(), showerStart)]->Fill(partialScifiQDCSum);
  }
}

void runAnalysis(int runNumber, int nFiles, bool isTB, bool isMulticore = false) //(int runN, int partN)
{

  auto start = std::chrono::system_clock::now();
  auto now = std::chrono::system_clock::to_time_t(start);
  std::cout << "Start: " << std::ctime(&now)  << "\n" <<std::flush;

  // ##################### Set right parameters for data type (TB/TI18) #####################
  cfg configuration = setCfg(isTB);
  
  // ##################### Read file #####################
  //int runNumber{100633};  
  // 100673: pion 100 GeV 3 walls 16 files 
  // 100633: pion 140 GeV 3 walls 22 files
  // 100671: pion 180 GeV 3 walls 33 files
  // 100648: pion 240 GeV 3 walls 59 files
  // 100639: pion 300 GeV 3 walls 58 files

  auto *fEventTree = new TChain("rawConv");
  TFile* outputFile;
  if (isMulticore){
    fEventTree->Add(Form("%srun_%06d/sndsw_raw-%04d.root", configuration.INFILENAME.c_str(), runNumber, nFiles));
    outputFile = new TFile(Form("%sRun_%d_%d.root", configuration.OUTFILENAME.c_str(), runNumber, nFiles), "RECREATE");
  }
  else {
    for (int i = 0; i<nFiles; ++i){
      fEventTree->Add(Form("%srun_%06d/sndsw_raw-%04d.root", configuration.INFILENAME.c_str(), runNumber, i)); 
    }
    outputFile = new TFile(Form("%sRun_%d.root", configuration.OUTFILENAME.c_str(), runNumber), "RECREATE");
  }
 
  outputFile->cd();

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
    //if (m % 100 == 0) std::cout << "Processing event: " << m << '\r' << std::flush;
    //if (m >10000) break;
    fEventTree->GetEntry(m);


    int sf_max=sf_hits->GetEntries();
    int mu_max=mu_hits->GetEntries();

    if (sf_max < 15) continue;
    
    auto scifi_planes = fillSciFi(configuration, sf_hits);
    auto us_planes = fillUS(configuration, mu_hits);

    //Before cut
    int showerStart = checkShower(scifi_planes);
    plots[Form("%s_ShowerStart", tags[0].c_str())]->Fill(showerStart);
    for (auto &plane : scifi_planes)  plane.findCentroid(6);
    fillPlots(scifi_planes, us_planes, plots, tags[0], showerStart);

    //After cut
    if ( !hitCut(scifi_planes) ) continue;
    timeCut(scifi_planes, us_planes);
    showerStart = checkShower(scifi_planes);

    plots[Form("%s_ShowerStart", tags[1].c_str())]->Fill(showerStart);
    for (auto &plane : scifi_planes) plane.findCentroid(6);
    fillPlots(scifi_planes, us_planes, plots, tags[1], showerStart);
    
    //After cluster
    for (auto &plane : scifi_planes){
      plane.findCluster();
      plane.findCentroid(6);
    }
    showerStart = checkShower(scifi_planes);
    plots[Form("%s_ShowerStart", tags[2].c_str())]->Fill(showerStart);
    fillPlots(scifi_planes, us_planes, plots, tags[2], showerStart);
    
  }
 
  auto stop = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = stop-start;
  auto end = std::chrono::system_clock::to_time_t(stop);
  std::cout << "\nDone: " << std::ctime(&end)  << std::endl;

  // ##################### Write results to file #####################
  outputFile->Write();
  outputFile->Close();
}
