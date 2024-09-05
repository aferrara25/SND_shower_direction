#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "Inclusion.h"
#include "SciFiPlaneView.h"
#include "USPlaneView.h"


cfg setCfg( bool istb, int runN ) {
  cfg config;
  if (istb) {
    config.SCIFI_STATIONS = 4;
    config.SCIFI_BOARDPERPLANE = 1;
    config.SCIFI_NCHANNELS = 512;
    config.SCIFI_TIMECUT = 1;  
    config.SCIFI_DIMCLUSTER = 35;
    config.SCIFI_GAPCLUSTER = 2;
    config.SCIFI_DENSITYWINDOW = 128;
    config.SCIFI_DENSITYHITS = 36;
    config.SCIFI_F = 0.17; 
    // config.SCIFI_maxX = ;
    // config.SCIFI_maxY = ;
    // config.SCIFI_minX = ;
    // config.SCIFI_minY = ;

    config.US_STATIONS = 5;
    config.US_TIMECUT = 3;

    config.SCIFI_DIM = 13;

    config.INFILENAME = "root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/";
    config.OUTFILENAME = "output/TB_output";
  }
  else {
    config.SCIFI_STATIONS = 5;
    config.SCIFI_BOARDPERPLANE = 3;
    config.SCIFI_NCHANNELS = 512*3;
    config.SCIFI_TIMECUT = 1;  
    config.SCIFI_DIMCLUSTER = 18;
    config.SCIFI_GAPCLUSTER = 2;
    config.SCIFI_DENSITYWINDOW = 128;
    config.SCIFI_DENSITYHITS = 18;
    config.SCIFI_F = 0.17;

    config.US_STATIONS = 5;
    config.US_TIMECUT = 3;

    config.SCIFI_DIM = 13*3;
    if (runN < 5422) {
      config.INFILENAME = "root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/physics/2022/";
    }
    else {
      config.INFILENAME = "root://eospublic.cern.ch//eos/experiment/sndlhc/convertedData/physics/2023/";
    }
    config.OUTFILENAME = "output/TI18_output";
  }
  return config;
}


R__LOAD_LIBRARY(/cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/sw/slc9_x86-64/ROOT/v6-28-04-local1/lib/libROOTTPython.so)

void definePlots( cfg configuration, std::map<std::string, TH1*> &plots, std::map<std::string, double> &m_counters, std::vector<std::string> &tags, std::vector<std::string> &shower_tags) {
  int nChannels = configuration.SCIFI_NCHANNELS;
  plots["Delta_timestamp"] = new TH1D("Delta_timestamp", "Delta_timestamp; delta_t (ns); entries", 300, 0, 3000);

  for (auto stag : shower_tags) {
    auto name = Form("%s_%s_vs_%s_%s", tags[1].c_str(), stag.c_str(), tags[2].c_str(), stag.c_str());
    plots[name] = new TH2D(name, Form("%s; station; station;", name), 8, -2.5, 5.5, 8, -2.5, 5.5);
  }

  auto name = Form("%s_%s_vs_%s_%s", tags[2].c_str(), shower_tags[1].c_str(), tags[2].c_str(), shower_tags[2].c_str());
  plots[name] = new TH2D(name, Form("%s; station; station;", name), 8, -2.5, 5.5, 8, -2.5, 5.5);
  
  for (auto tag : tags) {

    const auto t{tag.c_str()};
    //plot per event
    plots[Form("%s_ShowerStart_with_clusters", t)] = new TH1D(Form("%s_ShowerStart_with_clusters", t), Form("%s_ShowerStart_with_clusters; station; entries", t), 7, -1.5, 5.5);
    plots[Form("%s_ShowerStart_with_density", t)] = new TH1D(Form("%s_ShowerStart_with_density", t), Form("%s_ShowerStart_with_density; station; entries", t), 7, -1.5, 5.5);
    plots[Form("%s_ShowerStart_with_F", t)] = new TH1D(Form("%s_ShowerStart_with_F", t), Form("%s_ShowerStart_with_F; station; entries", t), 7, -1.5, 5.5);
    plots[Form("%s_Times", t)] = new TH1D(Form("%s_Times", t), Form("%s_Times; time (clk cycles) ; entries", t), 600, -5, 25);
    plots[Form("%s_Station", t)] = new TH1D(Form("%s_Station", t), Form("%s_Station; station ; entries", t), 6, -0.5, 5.5);
    plots[Form("%s_QDCUS_vs_QDCScifi", t)] = new TH2D(Form("%s_QDCUS_vs_QDCScifi", t), Form("%s_QDCUS_vs_QDCScifi; US qdc; SciFi qdc;", t), 1500, -1000, 20000, 1500, -1000, 8000);
    plots[Form("%s_Shower_development_X", t)] = new TH2D(Form("%s_Shower_development_X", t), Form("%s_Shower_development_X; x (cm); SciFi plane;", t), nChannels+1, -0.5*.025, (nChannels+0.5)*.025, configuration.SCIFI_STATIONS*5, 1, configuration.SCIFI_STATIONS + 1);
    plots[Form("%s_Shower_development_Y", t)] = new TH2D(Form("%s_Shower_development_Y", t), Form("%s_Shower_development_Y; y (cm); SciFi plane;", t), nChannels+1, -0.5*.025, (nChannels+0.5)*.025, configuration.SCIFI_STATIONS*5, 1, configuration.SCIFI_STATIONS + 1);

    //plot per station
    for (int st = 1; st < configuration.SCIFI_STATIONS+1; ++st){
      plots[Form("%s_ClusterSize_st%dX", t, st)] = new TH1D (Form("%s_ClusterSize_st%dX", t, st), Form("%s_ClusterSize_st%dX; cluster size; entries", t, st), nChannels, 0, nChannels);
      plots[Form("%s_ClusterSize_st%dY", t, st)] = new TH1D (Form("%s_ClusterSize_st%dY", t, st), Form("%s_ClusterSize_st%dY; cluster size; entries", t, st), nChannels, 0, nChannels);
      plots[Form("%s_HitsperStation_st%dX", t, st)] = new TH1D (Form("%s_HitsperStation_%dX", t, st), Form("%s_HitsperStation_%dX;n hit in event;entries", t, st), nChannels, 0, nChannels);
      plots[Form("%s_HitsperStation_st%dY", t, st)] = new TH1D (Form("%s_HitsperStation_%dY", t, st), Form("%s_HitsperStation_%dY;n hit in event;entries", t, st), nChannels, 0, nChannels);
      plots[Form("%s_Position_st%dX", t, st)] = new TH1D(Form("%s_Position_st%dX", t, st), Form("%s_Position_st%dX; x (cm); entries", t, st), nChannels+1, -0.5*.025, (nChannels+0.5)*.025);
      plots[Form("%s_Position_st%dY", t, st)] = new TH1D(Form("%s_Position_st%dY", t, st), Form("%s_Position_st%dY; y (cm); entries", t, st), nChannels+1, -0.5*.025, (nChannels+0.5)*.025);
      plots[Form("%s_Signals_st%dX", t, st)] = new TH1D(Form("%s_Signals_st%dX", t, st), Form("%s_Signals_st%dX; qdc (a.u.) ; entries", t, st), 100, -30, 80);
      plots[Form("%s_Signals_st%dY", t, st)] = new TH1D(Form("%s_Signals_st%dY", t, st), Form("%s_Signals_st%dY; qdc (a.u.) ; entries", t, st), 100, -30, 80);
      plots[Form("%s_Tofpet_st%dX", t, st)] = new TH1D(Form("%s_Tofpet_st%dX", t, st), Form("%s_Tofpet_st%dX; tofpet number; entries", t, st), TOFPETperBOARD*configuration.SCIFI_BOARDPERPLANE, 0, TOFPETperBOARD*configuration.SCIFI_BOARDPERPLANE);
      plots[Form("%s_Tofpet_st%dY", t, st)] = new TH1D(Form("%s_Tofpet_st%dY", t, st), Form("%s_Tofpet_st%dY; tofpet number; entries", t, st), TOFPETperBOARD*configuration.SCIFI_BOARDPERPLANE, 0, TOFPETperBOARD*configuration.SCIFI_BOARDPERPLANE);


      plots[Form("%s_Hits_Position_st%dX", t, st)] = new TH1D (Form("%s_Hits_Position_%dX", t, st), Form("%s_Hits_Position_%dX; x (cm); entries", t, st), 1000, -100, 100);
      plots[Form("%s_Hits_Position_st%dY", t, st)] = new TH1D (Form("%s_Hits_Position_%dY", t, st), Form("%s_Hits_Position_%dY; y (cm); entries", t, st), 1000, -100, 100);


      //plots[Form("%s_Centroid_Position_st%d", t, st)] = new TH2D(Form("%s_Centroid_Position_st%d", t, st), Form("%s_Centroid_Position_st%d; x (cm); y (cm)", t, st), nChannels+1, -0.5*.025, (nChannels+0.5)*.025, nChannels+1, -0.5*.025, (nChannels+0.5)*.025);
      plots[Form("%s_HitDistribution_st%d", t, st)] = new TH2D (Form("%s_HitDistribution_st%d", t, st), Form("%s_HitDistribution_st%d; n hit %dX; n hit %dY", t,  st, st, st), nChannels, 0, nChannels, nChannels, 0, nChannels);
      plots[Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d", t, st)] = new TH2D(Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d", t, st), Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d; US qdc; SciFi qdc;", t, st), 1500, -1000, 20000, 1500, -1000, 8000);  
    
      plots[Form("%s_Muon_Position_X_st%d", t, st)] = new TH1D(Form("%s_Muon_Position_X_st%d", t, st), Form("%s_Muon_Position_X_st%d; x (cm); entries", t, st), 1000, -44.75, -44.6);
      plots[Form("%s_Muon_Position_Y_st%d", t, st)] = new TH1D(Form("%s_Muon_Position_Y_st%d", t, st), Form("%s_Muon_Position_Y_st%d; y (cm); entries", t, st), 1000, 37.77, 37.93);
      plots[Form("%s_Muon_Position_Z.x_st%d", t, st)] = new TH1D(Form("%s_Muon_Position_Z.x_st%d", t, st), Form("%s_Muon_Position_Z.x_st%d; z (cm); entries", t, st), 1000, 315, 360);
      plots[Form("%s_Muon_Position_Z.y_st%d", t, st)] = new TH1D(Form("%s_Muon_Position_Z.y_st%d", t, st), Form("%s_Muon_Position_Z.y_st%d; z (cm); entries", t, st), 1000, 315, 360);
 
    }

    plots[Form("%s_Slope_Muon_X", t)] = new TH1D(Form("%s_Slope_Muon_X", t), Form("%s_Slope_Muon_X; slope; entries", t), 1000, -0.005, 0.005);
    plots[Form("%s_Slope_Muon_Y", t)] = new TH1D(Form("%s_Slope_Muon_Y", t), Form("%s_Slope_Muon_Y; slope; entries", t), 1000, -0.005, 0.005);

    plots[Form("%s_Intercept_Muon_X", t)] = new TH1D(Form("%s_Intercept_Muon_X", t), Form("%s_Intercept_Muon_X; Intercept; entries", t), 1000, -47, -42);
    plots[Form("%s_Intercept_Muon_Y", t)] = new TH1D(Form("%s_Intercept_Muon_Y", t), Form("%s_Intercept_Muon_Y; Intercept; entries", t), 1000, 35, 40);

    for (int start_st = 1; start_st < configuration.SCIFI_STATIONS+1; ++start_st){
      for (int st = start_st; st < configuration.SCIFI_STATIONS+1; ++st){
        plots[Form("%s_Centroid_Position_st%d_start%d", t, st, start_st)] = new TH2D(Form("%s_Centroid_Position_st%d_start%d", t, st, start_st), Form("%s_Centroid_Position_st%d_start%d; x (cm); y (cm)", t, st, start_st), 1000, -44.69, -44.66, 1000, 37.84, 37.9);
        plots[Form("%s_Centroid_Position_X_st%d_start%d", t, st, start_st)] = new TH1D(Form("%s_Centroid_Position_X_st%d_start%d", t, st, start_st), Form("%s_Centroid_Position_X_st%d_start%d; x (cm); entries", t, st, start_st), 1000, -44.69, -44.66);
        plots[Form("%s_Centroid_Position_Y_st%d_start%d", t, st, start_st)] = new TH1D(Form("%s_Centroid_Position_Y_st%d_start%d", t, st, start_st), Form("%s_Centroid_Position_Y_st%d_start%d; y (cm); entries", t, st, start_st), 1000, 37.84, 37.9);
        plots[Form("%s_Centroid_Position_Z.x_st%d_start%d", t, st, start_st)] = new TH1D(Form("%s_Centroid_Position_Z.x_st%d_start%d", t, st, start_st), Form("%s_Centroid_Position_Z.x_st%d_start%d; z (cm); entries", t, st, start_st), 1000, 315, 360);
        plots[Form("%s_Centroid_Position_Z.y_st%d_start%d", t, st, start_st)] = new TH1D(Form("%s_Centroid_Position_Z.y_st%d_start%d", t, st, start_st), Form("%s_Centroid_Position_Z.y_st%d_start%d; z (cm); entries", t, st, start_st), 1000, 315, 360);
      }

      plots[Form("%s_Slope_X_start%d", tag.c_str(), start_st)] = new TH1D(Form("%s_Slope_X_start%d", tag.c_str(), start_st), Form("%s_Slope_X_start%d; slope; entries", tag.c_str(), start_st), 80, -0.00004, 0.00004);
      plots[Form("%s_Slope_Y_start%d", tag.c_str(), start_st)] = new TH1D(Form("%s_Slope_Y_start%d", tag.c_str(), start_st), Form("%s_Slope_Y_start%d; slope; entries", tag.c_str(), start_st), 80, 0, 0.0016);

      plots[Form("%s_Intercept_X_start%d", tag.c_str(), start_st)] = new TH1D(Form("%s_Intercept_X_start%d", tag.c_str(), start_st), Form("%s_intercept_X_start%d; intercept; entries", tag.c_str(), start_st), 80, -44.69, -44.655);
      plots[Form("%s_Intercept_Y_start%d", tag.c_str(), start_st)] = new TH1D(Form("%s_Intercept_Y_start%d", tag.c_str(), start_st), Form("%s_intercept_Y_start%d; intercept; entries", tag.c_str(), start_st), 80, 37.3, 37.85);

    }

    for (int st = 1; st < configuration.SCIFI_STATIONS+1; ++st){
      plots[Form("%s_Centroid_Residuals_st%dX", t, st)] = new TH1D (Form("%s_Centroid_Residuals_st%dX", t, st), Form("%s_Centroid_Residuals_st%dX; x-x_ref (cm);entries", t, st), 2*nChannels + 1, -(nChannels+0.5)*.025, (nChannels+0.5)*.025);
      plots[Form("%s_Centroid_Residuals_st%dY", t, st)] = new TH1D (Form("%s_Centroid_Residuals_st%dY", t, st), Form("%s_Centroid_Residuals_st%dY; y-y_ref (cm);entries", t, st), 2*nChannels + 1, -(nChannels+0.5)*.025, (nChannels+0.5)*.025);
    }
    for (int start_st = 1; start_st < configuration.SCIFI_STATIONS+1; ++start_st){
      plots[Form("%s_Shower_SciFi_QDC_shStart%d", t, start_st)] = new TH1D (Form("%s_Shower_SciFi_QDC_shStart%d", t, start_st), Form("%s_Shower_SciFi_QDC_shStart%d; qdc sum (a.u.);entries", t, start_st), 1000, -100, 8000);
      for (int st = start_st; st < configuration.SCIFI_STATIONS+1; ++st){
        plots[Form("%s_Hits_st%d_start%d", t, st, start_st)] = new TH1D (Form("%s_Hits_st%d_start%d", t, st, start_st), Form("%s_Hits_st%d_start%d;n hit in event;entries", t, st, start_st), nChannels, 0, 2*nChannels);
      }
      for (int st = start_st+1; st < configuration.SCIFI_STATIONS+1; ++st){
        plots[Form("%s_Hits_st%d_vs_start%d", t, st, start_st)] = new TH2D (Form("%s_Hits_st%d_vs_start%d", t, st, start_st), Form("%s_Hits_st%d_vs_start%d;n hit in st%d;n hit in st%d", t, st, start_st, st, start_st), nChannels, 0, 2*nChannels, nChannels, 0, 2*nChannels);
      }
    }
    // US planes
    for (int pl{1}; pl < configuration.US_STATIONS + 1; ++pl) {
      plots[Form("%s_US_Timestamps_pl%ds", t, pl)] = new TH1D (Form("%s_US_Timestamps_pl%ds", t, pl), Form("%s_US_Timestamps_pl%ds;timestamps - tref (clk cycles);entries", t, pl), 100, -5, 15);
      plots[Form("%s_US_Timestamps_pl%dl", t, pl)] = new TH1D (Form("%s_US_Timestamps_pl%dl", t, pl), Form("%s_US_Timestamps_pl%dl;timestamps - tref (clk cycles);entries", t, pl), 100, -5, 15);
      plots[Form("%s_US_Timestamps_vs_QDC_pl%ds", t, pl)] = new TH2D (Form("%s_US_Timestamps_vs_QDC_pl%ds", t, pl), Form("%s_US_Timestamps_vs_QDC_pl%ds;timestamps - tref (clk cycles);qdc (a.u.)", t, pl), 1000, -5, 15, 1000, -20, 300);   
      plots[Form("%s_US_Timestamps_vs_QDC_pl%dl", t, pl)] = new TH2D (Form("%s_US_Timestamps_vs_QDC_pl%dl", t, pl), Form("%s_US_Timestamps_vs_QDC_pl%dl;timestamps - tref (clk cycles);qdc (a.u.)", t, pl), 1000, -5, 15, 1000, -20, 300);   
    }
  }
}

std::vector<SciFiPlaneView> fillSciFi(cfg configuration, TClonesArray *sf_hits, Scifi* ScifiDet){

  std::vector<SciFiPlaneView> scifi_planes;

  int begin{0};
  int count{0};

  int n_sf_hits{sf_hits->GetEntries()};

  for (int st{1}; st <= configuration.SCIFI_STATIONS; ++st) {
    begin = count;
    while (count < n_sf_hits &&
           st == static_cast<sndScifiHit *>(sf_hits->At(count))->GetStation()) {
      ++count;
    }

    auto plane = SciFiPlaneView(configuration, sf_hits, ScifiDet, begin, count, st);
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
  for (int pl{0}; pl < configuration.US_STATIONS; ++pl) {
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

int checkShower_with_clusters(std::vector<SciFiPlaneView> scifi_planes ) {
  //find start of shower
  for (auto &plane : scifi_planes) {
    if (plane.infoCluster()) return plane.getStation(); 
  }
  return -1;
}

int checkShower_with_density(std::vector<SciFiPlaneView> scifi_planes ) {
  //find start of shower
  for (auto &plane : scifi_planes) {
    if (plane.infoDensity(plane.getConfig().SCIFI_DENSITYWINDOW, plane.getConfig().SCIFI_DENSITYHITS)) return plane.getStation();
  }
  return -1;
}

int checkShower_with_F(std::vector<SciFiPlaneView> scifi_planes) {
  SciFiPlaneView::xy_pair<int> previous_hits;
  int k{0};
  //find start of shower
  for (auto &plane : scifi_planes) {
    auto f = plane.getConfig().SCIFI_F;
    if (k>0 && plane.sizes().x>0 && plane.sizes().y>0) { // && previous_hits.x>0 && previous_hits.y>0
      if (((static_cast<float>(previous_hits.x / (plane.sizes().x + previous_hits.x)) < f) && plane.sizes().x > 15) || ((static_cast<float>(previous_hits.y / (plane.sizes().y + previous_hits.y)) < f) && plane.sizes().y > 15)) {
        return plane.getStation();
      }
    }
    previous_hits = plane.sizes();
    k++;
  }
  return -1;
}

// timecut -> vector scifiplaneview time cut 
// nel file skimmato ho comunque sempre solo un evento in stazione 1-> leggo tempo di quello (parametro è vector scifiplaneview con tutti hit) e butto via i fuori tempo

bool hitCut (std::vector<SciFiPlaneView> &detector){
  for (auto &plane : detector){
    if (plane.getStation() == 1 && plane.sizes().x == 1 && plane.sizes().y == 1 ) return true;
    /*else if (plane.getStation() > 1){
      int thr = plane.getConfig().SCIFI_DIMCLUSTER;
      if (plane.sizes().x > thr && plane.sizes().y > thr) return true;
    }*/
  }
  return false;
}

void evaluateAveragePosition(const std::vector<SciFiPlaneView>& scifi, int clustermaxsize, int max_miss) {
  for (const auto& view : scifi) {
    if (view.evaluateNeighboringHits(clustermaxsize, max_miss) == true) {
      //std::cout << "The condition is satisfied " << std::endl;

      //std::vector<int> valid_positionsx = view.calculateValidPositionsx(clustermaxsize, max_miss);
      //std::vector<int> valid_positionsy = view.calculateValidPositionsy(clustermaxsize, max_miss);
      // if (*std::max_element(valid_positions.begin(),valid_positions.end())-*std::min_element(valid_positions.begin(),valid_positions.end())>=valid_positions.size()) {
      //  std::cout << "Posizioni valide: ";
      // for (std::size_t i = 0; i < valid_positions.size(); i++) {
      //   std::cout << valid_positions[i];
      //   if (i != valid_positions.size() - 1) {
      //     std::cout << ", ";
      //   }       

  //}
    std::vector<int> valid_positionsx = view.calculateValidPositionsx(clustermaxsize, max_miss);
    /*std::cout << "Positions in x: ";
    for (std::size_t i = 0; i < valid_positionsx.size(); i++) {
      std::cout << valid_positionsx[i];
      if (i != valid_positionsx.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << std::endl;
    */

    std::vector<int> valid_positionsy = view.calculateValidPositionsy(clustermaxsize, max_miss);
    /*std::cout << "Positions in in y: ";
    for (std::size_t i = 0; i < valid_positionsy.size(); i++) {
      std::cout << valid_positionsy[i];
      if (i != valid_positionsy.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << std::endl; */
  }
  }
}

void timeCut (std::vector<SciFiPlaneView> &Scifi, std::vector<USPlaneView> &US) {
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
      plane.timeCut(referenceTime, referenceTime + plane.getConfig().SCIFI_TIMECUT);
    }
  }
  for (auto &plane : US) {
    if (referenceTime == -1) continue;
    plane.timeCut(referenceTime);
  }
}

double timeCutGuil (std::vector<SciFiPlaneView> &Scifi, std::vector<USPlaneView> &US) {
  TH1D* times = new TH1D ("times", "times; clk cycles; entries", 1000, 0, 50);

  for (auto &plane : Scifi) {
    auto time = plane.getTime();
    for (int i{0}; i<time.x.size(); ++i) {
      if (time.x[i]>DEFAULT) {
        times->Fill(time.x[i]);
      }
      if (time.y[i]>DEFAULT) {
        times->Fill(time.y[i]);
      }
    }
  }

  double referenceTime = times->GetXaxis()->GetBinCenter(times->GetMaximumBin()); // in clk cycles

  for (auto &plane : Scifi) {
    plane.timeCut(referenceTime - 0.5, referenceTime + 0.5);
  }

  delete times;

  for (auto &plane : US) {
    plane.timeCut(referenceTime);
  }

  return referenceTime;
}

void fitAndStoreSlopes(const std::vector<double> &x_positions, const std::vector<double> &z_positions_x, const std::vector<double> &y_positions, const std::vector<double> &z_positions_y, std::map<std::string, TH1*> &plots, const std::string &t, int shStart);

void fitMuonSlopes(const std::vector<double> &mux_positions, const std::vector<double> &muz_positions_x, const std::vector<double> &muy_positions, const std::vector<double> &muz_positions_y, std::map<std::string, TH1*> &plots, const std::string &t);


void fillPlots (std::vector<SciFiPlaneView> &Scifi_detector, std::vector<USPlaneView> US, std::map<std::string, TH1*> &plots, std::string &t, int shStart, double referenceTime=0) {
  int showerHits{0};
  double ScifiQDCSum{0}, partialScifiQDCSum{0};
  double USQDCSum{0};
  double Small_USQDCSum{0}, Large_USQDCSum{0};
  auto refCentroid{Scifi_detector[0].getCentroid()};

  std::vector<double> x_positions, y_positions, z_positions_x, z_positions_y;

  std::vector<double> mux_positions, muy_positions, muz_positions_x, muz_positions_y;
  
  for (auto plane : Scifi_detector){
    auto centroid{plane.getCentroid()};
    auto centroid_depth{plane.getCentroidDepth()};
    auto muonpos{plane.getMuon()};
    auto muondepth{plane.getMuonDepth()};
    auto station{plane.getStation()};
    auto pos{plane.getGeometry()};


    bool xValid = std::any_of(pos.x.begin(), pos.x.end(), [](double value) { return value > DEFAULT; });
    bool yValid = std::any_of(pos.y.begin(), pos.y.end(), [](double value) { return value > DEFAULT; });

    if (shStart != -1 && xValid && yValid) {
        for (const auto& x : pos.x) {
            if (x > DEFAULT) {
                //std::cout <<"X: " <<x <<" nel piano " <<plane.getStation() <<std::endl;
                plots[Form("%s_Hits_Position_st%dX", t.c_str(), plane.getStation())]->Fill(x);
            }
        }
        for (const auto& y : pos.y) {
            if (y > DEFAULT) {
                //std::cout <<"Y: " <<y <<" nel piano " <<plane.getStation() <<std::endl;
                plots[Form("%s_Hits_Position_st%dY", t.c_str(), plane.getStation())]->Fill(y);
            }
        }
    }


    if (shStart != -1 && station >= shStart && centroid.x > DEFAULT && centroid.y > DEFAULT) {
        plots[Form("%s_Centroid_Position_st%d_start%d", t.c_str(), station, shStart)]->Fill(centroid.x, centroid.y);
        plots[Form("%s_Centroid_Position_X_st%d_start%d", t.c_str(), station, shStart)]->Fill(centroid.x);
        plots[Form("%s_Centroid_Position_Y_st%d_start%d", t.c_str(), station, shStart)]->Fill(centroid.y);
        plots[Form("%s_Centroid_Position_Z.x_st%d_start%d", t.c_str(), station, shStart)]->Fill(centroid_depth.x);
        plots[Form("%s_Centroid_Position_Z.y_st%d_start%d", t.c_str(), station, shStart)]->Fill(centroid_depth.y);
    
        x_positions.push_back(centroid.x);
        y_positions.push_back(centroid.y);
        z_positions_x.push_back(centroid_depth.x);
        z_positions_y.push_back(centroid_depth.y);
    }

    if (muonpos.x > DEFAULT && muonpos.y > DEFAULT) {
        plots[Form("%s_Muon_Position_X_st%d", t.c_str(), station)]->Fill(muonpos.x);
        plots[Form("%s_Muon_Position_Y_st%d", t.c_str(), station)]->Fill(muonpos.y);
        plots[Form("%s_Muon_Position_Z.x_st%d", t.c_str(), station)]->Fill(muondepth.x);
        plots[Form("%s_Muon_Position_Z.y_st%d", t.c_str(), station)]->Fill(muondepth.y);
    
        mux_positions.push_back(muonpos.x);
        muy_positions.push_back(muonpos.y);
        muz_positions_x.push_back(muondepth.x);
        muz_positions_y.push_back(muondepth.y);
    }

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
    auto timeX{plane.getTime().x};
    auto timeY{plane.getTime().y};
    auto qdcX{plane.getQDC().x};
    auto qdcY{plane.getQDC().y};
    auto sumScifi = plane.getTotQDC();
    auto clusterSize{plane.getClusterSize()};

    ScifiQDCSum += sumScifi.x;
    ScifiQDCSum += sumScifi.y;
    if (station >= shStart) {
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
    for (int i{0}; i<plane.getConfig().SCIFI_NCHANNELS; ++i) {
      if (qdcX[i] != DEFAULT) {
        plots[Form("%s_Signals_st%dX", t.c_str(), station)]->Fill(qdcX[i]);
        plots[Form("%s_Times", t.c_str())]->Fill(timeX[i]);
        plots[Form("%s_Position_st%dX", t.c_str(), station)]->Fill(i*0.025);
        ((TH2D*)plots[Form("%s_Shower_development_X", t.c_str())])->Fill(i*0.025, station + std::min(timeX[i],4.99)/5, qdcX[i]); 
        plots[Form("%s_Tofpet_st%dX", t.c_str(), station)]->Fill(static_cast<int>(i/64));
        if (station == 1 && nhitsX == 1 && nhitsY == 1){
          refCentroid.x = i*0.025;
        }
      }
      if (qdcY[i] != DEFAULT) {
        plots[Form("%s_Signals_st%dY", t.c_str(), station)]->Fill(qdcY[i]);
        plots[Form("%s_Times", t.c_str())]->Fill(timeY[i]);
        plots[Form("%s_Position_st%dY", t.c_str(), station)]->Fill(i*0.025);
        ((TH2D*)plots[Form("%s_Shower_development_Y", t.c_str())])->Fill(i*0.025, station + std::min(timeY[i],4.99)/5, qdcY[i]);
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
    auto qdc = plane.getQDC();
    int pl = plane.getStation();
    Small_USQDCSum += sumQDC.s;
    Large_USQDCSum += sumQDC.l;
    USQDCSum += (sumQDC.s + sumQDC.l);
    for (int i{0}; i<plane.getConfig().US_NCHANNELS; ++i) {
      if (timesUS[i] != DEFAULT) {
        if ((i%16)%8==2 || (i%16)%8==5) {
          plots[Form("%s_US_Timestamps_pl%ds", t.c_str(), pl)]->Fill(timesUS[i] - referenceTime);
          plots[Form("%s_US_Timestamps_vs_QDC_pl%ds", t.c_str(), pl)]->Fill(timesUS[i] - referenceTime, qdc[i]);
        }
        else {
          plots[Form("%s_US_Timestamps_pl%dl", t.c_str(), pl)]->Fill(timesUS[i] - referenceTime);    ////// TO BE CHANGED
          plots[Form("%s_US_Timestamps_vs_QDC_pl%dl", t.c_str(), pl)]->Fill(timesUS[i] - referenceTime, qdc[i]);
        }
      }
    }    
  }
  plots[Form("%s_QDCUS_vs_QDCScifi", t.c_str())]->Fill(Large_USQDCSum, ScifiQDCSum);  // only large?
  if (shStart > 0) {
    plots[Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d", t.c_str(), shStart)]->Fill(Large_USQDCSum, partialScifiQDCSum); // only large?
    plots[Form("%s_Shower_SciFi_QDC_shStart%d", t.c_str(), shStart)]->Fill(partialScifiQDCSum);
  }

  fitAndStoreSlopes(x_positions, z_positions_x, y_positions, z_positions_y, plots, t, shStart);

  fitMuonSlopes(mux_positions, muz_positions_x, muy_positions, muz_positions_y, plots, t);
  
  //std::cout<<"SciFi:\t"<<partialScifiQDCSum*0.063194<<"\t US:\t"<<USQDCSum*0.0130885<<"\t Tot Energy:\t"<<partialScifiQDCSum*0.063194 + USQDCSum*0.0130885<<"\n";
}

void fitAndStoreSlopes(const std::vector<double> &x_positions, const std::vector<double> &z_positions_x, const std::vector<double> &y_positions, const std::vector<double> &z_positions_y, std::map<std::string, TH1*> &plots, const std::string &t, int shStart) {
    if (x_positions.size() < 2 || y_positions.size() < 2) return; 

    TGraph *graph_x = new TGraph(z_positions_x.size(), &z_positions_x[0], &x_positions[0]);
    TF1 *fit_x = new TF1("fit_x", "pol1", z_positions_x.front(), z_positions_x.back());
    graph_x->Fit(fit_x, "Q");
    double slope_x = fit_x->GetParameter(1);
    double intercept_x = fit_x->GetParameter(0);
    plots[Form("%s_Slope_X_start%d", t.c_str(), shStart)]->Fill(slope_x);
    plots[Form("%s_Intercept_X_start%d", t.c_str(), shStart)]->Fill(intercept_x);

    TGraph *graph_y = new TGraph(z_positions_y.size(), &z_positions_y[0], &y_positions[0]);
    TF1 *fit_y = new TF1("fit_y", "pol1", z_positions_y.front(), z_positions_y.back());
    graph_y->Fit(fit_y, "Q");
    double slope_y = fit_y->GetParameter(1);
    double intercept_y = fit_y->GetParameter(0);
    plots[Form("%s_Slope_Y_start%d", t.c_str(), shStart)]->Fill(slope_y);
    plots[Form("%s_Intercept_Y_start%d", t.c_str(), shStart)]->Fill(intercept_y);

    // Clean up
    delete graph_x;
    delete graph_y;
    delete fit_x;
    delete fit_y;
}

void fitMuonSlopes(const std::vector<double> &mux_positions, const std::vector<double> &muz_positions_x, const std::vector<double> &muy_positions, const std::vector<double> &muz_positions_y, std::map<std::string, TH1*> &plots, const std::string &t) {
    if (mux_positions.size() == 4 && muy_positions.size() == 4) { 

      TGraph *graph_x = new TGraph(muz_positions_x.size(), &muz_positions_x[0], &mux_positions[0]);
      TF1 *fit_x = new TF1("fit_x", "pol1", muz_positions_x.front(), muz_positions_x.back());
      graph_x->Fit(fit_x, "Q");
      double muslope_x = fit_x->GetParameter(1);
      double muintercept_x = fit_x->GetParameter(0);      
      plots[Form("%s_Slope_Muon_X", t.c_str())]->Fill(muslope_x);
      plots[Form("%s_Intercept_Muon_X", t.c_str())]->Fill(muintercept_x);

      TGraph *graph_y = new TGraph(muz_positions_y.size(), &muz_positions_y[0], &muy_positions[0]);
      TF1 *fit_y = new TF1("fit_y", "pol1", muz_positions_y.front(), muz_positions_y.back());
      graph_y->Fit(fit_y, "Q");
      double muslope_y = fit_y->GetParameter(1);
      double muintercept_y = fit_y->GetParameter(0);  
      plots[Form("%s_Slope_Muon_Y", t.c_str())]->Fill(muslope_y);
      plots[Form("%s_Intercept_Muon_Y", t.c_str())]->Fill(muintercept_y);

      delete graph_x;
      delete graph_y;
      delete fit_x;
      delete fit_y;

    }
    else return;
}

void fitGaussianAndStoreSlopesIntercepts(std::map<std::string, TH1*> &plots) {
    for (auto &plot : plots) {
        const std::string& histName = plot.first;

        // Check if the histogram name corresponds to slope or intercept
        if (histName.find("Slope") != std::string::npos || histName.find("Intercept") != std::string::npos) {
            if (plot.second->GetEntries() > 0) {
                // Perform Gaussian fit
                TF1 *gausFit = new TF1("gausFit", "gaus");
                plot.second->Fit(gausFit, "Q");

                // Extract fit parameters
                double mean = gausFit->GetParameter(1);
                double sigma = gausFit->GetParameter(2);

                // Print fit results to the console
                std::cout << "Histogram: " << histName 
                          << " | Mean: " << mean 
                          << " | Sigma: " << sigma << std::endl;
            }
        }
    }
}

void read_geo(Scifi*& ScifiDet, MuFilter*& MufiDet, const std::string geoFilePath) {
    TPython::Exec("import SndlhcGeo");
    TPython::Exec(("SndlhcGeo.GeoInterface('" + geoFilePath + "')").c_str());

    // Init detectors
    ScifiDet = new Scifi("Scifi", kTRUE);
    MufiDet = new MuFilter("MuFilter", kTRUE);

    // Retrieve the detectors from ROOT's global list
    ScifiDet = (Scifi*)gROOT->GetListOfGlobals()->FindObject("Scifi");
    MufiDet = (MuFilter*)gROOT->GetListOfGlobals()->FindObject("MuFilter");

    // Print some configuration parameters as checks
    std::cout << ScifiDet->GetConfParF("Scifi/station2t") << std::endl;
    std::cout << ScifiDet->GetConfParI("Scifi/nscifi") << std::endl;
}

void runAnalysis(int runNumber, int nFiles, bool isTB, bool isMulticore = false, int target = -1) //(int runN, int partN)
{

  auto start = std::chrono::system_clock::now();
  auto now = std::chrono::system_clock::to_time_t(start);
  if (target == -1){
    std::cout << "Start: " << std::ctime(&now)  << "\n" <<std::flush;
  }
  else{
    std::cout << "Run: " << runNumber << "\t File: "  << nFiles << "\t Evt: " << target << "\n" <<std::flush;
  }

  // ##################### Set right parameters for data type (TB/TI18) #####################
  cfg configuration = setCfg(isTB, runNumber);
  
  // ##################### Read file #####################

  auto *fEventTree = new TChain("rawConv");
  TFile* outputFile;
  //TFile *treeFile; = new TFile("tree.root", "RECREATE");
  if (isMulticore){
    fEventTree->Add(Form("%srun_%06d/sndsw_raw-%04d.root", configuration.INFILENAME.c_str(), runNumber, nFiles));
    if (target != -1){
      outputFile = new TFile(Form("%sRun_%d_%d_%d.root", configuration.OUTFILENAME.c_str(), runNumber, nFiles, target), "RECREATE");
    }
    else {
      outputFile = new TFile(Form("%sRun_%d_%d.root", configuration.OUTFILENAME.c_str(), runNumber, nFiles), "RECREATE");
    }
    //treeFile = new TFile(Form("%s_tree_Run_%d_%d.root", configuration.OUTFILENAME.c_str(), runNumber, nFiles), "RECREATE");
  }
  else {
    for (int i = 0; i<nFiles; ++i){
      fEventTree->Add(Form("%srun_%06d/sndsw_raw-%04d.root", configuration.INFILENAME.c_str(), runNumber, i)); 
    }
    outputFile = new TFile(Form("%sRun_%d.root", configuration.OUTFILENAME.c_str(), runNumber), "RECREATE");
    //treeFile = new TFile(Form("%s_tree_Run_%d.root", configuration.OUTFILENAME.c_str(), runNumber), "RECREATE");
  }
 
  outputFile->cd();

  std::map<std::string, double> counters;
  std::map<std::string, TH1*> plots;
  std::vector<std::string> tags;
  tags.push_back("NoCut");
  tags.push_back("Cut");
  tags.push_back("GuilCut");
  //tags.push_back("Cluster");
  std::vector<std::string> shower_tags;
  shower_tags.push_back("clusters");
  shower_tags.push_back("density");
  shower_tags.push_back("F");


  definePlots(configuration, plots, counters, tags, shower_tags);
  
  // ##################### Read hits from Scifi and Mufilter  #####################

  auto mu_hits = new TClonesArray("MuFilterHit");
  fEventTree->SetBranchAddress("Digi_MuFilterHits", &mu_hits);
  auto sf_hits = new TClonesArray("sndScifiHit");
  fEventTree->SetBranchAddress("Digi_ScifiHits", &sf_hits);
  auto header = new SNDLHCEventHeader();
  if (!isTB && runNumber<5422){
    fEventTree->SetBranchAddress("EventHeader", &header);
  }
  else {
    fEventTree->SetBranchAddress("EventHeader.", &header);
  }

  //geom for TI18
  std::string geoFilePath = "/eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V3_2023.root";

  //geom for Test Beam
  // std::string geoFilePath = "/eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/geofile_sndlhc_H8_2023_3walls.root";
  
  Scifi* ScifiDet = nullptr;
  MuFilter* MufiDet = nullptr;
  
  read_geo(ScifiDet, MufiDet, geoFilePath);
  ScifiDet->InitEvent(header);

  // ###################### Create tree to store shower tagging info   ####################

  // TTree *tree = new TTree("ShowerTags", "ShowerTags");
  // Int_t run_number, event_number, wall;
  // tree->Branch("run_number", &run_number);
  // tree->Branch("event_number", &event_number);
  // tree->Branch("wall", &wall);
  // run_number = runNumber;

  const float TDC2ns{1000/160.316};
  double last_timestamp{-1};
  bool is_one_hit = true;
  bool is_apart = true;
  int first{target};
  int last{target+1};

  // Loop over events
  
  int iMax = fEventTree->GetEntries();
  if (target == -1) {
    first = 0;
    last = iMax;
  }
  //for ( int m = 0; m < iMax; m++ ){ 
  for ( int m = first; m < last; m++ ){ 
    if (m % 100 == 0) std::cout << "Processing event: " << m << '\r' << std::flush;
    //if (m >1000000) break;
    fEventTree->GetEntry(m);

    // Check timestamp difference with previous event
    if (last_timestamp != -1) {
      plots["Delta_timestamp"]->Fill((header->GetEventTime()-last_timestamp)*TDC2ns);
      if (((header->GetEventTime()-last_timestamp)*TDC2ns) > 150) {
        is_apart = true;
      }
      else {
        is_apart = false;
      }
    }
    last_timestamp = header->GetEventTime();
    

    int sf_max=sf_hits->GetEntries();
    int mu_max=mu_hits->GetEntries();

    // if analyze showers -> 15
    // if analyze Muons -> 8
    if (sf_max < 15) continue;
    // if (sf_max < 8) continue;
    
    auto scifi_planes = fillSciFi(configuration, sf_hits, ScifiDet);
    auto scifi_planes_guil = scifi_planes;
    auto us_planes = fillUS(configuration, mu_hits);
    auto us_planes_guil = us_planes;

    //Before cut
    // int showerStart = checkShower_with_clusters(scifi_planes);
    // plots[Form("%s_ShowerStart_with_clusters", tags[0].c_str())]->Fill(showerStart);
    // plots[Form("%s_ShowerStart_with_density", tags[0].c_str())]->Fill(checkShower_with_density(scifi_planes));
    // plots[Form("%s_ShowerStart_with_F", tags[0].c_str())]->Fill(checkShower_with_F(scifi_planes));
    // for (auto &plane : scifi_planes)  plane.findCentroid();
    // fillPlots(scifi_planes, us_planes, plots, tags[0], showerStart);

    std::vector<int> sh_start(6, -2);
    //After cut
    // is_one_hit = hitCut(scifi_planes);
    // if (is_one_hit) {
    //   timeCut(scifi_planes, us_planes);
    //   sh_start[0] = checkShower_with_clusters(scifi_planes);
    //   sh_start[1] = checkShower_with_density(scifi_planes);
    //   sh_start[2] = checkShower_with_F(scifi_planes);

    //   plots[Form("%s_ShowerStart_with_clusters", tags[1].c_str())]->Fill(sh_start[0]);
    //   plots[Form("%s_ShowerStart_with_density", tags[1].c_str())]->Fill(sh_start[1]);
    //   plots[Form("%s_ShowerStart_with_F", tags[1].c_str())]->Fill(sh_start[2]);
    //   for (auto &plane : scifi_planes) plane.findCentroid();
    //   fillPlots(scifi_planes, us_planes, plots, tags[1], sh_start[0]);

    // }
    // else {
    //   evaluateAveragePosition(scifi_planes, 4, 1);
    // }

    if (is_apart) {
      double refT = timeCutGuil(scifi_planes_guil, us_planes_guil);
      sh_start[3] = checkShower_with_clusters(scifi_planes_guil);
      sh_start[4] = checkShower_with_density(scifi_planes_guil);
      sh_start[5] = checkShower_with_F(scifi_planes_guil);

      // if (target != -1){
      //   std::cout<<"RUN "<<runNumber<<"\t ev:\t"<<m<<"\t clusters:\t"<<sh_start[3]<<"\t density:\t"<<sh_start[4]<<"\t F:\t"<<sh_start[5]<<std::endl;
      // }
      // if (sh_start[4] > 0 && isTB) {
      //   event_number = header->GetEventNumber();
      //   wall = sh_start[4] - 1;
      //   tree->Fill();
      // }
      // if (sh_start[4] > 0 && !isTB) {
      //   event_number = header->GetEventNumber();
      //   wall = sh_start[4];
      //   tree->Fill();
      // }

      // plots[Form("%s_ShowerStart_with_clusters", tags[2].c_str())]->Fill(sh_start[3]);
      plots[Form("%s_ShowerStart_with_density", tags[2].c_str())]->Fill(sh_start[4]);
      // plots[Form("%s_ShowerStart_with_F", tags[2].c_str())]->Fill(sh_start[5]);
      for (auto &plane : scifi_planes_guil) plane.findCentroid();
      for (auto &plane : scifi_planes_guil) plane.findMuon();
      fillPlots(scifi_planes_guil, us_planes_guil, plots, tags[2], sh_start[4], refT);  
      }

    // for (int k{0}; k < shower_tags.size(); ++k) {
    //   plots[Form("%s_%s_vs_%s_%s", tags[1].c_str(), shower_tags[k].c_str(), tags[2].c_str(), shower_tags[k].c_str())]->Fill(sh_start[k],sh_start[k+3]);
    // }

    //Guil cut density vs F
    // plots[Form("%s_%s_vs_%s_%s", tags[2].c_str(), shower_tags[1].c_str(), tags[2].c_str(), shower_tags[2].c_str())]->Fill(sh_start[4],sh_start[5]);

    
    //After cluster
    /*
    for (auto &plane : scifi_planes){
      plane.findCluster();
      plane.findCentroid(6);
    }
    showerStart = checkShower_with_clusters(scifi_planes);
    plots[Form("%s_ShowerStart_with_clusters", tags[3].c_str())]->Fill(showerStart);
    plots[Form("%s_ShowerStart_with_density", tags[3].c_str())]->Fill(checkShower_with_density(scifi_planes));
    plots[Form("%s_ShowerStart_with_F", tags[3].c_str())]->Fill(checkShower_with_F(scifi_planes));
    fillPlots(scifi_planes, us_planes, plots, tags[3], showerStart);
    */
  }

  auto stop = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = stop-start;
  auto end = std::chrono::system_clock::to_time_t(stop);
  if (target == -1){
    std::cout << "\nDone: " << std::ctime(&end)  << std::endl;
  }

  // ##################### Write results to file #####################
  outputFile->Write();
  outputFile->Close();
}
