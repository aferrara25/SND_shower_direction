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
      plots[Form("%s_Centroid_Position_st%d", t, st)] = new TH2D(Form("%s_Centroid_Position_st%d", t, st), Form("%s_Centroid_Position_st%d; x (cm); y (cm)", t, st), nChannels+1, -0.5*.025, (nChannels+0.5)*.025, nChannels+1, -0.5*.025, (nChannels+0.5)*.025);
      plots[Form("%s_HitDistribution_st%d", t, st)] = new TH2D (Form("%s_HitDistribution_st%d", t, st), Form("%s_HitDistribution_st%d; n hit %dX; n hit %dY", t,  st, st, st), nChannels, 0, nChannels, nChannels, 0, nChannels);
      plots[Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d", t, st)] = new TH2D(Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d", t, st), Form("%s_QDCUS_vs_QDCScifi_ShStart_st%d; US qdc; SciFi qdc;", t, st), 1500, -1000, 20000, 1500, -1000, 8000);
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


std::vector<SciFiPlaneView> fillSciFi(cfg configuration, TClonesArray *sf_hits){

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

double evaluateNeighboringHits (std::vector<SciFiPlaneView> &Scifi, std::vector<SciFiPlaneView>) {
  double position = 0.0;
  if (position != DEFAULT) {
    std::cout <<"Posizione media dei vicini: " <<position <<std::endl;
  } else {
    std::cout <<"Non ci sono abbastanza hit vicini" <<std::endl;
  }
std::cout<<"SciFi:\t"<<position <<std::endl
return 0;
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

void fillPlots (std::vector<SciFiPlaneView> &Scifi_detector, std::vector<USPlaneView> US, std::map<std::string, TH1*> &plots, std::string &t, int shStart, double referenceTime=0) {
  int showerHits{0};
  double ScifiQDCSum{0}, partialScifiQDCSum{0};
  double USQDCSum{0};
  double Small_USQDCSum{0}, Large_USQDCSum{0};
  auto refCentroid{Scifi_detector[0].getCentroid()};
  
  for (auto plane : Scifi_detector){
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
  //std::cout<<"SciFi:\t"<<partialScifiQDCSum*0.063194<<"\t US:\t"<<USQDCSum*0.0130885<<"\t Tot Energy:\t"<<partialScifiQDCSum*0.063194 + USQDCSum*0.0130885<<"\n";
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
    //if (m % 100 == 0) std::cout << "Processing event: " << m << '\r' << std::flush;
    //if (m >1000) break;
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

    if (sf_max < 15) continue;
    
    auto scifi_planes = fillSciFi(configuration, sf_hits);
    auto scifi_planes_guil = scifi_planes;
    auto us_planes = fillUS(configuration, mu_hits);
    auto us_planes_guil = us_planes;

    //Before cut
    int showerStart = checkShower_with_clusters(scifi_planes);
    plots[Form("%s_ShowerStart_with_clusters", tags[0].c_str())]->Fill(showerStart);
    plots[Form("%s_ShowerStart_with_density", tags[0].c_str())]->Fill(checkShower_with_density(scifi_planes));
    plots[Form("%s_ShowerStart_with_F", tags[0].c_str())]->Fill(checkShower_with_F(scifi_planes));
    for (auto &plane : scifi_planes)  plane.findCentroid(6);
    fillPlots(scifi_planes, us_planes, plots, tags[0], showerStart);

    std::vector<int> sh_start(6, -2);
    //After cut
    is_one_hit = hitCut(scifi_planes);
    if (is_one_hit) {
      timeCut(scifi_planes, us_planes);
      sh_start[0] = checkShower_with_clusters(scifi_planes);
      sh_start[1] = checkShower_with_density(scifi_planes);
      sh_start[2] = checkShower_with_F(scifi_planes);

      plots[Form("%s_ShowerStart_with_clusters", tags[1].c_str())]->Fill(sh_start[0]);
      plots[Form("%s_ShowerStart_with_density", tags[1].c_str())]->Fill(sh_start[1]);
      plots[Form("%s_ShowerStart_with_F", tags[1].c_str())]->Fill(sh_start[2]);
      for (auto &plane : scifi_planes) plane.findCentroid(6);
      fillPlots(scifi_planes, us_planes, plots, tags[1], sh_start[0]);

    }

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

      plots[Form("%s_ShowerStart_with_clusters", tags[2].c_str())]->Fill(sh_start[3]);
      plots[Form("%s_ShowerStart_with_density", tags[2].c_str())]->Fill(sh_start[4]);
      plots[Form("%s_ShowerStart_with_F", tags[2].c_str())]->Fill(sh_start[5]);
      for (auto &plane : scifi_planes_guil) plane.findCentroid(6);
      fillPlots(scifi_planes_guil, us_planes_guil, plots, tags[2], sh_start[4], refT);
    }

    for (int k{0}; k < shower_tags.size(); ++k) {
      plots[Form("%s_%s_vs_%s_%s", tags[1].c_str(), shower_tags[k].c_str(), tags[2].c_str(), shower_tags[k].c_str())]->Fill(sh_start[k],sh_start[k+3]);
    }

    //Guil cut density vs F
    plots[Form("%s_%s_vs_%s_%s", tags[2].c_str(), shower_tags[1].c_str(), tags[2].c_str(), shower_tags[2].c_str())]->Fill(sh_start[4],sh_start[5]);

    
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
