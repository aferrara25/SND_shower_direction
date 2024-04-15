#ifndef SciFiPlaneView_h
#define SciFiPlaneView_h

#include<vector>
#include<array>
#include<string>
#include "TChain.h"
#include "TClonesArray.h"



struct cfg
{
  int SCIFI_STATIONS{-1};
  int SCIFI_BOARDPERPLANE{-1};
  int SCIFI_NCHANNELS{-1};
  double SCIFI_TIMECUT{-1};  
  int SCIFI_DIMCLUSTER{-1};
  int SCIFI_GAPCLUSTER{-1};
  int SCIFI_DENSITYWINDOW{-1};
  int SCIFI_DENSITYHITS{-1};
  double SCIFI_F{-1}; 

  int US_STATIONS{-1};
  const int US_NCHANNELS{160};
  const int US_NSIPM{8};
  const int US_NSIDES{2};
  double US_TIMECUT{-1};

  //geometry parameters
  double SCIFI_DIM{-1};
  
  std::string INFILENAME;
  std::string OUTFILENAME;
};



class SciFiPlaneView {

public:
  template <class T> struct xy_pair {
    T x{};
    T y{};
  };

  SciFiPlaneView(cfg c, TClonesArray *h, int b, int e, int s);
  
  const xy_pair<int> sizes() const;
  const int getStation() const;
  const cfg getConfig() const;
  const int getBegin() const;
  const int getEnd() const;
  const xy_pair<double> getCentroid() const;
  const xy_pair<std::vector<double>> getTime() const;
  const xy_pair<std::vector<double>> getQDC() const;
  const xy_pair<double> getTotQDC() const;
  const xy_pair<int>  getClusterSize() const;

  void findCluster();
  bool infoCluster();
  bool infoDensity(int window, int min_hits);
  void findCentroid(int windowSize);
  void timeCut(double minTime, double maxTime);


  double evaluateNeighboringHits(int windowSize, int min_hits);

private:
  xy_pair<std::vector<double>> qdc;
  xy_pair<std::vector<double>> hitTimestamps;
  xy_pair<int> clusterBegin;
  xy_pair<int> clusterEnd;
  xy_pair<double> centroid;

  int begin{};
  int end{};
  int station{};
  cfg config;

  TClonesArray *sf_hits{nullptr};

  void fillQDC();
  void fillTimestamps();
  void resetHit(bool isVertical, int index);

};
#endif
