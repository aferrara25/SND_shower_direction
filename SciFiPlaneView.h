#ifndef SciFiPlaneView_h
#define SciFiPlaneView_h

#include<vector>
#include<array>
#include<string>
#include "TChain.h"
#include "TClonesArray.h"



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

  double TIMECUT{-1};
  double MUTIMECUT{-1};

  //geometry parameters
  double SCIFIDIM{-1};
  
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
  const xy_pair<std::array<double, 512>> getTime() const;
  const xy_pair<std::array<double, 512>> getQDC() const;
  const xy_pair<double> getTotQDC() const;
  const xy_pair<int>  getClusterSize() const;

  void findCluster();
  bool infoCluster();
  bool infoDensity(int radius, int min_hits);
  void findCentroid(int windowSize);
  void timeCut(double referenceTime);

private:
  xy_pair<std::array<double,512>> qdc;
  xy_pair<std::array<double,512>> hitTimestamps;
  xy_pair<int> clusterBegin;
  xy_pair<int> clusterEnd;
  xy_pair<double> centroid;

  int begin{};
  int end{};
  int station{};

  TClonesArray *sf_hits{nullptr};
  cfg config;

  void fillQDC();
  void fillTimestamps();
  void resetHit(bool isVertical, int index);

};
#endif
