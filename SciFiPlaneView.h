#ifndef SciFiPlaneView_h
#define SciFiPlaneView_h

#include<vector>
#include<array>
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

  //geometry parameters
  double SCIFIDIM{-1};
  
  const char *INFILENAME;
  const char *OUTFILENAME;
};



class SciFiPlaneView {

  int begin{};
  int end{};
  int station{};

  TClonesArray *sf_hits{nullptr};
  cfg config;

public:
  template <class T> struct xy_pair {
    T x{};
    T y{};
  };

  xy_pair<std::array<double,512>> qdc;

  SciFiPlaneView(cfg c, TClonesArray *h, int b, int e, int s);
  auto sizes() const;
  void fillQDC();

  const int getStation() const;
  const cfg getConfig() const;
  const int getBegin() const;
  const int getEnd() const;

};
#endif