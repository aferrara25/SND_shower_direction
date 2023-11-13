#ifndef SciFiPlaneView_h
#define SciFiPlaneView_h

#include<vector>
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

  SciFiPlaneView(cfg c, TClonesArray *h, int b, int e, int s);
  auto sizes() const ;

  const int getStation() const;
  const cfg getConfig() const;
  const int getBegin() const;
  const int getEnd() const;

};
#endif