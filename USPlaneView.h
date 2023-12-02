#ifndef USPlaneView_h
#define USPlaneView_h

#include<vector>
#include<array>
#include<string>
#include "TChain.h"
#include "TClonesArray.h"

const int NCHANNELS = 160;
/*
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
  
  std::string INFILENAME;
  std::string OUTFILENAME;
};
*/


class USPlaneView {

public:
  template <class T> struct sl_pair {
    T s{};
    T l{};
  };
  USPlaneView(cfg c, TClonesArray *h, int b, int e, int s);
  
  const sl_pair<int> sizes() const;
  const int getStation() const;
  const cfg getConfig() const;
  const int getBegin() const;
  const int getEnd() const;
  const std::array<double, NCHANNELS> getTime() const;
  const std::array<double, NCHANNELS> getQDC() const;
  const sl_pair<double> getTotQDC() const;

  void timeCut(double referenceTime);

private:
  std::array<double,NCHANNELS> qdc;
  std::array<double,NCHANNELS> hitTimestamps;

  int begin{};
  int end{};
  int station{};

  TClonesArray *mufi_hits{nullptr};
  cfg config;

  void fillQDC();
  void fillTimestamps();
  void resetHit(int index);

};
#endif