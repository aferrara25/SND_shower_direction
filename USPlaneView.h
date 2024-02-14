#ifndef USPlaneView_h
#define USPlaneView_h

#include<vector>
#include<array>
#include<string>
#include "TChain.h"
#include "TClonesArray.h"

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
  const std::vector<double> getTime() const;
  const std::vector<double> getQDC() const;
  const sl_pair<double> getTotQDC() const;

  void timeCut(double referenceTime);

private:
  std::vector<double> qdc;
  std::vector<double> hitTimestamps;

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