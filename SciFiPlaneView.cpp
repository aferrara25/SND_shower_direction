#include "SciFiPlaneView.h"
#include <stdexcept>



SciFiPlaneView::SciFiPlaneView(cfg c, TClonesArray *h, int b, int e,
                int s) : config(c), sf_hits(h), begin(b), end(e), station(s) {
    if (b > e) {
        throw std::runtime_error{"Begin index > end index"};
    }
}


auto SciFiPlaneView::sizes() const{
    xy_pair<int> counts{0, 0};

    for (int i{begin}; i < end; ++i) {
      if (static_cast<sndScifiHit *>(sf_hits->At(i))->isVertical()) {
        ++counts.y;
      } else {
        ++counts.x;
      }
    }
    return counts;
}

const int SciFiPlaneView::getStation() const {
    return station;
}

const cfg SciFiPlaneView::getConfig() const {
    return config;
}

const int SciFiPlaneView::getBegin() const {
  return begin;
}

const int SciFiPlaneView::getEnd() const {
  return end;
}
