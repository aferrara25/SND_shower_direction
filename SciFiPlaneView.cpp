#include "SciFiPlaneView.h"
#include <stdexcept>



SciFiPlaneView::SciFiPlaneView(cfg c, TClonesArray *h, int b, int e,
                int s) : config(c), sf_hits(h), begin(b), end(e), station(s) {
    if (b > e) {
        throw std::runtime_error{"Begin index > end index"};
    }
    std::fill(qdc.x.begin(), qdc.x.end(), -999);
    std::fill(qdc.y.begin(), qdc.y.end(), -999);
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

void SciFiPlaneView::fillQDC() {
    for (int i{begin}; i<end; ++i) {
        int position = 64*static_cast<sndScifiHit *>(sf_hits->At(i))->GetTofpetID(0) + 63 - static_cast<sndScifiHit *>(sf_hits->At(i))->Getchannel(0);
        if (static_cast<sndScifiHit *>(sf_hits->At(i))->isVertical()) {
            qdc.y[position] = static_cast<sndScifiHit *>(sf_hits->At(i))->GetSignal(0);
        }
        else {
            qdc.x[position] = static_cast<sndScifiHit *>(sf_hits->At(i))->GetSignal(0);
        }
    }
}

const int SciFiPlaneView::getStation() const {
    return station;
}

const cfg SciFiPlaneView::getConfig() const {
    return config;
}
