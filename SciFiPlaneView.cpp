#include "SciFiPlaneView.h"
#include <stdexcept>



SciFiPlaneView::SciFiPlaneView(cfg c, TClonesArray *h, int b, int e,
                int s) : config(c), sf_hits(h), begin(b), end(e), station(s) {
    if (b > e) {
        throw std::runtime_error{"Begin index > end index"};
    }
    std::fill(qdc.x.begin(), qdc.x.end(), -999);
    std::fill(qdc.y.begin(), qdc.y.end(), -999);
    clusterBegin.x = -1;
    clusterBegin.y = -1;
    clusterEnd.x = -1;
    clusterEnd.y = -1;
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

void SciFiPlaneView::findCluster() {
    int currentStart{-1};
    int currentEnd{-1};
    int longestStart{-1};
    int longestEnd{-1};
    int currentLength{0};
    int maxLength{0};
    int consecutiveGaps{0};

    // Loop on x
    for (int i = 0; i < qdc.x.size(); ++i) {
        if (qdc.x[i] != -999) {
            // Inside a cluster
            if (currentStart == -1) {
                currentStart = i;
            }
            currentEnd = i;
            currentLength = currentEnd - currentStart + 1;

            // Update longest cluster if needed
            if (currentLength > maxLength) {
                maxLength = currentLength;
                longestStart = currentStart;
                longestEnd = currentEnd;
            }

            consecutiveGaps = 0; // Reset consecutive gaps counter
        } else {
            // Outside a cluster
            consecutiveGaps++;

            if (consecutiveGaps > config.SCIFIMAXGAP) {
                // Too many consecutive gaps, reset cluster start
                currentStart = -1;
                currentEnd = -1;
            }
        }
    }

    if (maxLength>=config.SCIFITHRESHOLD) {
        clusterBegin.x = longestStart;
        clusterEnd.x = longestEnd;
    }

    // Repeat for Y
    currentStart = -1;
    currentEnd = -1;
    longestStart = -1;
    longestEnd = -1;
    currentLength = 0;
    maxLength = 0;
    consecutiveGaps = 0;

    for (int i = 0; i < qdc.y.size(); ++i) {
        if (qdc.y[i] != -999) {
            // Inside a cluster
            if (currentStart == -1) {
                currentStart = i;
            }
            currentEnd = i;
            currentLength = currentEnd - currentStart + 1;

            // Update longest cluster if needed
            if (currentLength > maxLength) {
                maxLength = currentLength;
                longestStart = currentStart;
                longestEnd = currentEnd;
            }

            consecutiveGaps = 0; // Reset consecutive gaps counter
        } else {
            // Outside a cluster
            consecutiveGaps++;

            if (consecutiveGaps > config.SCIFIMAXGAP) {
                // Too many consecutive gaps, reset cluster start
                currentStart = -1;
                currentEnd = -1;
            }
        }
    }

    if (maxLength>=config.SCIFITHRESHOLD) {
        clusterBegin.y = longestStart;
        clusterEnd.y = longestEnd;
    }    
}

const int SciFiPlaneView::getStation() const {
    return station;
}

const cfg SciFiPlaneView::getConfig() const {
    return config;
}
