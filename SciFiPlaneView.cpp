#include "SciFiPlaneView.h"
#include <stdexcept>

const double DEFAULT = -999.0;
const double FREQ{160.316E6};
const double TDC2ns = 1E9/FREQ;
const int NSIPM{8};
const int NSIDE{2};
const int NsidesNch{16};
const int TOFPETperBOARD{8};
const int TOFPETCHANNELS{64};


SciFiPlaneView::SciFiPlaneView(cfg c, TClonesArray *h, int b, int e,
                int s) : config(c), sf_hits(h), begin(b), end(e), station(s) {
    if (b > e) {
        throw std::runtime_error{"Begin index > end index"};
    }
    std::fill(qdc.x.begin(), qdc.x.end(), DEFAULT);
    std::fill(qdc.y.begin(), qdc.y.end(), DEFAULT);
    std::fill(hitTimestamps.x.begin(), hitTimestamps.x.end(), DEFAULT);
    std::fill(hitTimestamps.y.begin(), hitTimestamps.y.end(), DEFAULT);
    clusterBegin.x = -1;
    clusterBegin.y = -1;
    clusterEnd.x = -1;
    clusterEnd.y = -1;
    centroid.x = -1;
    centroid.y = -1;
    fillQDC();
    fillTimestamps();
}


auto SciFiPlaneView::sizes() const{
    xy_pair<int> counts{0, 0};

    counts.x = std::count_if(qdc.x.begin(), qdc.x.end(), [] (double t) {return t > DEFAULT;});
    counts.y = std::count_if(qdc.y.begin(), qdc.y.end(), [] (double t) {return t > DEFAULT;});

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

void SciFiPlaneView::fillTimestamps() {
    for (int i{begin}; i<end; ++i) {
        int position = 64*static_cast<sndScifiHit *>(sf_hits->At(i))->GetTofpetID(0) + 63 - static_cast<sndScifiHit *>(sf_hits->At(i))->Getchannel(0);
        if (static_cast<sndScifiHit *>(sf_hits->At(i))->isVertical()) {
            hitTimestamps.y[position] = static_cast<sndScifiHit *>(sf_hits->At(i))->GetTime(0);
        }
        else {
            hitTimestamps.x[position] = static_cast<sndScifiHit *>(sf_hits->At(i))->GetTime(0);
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
        if (qdc.x[i] != DEFAULT) {
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
        if (qdc.y[i] != DEFAULT) {
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

    // Set values outside cluster to default
    if (clusterBegin.x != -1 && clusterBegin.y != -1) {
        std::fill(hitTimestamps.x.begin(), hitTimestamps.x.begin() + clusterBegin.x, DEFAULT);   
        std::fill(hitTimestamps.x.begin() + clusterEnd.x + 1, hitTimestamps.x.end(), DEFAULT);
        std::fill(hitTimestamps.y.begin(), hitTimestamps.y.begin() + clusterBegin.y, DEFAULT);
        std::fill(hitTimestamps.y.begin() + clusterEnd.y + 1, hitTimestamps.y.end(), DEFAULT);

        std::fill(qdc.x.begin(), qdc.x.begin() + clusterBegin.x, DEFAULT);   
        std::fill(qdc.x.begin() + clusterEnd.x + 1, qdc.x.end(), DEFAULT);
        std::fill(qdc.y.begin(), qdc.y.begin() + clusterBegin.y, DEFAULT);
        std::fill(qdc.y.begin() + clusterEnd.y + 1, qdc.y.end(), DEFAULT);
    }
}

void SciFiPlaneView::findCentroid(int windowSize) {
  std::array<double, 2> centroidCoordinates;
  // find shower centroid position in cm   
  double maxSignal{0};
  for (int index{0}; index <2; ++index) {
    for (int i{0}; i < (config.BOARDPERSTATION*TOFPETperBOARD*TOFPETCHANNELS -windowSize); ++i) {
  
      double signalSum{0};

      for (int j{i}; j < (i+windowSize); ++j) {
        int den = windowSize;
        double signal;
        if (index == 0) signal = qdc.x[j];
        else signal = qdc.y[j];
        
        if ( signal != -999 ){
          signalSum += signal;
        } else {
          den -=1;
        }

        if (j == i+windowSize-1) {
          if (den < 2) continue;
          double ratio =  signalSum/den;
          if ( maxSignal < ratio ) {
            maxSignal = ratio;
            centroidCoordinates[index] = (i+(windowSize*.5)) *.025;    // conversion in cm
          }
        }
      }
    }
  }
  centroid.x = centroidCoordinates[0];
  centroid.y = centroidCoordinates[1];
}


void SciFiPlaneView::resetHit( bool isVertical, int index){
    if (isVertical) {
        qdc.y[index] = DEFAULT;
        hitTimestamps.y[index] = DEFAULT;
    }
    else {
        qdc.x[index] = DEFAULT;
        hitTimestamps.x[index] = DEFAULT;
    }
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

const SciFiPlaneView::xy_pair<double> SciFiPlaneView::getCentroid() const{
    xy_pair<double> c{centroid.x, centroid.y};
    return c;
}

const SciFiPlaneView::xy_pair<std::array<double, 512>> SciFiPlaneView::getTime() const{
    xy_pair<std::array<double, 512>> time{hitTimestamps.x, hitTimestamps.y};
    return time;
}

const SciFiPlaneView::xy_pair<std::array<double, 512>> SciFiPlaneView::getQDC() const{
    xy_pair<std::array<double, 512>> charge{qdc.x, qdc.y};
    return charge;
}