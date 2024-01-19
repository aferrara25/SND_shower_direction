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


SciFiPlaneView::SciFiPlaneView(cfg c, TClonesArray *h, int b, int e, int s) : 
            config(c), sf_hits(h), begin(b), end(e), station(s), clusterBegin({-1,-1}), clusterEnd({-1,-1}), centroid({-1,-1}) 
                {
    if (b > e) {
        throw std::runtime_error{"Begin index > end index"};
    }
    std::fill(qdc.x.begin(), qdc.x.end(), DEFAULT);
    std::fill(qdc.y.begin(), qdc.y.end(), DEFAULT);
    std::fill(hitTimestamps.x.begin(), hitTimestamps.x.end(), DEFAULT);
    std::fill(hitTimestamps.y.begin(), hitTimestamps.y.end(), DEFAULT);
    fillQDC();
    fillTimestamps();
}


const SciFiPlaneView::xy_pair<int> SciFiPlaneView::sizes() const{
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

    if (station == 1 && sizes().x == 1 && sizes().y == 1) {return;}
    auto clusterize = [&] (std::array<double, 512> &qdcarr, std::array<double, 512> &timearr, int &clusterB, int &clusterE) {
        int currentStart{-1};
        int currentEnd{-1};
        int longestStart{-1};
        int longestEnd{-1};
        int currentLength{0};
        int maxLength{0};
        int consecutiveGaps{0};
        // Loop on x
        for (int i = 0; i < qdcarr.size(); ++i) {
            if (qdcarr[i] != DEFAULT) {
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
            clusterB = longestStart;
            clusterE = longestEnd;
        }
        if (clusterB != -1) {
            std::fill(timearr.begin(), timearr.begin() + clusterB, DEFAULT);   
            std::fill(timearr.begin() + clusterE + 1, timearr.end(), DEFAULT);

            std::fill(qdcarr.begin(), qdcarr.begin() + clusterB, DEFAULT);   
            std::fill(qdcarr.begin() + clusterE + 1, qdcarr.end(), DEFAULT);
        }
        else {
            std::fill(timearr.begin(), timearr.end(), DEFAULT);   

            std::fill(qdcarr.begin(), qdcarr.end(), DEFAULT);   ;
        }        
    };
    clusterize(qdc.x, hitTimestamps.x, clusterBegin.x, clusterEnd.x);
    clusterize(qdc.y, hitTimestamps.y, clusterBegin.y, clusterEnd.y);
}

bool SciFiPlaneView::infoCluster() {

    if (station == 1 && sizes().x == 1 && sizes().y == 1) {return false;}
    auto clusterize = [&] (std::array<double, 512> &qdcarr, std::array<double, 512> &timearr, int &clusterB, int &clusterE) {
        int currentStart{-1};
        int currentEnd{-1};
        int longestStart{-1};
        int longestEnd{-1};
        int currentLength{0};
        int maxLength{0};
        int consecutiveGaps{0};
        // Loop on x
        for (int i = 0; i < qdcarr.size(); ++i) {
            if (qdcarr[i] != DEFAULT) {
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
        clusterB = longestStart;
        clusterE = longestEnd;
	    return true;
        }
        else return false;
       
    };
    if (clusterize(qdc.x, hitTimestamps.x, clusterBegin.x, clusterEnd.x) && clusterize(qdc.y, hitTimestamps.y, clusterBegin.y, clusterEnd.y)) {
	return true;
    }
    else {
	return false;
    }
}

bool SciFiPlaneView::infoDensity(int window, int min_hits) {

    if (station == 1 && sizes().x == 1 && sizes().y == 1) {return false;}
    if (min_hits>window) {throw std::runtime_error{"min_hits > radius"};}
    auto density = [&] (std::array<double, 512> &qdcarr) {
        for (int i{0}; i < 512-window+1; ++i) {
            if (std::count_if(qdcarr.begin()+i, qdcarr.begin()+i+window, [] (double t) {return t > DEFAULT;}) >= min_hits) {
                return true;
            }
        }
        return false;
    };
    if (density(qdc.x) && density(qdc.y)) {
        return true;
    }
    else {
        return false;
    }
}

void SciFiPlaneView::findCentroid(int windowSize) {
  // find shower centroid position in cm   
  auto slide = [&] (std::array<double, 512> &qdcarr, double &centroidCoordinate) {
    double maxSignal{0};    
    for (int i{0}; i < (config.BOARDPERSTATION*TOFPETperBOARD*TOFPETCHANNELS -windowSize); ++i) {      
      int off = std::count(qdcarr.begin() + i, qdcarr.begin() + i + windowSize, DEFAULT);
      if (off >= static_cast<int>(2*windowSize/3)) continue;

      double signalSum = std::accumulate(qdcarr.begin() + i, qdcarr.begin() + i + windowSize, 0, 
        [](int current_sum, int value) {
            return (value > DEFAULT) ? (current_sum + value) : current_sum;
        });

      double ratio =  signalSum/(windowSize-off);
      if ( maxSignal < ratio ) {
        maxSignal = ratio;
        centroidCoordinate = (i+(windowSize*.5)) *.025;    // conversion in cm
      }
    }
  };
  slide(qdc.x, centroid.x);
  slide(qdc.y, centroid.y);
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

void SciFiPlaneView::timeCut (double referenceTime) {

    if (station == 1) return;
    else {
        for (int i{0}; i < 512; ++i) {
        if ( (hitTimestamps.x[i] - referenceTime) > config.TIMECUT) resetHit(false, i); 
        if ( (hitTimestamps.y[i] - referenceTime) > config.TIMECUT) resetHit(true, i);  
        }
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

const SciFiPlaneView::xy_pair<double> SciFiPlaneView::getTotQDC() const{
    xy_pair<double> qdcSum{0,0};
    qdcSum.x = std::accumulate(qdc.x.begin(), qdc.x.end(), 0, [](int current_sum, int value) {return (value > DEFAULT) ? (current_sum + value) : current_sum;});
    qdcSum.y = std::accumulate(qdc.y.begin(), qdc.y.end(), 0, [](int current_sum, int value) {return (value > DEFAULT) ? (current_sum + value) : current_sum;});
    
    return qdcSum;
}

  const SciFiPlaneView::xy_pair<int> SciFiPlaneView::getClusterSize() const{
    xy_pair<int> clSize{-1, -1};
    clSize.x = clusterEnd.x - clusterBegin.x;
    clSize.y = clusterEnd.y - clusterBegin.x;

    return clSize;
  }
