#include "SciFiPlaneView.h"
#include <stdexcept>

const double DEFAULT{-999.0};
const double FREQ{160.316E6};
const double TDC2ns = 1E9/FREQ;
const int NSIPM{8};
const int NSIDE{2};
const int TOFPETperBOARD{8};
const int TOFPETCHANNELS{64};


SciFiPlaneView::SciFiPlaneView(cfg c, TClonesArray *h, Scifi *det, int b, int e, int s) : 
            config(c), sf_hits(h), ScifiDet(det), begin(b), end(e), station(s), clusterBegin({-1,-1}), clusterEnd({-1,-1}), centroid({-1,-1}) 
                {
    if (b > e) {
        throw std::runtime_error{"Begin index > end index"};
    }
    qdc.x.resize(config.SCIFI_NCHANNELS, DEFAULT);
    qdc.y.resize(config.SCIFI_NCHANNELS, DEFAULT);
    hitTimestamps.x.resize(config.SCIFI_NCHANNELS, DEFAULT);
    hitTimestamps.y.resize(config.SCIFI_NCHANNELS, DEFAULT);
    geom.x.resize(config.SCIFI_NCHANNELS, DEFAULT);
    geom.y.resize(config.SCIFI_NCHANNELS, DEFAULT);
    depth.x.resize(config.SCIFI_NCHANNELS, DEFAULT);
    depth.y.resize(config.SCIFI_NCHANNELS, DEFAULT);
    fillQDC();
    fillTimestamps();
    fillGeometry();
}

const SciFiPlaneView::xy_pair<int> SciFiPlaneView::sizes() const{
    xy_pair<int> counts{0, 0};

    counts.x = std::count_if(qdc.x.begin(), qdc.x.end(), [] (double t) {return t > DEFAULT;});
    counts.y = std::count_if(qdc.y.begin(), qdc.y.end(), [] (double t) {return t > DEFAULT;});

    return counts;
}

void SciFiPlaneView::fillQDC() {
    for (int i{begin}; i<end; ++i) {
        //std::cout<<"Getchannel(0):\t"<<static_cast<sndScifiHit *>(sf_hits->At(i))->Getchannel(0)<<"\t GetSiPM():\t"<<static_cast<sndScifiHit *>(sf_hits->At(i))->GetSiPM()<<"\t GetSiPMChan()\t"<<static_cast<sndScifiHit *>(sf_hits->At(i))->GetSiPMChan()
        //<<"\t GetChannelID():\t"<<static_cast<sndScifiHit *>(sf_hits->At(i))->GetChannelID()<<"\t GetMat():\t"<<static_cast<sndScifiHit *>(sf_hits->At(i))->GetMat()<<std::endl;
        auto hit = static_cast<sndScifiHit *>(sf_hits->At(i));
        int mat = hit->GetMat();
        int sipm = hit->GetSiPM();
        int channel = hit->GetSiPMChan();
        int new_pos = channel + sipm*128 + mat*512;
        int position = 512*static_cast<sndScifiHit *>(sf_hits->At(i))->GetMat() + 64*static_cast<sndScifiHit *>(sf_hits->At(i))->GetTofpetID(0) + 63 - static_cast<sndScifiHit *>(sf_hits->At(i))->Getchannel(0);
        if (position != new_pos) {
            std::cout<<"mat:\t"<<mat<<"\t sipm:\t"<<sipm<<"\t channel:\t"<<channel<<"\t new:\t"<< new_pos <<"\t old:\t"<<position<<"\n";
        }
        
        if (static_cast<sndScifiHit *>(sf_hits->At(i))->isVertical()) {
            qdc.x[position] = static_cast<sndScifiHit *>(sf_hits->At(i))->GetSignal(0);
        }
        else {
            qdc.y[position] = static_cast<sndScifiHit *>(sf_hits->At(i))->GetSignal(0);
        }
    }
}

void SciFiPlaneView::fillTimestamps() {
    for (int i{begin}; i<end; ++i) {
        int position = 512*static_cast<sndScifiHit *>(sf_hits->At(i))->GetMat() + 64*static_cast<sndScifiHit *>(sf_hits->At(i))->GetTofpetID(0) + 63 - static_cast<sndScifiHit *>(sf_hits->At(i))->Getchannel(0);
        if (static_cast<sndScifiHit *>(sf_hits->At(i))->isVertical()) {
            hitTimestamps.x[position] = static_cast<sndScifiHit *>(sf_hits->At(i))->GetTime(0);
        }
        else {
            hitTimestamps.y[position] = static_cast<sndScifiHit *>(sf_hits->At(i))->GetTime(0);
        }
    }
}

// This function fills the geometry and depth arrays for a set of SciFi detector hits. 
// For each hit, it calculates the position based on the sensor (SiPM) channel and 
// determines if the hit is vertical or horizontal. Depending on the orientation, it 
// assigns the appropriate X or Y coordinates for the geometry, along with the depth (Z).

void SciFiPlaneView::fillGeometry() {
    TVector3 A, B;
    for (int i{begin}; i<end; ++i) {
        auto hit = static_cast<sndScifiHit *>(sf_hits->At(i));
        int position = 512*static_cast<sndScifiHit *>(sf_hits->At(i))->GetMat() + 64*static_cast<sndScifiHit *>(sf_hits->At(i))->GetTofpetID(0) + 63 - static_cast<sndScifiHit *>(sf_hits->At(i))->Getchannel(0);
        
        //hit->GetPosition(A,B);

        //int SiPMChan = hit->GetSiPMChan();
        int detectorID = hit->GetDetectorID();
        ScifiDet->GetSiPMPosition(detectorID, A, B);
        //ScifiDet->GetPosition(detectorID, A, B);     //error in Path of sndsw

        //std::cout<<"X= " <<A.X() <<", Y= " <<A.Y() <<std::endl;
        
        if (hit->isVertical()) {
            geom.x[position] = A.X();
            depth.x[position] = A.Z();
            // std::cout<<"A.X= " <<A.X() <<" B.X= " <<B.X() <<" Difference in X " <<B.Y()-A.Y() <<std::endl;  //check if the lenght of SciFi planes is 13cm
        }
        else {
            geom.y[position] = A.Y();
            depth.y[position] = A.Z();
            // std::cout<<"A.Y= " <<A.Y() <<" B.Y= " <<B.Y() <<" Difference in Y " <<B.X()-A.X() <<std::endl;  //check if the lenght of SciFi planes is 13cm
        }
    }
}

void SciFiPlaneView::findCluster() {

    if (station == 1 && sizes().x == 1 && sizes().y == 1) {return;}
    auto clusterize = [&] (std::vector<double> &qdcarr, std::vector<double> &timearr, int &clusterB, int &clusterE) {
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

                if (consecutiveGaps > config.SCIFI_GAPCLUSTER) {
                    // Too many consecutive gaps, reset cluster start
                    currentStart = -1;
                    currentEnd = -1;
                }
            }
        }

        if (maxLength>=config.SCIFI_DIMCLUSTER) {
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
    auto clusterize = [&] (std::vector<double> &qdcarr, std::vector<double> &timearr, int &clusterB, int &clusterE) {
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

                if (consecutiveGaps > config.SCIFI_GAPCLUSTER) {
                    // Too many consecutive gaps, reset cluster start
                    currentStart = -1;
                    currentEnd = -1;
                }
            }
        }


        if (maxLength>=config.SCIFI_DIMCLUSTER) {
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
    auto density = [&] (std::vector<double> &qdcarr) {
        for (int i{0}; i < config.SCIFI_NCHANNELS-window+1; ++i) {
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

// Initially used for calculating the centroid of the shower, then it was modified.

// void SciFiPlaneView::findCentroid(int windowSize) {
//   // find shower centroid position in cm   
//   auto slide = [&] (std::vector<double> &qdcarr, double &centroidCoordinate) {
//     double maxSignal{0};    
//     for (int i{0}; i < (config.SCIFI_NCHANNELS - windowSize); ++i) {      
//       int off = std::count(qdcarr.begin() + i, qdcarr.begin() + i + windowSize, DEFAULT);
//       if (off >= static_cast<int>(2*windowSize/3)) continue;

//       double signalSum = std::accumulate(qdcarr.begin() + i, qdcarr.begin() + i + windowSize, 0, 
//         [](int current_sum, int value) {
//             return (value > DEFAULT) ? (current_sum + value) : current_sum;
//         });

//       double ratio =  signalSum/(windowSize-off);
//       if ( maxSignal < ratio ) {
//         maxSignal = ratio;
//         centroidCoordinate = (i+(windowSize*.5)) *.025;    // conversion in cm
//       }
//     }
//   };
//   slide(qdc.x, centroid.x);
//   slide(qdc.y, centroid.y);
// }

// This function calculates the centroid position and depth for hits in both the X and Y directions
// by using QDC signal weights. It sums up the weighted positions and depths and computes the 
// centroids accordingly, avoiding negative signals. The results are stored in `centroid` and `centroid_depth`.

void SciFiPlaneView::findCentroid() {
    auto calculateCentroid = [](const std::vector<double>& signals, const std::vector<double>& positions, const std::vector<double>& depths, double& centroid, double& centroid_depth) {
        double sumSignal = 0.0;
        double sumWeighted = 0.0;
        double sumZ = 0.0;

        for (size_t position = 0; position < signals.size(); ++position) {
            double signal = signals[position];

            if (signal >= 0) { // Process valid signals only (non-negative)
                double pos = positions[position];
                double d = depths[position];
                // Accumulate signal, weighted position, and depth
                // std::cout <<"pos: " <<pos <<" signal: " <<signal <<std::endl;
                sumSignal += signal;
                sumWeighted += signal * pos;
                sumZ += signal * d;
            }
        }
        if (signals.size()<1) { // If no signals, assign default values
            centroid = DEFAULT;
            centroid_depth = DEFAULT;
        }
        else { // Otherwise, calculate the centroid and depth
            centroid = sumWeighted / sumSignal;
            centroid_depth = sumZ / sumSignal;
        }
    };

    calculateCentroid(qdc.x, geom.x, depth.x, centroid.x, centroid_depth.x);
    calculateCentroid(qdc.y, geom.y, depth.y, centroid.y, centroid_depth.y);

    //std::cout<<centroid.x <<", " <<centroid.y <<std::endl;
    //std::cout<<centroid_depth.x <<", " <<centroid_depth.y <<std::endl;
}

// Similar analysis performed for muons (MIP), QDC was omitted as it averages to zero overall.

void SciFiPlaneView::findMuon() {
    auto positionMuon = [](const std::vector<double>& positions, const std::vector<double>& depths, double& muonpos, double& muondepth) {
        double sumPositions = 0.0;
        double sumDepths = 0.0;
        size_t count = 0;

        for (size_t i = 0; i < positions.size(); ++i) {
            if (positions[i] > DEFAULT) {
                sumPositions += positions[i];
                sumDepths += depths[i];
                ++count;
            }
        }

        if (count == 0) {
            muonpos = DEFAULT;
            muondepth = DEFAULT;
        } else {
            muonpos = sumPositions / count;
            muondepth = sumDepths / count;
        }
    };

    positionMuon(geom.x, depth.x, muonpos.x, muondepth.x);
    positionMuon(geom.y, depth.y, muonpos.y, muondepth.y);

    // std::cout << muonpos.x << ", " << muonpos.y << std::endl;
    // std::cout << muondepth.x << ", " << muondepth.y << std::endl;
}

// Retrieves the geometry positions for x and y, then prints them to the console.

void SciFiPlaneView::findPosition() {
    auto positions = getGeometry();
    
    std::cout << "Positions x: ";
    for (const auto& p : positions.x) std::cout << p << " ";
    std::cout << std::endl;

    std::cout << "Positions y: ";
    for (const auto& p : positions.y) std::cout << p << " ";
    std::cout << std::endl;
}

void SciFiPlaneView::resetHit( bool isVertical, int index){
    if (isVertical) {
        qdc.y[index] = DEFAULT;
        hitTimestamps.y[index] = DEFAULT;
        geom.y[index] = DEFAULT;
    }
    else {
        qdc.x[index] = DEFAULT;
        hitTimestamps.x[index] = DEFAULT;
        geom.x[index] = DEFAULT;
    }
}

void SciFiPlaneView::timeCut (double minTime, double maxTime) {
    if (maxTime < minTime) {
        throw std::runtime_error{"maxTime < minTime"};
    }
    for (int i{0}; i < config.SCIFI_NCHANNELS; ++i) {
        if (hitTimestamps.x[i] < minTime || hitTimestamps.x[i] > maxTime) resetHit(false, i); 
        if (hitTimestamps.y[i] < minTime || hitTimestamps.y[i] > maxTime) resetHit(true, i);  
    }
}

bool SciFiPlaneView::evaluateNeighboringHits(int clustermaxsize, int max_miss) const {
    auto slide = [&](const std::vector<double> &qdcarr, int n) {
        if (n > clustermaxsize || n == 0) {
            return false;
        }

        int hits{0};
        int miss{0};

        for (int i{0}; i < qdcarr.size(); ++i) {
            if (qdcarr[i] == DEFAULT) {
                if (hits > 0) {
                    miss++;
                    if (miss > max_miss) {
                        return false;
                    }
                }
            }
            else {
                hits++;
                if (hits + miss > clustermaxsize) {
                    return false;
                }
                if (n == hits) {
                    return true;
                }
            }
        }
        return true;
    };

    int sizeX = sizes().x;
    int sizeY = sizes().y;

    return slide(qdc.x, sizeX) && slide(qdc.y, sizeY);
}

std::vector<int> SciFiPlaneView::calculateValidPositionsx(int clustermaxsize, int max_miss) const {
    std::vector<int> posx;
    
    int i=0;
    for (int value : qdc.x) {
        if (SciFiPlaneView::evaluateNeighboringHits(clustermaxsize, max_miss) == true) {
            if (value != DEFAULT) { 
                posx.push_back(i);
            }
        }
        i++;
    }

    return posx;
}

std::vector<int> SciFiPlaneView::calculateValidPositionsy(int clustermaxsize, int max_miss) const {
    std::vector<int> posy;
    
    int i=0;
    for (int value : qdc.y) {
        if (SciFiPlaneView::evaluateNeighboringHits(clustermaxsize, max_miss) == true) {
            if (value != DEFAULT) { 
                posy.push_back(i);
            }
        }
        i++;
    }

    return posy;
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

// Returns the centroid (x, y) position
const SciFiPlaneView::xy_pair<double> SciFiPlaneView::getCentroid() const {
    xy_pair<double> c{centroid.x, centroid.y};
    return c;
}

// Returns the centroid depth (x, y)
const SciFiPlaneView::xy_pair<double> SciFiPlaneView::getCentroidDepth() const {
    xy_pair<double> c{centroid_depth.x, centroid_depth.y};
    return c;
}

// Returns the muon position (x, y) 
const SciFiPlaneView::xy_pair<double> SciFiPlaneView::getMuon() const {
    xy_pair<double> c{muonpos.x, muonpos.y};
    return c;
}

// Returns the muon depth (x, y) 
const SciFiPlaneView::xy_pair<double> SciFiPlaneView::getMuonDepth() const {
    xy_pair<double> c{muondepth.x, muondepth.y};
    return c;
}

// Returns the geometry (x, y) as a pair of vectors of doubles.
const SciFiPlaneView::xy_pair<std::vector<double>> SciFiPlaneView::getGeometry() const {
    return geom;
}

const SciFiPlaneView::xy_pair<std::vector<double>> SciFiPlaneView::getTime() const{
    xy_pair<std::vector<double>> time{hitTimestamps.x, hitTimestamps.y};
    return time;
}

const SciFiPlaneView::xy_pair<std::vector<double>> SciFiPlaneView::getQDC() const{
    xy_pair<std::vector<double>> charge{qdc.x, qdc.y};
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

