#include "USPlaneView.h"
#include <stdexcept>

/*
const double DEFAULT = -999.0;
const double FREQ{160.316E6};
const double TDC2ns = 1E9/FREQ;
const int NSIPM{8};
const int NSIDE{2};
const int NsidesNch{16};
const int TOFPETperBOARD{8};
const int TOFPETCHANNELS{64};
*/

USPlaneView::USPlaneView(cfg c, TClonesArray *h, int b, int e, int s) : 
            config(c), mufi_hits(h), begin(b), end(e), station(s) 
                {
    if (b > e) {
        throw std::runtime_error{"Begin index > end index"};
    }
    std::fill(qdc.begin(), qdc.end(), DEFAULT);
    std::fill(hitTimestamps.begin(), hitTimestamps.end(), DEFAULT);
    fillQDC();
    fillTimestamps();
}


const USPlaneView::sl_pair<int> USPlaneView::sizes() const{
    sl_pair<int> counts{0, 0};
    for (int i{0}; i<NCHANNELS; ++i) {
        if (qdc[i]>DEFAULT) {
            if (i%8==2 || i%8==5) {
                counts.s++;
            }
            else {
                counts.l++;
            }
        }
    }
    return counts;
}

void USPlaneView::fillQDC() {
    for (int j{begin}; j<end; ++j) {
        auto hit = static_cast<MuFilterHit *>(mufi_hits->At(j)); 
        int bar = static_cast<int>(hit->GetDetectorID()%1000);
        for (int i{0}; i<16; ++i) {
            if (!hit->isMasked(i) && hit->GetSignal(i) > -999) {
                if (16*bar + i>159) {
                    std::cout<<"Out of range\t"<<bar<<"\t"<<i<<"\n";
                    continue;
                }
                qdc[16*bar + i] = hit->GetSignal(i);
                //cout<<i<<"\t"<<hit->isShort(i)<<"\t"<<(i%8==2 || i%8==5)<<"\n";
            }
        }
    }
}

void USPlaneView::fillTimestamps() {
    for (int j{begin}; j<end; ++j) {
        auto hit = static_cast<MuFilterHit *>(mufi_hits->At(j)); 
        int bar = static_cast<int>(hit->GetDetectorID()%1000);
        for (int i{0}; i<16; ++i) {
            if (!hit->isMasked(i) && hit->GetSignal(i) > -999) {
                if (16*bar + i>159) {
                    continue;
                }
                hitTimestamps[16*bar + i] = hit->GetTime(i);
            }
        }
    }
}

void USPlaneView::resetHit(int index){
    qdc[index] = DEFAULT;
    hitTimestamps[index] = DEFAULT;

}

void USPlaneView::timeCut (double referenceTime) {
    for (int i{0}; i < NCHANNELS; ++i) {
        if ( (hitTimestamps[i] - referenceTime) > config.MUTIMECUT) resetHit(i);   
    }
}

const int USPlaneView::getStation() const {
    return station;
}

const cfg USPlaneView::getConfig() const {
    return config;
}

const int USPlaneView::getBegin() const {
    return begin;
}

const int USPlaneView::getEnd() const {
    return end;
}

const std::array<double, NCHANNELS> USPlaneView::getTime() const{
    return hitTimestamps;
}

const std::array<double, NCHANNELS> USPlaneView::getQDC() const{
    return qdc;
}

const USPlaneView::sl_pair<double> USPlaneView::getTotQDC() const{
    sl_pair<double> totQDC{0, 0};
    for (int i{0}; i<NCHANNELS; ++i) {
        if (qdc[i]>DEFAULT) {
            if (i%8==2 || i%8==5) {
                totQDC.s += qdc[i];
            }
            else {
                totQDC.l += qdc[i];
            }
        }
    }
    return totQDC;
}