#include "USPlaneView.h"
#include <stdexcept>

USPlaneView::USPlaneView(cfg c, TClonesArray *h, int b, int e, int s) : 
            config(c), mufi_hits(h), begin(b), end(e), station(s) 
                {
    if (b > e) {
        throw std::runtime_error{"Begin index > end index"};
    }
    qdc.resize(config.US_NCHANNELS, DEFAULT);
    hitTimestamps.resize(config.US_NCHANNELS, DEFAULT);
    fillQDC();
    fillTimestamps();
}


const USPlaneView::sl_pair<int> USPlaneView::sizes() const{
    sl_pair<int> counts{0, 0};
    for (int i{0}; i<config.US_NCHANNELS; ++i) {
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
            if (!hit->isMasked(i) && hit->GetSignal(i) > DEFAULT) {
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
    for (int i{0}; i < config.US_NCHANNELS; ++i) {
        if ( (hitTimestamps[i] - referenceTime) >= config.US_TIMECUT) resetHit(i);   
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

const std::vector<double> USPlaneView::getTime() const{
    return hitTimestamps;
}

const std::vector<double> USPlaneView::getQDC() const{
    return qdc;
}

const USPlaneView::sl_pair<double> USPlaneView::getTotQDC() const{
    sl_pair<double> totQDC{0, 0};
    for (int i{0}; i<config.US_NCHANNELS; ++i) {
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