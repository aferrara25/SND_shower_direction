#include "Scifihits.h"

#include <algorithm>
#include <iostream>

Scifihits::Scifihits(int dim): size{dim} {
    
    for (int i = 0; i < size; ++i) {
      xHits.emplace_back(0);
      yHits.emplace_back(0);
    }
    //double position = (64*tpID+63-ch)*0.025;
    
}

void Scifihits::addHit(int tofpet, int channel, bool isVertical) {
    if (isVertical) {
        yHits[64*tofpet+63-channel] +=1;
        ySize +=1;
    }
    else {
        xHits[64*tofpet+63-channel] +=1;
        xSize +=1;
    }
}

bool Scifihits::checkShower(int thr, int ngaps, int sidecut){
    // if there's no hits there's for sure no shower
    if (std::all_of(yHits.begin(), yHits.end(), [](int i) { return i==0; }) || std::all_of(xHits.begin(), xHits.end(), [](int i) { return i==0; })) return false;
    //std::cout << " X hits: " << xSize << " Y hits: " << ySize << std::endl;

    int xCounter{0};
    int xMax{0};
    int xGaps{0};

    int yCounter{0};
    int yMax{0};
    int yGaps{0};

    for (int index = sidecut; index < (size-sidecut); ++index ) {
        //VERTICAL
        //if the bin is full add to counter else add to gaps
        if (yHits[index]>0 ) yCounter +=1;
        else yGaps +=1;

        
        // if the limit of ngaps is reached reset the counter
        if (yGaps > ngaps ) {
            //save the max number of closeby channels on I find
            if (yCounter > yMax) yMax = yCounter;
            //reset counter
            yCounter = 0;
            yGaps = 0;
        }

        //HORIZONTAL
        //if the bin is full add to counter else add to gaps
        if (xHits[index]>0 ) xCounter +=1;
        else xGaps +=1;
        // if the limit of ngaps is reached reset the counter
        if (xGaps > ngaps ) {
            //save the max number of closeby channels on I find
            if (xCounter > xMax) xMax = xCounter;
            //reset counter
            xCounter = 0;
            xGaps = 0;
        }
    }
    // in case counter never > gaps 
    if (yCounter > yMax) yMax = yCounter;
    if (xCounter > xMax) xMax = xCounter;

    if (xMax > thr && yMax > thr) {
        isShower = true;
        return true;
    }
    else return false;

}

bool Scifihits::checkMultipleHits(){

    int xCounter = std::count_if(xHits.begin(), xHits.end(), [&](int val){return val > 1;});
    int yCounter = std::count_if(yHits.begin(), yHits.end(), [&](int val){return val > 1;});

    if (xCounter > 0 || yCounter > 0) return true;
    return false;
}

const int Scifihits::getYhits() const { return ySize; } 

const int Scifihits::getXhits() const { return xSize; }


