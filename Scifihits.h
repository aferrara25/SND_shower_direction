#ifndef Scifihits_h
#define Scifihits_h

#include <vector>

class Scifihits{
    private:
    int size{-1};
        
    std::vector<int> xHits;
    std::vector<int> yHits;

    int xSize{0};
    int ySize{0};

    public:
    bool isShower{false};    


    Scifihits(){};
    Scifihits(int dim);

    void addHit(int tofpet, int channel, bool isVertical);

    bool checkShower(int thr, int ngaps, int sidecut);
    bool checkMultipleHits();

    int getYhits();
    int getXhits();
    
};
#endif