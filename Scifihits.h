#ifndef Scifihits_h
#define Scifihits_h

#include <vector>

class Scifihits{
    private:
    int size{-1};
        


    int xSize{0};
    int ySize{0};

    public:
    std::vector<int> xHits;
    std::vector<int> yHits;

    bool isShower{false};    


    Scifihits(){};
    Scifihits(int dim);

    void addHit(int tofpet, int channel, bool isVertical);

    bool checkShower(int thr, int ngaps, int sidecut);
    bool checkMultipleHits();

    const int getYhits() const;
    const int getXhits() const;

    int getSize(){return size;};
    
};

inline std::ostream& operator<<(std::ostream& os, Scifihits const hit) {
  os << "This plane has " << hit.getXhits() << " hits in the X module and " << hit.getYhits() << " in the Y module" << std::endl;
  return os;
}
#endif