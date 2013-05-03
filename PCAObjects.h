#ifndef PCA_OBJECTS_JH
#define PCA_OBJECTS_JH 1

#include "Bounds.h"
#include <vector>
#include <map>
#include <algorithm>
#include "TMV.h"
// Generic detection object that can hold any type of data
// for the PCA

namespace PCA {

  template <class T>
  T percentile(std::vector<T>& v, float pct) {
    if (v.empty()) return 0;
    int pctileIndex= static_cast<int> ((v.size()-1)*pct);
    if (pctileIndex>=v.size()) pctileIndex=v.size()-1;
    std::nth_element(v.begin(), v.begin()+pctileIndex, v.end());
    return v[pctileIndex];
  }
  
  template <class T>
  T median(std::vector<T>& v) {return percentile(v,0.5);}





class Detection {

 public:

  Detection(float x,float y,int nobs) :pos(x,y),vals(nobs,-9999.0),nval(nobs) {}
  
  Position<float> getPos() {return pos;}
  std::vector<float> getVals() {return vals;}
  float getVal(int i) {
    if(i>=0 && i<nval) return vals[i];
    return -9999.0;
  }
  void setVal(int i,float val) {
    if(i>=0 && i<nval) vals[i]=val;
  }
  // need to think how we will get these values for vector settings
  void setVals(std::vector<float> &_vals) {
    vals=_vals;
   }

  
protected:
  int nval;
  Position<float> pos;
  std::vector<float> vals;

};

// Cell contains potentially many detections, the ares is rectangular (for now)
// May want outlier rejection and estimation mean/median
// done in this class
class Cell {

public:

  Cell(int _nvar,float xmin, float xmax,float ymin, float ymax): nvar(_nvar),bounds(xmin,xmax,ymin,ymax) {}
  Cell(int _nvar,Bounds<float> b): nvar(_nvar),bounds(b) {}
  void addDet(Detection* _det) {dets.push_back(_det);}
  std::vector<float> getMeanVals();
  std::vector<float> getMedianVals();
  //tmv::Vector<float> getTMVVals();

private:
  int nvar;
  Bounds<float> bounds;
  std::vector<Detection*> dets;
};



class Chip {

public:
  // assume chip starts at zero
  Chip(int _nvar,int _label,float xmax,float ymax): nvar(_nvar),label(_label),bounds(0,xmax,0,ymax) {}
  void addDet(Detection* det);
  void divide(int _nx,int _ny); // setup the cell sizes
  std::vector<float> getMeanVals();
  std::vector<float> getMedianVals();
  std::vector<Bounds<float> > getCellBounds() {return cbounds;}

private:
  int nvar;
  int label;
  Bounds<float> bounds;
  std::vector<Bounds<float> > cbounds;
  std::vector<Cell> cells;
  int nx;
  int ny;
  
};


class Exposure {

public:
  Exposure(std::string _label,int _nchip, int nvar,double ra=-9999.0,double dec=-9999.0,float airmass=-9999.0);
  
  void setChipDivide(int nx,int ny) {
    ny_chip=ny;
    nx_chip=nx;
  }
  void setChipMax(int xmax,int ymax) {
    xmax_chip=xmax; 
    ymax_chip=ymax;
  }
  
  Chip* operator[](int i) {return chips[i];}
  
  void setShapeStart(int start) {shapeStart=start;}
  void addSkip(int ichip) { skip.push_back(ichip);}
  void addChip(int ichip,Chip *chip) { chips[ichip]=chip;}
  int nSkip() {return skip.size();}
  bool readShapelet(std::string dir,std::string exposure="");
  // std::vector<float> getMeanVals();
//   std::vector<float> getMedianVals();
  tmv::Vector<float> getMeanVals();
  tmv::Vector<float> getMedianVals();
private:
  float xmax_chip;
  float ymax_chip;
  int nx_chip;
  int ny_chip;
  std::string label;
  std::map<int,Chip*> chips;
  std::vector<int> skip;
  int nchip;
  int nvar;
  int shapeStart;
  double ra;
  double dec;
  float airmass;

};

};
#endif
