#ifndef PCA_OBJECTS_H
#define PCA_OBJECTS_H 1

#include "Bounds.h"
#include <vector>
#include <valarray>
#include <map>
#include <algorithm>
#include "TMV.h"
#include <cassert>
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

  template <class T>
  double median_mad(T& v,double &mad) {
    double med=median(v);
    T medr(v.size());
    for(int i=0;i<v.size();++i) medr[i]=std::abs(v[i]-med);
    mad=1.4826*median(medr);
    return med;
  }
  

  static const float defaultVal  = -9999.0;

template<class T=double>
class Detection {

 public:

  Detection(float x,float y,int nobs,double ra=-1,double dec=-1) :
    pos(x,y),vals(nobs,defaultVal),nval(nobs),clip(false) {}
  
  Position<float> getPos() {return pos;}
  Position<double> getSky() {return sky;}
  void setSky(Position<double> &_sky) {sky=_sky;}

  std::vector<T> getVals() {return vals;}
  std::valarray<T> getVVals() {
    std::valarray<T> v(vals.size());
    for(int i=0;i<vals.size();++i) v[i]=vals[i];
    return v;
  }
  float getVal(int i) {
    if(i>=0 && i<nval) return vals[i];
    return -9999.0;
  }
  void setVal(int i,T val) {
    if(i>=0 && i<nval) vals[i]=val;
  }
  // need to think how we will get these values for vector settings
  void setVals(std::vector<T> &_vals) {
    vals=_vals;
   }
  void setClip(bool val) {clip=val;}
  bool isClipped() {return clip;}
  
protected:
  int nval;
  Position<float> pos;
  Position<double> sky;
  std::vector<T> vals;
  bool clip;
};



// Cell contains potentially many detections, the ares is rectangular (for now)
// May want outlier rejection and estimation mean/median
// done in this class
template<class T=double>
class Cell {


public:

  Cell(int _nvar,float xmin, float xmax,float ymin, float ymax): nvar(_nvar),missing(false),
								 bounds(xmin,xmax,ymin,ymax),
								 clipped(false) {missing=false;}
  Cell(int _nvar,Bounds<float> b): nvar(_nvar),bounds(b),clipped(false),missing(false) {}
  void addDet(Detection<T>* _det) {dets.push_back(_det);}
  int getNDet() {return dets.size();}
  int getNClip();
  int getNGood();
  bool isClipped() {return clipped;}
  void setClipped(bool miss) {clipped=clipped;}
  bool isMissing() {return missing || getNGood()==0;}
 
  void setMissing(bool miss) {missing=miss;}

  std::vector<T> getVals(std::string type,std::vector<float> &params);

  std::vector<T> getMeanVals();
  std::vector<std::vector<T> > getDiff(tmv::ConstVectorView<T> &vals,std::string type,
				       std::vector<float> params,
				       const std::vector<double> &mean,
				       const std::vector<double> &sigma,
				       bool clip=false,
				       double nclip=-1);
                                     
  std::vector<std::valarray<T> >getDetVals(tmv::ConstVectorView<T> &vals,
					   std::string type,std::vector<float> params);
  std::vector<T> getMeanClipVals(float clip);
  std::vector<T> getMedianVals();
  std::vector<T> getFitVals(int order=1,float clip=5);
 

  int getNVal(std::string type,std::vector<float> &params);
  int getNVar() {return nvar;}

  Detection<T> * getDet(int i) {
    if (i>=0 && i<dets.size()) return dets[i];
    return 0;
  }

  enum ValType {Mean=1,MeanClip=2,Median=3,Fit=4};

  ValType getTypeFromString(std::string type) {
    if(type=="mean") return Mean;
    if(type=="mean_clip") return MeanClip;
    if(type=="median") return Median;
    if(type=="fit") return Fit;
    std::cout<<"Not a valid type: "<<type<<std::endl;
    assert(0);
    return Mean;
  };


private:
  int nvar;
  int fitorder;
  Bounds<float> bounds;
  std::vector<Detection<T>* > dets;
  bool clipped;
  bool missing;

};


template<class T=double>
class Chip {

public:
  // assume chip starts at zero,zero
  Chip(int _label,float xmax,float ymax): label(_label),bounds(0,xmax,0,ymax) {}
  void addDet(Detection<T>* det);
  void divide(int nvar, int _nx,int _ny); // setup the cell sizes
  std::vector<T> getVals(std::string type,std::vector<float> &params);
  std::vector<bool> getMissing();
  void setMissing(float prob);
  std::vector<Bounds<float> > getCellBounds() {return cbounds;}
  Cell<T>* getCell(int i) {return cells[i];}
  int getNCell() {return cells.size();}
  Cell<T>* operator[](int i) {return cells[i];}
  int getNClip();
  int getNGood();
  int getNDet();
  int getLabel() {return label;}
  //int getNVal(std::string type,std::vector<float> params);

private:
  int label;
  int nvar;
  Bounds<float> bounds;
  std::vector<Bounds<float> > cbounds;
  std::vector<Cell<T>*> cells;
  int nx;
  int ny;
  
};

template<class T=double>
class Exposure {

public:
  Exposure(std::string _label,int _nchip, int _shapestart=3,double ra=-9999.0,double dec=-9999.0,float airmass=-9999.0);
  
  void setChipDivide(int nx,int ny) {
    ny_chip=ny;
    nx_chip=nx;
  }
  void setChipMax(int xmax,int ymax) {
    xmax_chip=xmax; 
    ymax_chip=ymax;
  }
  
  Chip<T>* operator[](int i) {return chips[i];}
  Chip<T>* getChip(int i) {return chips[i];}
  int getCells() {return nx_chip*ny_chip*chips.size();}
  //int getNVal(std::string type,std::vector<float> params);
  int getNChip() {return nchip;}
  void setShapeStart(int start) {shapeStart=start;}
  void addSkip(int ichip) { skip.push_back(ichip);}
  void addChip(int ichip,Chip<T> *chip) { chips[ichip]=chip;}
  int nSkip() {return skip.size();}
  bool readShapelet(std::string dir,int nvar,bool add_size=false,bool 
		    include_miss=false,bool use_dash=false,std::string suffix="psf.fits",
		    std::string exposure="",float max=1.0,
		    std::string use_dir=".",std::string cdir="");
  bool readPixels(std::string dir,int npix,int nvar,std::string sdir, 
		  bool use_dash=false,std::string exposure="");
  tmv::Vector<T> getVals(std::string type,std::vector<float> &params);
  std::vector<bool> getMissing();
  void setMissing(float prob);
  std::string getLabel() {return label;}
  bool isOutlier() {return outlier;}
  void setOutlier(bool val) {outlier=val;}
  int getNClip();
  int getNGood();
  int getNDet();
  std::vector<double> outlierReject(const tmv::Vector<T> &data_r,
				    double sigma,std::string type,std::vector<float> param);

  std::map<int,Chip<T>*> chips; // do I want to put this in public?  way to iterate through chips
private:
  float xmax_chip;
  float ymax_chip;
  int nx_chip;
  int ny_chip;
  std::string label;

  std::vector<int> skip;
  int nchip;
  int nvar;
  int shapeStart;
  double ra;
  double dec;
  float airmass;
  bool outlier;

};

};
#endif
