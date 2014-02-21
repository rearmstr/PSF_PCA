#ifndef DETECTION_H
#define DETECTION_H 
#include "dbg.h"
#include "UseTMV.h"
#include "Bounds.h"
#include <valarray>
#include <map>
#include "Utils.h"


// This class represents a single detection.
// It inherits from the DVector class

namespace PCA {
  template<class T=float>
  class Detection :public tmv::Vector<T>
  {
    
  public:
    
    Detection(float x,float y,int nobs) :
      tmv::Vector<T>(nobs,defaultVal),pos(x,y),nvar(nobs),clip(false) {}
  
    int getNVar() const {return nvar;}
    Position<float> getPos() const {return pos;}
    float getX() const {return pos.x;}
    float getY() const {return pos.y;}

    void setClip(bool val) {clip=val;}
    bool isClipped() const {return clip;}
    
    // This is to write to fits files using ccfits
    std::valarray<T> const getValarray() {
      std::valarray<T> v(nvar);
      for(int i=0;i<nvar;++i) v[i]=(*this)(i);
      return v;
    }
    
    
  protected:
    int nvar;
    Position<float> pos;
    bool clip;
  };



  
template<class T=float>
class Chip {

public:
  // assume chip starts at zero,zero
  Chip(int _label,float xmin,float xmax,float ymin,float ymax): 
    label(_label),bounds(ymin,xmax,ymin,ymax) {}

  Chip(int _label,Bounds<float> _bounds): 
    label(_label),bounds(_bounds) {}

  ~Chip() {
    DetIter it=DetBegin();
    for(; it!=DetEnd();++it) {
      delete *it;
    }
    det.clear();
  }

  void addDet(Detection<T>* _det) {det.push_back(_det);}
  int getLabel() {return label;}
  int getNDet() {return det.size();}

  typedef typename std::vector<Detection<T>*>::iterator DetIter;
  DetIter DetBegin() {return det.begin();}
  DetIter DetEnd() {return det.end();}
  std::vector<Detection<T>* > det;
  Bounds<float> getBounds() const {return bounds;}

private:
  int label;
  Bounds<float> bounds;

  
};

  // Probably need to generalize this because there are shapelet
  // specific items in the class
template<class T=float>
class Exposure {

public:


  Exposure(std::string _label="",int nvar=0,
	   float x_min=0, float x_max=2047,
	   float y_min=0, float y_max=4095,
	   std::vector<int> skip=std::vector<int>(),
	   int _shapeStart=3,
	   bool _addSize=true,
	   bool _include_miss=true);
  Exposure(std::string _label="",int nvar=0,
	   Bounds<float> _bounds=Bounds<float>(0,2047,0,4095),
	   std::vector<int> skip=std::vector<int>(),
	   int _shapeStart=3,
	   bool _addSize=true,
	   bool _include_miss=true);

  ~Exposure() {
    ChipIter it=ChipBegin();
    for(; it!=ChipEnd();++it) {
      delete (*it).second;
    }
    chips.clear();
  }
  
  Chip<T>* getChip(int i) {return chips[i];}

  int getNChip() const {return chips.size();}
  void addChip(int ichip,Chip<T> *chip) { 
    Assert(chip->getBounds()==bounds);
    chips[ichip]=chip;
  }

  int nSkip() const {return skip.size();}
  bool readShapeletDir(std::string dir,int nchip,
		       std::string suffix="psf.fits",
		       std::string exposure="");

  bool readShapelet(std::string file,int _nvar=-1,int _nccd=-1);


  bool readPixels(std::string dir,int npix,int nvar,std::string sdir,
		  std::string exposure="");
  
  // This writes out all ccds to single file
  // maybe this will speed reading time
  void writeFits(std::string dir,std::string suffix="_det.fits");

  std::string getLabel() {return label;}
  bool isOutlier() {return outlier;}
  void setOutlier(bool val) {outlier=val;}
  int getNDet();

  typedef typename std::map<int,Chip<T>*>::iterator ChipIter;
  ChipIter ChipBegin() {return chips.begin();}
  ChipIter ChipEnd() {return chips.end();}

  void detectOutlier(float sigma);

private:
  std::map<int,Chip<T> *> chips; 
  std::vector<int> skip;
  Bounds<float> bounds;
  std::string label;
  bool outlier;
  int shapeStart;
  bool addSize;
  bool includeMiss;
  int nvar;
  int max_ccd;

};
  
    


}




#endif
