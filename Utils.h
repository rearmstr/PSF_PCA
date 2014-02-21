#ifndef UTILS_H
#define UTILS_H
#include <string>
#include <vector>
#include <algorithm>
#include "Bounds.h"
#include "UseTMV.h"
#include <CCfits/CCfits>
using std::string;
using namespace CCfits;

namespace PCA {

  static const float defaultVal  = -9999.0;
  
  // Useful stats functions
  float percentile(std::vector<float>& v, float pct);
  
  float median(std::vector<float>& v);
  
  // compute median and median absolute deviation
  float median_mad(std::vector<float> & v,float &mad);

  // flag objects that are greater than n sigma from median
  // compute mad, median and mean of good objects
  float median_mad_flag(std::vector<float> & v,float &mad,
			 float &mean,float sigma,
			 std::vector<int> &ind);

  // split string into vector separated by the delimiters
  void Tokenize(const string& str,
		std::vector<string>& tokens,
		const string& delimiters);
  

  // for polynomiial fitting
  DVector definePXY(int order, float x, float xmin, float xmax);
  void setPRow(int fitorder, Position<float> pos, 
	       const Bounds<float>& bounds, DVectorView prow);
  
  
  void meanRemove(FMatrix &m,FVector &mean,
		     const std::vector<std::vector<bool> >  &missing);


  // methods to read/write TMV matrices/vectors to fits files
  
  template<class T>
  ExtHDU* writeMatrixToFits(FITS* file,const T &mat,string name)
  {
    std::vector<long> ax(2,0);
    ax[0]=mat.nrows();
    ax[1]=mat.ncols();
    int n=mat.nrows()*mat.ncols();

    std::valarray<float> data(n);
    std::copy(mat.cptr(),mat.cptr()+n,&data[0]);
    ExtHDU* ext = file->addImage(name,FLOAT_IMG,ax);
    long fpixel(1);
    ext->write(fpixel,n,data);
    return ext;
  }

  template< class T>
  ExtHDU* writeVectorToFits(FITS* file,const T &vec,string name)
  {
    std::vector<long> ax(2,0);
    ax[0]=vec.size();
    ax[1]=1;
    int n=vec.size();

    std::valarray<float> data(n);
    std::copy(vec.cptr(),vec.cptr()+n,&data[0]);
    ExtHDU* ext = file->addImage(name,FLOAT_IMG,ax);
    long fpixel(1);
    ext->write(fpixel,n,data);
    return ext;
  }

  template<class T>
  void readMatrixFromFits(FITS* file,string name,T &mat)
  {
    ExtHDU& e=file->extension(name);
    e.readAllKeys();
    int nrows=e.axis(0);
    int ncols=e.axis(1);
    mat.resize(nrows,ncols);

    valarray <float > data;
    e.read(data);

    std::copy(&data[0],&data[0]+nrows*ncols,mat.ptr());
  }

  template<class T>
  void readVecFromFits(FITS* file,string name,T &mat)
  {
    ExtHDU& e=file->extension(name);
    e.readAllKeys();
    int nrows=e.axis(0);
    int ncols=e.axis(1);
    mat.resize(nrows);

    valarray <float > data;
    e.read(data);

    std::copy(&data[0],&data[0]+nrows*ncols,mat.ptr());
  }

  
}








#endif
