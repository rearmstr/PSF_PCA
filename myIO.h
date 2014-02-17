#ifndef myIO_H
#define myIO_H
#include "myTypeDef.h"
#include <string>
#include <vector>
#include <fstream>
#include <CCfits/CCfits>
using namespace CCfits;

template<class T>
void writeMatrix(const T & mat,std::string fileName)
{                  
  std::ofstream myfile(fileName.c_str());                     

  myfile.precision(7);

  for (size_t i=0; i < mat.nrows(); i++) { 
    for (size_t j=0; j < mat.ncols(); j++) {
      myfile << mat(i,j) << "  "; 
    }
    myfile<<std::endl;
  }
  
}

template<class T>
void writeVector(const T & vec,std::string fileName)
{                  
  std::ofstream myfile(fileName.c_str());                     

  myfile.precision(7);

  for (size_t i=0; i < vec.size(); i++) { 
      myfile << vec(i) << "  "; 
  }
  myfile<<std::endl;
  
}

template<class T>
ExtHDU* writeMatrixToFits(FITS* file,const T &mat,string name)
{
  std::vector<long> ax(2,0);
  ax[0]=mat.nrows();
  ax[1]=mat.ncols();
  int n=mat.nrows()*mat.ncols();

  std::valarray<double> data(n);
  std::copy(mat.cptr(),mat.cptr()+n,&data[0]);
  ExtHDU* ext = file->addImage(name,DOUBLE_IMG,ax);
  long fpixel(1);
  ext->write(fpixel,n,data);
  return ext;
}

template<class T>
ExtHDU* writeVVectorToFits(FITS* file,const std::vector<std::vector<T> > &vec,
			   string name)
{
  std::vector<long> ax(2,0);
  
  ax[0]=vec.size();
  ax[1]=vec[0].size();
  int n=ax[0]*ax[1];

  std::valarray<float> data(n);
  int c=0;
  for(int i=0;i<ax[0];++i) {
    for(int j=0;j<ax[1];++j) {

      int in=ax[0]*j+i;
      data[in]=int(vec[i][j]);
    }
  }
  
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

  std::valarray<double> data(n);
  std::copy(vec.cptr(),vec.cptr()+n,&data[0]);
  ExtHDU* ext = file->addImage(name,DOUBLE_IMG,ax);
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

  valarray <double > data;
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

  valarray <double > data;
  e.read(data);

  std::copy(&data[0],&data[0]+nrows*ncols,mat.ptr());
}

#endif
