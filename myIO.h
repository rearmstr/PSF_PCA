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


ExtHDU* writeToFits(FITS* file,const DMatrix &mat,string name);

#endif
