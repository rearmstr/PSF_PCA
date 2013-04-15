#ifndef myIO_H
#define myIO_H

#include <string>
#include <vector>

void inputFromFile(std::string,c_ControlParam &,c_Data &);  // Xmat; all data in one file
                                                            // vector is a line in file
void readVector(int, std::string, DMatrix &, int);          // Xmat; read a single vector
void readVectorDefocus(int, std::string, DMatrix &, int);   // zTab; read a single vector

std::vector< std::vector<double> > readIn2dData(const char* filename);

void outputToFile(const DMatrix &, std::string);                  // output a double TMV matrix to file
                                                           // write every column to a line
void outputToFitsFile(c_ControlParam &,const DMatrix & matrix, 
		      std::string filename);                  // output a double TMV matrix to file
                                                           // write every column to a line
void outputToFileTMVint(tmv::Matrix<int>, std::string);   // output an int TMV matrix to file
void outputVector(int, DMatrix&, std::string);             // output a vector to file
void outputVectorAll(DMatrix&, std::string);               // all vectors, each to a file
void outputIntArr(int,std::vector<int>&,std::string);      // output an int array to file
void outputArr(int,std::vector<double>&,std::string);      // output a double array to file

std::string convertInt(int);                               // convert an integer to string
bool FileExists(std::string);                              // check file or dir existance

#endif
