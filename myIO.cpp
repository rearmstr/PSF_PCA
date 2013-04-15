#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "TMV.h"
#include "TMV_Sym.h"

#include <vector>

#include <sys/stat.h>  // needed by FileExists()

#include "myClass.h"
#include "myTypeDef.h"

using namespace std;


/* --------------------------------------------------------------------- */
void inputFromFile (string fileName,c_ControlParam &contParam,c_Data &myData)
{
  int i,j;
  string word;                            // Xmat; input from file
  fstream myfile;

  myfile.open(fileName.c_str(), ios::in);
  if (myfile.is_open())
  {
    cout << "\t opened file " << fileName << " for input" << endl;
    i=0;
    j=0;
    while (! myfile.eof())
    {
      while ( getline (myfile, word, ' '))   // whitespace delimited
      {            // requires a whitespace at the end of every line
                   // except the last line
                   // can't have space at the beginning of the line
                   // can't have more than one space together
         istringstream ss( word );
         ss >> myData.Xmat(j,i);
         // cout << i << "\t" << j << "\t" << myData.Xmat(j,i) << "\n";
         j=j+1;
         if (j == contParam.nrows) {
           j=0;
           ++i;
         }
         if (i >= contParam.ncols) { break; }         // only read in ncols vectors
      }

    }
    myfile.close();
  }

  else
  {
    cout << "Unable to open file " << fileName << endl;
  }

}


/* --------------------------------------------------------------------- */
void readVector (int iCol, string fileName, DMatrix &Xmat, int nrows)
{
  fstream myfile;                          // Xmat; read a single vector
  int i;

  myfile.open(fileName.c_str(), ios::in);
  if (myfile.is_open())
  {
    i=0;
    while (! myfile.eof())
    {
       myfile >> Xmat(i,iCol);
       i++;
       if (i >= nrows) { break; }         // only read in nrows dimmension
    }
    myfile.close();
  }

  else { cout << "Unable to open file " << fileName << endl; }
}


/* --------------------------------------------------------------------- */
void readVectorDefocus (int iCol, string fileName, DMatrix &zTab, int nrows)
{
  fstream myfile;                          // zTab; read a single vector
  int i;

  myfile.open(fileName.c_str(), ios::in);
  if (myfile.is_open())
  {
    i=0;
    while (! myfile.eof())
    {
       myfile >> zTab(i,iCol);
       i++;
       if (i >= nrows) { break; }         // only read in nrows dimmension
    }
    myfile.close();
  }

  else { cout << "Unable to open file " << fileName << endl; }
}


/* --------------------------------------------------------------------- */
std::vector< std::vector<double> > readIn2dData(const char* filename)
{
    /* Function takes a char* filename argument and returns a
     * 2d dynamic array containing the data
     */

    std::vector< std::vector<double> > table;
    std::fstream ifs;

    /*  open file  */
    ifs.open(filename);

    if (ifs.is_open()) {

      while (true)
      {
          std::string line;
          double buf;

          getline(ifs, line);

          std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);

          if (!ifs)                  // mainly catch EOF
              break;

          if (line[0] == '#' || line.empty())   // catch empty lines or comment lines
              continue;

          std::vector<double> row;

          while (ss >> buf)
              row.push_back(buf);

          table.push_back(row);

      }
      ifs.close();

    }
    else {
      cout << "Unable to open file " << filename << endl;
      exit(1);
    }

    return table;
}


/* --------------------------------------------------------------------- */
void outputArr (int nArr, vector<double>& arr, string fileName)
{
  fstream myfile;                     // output a int array to file

  myfile.open(fileName.c_str(), ios::out | ios::binary);
  if (myfile.is_open())
  {
    myfile.width(20);            // format
    myfile.precision(5);
    myfile.setf( ios_base::scientific | ios_base::uppercase |
                 ios::showpos | ios_base::right);

    for (size_t i=0; i < nArr; i++) { myfile << arr[i] << "\n"; }

    myfile.close();
  }
  else { cout << "Failed to open file " << fileName << endl; }

  // cout << "\t" << arr[25] << endl;
}


/* --------------------------------------------------------------------- */
void outputIntArr (int nArr, vector<int>& intArr, string fileName)
{
  fstream myfile;                     // output a int array to file

  myfile.open(fileName.c_str(), ios::out | ios::binary);
  if (myfile.is_open())
  {
    for (size_t i=0; i < nArr; i++) { myfile << intArr[i] << "\n"; }

    myfile.close();
  }
  else { cout << "Failed to open file " << fileName << endl; }

  // cout << "\t" << intArr[25] << endl;
}


/* --------------------------------------------------------------------- */
void outputVector (int iCol, DMatrix &mat, string fileName)
{
  fstream myfile;                     // output a vector to file

  // myfile.open(fileName.c_str(), ios::out | ios::app | ios::binary);
  myfile.open(fileName.c_str(), ios::out | ios::binary);
  if (myfile.is_open())
  {
    // cout << "\t opened file " << fileName << " for output" << endl;

    myfile.width(20);            // format
    myfile.precision(5);
    myfile.setf( ios_base::scientific | ios_base::uppercase |
                 ios::showpos | ios_base::right);

    for (size_t i=0; i < mat.nrows(); i++) { myfile << mat(i,iCol) << "\n"; }

    myfile.close();
  }
  else { cout << "Failed to open file " << fileName << endl; }
}


/* --------------------------------------------------------------------- */
string convertInt(int number)
{                        // convert an integer to string
   stringstream ss;      // create a stringstream
   ss << number;         // add number to the stream
   return ss.str();      // return a string with the contents of the stream
}


/* --------------------------------------------------------------------- */
void outputVectorAll (DMatrix &mat, string nameBase)
{
  int i;                                  // output all vectors, each one
  string s,fileName;                      // a separate file. Format is
                                          // ${nameBase}${colNum}

  for (i=0; i<mat.ncols(); i++)
    {
       s = convertInt(i+1);
       fileName=nameBase+s;
       outputVector(i,mat,fileName);
    }
}


/* --------------------------------------------------------------------- */
void outputToFile (const DMatrix& mat, string fileName)
{                                     // write every column to a line
  fstream myfile;                     // output matrix to file

  // myfile.open("example.txt", ios::out | ios::app | ios::binary);
  myfile.open(fileName.c_str(), ios::out | ios::binary);
       // single or double quote makes difference
       // single quote for a single character, e.g. 'a', it takes 1 byte;
       // double quote for a string, e.g. "abc", there is a termination byte
       // in the end.
       // 'a' and "a" differs by the termination byte.
  if (myfile.is_open())
  {
    cout << "\t opened file " << fileName << " for output" << endl;

    myfile.width(20);            // format
    myfile.precision(5);
    myfile.setf( ios_base::scientific | ios_base::uppercase |
                 ios::showpos | ios_base::right);

    for (size_t j=0; j < mat.ncols(); j++)
      {
        for (size_t i=0; i < mat.nrows(); i++) { myfile << mat(i,j) << "  "; }
        myfile << "\n";
      }

    myfile.close();
  }
  else { cout << "Failed to open file " << fileName << endl; }
}




/* --------------------------------------------------------------------- */
void outputToFileTMVint (tmv::Matrix<int> & mat, string fileName)
{                                     // write every column to a line
  fstream myfile;                     // output matrix to file

  myfile.open(fileName.c_str(), ios::out | ios::binary);

  if (myfile.is_open())
  {
    cout << "\t opened file " << fileName << " for output" << endl;

    myfile.width(5);            // format

    for (size_t j=0; j < mat.ncols(); j++)
      {
        for (size_t i=0; i < mat.nrows(); i++) { myfile << mat(i,j) << "  "; }
        myfile << "\n";
      }

    myfile.close();
  }
  else { cout << "Failed to open file " << fileName << endl; }
}


/* --------------------------------------------------------------------- */
bool FileExists(string strFilename) { 
  /* source: http://www.techbytes.ca/techbyte103.html
     need to include sys/stat.h
   */
  struct stat stFileInfo; 
  bool blnReturn; 
  int intStat; 

  // Attempt to get the file attributes 
  intStat = stat(strFilename.c_str(),&stFileInfo); 
  if(intStat == 0) { 
    // We were able to get the file attributes 
    // so the file obviously exists. 
    blnReturn = true; 
  } else { 
    // We were not able to get the file attributes. 
    // This may mean that we don't have permission to 
    // access the folder which contains this file. If you 
    // need to do that level of checking, lookup the 
    // return values of stat which will give you 
    // more details on why stat failed. 
    blnReturn = false; 
  } 
   
  return(blnReturn); 
}
